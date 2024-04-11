program psb_dbf_sample

   use psb_base_mod
   use psb_prec_mod
   use psb_krylov_mod
   use psb_util_mod
   use psb_cuda_mod
   use getp

   implicit none

   ! input parameters
   character(len=40) :: mtrx_file = "pde80.mtx"
   character(len=2)  :: filefmt   = "MM"
   character(len=40) :: kmethd 	  = "GMRES"
   character(len=40) :: ptype     = "NONE"
   character(len=20) :: part	  = "GRAPH"
   character(len=5)  :: afmt      = "CSR"
   character(len=5)  :: agfmt     = "ELG"
   integer(psb_ipk_) :: nrhs      = 6
   integer(psb_ipk_) :: istopbg   = 1
   integer(psb_ipk_) :: istoprg   = 2
   integer(psb_ipk_) :: itmax     = 500
   integer(psb_ipk_) :: itrace	  = -1
   integer(psb_ipk_) :: itrs	  = 30
   real(psb_dpk_)    :: eps	      = 1.d-7

   ! sparse matrices
   type(psb_dspmat_type)  :: a, agpu
   type(psb_ldspmat_type) :: aux_a

   ! preconditioner data
   type(psb_dprec_type) :: prec

   ! dense matrices
   real(psb_dpk_), allocatable, target :: aux_b(:,:)
   real(psb_dpk_), allocatable, save   :: x_mv_glob(:,:), r_mv_glob(:,:)
   real(psb_dpk_), allocatable, save   :: x_col_glob(:), r_col_glob(:)
   real(psb_dpk_), pointer             :: b_mv_glob(:,:)
   real(psb_dpk_), pointer  		   :: b_col_glob(:)
   type(psb_d_multivect_type)          :: b_mv, x_mv, r_mv
   type(psb_d_vect_type)               :: b_col, x_col, r_col
   type(psb_d_multivect_cuda)          :: gpumold_mv
   type(psb_d_vect_cuda)               :: gpumold_col
   type(psb_i_vect_cuda)               :: imold
   integer(psb_ipk_)                   :: m
   real(psb_dpk_)                      :: random_value

   ! molds
   type(psb_d_cuda_csrg_sparse_mat), target  :: acsrg
   type(psb_d_cuda_hlg_sparse_mat), target   :: ahlg
   type(psb_d_cuda_elg_sparse_mat), target   :: aelg
   class(psb_d_base_sparse_mat), pointer     :: agmold

   ! communications data structure
   type(psb_desc_type) :: desc_a
   type(psb_ctxt_type) :: ctxt
   integer(psb_ipk_)   :: iam, np
   integer(psb_lpk_)   :: lnp

   ! solver paramters
   integer(psb_ipk_) :: iter, ierr, ircode, methd, istopc
   integer(psb_epk_) :: amatsize, precsize, descsize
   real(psb_dpk_)    :: err, cond

   ! other variables
   character(len=20)              :: name
   integer(psb_ipk_)              :: reps = 50
   integer(psb_ipk_)              :: i, j, info, rep
   real(psb_dpk_)                 :: tb1, tb2, tr1, tr2
   real(psb_dpk_), allocatable    :: resmx(:), tb(:), tr(:)
   real(psb_dpk_)                 :: resmxp
   integer(psb_ipk_), allocatable :: ivg(:)
   integer(psb_ipk_), parameter   :: iunit=12
   logical                        :: print_matrix = .true.

   call psb_init(ctxt)
   call psb_info(ctxt,iam,np)
   call psb_cuda_init(ctxt)
   if (iam < 0) then
      ! This should not happen, but just in case
      call psb_exit(ctxt)
      stop
   endif

   name='psb_perf_test'
   if(psb_errstatus_fatal()) goto 9999
   info=psb_success_
   call psb_set_errverbosity(itwo)
   !
   ! Hello world
   !
   if (iam == psb_root_) then
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Welcome to PSBLAS version: ",a)') psb_version_string_
      write(psb_out_unit,'("This is the ",a," sample program")') trim(name)
      write(psb_out_unit,'("Number of processors: ",i8)')np
      write(psb_out_unit,'(" ")')
   end if
   write(*,*) 'Process ',iam,' running on device: ', psb_cuda_getDevice(),' out of', psb_cuda_getDeviceCount()
   write(*,*) 'Process ',iam,' device ', psb_cuda_getDevice(),' is a: ', trim(psb_cuda_DeviceName())

   select case(psb_toupper(agfmt))
    case('ELG')
      agmold => aelg
    case('HLG')
      agmold => ahlg
    case('CSRG')
      agmold => acsrg
    case default
      write(*,*) 'Unknown format defaulting to CSRG'
      agmold => acsrg
   end select

   if (iam == psb_root_) then
      select case(psb_toupper(filefmt))
       case('MM')
         call mm_mat_read(aux_a,info,iunit=iunit,filename=mtrx_file)
       case ('HB')
         call hb_read(aux_a,info,iunit=iunit,b=aux_b,filename=mtrx_file)
       case default
         info = -1
         write(psb_err_unit,*) 'Wrong choice for fileformat ', filefmt
      end select
      if (info /= psb_success_) then
         write(psb_err_unit,*) 'Error while reading input matrix '
         call psb_abort(ctxt)
      end if

      m = aux_a%get_nrows()
      call psb_bcast(ctxt,m)

      ! At this point aux_b may still be unallocated
      if (size(aux_b) == m*nrhs) then
         ! if any rhs were present, broadcast the first one
         write(psb_err_unit,'("Ok, got an rhs ")')
         b_mv_glob =>aux_b(:,:)
      else
         write(psb_out_unit,'("Generating an rhs...")')
         write(psb_out_unit,'("Number of RHS: ",i3)') nrhs
         write(psb_out_unit,'(" ")')
         call psb_realloc(m,nrhs,aux_b,ircode)
         if (ircode /= 0) then
            call psb_errpush(psb_err_alloc_dealloc_,name)
            goto 9999
         endif

         b_mv_glob => aux_b(:,:)
         do i=1, m
            do j=1, nrhs
               b_mv_glob(i,j) = done
               !call random_number(random_value)
               !b_mv_glob(i,j) = random_value
            enddo
         enddo
      endif
   else
      call psb_bcast(ctxt,m)
   end if

   ! switch over different partition types
   select case(psb_toupper(part))
    case('BLOCK')
      if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
      call psb_matdist(aux_a,a,ctxt,desc_a,info,fmt=afmt,parts=part_block)
    case('GRAPH')
      if (iam == psb_root_) then
         write(psb_out_unit,'("Partition type: graph vector")')
         write(psb_out_unit,'(" ")')
         call aux_a%cscnv(info,type='csr')
         lnp = np
         call build_mtpart(aux_a,lnp)
      endif
      call psb_barrier(ctxt)
      call distr_mtpart(psb_root_,ctxt)
      call getv_mtpart(ivg)
      call psb_matdist(aux_a,a,ctxt,desc_a,info,fmt=afmt,vg=ivg)
    case default
      if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
      call psb_matdist(aux_a,a,ctxt,desc_a,info,fmt=afmt,parts=part_block)
   end select

   ! building the preconditioner
   call prec%init(ctxt,ptype,info)
   call prec%build(a,desc_a,info)
   if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_precbld')
      goto 9999
   end if

   ! Set GPU matrix
   call a%cscnv(agpu,info,mold=agmold)
   call desc_a%cnv(mold=imold)
   if ((info /= 0).or.(psb_get_errstatus()/=0)) then
      write(0,*) 'From cscnv ',info
      call psb_error()
      stop
   end if

   if (iam == psb_root_) then
      allocate(tb(reps),tr(reps))
      tb = dzero
      tr = dzero
   end if

   if(iam == psb_root_) then
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Starting BGMRES")')
      write(psb_out_unit,'(" ")')
   end if

   do rep=1,reps
      call psb_scatter(b_mv_glob,b_mv,desc_a,info,root=psb_root_,mold=gpumold_mv)
      call psb_geall(x_mv,desc_a,info,nrhs)
      call x_mv%zero()
      call psb_geasb(x_mv,desc_a,info,mold=gpumold_mv)
      call psb_geall(r_mv,desc_a,info,nrhs)
      call r_mv%zero()
      call psb_geasb(r_mv,desc_a,info)

      call psb_barrier(ctxt)
      tb1 = psb_wtime()

      call psb_krylov(kmethd,a,prec,b_mv,x_mv,eps,desc_a,info,&
      & itmax=itmax,iter=iter,err=err,itrace=itrace,&
      & itrs=itrs,istop=istopbg)

      call psb_barrier(ctxt)
      tb2 = psb_wtime() - tb1
      call psb_amx(ctxt,tb2)

      if (iam == psb_root_) then
         tb(rep) = tb2
         write(*,*) 'Time', rep, tb(rep)
      end if
      call psb_barrier(ctxt)

      call psb_gefree(b_mv,desc_a,info)
      call psb_gefree(x_mv,desc_a,info)
   end do

   if(iam == psb_root_) then
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Finished BGMRES")')
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Starting sGMRES")')
      write(psb_out_unit,'(" ")')
   end if

   call psb_barrier(ctxt)

   do rep=1,reps
      do i=1,nrhs

         b_col_glob => aux_b(:,i)
         call psb_scatter(b_col_glob,b_col,desc_a,info,root=psb_root_,mold=gpumold_col)
         call psb_geall(x_col,desc_a,info)
         call x_col%zero()
         call psb_geasb(x_col,desc_a,info,mold=gpumold_col)
         call psb_geall(r_col,desc_a,info)
         call r_col%zero()
         call psb_geasb(r_col,desc_a,info)

         cond = dzero

         call psb_barrier(ctxt)
         tr1 = psb_wtime()

         call psb_krylov(kmethd,a,prec,b_col,x_col,eps,desc_a,info,&
         & itmax=itmax,iter=iter,err=err,itrace=itrace,&
         & istop=istoprg,irst=itrs,cond=cond)

         call psb_barrier(ctxt)
         tr2 = psb_wtime() - tr1
         call psb_amx(ctxt,tr2)

         if (iam == psb_root_) then
            tr(rep) = tr(rep) + tr2
         end if
         call psb_barrier(ctxt)

         call psb_gefree(b_col, desc_a,info)
         call psb_gefree(x_col, desc_a,info)
      end do

      if (iam == psb_root_) then
        write(*,*) 'Time', rep, tr(rep)
      end if
   end do

   call psb_barrier(ctxt)

   if(iam == psb_root_) then
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Finished sGMRES")')
      write(psb_out_unit,'(" ")')
   end if

   call psb_spfree(a, desc_a,info)
   call prec%free(info)
   call psb_cdfree(desc_a,info)
   call psb_cuda_exit()
   call psb_exit(ctxt)

   return

   call psb_geaxpby(done,b_mv,dzero,r_mv,desc_a,info)
   call psb_spmm(-done,a,x_mv,done,r_mv,desc_a,info)

   resmx  = psb_genrm2(r_mv,desc_a,info)
   resmxp = psb_geamax(r_mv,desc_a,info)

   amatsize = a%sizeof()
   descsize = desc_a%sizeof()
   precsize = prec%sizeof()

   call psb_sum(ctxt,amatsize)
   call psb_sum(ctxt,descsize)
   call psb_sum(ctxt,precsize)

   call psb_gather(x_mv_glob,x_mv,desc_a,info,root=psb_root_)
   if (info == psb_success_) call psb_gather(r_mv_glob,r_mv,desc_a,info,root=psb_root_)
   if (info /= psb_success_) goto 9999

   if (iam == psb_root_) then
      call prec%descr(info)
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Computed solution on:         ",i8," processors")')np
      write(psb_out_unit,'("Matrix:                              ",a)')mtrx_file
      write(psb_out_unit,'("Storage format for A:                ",a)')a%get_fmt()
      write(psb_out_unit,'("Storage format for DESC_A:           ",a)')desc_a%get_fmt()
      write(psb_out_unit,'("Total memory occupation for A:      ",i12)')amatsize
      write(psb_out_unit,'("Total memory occupation for PREC:   ",i12)')precsize
      write(psb_out_unit,'("Total memory occupation for DESC_A: ",i12)')descsize
      write(psb_out_unit,'("Iterations to convergence:          ",i12)')iter
      write(psb_out_unit,'("Error estimate on exit:             ",es12.5)')err
      !write(psb_out_unit,'("Time to buil prec.:                 ",es12.5)')tprec
      ! write(psb_out_unit,'("Time to solve system:               ",es12.5)')t2
      !write(psb_out_unit,'("Time per iteration:                 ",es12.5)')t2/(iter)
      !write(psb_out_unit,'("Total time:                         ",es12.5)')t2+tprec
      write(psb_out_unit,'("Residual norm 2:                    ",es12.5)')maxval(resmx)
      write(psb_out_unit,'("Residual norm inf:                  ",es12.5)')resmxp
      write(psb_out_unit,'(a8,4(2x,a20))') 'I','X(I)','R(I)','B(I)'
      if (print_matrix) then
         do i=1,m
            write(psb_out_unit,993) i, x_mv_glob(i,:), r_mv_glob(i,:), b_mv_glob(i,:)
         end do
      end if
   end if

998 format(i8,4(2x,g20.14))
993 format(i6,4(1x,e12.6))

   call psb_gefree(b_mv,desc_a,info)
   call psb_gefree(x_mv,desc_a,info)
   call psb_spfree(a,desc_a,info)
   call prec%free(info)
   call psb_cdfree(desc_a,info)
   call psb_cuda_exit()
   call psb_exit(ctxt)

   return

9999 call psb_error(ctxt)

   return

end program psb_dbf_sample
