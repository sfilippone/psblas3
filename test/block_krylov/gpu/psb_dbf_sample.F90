program psb_dbf_sample

   use psb_base_mod
   use psb_prec_mod
   use psb_krylov_mod
   use psb_util_mod
   use psb_cuda_mod

   implicit none

   ! input parameters
   character(len=40) :: kmethd, ptype, mtrx_file, rhs_file

   ! sparse matrices
   type(psb_dspmat_type)  :: a, agpu
   type(psb_ldspmat_type) :: aux_a

   ! preconditioner data
   type(psb_dprec_type) :: prec

   ! dense matrices
   real(psb_dpk_), allocatable, target :: aux_b(:,:)
   real(psb_dpk_), allocatable, save   :: x_mv_glob(:,:), r_mv_glob(:,:)
   real(psb_dpk_), pointer             :: b_mv_glob(:,:)
   type(psb_d_multivect_type)          :: b_mv, x_mv, r_mv
   type(psb_d_multivect_cuda)          :: gpumold
   type(psb_i_vect_cuda)               :: imold
   integer(psb_ipk_)                   :: m, nrhs
   real(psb_dpk_)                      :: random_value

   ! molds
   type(psb_d_cuda_csrg_sparse_mat), target  :: acsrg
   type(psb_d_cuda_hlg_sparse_mat), target   :: ahlg
   type(psb_d_cuda_elg_sparse_mat), target   :: aelg
   class(psb_d_base_sparse_mat), pointer     :: amold

   ! communications data structure
   type(psb_desc_type) :: desc_a
   type(psb_ctxt_type) :: ctxt
   integer(psb_ipk_)   :: iam, np
   integer(psb_lpk_)   :: lnp

   ! solver paramters
   integer(psb_ipk_) :: iter, itmax, ierr, itrace, ircode, methd, istopc
   integer(psb_epk_) :: amatsize, precsize, descsize
   real(psb_dpk_)    :: err, eps

   ! input parameters
   character(len=5)             :: afmt
   character(len=20)            :: name, part
   character(len=2)             :: filefmt
   integer(psb_ipk_), parameter :: iunit=12

   ! other variables
   integer(psb_ipk_)              :: i, j, info
   real(psb_dpk_)                 :: t1, t2, tprec
   real(psb_dpk_), allocatable    :: resmx(:), res(:,:), temp(:,:), ttt(:,:)
   real(psb_dpk_)                 :: resmxp
   integer(psb_ipk_), allocatable :: ivg(:)
   logical                        :: print_matrix = .true.

   call psb_init(ctxt)
   call psb_info(ctxt,iam,np)
   call psb_cuda_init(ctxt)

   if (iam < 0) then
      ! This should not happen, but just in case
      call psb_exit(ctxt)
      stop
   endif

   name='psb_dbf_sample'
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
      write(psb_out_unit,'(" ")')
   end if
   write(*,*) 'Process ',iam,' running on device: ', psb_cuda_getDevice(),' out of', psb_cuda_getDeviceCount()
   write(*,*) 'Process ',iam,' device ', psb_cuda_getDevice(),' is a: ', trim(psb_cuda_DeviceName())
   !
   !  get parameters
   !
   call get_parms(ctxt,mtrx_file,rhs_file,filefmt,kmethd,ptype,&
   & part,afmt,nrhs,istopc,itmax,itrace,eps)

   select case(psb_toupper(afmt))
    case('ELG')
      amold => aelg
    case('HLG')
      amold => ahlg
    case('CSRG')
      amold => acsrg
    case default
      write(*,*) 'Unknown format defaulting to CSRG'
      amold => acsrg
   end select

   call psb_barrier(ctxt)
   t1 = psb_wtime()

   if (iam == psb_root_) then
      select case(psb_toupper(filefmt))
       case('MM')
         ! For Matrix Market we have an input file for the matrix
         ! and an (optional) second file for the RHS.
         call mm_mat_read(aux_a,info,iunit=iunit,filename=mtrx_file)
         if (info == psb_success_) then
            if (rhs_file /= 'NONE') then
               call mm_array_read(aux_b,info,iunit=iunit,filename=rhs_file)
            end if
         end if

       case ('HB')
         ! For Harwell-Boeing we have a single file which may or may not
         ! contain an RHS.
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
      call psb_matdist(aux_a,a,ctxt,desc_a,info,fmt='CSR  ',parts=part_block)

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
      call psb_matdist(aux_a,a,ctxt,desc_a,info,fmt='CSR  ',vg=ivg)

    case default
      if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
      call psb_matdist(aux_a,a,ctxt,desc_a,info,fmt='CSR  ',parts=part_block)
   end select

   call a%cscnv(agpu,info,mold=amold)
   if ((info /= 0).or.(psb_get_errstatus()/=0)) then
      write(0,*) 'From cscnv ',info
      call psb_error()
      stop
   end if
   call desc_a%cnv(mold=imold)

   call psb_scatter(b_mv_glob,b_mv,desc_a,info,root=psb_root_,mold=gpumold)
   call psb_geall(x_mv,desc_a,info,nrhs)
   call x_mv%zero()
   call psb_geasb(x_mv,desc_a,info,mold=gpumold)
   call psb_geall(r_mv,desc_a,info,nrhs)
   call r_mv%zero()
   call psb_geasb(r_mv,desc_a,info)

   t2 = psb_wtime() - t1
   call psb_amx(ctxt, t2)

   if (iam == psb_root_) then
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Time to read and partition matrix : ",es12.5)')t2
      write(psb_out_unit,'(" ")')
   end if

   ! building the preconditioner
   call prec%init(ctxt,ptype,info)
   t1 = psb_wtime()
   call prec%build(a,desc_a,info)
   tprec = psb_wtime()-t1
   if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_precbld')
      goto 9999
   end if

   call psb_amx(ctxt,tprec)

   if(iam == psb_root_) then
      write(psb_out_unit,'("Preconditioner time: ",es12.5)')tprec
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Starting algorithm")')
      write(psb_out_unit,'(" ")')
   end if

   call psb_barrier(ctxt)
   t1 = psb_wtime()

   call psb_krylov(kmethd,agpu,prec,b_mv,x_mv,eps,desc_a,info,&
   & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istopc)

   call psb_barrier(ctxt)
   t2 = psb_wtime() - t1
   call psb_amx(ctxt,t2)

   if(iam == psb_root_) then
      write(psb_out_unit,'("Finished algorithm")')
      write(psb_out_unit,'(" ")')
   end if

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
      write(psb_out_unit,'("Time to buil prec.:                 ",es12.5)')tprec
      write(psb_out_unit,'("Time to solve system:               ",es12.5)')t2
      write(psb_out_unit,'("Time per iteration:                 ",es12.5)')t2/(iter)
      write(psb_out_unit,'("Total time:                         ",es12.5)')t2+tprec
      write(psb_out_unit,'("Residual norm 2:                    ",es12.5)')sum(resmx)/size(resmx)
      write(psb_out_unit,'("Residual norm inf:                  ",es12.5)')resmxp
      if (print_matrix) then
         write(psb_out_unit,'(a8,4(2x,a20))') 'I','X(I)','R(I)','B(I)'
         do i=1,m
            write(psb_out_unit,993) i, x_mv_glob(i,:), r_mv_glob(i,:), b_mv_glob(i,:)
         end do
      end if
   end if

998 format(i8,10(2x,g20.14))
993 format(i6,10(1x,e12.6))

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

contains
   !
   ! get iteration parameters from standard input
   !
   subroutine get_parms(ctxt,mtrx_file,rhs_file,filefmt,kmethd,ptype,part,afmt,nrhs,istopc,itmax,itrace,eps)
      type(psb_ctxt_type) :: ctxt
      character(len=2)    :: filefmt
      character(len=40)   :: kmethd, mtrx_file, rhs_file, ptype
      character(len=20)   :: part
      character(len=5)    :: afmt
      integer(psb_ipk_)   :: nrhs, istopc, itmax, itrace
      real(psb_dpk_)      :: eps

      integer(psb_ipk_)   :: np, iam
      integer(psb_ipk_)   :: ip, inp_unit
      character(len=1024) :: filename

      call psb_info(ctxt, iam, np)

      if (iam == 0) then
         if (command_argument_count()>0) then
            call get_command_argument(1, filename)
            inp_unit = 30
            open(inp_unit, file=filename, action='read', iostat=info)
            if (info /= 0) then
               write(psb_err_unit,*) 'Could not open file ',filename,' for input'
               call psb_abort(ctxt)
               stop
            else
               write(psb_err_unit,*) 'Opened file ',trim(filename),' for input'
            end if
         else
            inp_unit = psb_inp_unit
         end if

         ! Read Input Parameters
         read(inp_unit,*) mtrx_file
         read(inp_unit,*) rhs_file
         read(inp_unit,*) filefmt
         read(inp_unit,*) kmethd
         read(inp_unit,*) ptype
         read(inp_unit,*) afmt
         read(inp_unit,*) part
         read(inp_unit,*) nrhs
         read(inp_unit,*) istopc
         read(inp_unit,*) itmax
         read(inp_unit,*) itrace
         read(inp_unit,*) eps

         call psb_bcast(ctxt,mtrx_file)
         call psb_bcast(ctxt,rhs_file)
         call psb_bcast(ctxt,filefmt)
         call psb_bcast(ctxt,kmethd)
         call psb_bcast(ctxt,ptype)
         call psb_bcast(ctxt,afmt)
         call psb_bcast(ctxt,part)
         call psb_bcast(ctxt,nrhs)
         call psb_bcast(ctxt,istopc)
         call psb_bcast(ctxt,itmax)
         call psb_bcast(ctxt,itrace)
         call psb_bcast(ctxt,eps)

         if (inp_unit /= psb_inp_unit) then
            close(inp_unit)
         end if
      else
         ! Receive Parameters
         call psb_bcast(ctxt,mtrx_file)
         call psb_bcast(ctxt,rhs_file)
         call psb_bcast(ctxt,filefmt)
         call psb_bcast(ctxt,kmethd)
         call psb_bcast(ctxt,ptype)
         call psb_bcast(ctxt,afmt)
         call psb_bcast(ctxt,part)
         call psb_bcast(ctxt,nrhs)
         call psb_bcast(ctxt,istopc)
         call psb_bcast(ctxt,itmax)
         call psb_bcast(ctxt,itrace)
         call psb_bcast(ctxt,eps)
      end if

      if (iam == 0) then
         write(psb_out_unit,'(" ")')
         write(psb_out_unit,'("Solving matrix           : ",a)') mtrx_file
         write(psb_out_unit,'("Data distribution        : ",a)') part
         write(psb_out_unit,'("Iterative method         : ",a)') kmethd
         write(psb_out_unit,'("Number of processors     : ",i0)') np
         write(psb_out_unit,'("Number of RHS            : ",i4)') nrhs
         write(psb_out_unit,'("Max number of iterations : ",i3)') itmax
         write(psb_out_unit,'("Storage format           : ",a)') afmt
         write(psb_out_unit,'(" ")')
      end if

      return

   end subroutine get_parms

end program psb_dbf_sample
