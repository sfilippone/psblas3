
!
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone
!        Alfredo Buttari
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!
!


program psb_df_sample
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use getp
  use mpi
  implicit none

  ! input parameters
  character(len=40) :: kmethd, ptype, mtrx_file, rhs_file,renum

  ! sparse matrices
  type(psb_dspmat_type) :: a, aux_a

  ! preconditioner data
  type(psb_dprec_type)  :: prec

  ! dense matrices
  real(psb_dpk_), allocatable, target ::  aux_b(:,:), d(:)
  real(psb_dpk_), allocatable , save  :: x_col_glob(:), r_col_glob(:)
  real(psb_dpk_), pointer  :: b_col_glob(:)
  type(psb_d_vect_type)    :: b_col, x_col, r_col

  ! communications data structure
  type(psb_desc_type):: desc_a

  integer(psb_ipk_) :: ictxt, me, np

  ! solver paramters
  integer(psb_ipk_) :: iter, itmax, ierr, itrace, ircode,&
       & methd, istopc, irst
  integer(psb_epk_) :: amatsize, precsize, descsize
  real(psb_dpk_)   :: err, eps, cond

  character(len=5)   :: afmt
  character(len=20)  :: name, part
  character(len=2)   :: filefmt
  integer(psb_ipk_), parameter :: iunit=12
  integer(psb_ipk_) :: iparm(20)

  ! other variables
  integer(psb_ipk_) :: i,info,j,m_problem, err_act
  integer(psb_ipk_) :: internal, m,ii,nnzero
  real(psb_dpk_) :: t1, t2, tprec, t_low, t_high, t_mean, t_sum
  real(psb_dpk_) :: r_amax, b_amax, scale,resmx,resmxp
  integer(psb_ipk_) :: nrhs, nrow, n_row, dim, ne, nv, ncol
  integer(psb_ipk_), allocatable :: ivg(:), perm(:)
  integer(psb_ipk_), allocatable :: ipv(:)
  character(len=40)  :: fname, fnout

  ! persistant communication
  integer(psb_lpk_), allocatable :: iv(:)
  real(psb_dpk_), allocatable :: xa(:)
  integer(psb_lpk_) :: num_iterations

  call psb_init(ictxt)
  call psb_info(ictxt,me,np)

  if (me < 0) then
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif


  name='psb_df_sample'
  if(psb_errstatus_fatal()) goto 9999
  info=psb_success_
  call psb_set_errverbosity(itwo)
  !
  ! Hello world
  !
  if (me == psb_root_) then
    write(*,*) 'Welcome to PSBLAS version: ',psb_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if
  !
  !  get parameters
  !
  call get_parms(ictxt,mtrx_file,rhs_file,filefmt,kmethd,ptype,&
       & part,afmt,istopc,itmax,itrace,irst,eps)

  call psb_barrier(ictxt)
  t1 = psb_wtime()
  ! read the input matrix to be processed and (possibly) the rhs
  nrhs = 1

  if (me == psb_root_) then
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
      call psb_abort(ictxt)
    end if

    m_problem = aux_a%get_nrows()
    call psb_bcast(ictxt,m_problem)
    call psb_mat_renum(psb_mat_renum_identity_,aux_a,info,perm)

    ! At this point aux_b may still be unallocated
    if (size(aux_b,dim=1) == m_problem) then
      ! if any rhs were present, broadcast the first one
      write(psb_err_unit,'("Ok, got an rhs ")')
      b_col_glob =>aux_b(:,1)
      call psb_gelp('N',perm(1:m_problem),&
           & b_col_glob(1:m_problem),info)
    else
      write(psb_out_unit,'("Generating an rhs...")')
      write(psb_out_unit,'(" ")')
      call psb_realloc(m_problem,1,aux_b,ircode)
      if (ircode /= 0) then
         call psb_errpush(psb_err_alloc_dealloc_,name)
         goto 9999
      endif

      b_col_glob => aux_b(:,1)
      do i=1, m_problem
         b_col_glob(i) = done
      enddo
    endif

  else
    call psb_bcast(ictxt,m_problem)

  end if

  ! switch over different partition types
  select case(psb_toupper(part))
  case('BLOCK')
    if (me == psb_root_) write(psb_out_unit,'("Partition type: block")')
    call psb_matdist(aux_a, a,  ictxt,desc_a,info,fmt=afmt,parts=part_block)
    ! call desc_a%cd_get_list
    ! call desc_a%cd_v_get_list
    ! print *, "desc_a%v_halo_index", desc_a%v_get_list()


  case('GRAPH')
    if (me == psb_root_) then
      write(psb_out_unit,'("Partition type: graph vector")')
      write(psb_out_unit,'(" ")')
      !      write(psb_err_unit,'("Build type: graph")')
      call build_mtpart(aux_a,np)

    endif
    call psb_barrier(ictxt)
    call distr_mtpart(psb_root_,ictxt)
    call getv_mtpart(ivg)
    call psb_matdist(aux_a, a, ictxt,desc_a,info,fmt=afmt,v=ivg)

  case default
    if (me == psb_root_) write(psb_out_unit,'("Partition type: block")')
    call psb_matdist(aux_a, a,  ictxt,desc_a,info,fmt=afmt,parts=part_block)

  end select

  ! call psb_scatter(b_col_glob,b_col,desc_a,info,root=psb_root_)
  call psb_geall(x_col,desc_a,info)

  call x_col%zero() ! change to mpi rank? Artless
  !
  nrow=desc_a%get_local_rows() ! number of rows
  ncol=desc_a%get_local_cols() ! number of columns
  iv = [ (i,i=1,ncol)]

!  call psb_loc_to_glob(iv,desc_a,info) !
  call desc_a%l2gv1(iv,info)
  xa = iv
  ! call desc_a%locgv1(iv,info)

  ! ---- test where halo regions are set to -1, sending 10+me

  ! xa(nrow+1:ncol) = -1 ! (10 + me) * -1
  ! ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
  ! ! do i=0,np
  ! !   if (me .eq. i) then
  ! !     print *, "======================="
  ! !     print *, "me = ", me, "iv ="
  ! !     print '(225 I4)', iv
  ! !     print *, "-----------------------"
  ! !   end if
  ! !   call MPI_Barrier(MPI_COMM_WORLD, ierr)
  ! ! end do
  ! ! call MPI_Barrier(MPI_COMM_WORLD, ierr)

  ! ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
  ! ! do i=0,np
  ! !   if (me .eq. i) then
  ! !     print *, "======================="
  ! !     print *, "me = ", me, "xa ="
  ! !     print '(225 F4.0)', xa
  ! !     print *, "-----------------------"
  ! !   end if
  ! !   call MPI_Barrier(MPI_COMM_WORLD, ierr)
  ! ! end do
  ! ! call MPI_Barrier(MPI_COMM_WORLD, ierr)

  ! call x_col%set_vect(xa)
  ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
  ! call psb_halo(x_col,desc_a,info,mode=psb_swap_persistent_)

  ! xa = x_col%get_vect()
  ! do i=nrow+1,ncol
  !   if (xa(i) /= iv(i)) then ! ha
  !     write(0,*) me, ': MISMATCH i=',i,"xa(i)=",xa(i),"iv(i)",iv(i)
  !     ! goto 9999
  !   end if
  ! end do
  num_iterations = ITERATIONS
  do iter = 1, num_iterations
    xa = iv
    xa(nrow+1:ncol) = -1
    call x_col%set_vect(xa)
    call psb_halo(x_col,desc_a,info,mode=psb_swap_persistent_)
    xa = x_col%get_vect()

    do i=nrow+1,ncol
      if (xa(i) /= iv(i)) then ! ha
        write(0,*) me, ': MISMATCH i=',i,"xa(i)=",xa(i),"iv(i)",iv(i)
        goto 9999
      end if
    end do
  end do

  ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
  ! do i=0,np
  !   if (me .eq. i) then
  !     print *, "======================="
  !     print *, "me = ", me, " xa after halo"
  !     print '(225 F4.0)', xa
      ! print *, "me = ", me, " iv after halo"
      ! print '(225 I4)', iv
  !     print *, "-----------------------"
  !   end if
  !   call MPI_Barrier(MPI_COMM_WORLD, ierr)
  ! end do
  ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
  ! goto 9999




  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call flush(output_unit)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (allocated(x_col%v)) then
    print *, me, "desc_a%p%comm_create_time = ", x_col%v%p%comm_create_time
    call MPI_Reduce(x_col%v%p%comm_create_time, t_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
        psb_root_, MPI_COMM_WORLD, ierr)
    if (me == psb_root_) then
      open(1,file='halo_persistant_output.txt', status='new')

      write(1,*) "flt     type max_buf[B] visits time[s] time[%] time/visit[us]  region"
      t_mean = t_sum / np
      write(1,*) "MPI",  t_sum, "  ", t_mean
    end if
  end if
  call flush(output_unit)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (me == psb_root_) then
    print *, "did halo communication", num_iterations, "times"
  end if

  ! artless: commenting these because geasb is breaking it for now
  ! call psb_geasb(x_col,desc_a,info)
  ! call psb_geall(r_col,desc_a,info)
  ! call r_col%zero()
  ! call psb_geasb(r_col,desc_a,info)
  ! t2 = psb_wtime() - t1


  ! call psb_amx(ictxt, t2)

  ! if (me == psb_root_) then
  !    write(psb_out_unit,'(" ")')
  !    write(psb_out_unit,'("Time to read and partition matrix : ",es12.5)')t2
  !    write(psb_out_unit,'(" ")')
  ! end if

  !

  ! call perc%init(ictxt,ptype,info)

  ! building the preconditioner
  ! t1 = psb_wtime()
  ! call prec%build(a,desc_a,info)
  ! tprec = psb_wtime()-t1
  ! if (info /= psb_success_) then
  !    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_precbld')
  !    goto 9999
  ! end if


  ! call psb_amx(ictxt,tprec)

  ! if(me == psb_root_) then
  !    write(psb_out_unit,'("Preconditioner time: ",es12.5)')tprec
  !    write(psb_out_unit,'(" ")')
  ! end if

  ! cond = dzero
  ! iparm = 0
  ! call psb_barrier(ictxt)
  ! t1 = psb_wtime()
  ! call psb_krylov(kmethd,a,prec,b_col,x_col,eps,desc_a,info,&
  !      & itmax=itmax,iter=iter,err=err,itrace=itrace,&
  !      & istop=istopc,irst=irst,cond=cond)
  ! call psb_barrier(ictxt)
  ! t2 = psb_wtime() - t1
  ! call psb_amx(ictxt,t2)
  ! call psb_geaxpby(done,b_col,dzero,r_col,desc_a,info)
  ! call psb_spmm(-done,a,x_col,done,r_col,desc_a,info)
  ! resmx  = psb_genrm2(r_col,desc_a,info)
  ! resmxp = psb_geamax(r_col,desc_a,info)

  ! amatsize = a%sizeof()
  ! descsize = desc_a%sizeof()
  ! precsize = prec%sizeof()
  ! call psb_sum(ictxt,amatsize)
  ! call psb_sum(ictxt,descsize)
  ! call psb_sum(ictxt,precsize)

  ! if (me == psb_root_) then
  !   call prec%descr()
  !   write(psb_out_unit,'("Matrix: ",a)')mtrx_file
  !   write(psb_out_unit,'("Computed solution on ",i8," processors")')np
  !   write(psb_out_unit,'("Iterations to convergence: ",i6)')iter
  !   write(psb_out_unit,'("Error estimate on exit   : ",es12.5)') err
  !   write(psb_out_unit,'("Time to buil prec.       : ",es12.5)')tprec
  !   write(psb_out_unit,'("Time to solve system     : ",es12.5)')t2
  !   write(psb_out_unit,'("Time per iteration       : ",es12.5)')t2/(iter)
  !   write(psb_out_unit,'("Total time               : ",es12.5)')t2+tprec
  !   write(psb_out_unit,'("Residual norm 2          : ",es12.5)')resmx
  !   write(psb_out_unit,'("Residual norm inf        : ",es12.5)')resmxp
  !   write(psb_out_unit,'("Condition number         : ",es12.5)')cond
  !   write(psb_out_unit,'("Total memory occupation for A:      ",i12)')amatsize
  !   write(psb_out_unit,'("Total memory occupation for PREC:   ",i12)')precsize
  !   write(psb_out_unit,'("Total memory occupation for DESC_A: ",i12)')descsize
  !   write(psb_out_unit,'("Storage format for A              : ",a)')&
  !        &  a%get_fmt()
  !   write(psb_out_unit,'("Storage format for DESC_A         : ",a)')&
  !        &  desc_a%get_fmt()
  ! end if

  ! call psb_gather(x_col_glob,x_col,desc_a,info,root=psb_root_)
  ! if (info == psb_success_) &
  !      & call psb_gather(r_col_glob,r_col,desc_a,info,root=psb_root_)
  ! if (info /= psb_success_) goto 9999

  ! if (me == psb_root_) then
  !   write(psb_err_unit,'(" ")')
  !   write(psb_err_unit,'("Saving x on file")')
  !   write(20,*) 'matrix: ',mtrx_file
  !   write(20,*) 'computed solution on ',np,' processors.'
  !   write(20,*) 'iterations to convergence: ',iter
  !   write(20,*) 'error estimate (infinity norm) on exit:', &
  !        & ' ||r||/(||a||||x||+||b||) = ',err
  !   write(20,'("Residual norm 2          : ",es12.5)')resmx
  !   write(20,'("Residual norm inf        : ",es12.5)')resmxp
  !   write(20,'(a8,4(2x,a20))') 'I','X(I)','R(I)','B(I)'
  !   do i=1,m_problem
  !     write(20,998) i,x_col_glob(i),r_col_glob(i),b_col_glob(i)
  !   enddo
  ! end if
998 format(i8,4(2x,g20.14))
993 format(i6,4(1x,e12.6))

  ! print *, "* calling free funcs"
  ! call psb_gefree(b_col, desc_a,info)
  call psb_gefree(x_col, desc_a,info)
  call psb_spfree(a, desc_a,info)
  call prec%free(info)
  call psb_cdfree(desc_a,info)
  call psb_exit(ictxt)
  print *, "* FIN for", me
  stop

9999 call psb_error(ictxt)

  stop
end program psb_df_sample
