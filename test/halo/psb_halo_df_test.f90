
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
  real(psb_dpk_) :: t1, t2, tprec, t_low, t_high, ave_request_create_t, request_create_t
  real(psb_dpk_) :: r_amax, b_amax, scale,resmx,resmxp
  integer(psb_ipk_) :: nrhs, nrow, n_row, dim, ne, nv, ncol
  integer(psb_ipk_), allocatable :: ivg(:), perm(:)
  integer(psb_ipk_), allocatable :: ipv(:)
  character(len=40)  :: fname, fnout

  ! persistent communication
  integer(psb_lpk_), allocatable :: iv(:)
  real(psb_dpk_), allocatable :: xa(:)
  integer(psb_ipk_) :: num_iterations, num_neighbors, &
      min_neighbors, max_neighbors, &
      sum_neighbors
  integer(psb_ipk_)  :: sum_snd, min_snd, max_snd, num_snd, tot_snd, sum_tot_snd
  integer(psb_ipk_)  :: sum_rcv, min_rcv, max_rcv, num_rcv, tot_rcv, sum_tot_rcv
  real(psb_spk_)     :: ave_neighbors, ave_snd_buf, ave_rcv_buf, ave_tot_snd, ave_tot_rcv
  real(psb_dpk_)     :: alltoall_comm_t, ave_alltoall_comm_t, total_time, &
      ave_time_pp, ave_time_pi
  real(psb_dpk_)     :: median_comm_t, ave_alltoall_median_comm_t, min_comm_t, max_comm_t, &
      ave_min_comm_t, ave_max_comm_t, sum_median_comm_t, sum_min_comm_t, sum_max_comm_t, &
      ave_alltoall_max_comm_t, ave_alltoall_min_comm_t
  real(psb_dpk_)     :: halo_t_start, halo_t_end, ave_halo_t_pi, &
      sum_halo_t, halo_time
  real(psb_spk_)     :: ave_snd, ave_rcv
  integer(psb_ipk_)  :: swap_mode
  logical            :: swap_persistent, swap_nonpersistent, file_exists

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
  swap_mode      = SWAPCHOICE
  swap_persistent     = iand(swap_mode,psb_swap_persistent_) /= 0
  swap_nonpersistent  = iand(swap_mode,psb_swap_nonpersistent_) /= 0
  halo_time = 0
  do iter = 1, num_iterations
    xa = iv
    xa(nrow+1:ncol) = -1
    call x_col%set_vect(xa)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    halo_t_start = MPI_Wtime()
    call psb_halo(x_col,desc_a,info,mode=swap_mode)
    halo_t_end = MPI_Wtime() - halo_t_start
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    halo_time = halo_time + halo_t_end

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
    if (allocated(x_col%v%p)) then
    ! -- collect times ---
    ! halo time
    call MPI_Reduce(halo_time, sum_halo_t, 1, &
        MPI_DOUBLE_PRECISION, MPI_SUM, psb_root_, MPI_COMM_WORLD, ierr)

    ! total time
    call MPI_Reduce(x_col%v%p%total_time, total_time, 1, &
        MPI_DOUBLE_PRECISION, MPI_SUM, psb_root_, MPI_COMM_WORLD, ierr)

    if (swap_persistent .or. swap_nonpersistent) then
      ! collect communication times
      call median(x_col%v%p%alltoall_comm_time_a(1:x_col%v%p%comm_count), &
          median_comm_t, min_comm_t, max_comm_t)

      call MPI_Reduce(median_comm_t, sum_median_comm_t, 1, &
          MPI_DOUBLE_PRECISION, MPI_SUM, psb_root_, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(min_comm_t, sum_min_comm_t, 1, &
          MPI_DOUBLE_PRECISION, MPI_SUM, psb_root_, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(max_comm_t, sum_max_comm_t, 1, &
          MPI_DOUBLE_PRECISION, MPI_SUM, psb_root_, MPI_COMM_WORLD, ierr)

      call MPI_Reduce(x_col%v%p%alltoall_comm_time, alltoall_comm_t, 1, &
          MPI_DOUBLE_PRECISION, MPI_SUM, psb_root_, MPI_COMM_WORLD, ierr)
    end if
    if (swap_persistent) then
      ! collect MPIX request initialization time
      call MPI_Reduce(x_col%v%p%request_create_time, request_create_t, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
          psb_root_, MPI_COMM_WORLD, ierr)
    end if
    ! collect min, max, and average snd/rcv elements
    tot_snd = sum(x_col%v%p%snd_counts)
    call MPI_Reduce(tot_snd, sum_tot_snd, 1, MPI_INTEGER, MPI_SUM, &
        psb_root_, MPI_COMM_WORLD, ierr)
    ave_tot_snd = sum_tot_snd / np
    ave_snd = real(tot_snd) / real(x_col%v%p%snd_count)
    call MPI_Reduce(ave_snd, ave_snd_buf, 1, MPI_REAL, &
        MPI_SUM, psb_root_, MPI_COMM_WORLD, ierr)
    ave_snd_buf = ave_snd_buf / np
    print *, "snd_counts = ", x_col%v%p%snd_counts
    call MPI_Reduce(minval(x_col%v%p%snd_counts), min_snd, 1, MPI_INTEGER, MPI_MIN, &
        psb_root_, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(maxval(x_col%v%p%snd_counts), max_snd, 1, MPI_INTEGER, MPI_MAX, &
        psb_root_, MPI_COMM_WORLD, ierr)
    ! collect min, max, and average rcv elements
    tot_rcv = sum(x_col%v%p%rcv_counts)
    call MPI_Reduce(tot_rcv, sum_tot_rcv, 1, MPI_INTEGER, MPI_SUM, &
        psb_root_, MPI_COMM_WORLD, ierr)
    ave_tot_rcv = sum_tot_rcv / np
    ave_rcv = real(tot_rcv) / real(x_col%v%p%rcv_count)
    call MPI_Reduce(ave_rcv, ave_rcv_buf, 1, MPI_REAL, &
        MPI_SUM, psb_root_, MPI_COMM_WORLD, ierr)
    ave_rcv_buf = ave_rcv_buf / np
    call MPI_Reduce(minval(x_col%v%p%rcv_counts), min_rcv, 1, MPI_INTEGER, MPI_MIN, &
        psb_root_, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(maxval(x_col%v%p%rcv_counts), max_rcv, 1, MPI_INTEGER, MPI_MAX, &
        psb_root_, MPI_COMM_WORLD, ierr)

    ! collect min, max, and average number of neighbors
    num_neighbors = size(x_col%v%p%snd_to, dim=1)
    call MPI_Reduce(num_neighbors, sum_neighbors, 1, MPI_INTEGER, MPI_SUM, &
        psb_root_, MPI_COMM_WORLD, ierr)
    ave_neighbors = REAL(sum_neighbors) / np
    call MPI_Reduce(num_neighbors, min_neighbors, 1, MPI_INTEGER, MPI_MIN, &
        psb_root_, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(num_neighbors, max_neighbors, 1, MPI_INTEGER, MPI_MAX, &
        psb_root_, MPI_COMM_WORLD, ierr)

    if (me == psb_root_) then
      inquire(file='halo_output.txt', exist=file_exists)
      if (.not. file_exists) then
        open(1,file='halo_output.txt')
        write(1,'(A)',advance='no') "partition;"
        write(1,'(A)',advance='no') "np;"
        write(1,'(A)',advance='no') "num_iterations;"
        write(1,'(A)',advance='no') "swap_mode;"
        write(1,'(A)',advance='no') "ave_halo_t_pi;"
        ! write(1,'(A)',advance='no') "total_time;"
        ! write(1,'(A)',advance='no') "ave_time_pp;"
        write(1,'(A)',advance='no') "ave_time_pi;"
        write(1,'(A)',advance='no') "alltoall_comm_t;"
        write(1,'(A)',advance='no') "ave_alltoall_comm_t;"
        write(1,'(A)',advance='no') "ave_alltoall_median_comm_t;"
        write(1,'(A)',advance='no') "ave_alltoall_max_comm_t;"
        write(1,'(A)',advance='no') "ave_alltoall_min_comm_t;"
        write(1,'(A)',advance='no') "request_create_t;"
        write(1,'(A)',advance='no') "ave_request_create_t;"
        write(1,'(A)',advance='no') "ave_neighbors;"

        write(1,'(A)',advance='no') "min_neighbors;"
        write(1,'(A)',advance='no') "max_neighbors;"
        write(1,'(A)',advance='no') "ave_tot_snd;"
        write(1,'(A)',advance='no') "ave_tot_rcv;"
        write(1,'(A)',advance='no') "ave_snd_buf;"
        write(1,'(A)',advance='no') "ave_rcv_buf;"
        write(1,'(A)',advance='no') "max_snd;"
        write(1,'(A)',advance='no') "max_rcv;"
        write(1,'(A)',advance='no') "min_snd;"
        write(1,'(A)',advance='yes')"min_rcv;"
      else
        open(1,file='halo_output.txt', access='APPEND')
      end if


      ! write(1,*) "flt     type max_buf[B] visits time[s] time[%] time/visit[us]  region"
      ! write(1,*) ""
      ! converting microseconds
      total_time = total_time * 1000000
      ave_time_pp = total_time / np
      ave_time_pi = total_time / (np * num_iterations)
      sum_halo_t  = sum_halo_t * 1000000
      ave_halo_t_pi = sum_halo_t / (np * num_iterations)
      alltoall_comm_t = alltoall_comm_t * 1000000
      ave_alltoall_comm_t = alltoall_comm_t / (np * num_iterations)
      sum_median_comm_t = sum_median_comm_t * 1000000
      ave_alltoall_median_comm_t = sum_median_comm_t / np
      sum_max_comm_t = sum_max_comm_t * 1000000
      ave_alltoall_max_comm_t = sum_max_comm_t / np
      sum_min_comm_t = sum_min_comm_t * 1000000
      ave_alltoall_min_comm_t = sum_min_comm_t / np

      ! print *, "alltoall", alltoall_comm_t, ave_alltoall_comm_t
      if (swap_persistent) then
        ave_request_create_t = (request_create_t / np) * 1000000
        request_create_t  = request_create_t        * 1000000
      else if (swap_nonpersistent) then
        ave_request_create_t = -1.0
        request_create_t  = -1.0
      end if

      ! write(1,*, advance='no') 37 65  73 64

      write(1,'(A5 A1)',advance='no')    psb_toupper(part), ';'
      write(1,'(I5 A1)',advance='no')    np, ';'
      write(1,'(I6 A1)',advance='no')    num_iterations, ';'
      write(1,'(I5 A1)'  ,advance='no')  swap_mode, ';'
      write(1,'(I5 A1)'  ,advance='no')  ave_halo_t_pi, ';'
      ! write(1,'(F20.4 A1)',advance='no') total_time, ';'
      ! write(1,'(F20.4 A1)',advance='no') ave_time_pp, ';'
      write(1,'(F20.4 A1)',advance='no') ave_time_pi, ';'
      write(1,'(F20.4 A1)',advance='no') alltoall_comm_t, ';'
      write(1,'(F20.4 A1)',advance='no') ave_alltoall_comm_t, ';'
      write(1,'(F20.4 A1)',advance='no') ave_alltoall_median_comm_t, ';'
      write(1,'(F20.4 A1)',advance='no') ave_alltoall_max_comm_t, ';'
      write(1,'(F20.4 A1)',advance='no') ave_alltoall_min_comm_t, ';'
      write(1,'(F10.2 A1)',advance='no') request_create_t,  ';'
      write(1,'(F10.2 A1)',advance='no') ave_request_create_t, ';'
      write(1,'(F10.2 A1)',advance='no') ave_neighbors, ';'

      write(1,'(I5 A1)'  ,advance='no') min_neighbors, ';'
      write(1,'(I5 A1)'  ,advance='no') max_neighbors, ';'
      write(1,'(F15.2 A1)',advance='no') ave_tot_snd, ';'
      write(1,'(F15.2 A1)',advance='no') ave_tot_rcv, ';'
      write(1,'(F15.2 A1)',advance='no') ave_snd_buf, ';'
      write(1,'(F15.2 A1)',advance='no') ave_rcv_buf, ';'
      write(1,'(I7 A1)',advance='no')   max_snd, ';'
      write(1,'(I7 A1)',advance='no')   max_rcv, ';'
      write(1,'(I7 A1)',advance='no')   min_snd, ';'
      write(1,'(I7 A1)',advance='no')   min_rcv

      ! write(1,'(F9.2 A1)',advance='no') , ';'
      ! write(1,'(I5 A1)',advance='no') , ';'
      write(1,'()')

      ! ---- need to write buffer size, min and average and max snd/rcv elements
      ! write(1,())
      ! ---- need to write buffer size, min and average and max neighbors
    end if
    end if
  end if
  call flush(output_unit)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (me == psb_root_) then
    print *, "did halo communication", num_iterations, "times"
  end if

998 format(i8,4(2x,g20.14))
993 format(i6,4(1x,e12.6))

  ! print *, "* calling free funcs"
  ! call psb_gefree(b_col, desc_a,info)
  call psb_gefree(x_col, desc_a,info)
  call psb_spfree(a, desc_a,info)
  call prec%free(info)
  call psb_cdfree(desc_a,info)
  call psb_exit(ictxt)
  if (me .eq. 0) then
    print *, "--HALO TEST FIN--"
  end if
  stop

9999 call psb_error(ictxt)

  stop


contains
  subroutine sort(n, a)
    implicit none
    integer :: n, i, j
    real(psb_dpk_) :: a(n), x

    do i = 2, n
      x = a(i)
      j = i - 1
      do while (j >= 1)
        if (a(j) <= x) exit
        a(j + 1) = a(j)
        j = j - 1
      end do
      a(j + 1) = x
    end do
  end subroutine sort

  subroutine median(a, median_val, min, max)
    real(psb_dpk_), dimension(:), intent(in) :: a
    real(psb_dpk_), intent(out) :: median_val, min, max

    integer :: l, u
    real(psb_dpk_), dimension(size(a,1)) :: ac

    ac = a
    ! this is not an intrinsic: peek a sort algo from
    ! Category:Sorting, fixing it to work with real if
    ! it uses integer instead.
    l = size(a,1)
    call sort(l, ac)

    if ( mod(l, 2) == 0 ) then
      median_val = (ac(l/2+1) + ac(l/2))/2.0
    else
      median_val = ac(l/2+1)
    end if
    min = ac(1)
    max = ac(ubound(ac, dim=1))
    u = ubound(ac, dim=1)
  end subroutine median
end program psb_df_sample
