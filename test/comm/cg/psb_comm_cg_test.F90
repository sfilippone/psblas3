program psb_comm_cg_test
  use psb_base_mod
  use psb_prec_mod
  use psb_linsolve_mod
  use psb_comm_factory_mod
  use psb_comm_neighbor_impl_mod, only: psb_comm_neighbor_handle
  use, intrinsic :: ieee_arithmetic

  implicit none

  type(psb_ctxt_type) :: ctxt
  type(psb_dspmat_type) :: a
  type(psb_desc_type) :: desc_a
  type(psb_d_vect_type) :: b, x
  type(psb_dprec_type)  :: prec

  integer(psb_ipk_) :: info, iam, np
  integer(psb_ipk_) :: idim, itmax, itrace, istop, iter
  integer(psb_ipk_) :: scheme_idx, prec_idx, rep, nrep, nwarm
  integer(psb_ipk_), parameter :: n_schemes=3, n_precs=2
  integer(psb_ipk_), allocatable :: iter_count(:,:,:), solve_info(:,:,:)
  integer(psb_ipk_) :: scheme_type(n_schemes)
  real(psb_dpk_) :: eps, err, t_start, t_elapsed
  real(psb_dpk_), allocatable :: prec_init_time(:,:,:), prec_bld_time(:,:,:)
  real(psb_dpk_), allocatable :: comm_set_time(:,:,:), krylov_time(:,:,:)
  real(psb_dpk_), allocatable :: setup_time(:,:,:), solve_time(:,:,:)
  real(psb_dpk_), allocatable :: total_time(:,:,:), final_error(:,:,:)
  real(psb_dpk_), allocatable :: krylov_it_time(:,:,:), total_it_time(:,:,:)
  real(psb_dpk_) :: iter_denom
  real(psb_dpk_) :: avg_t, std_t, med_t, p10_t, p90_t, min_t, max_t
  character(len=25) :: scheme_name(n_schemes)
  character(len=12) :: prec_type(n_precs)
  character(len=20) :: prec_name(n_precs)
  character(len=5) :: afmt
  character(len=256) :: arg
  logical :: setup_done

  info = psb_success_
  afmt = 'CSR'
  idim = 40
  itmax = 1000
  nrep = 5
  nwarm = 1
  ! Keep itrace positive to avoid modulo-by-zero paths in convergence logging.
  itrace = 1
  istop = 2
  eps = 1.d-6
  scheme_type = (/ psb_comm_isend_irecv_, psb_comm_ineighbor_alltoallv_, &
       & psb_comm_persistent_ineighbor_alltoallv_ /)
  scheme_name(1) = 'isend_irecv'
  scheme_name(2) = 'ineighbor_alltoallv'
  scheme_name(3) = 'persistent_ineighbor_a2av'
  prec_type(1) = 'NONE'
  prec_type(2) = 'DIAG'
  prec_name(1) = 'none'
  prec_name(2) = 'diag'

  call get_command_argument(1,arg)
  if (len_trim(arg) > 0) then
    read(arg,*,iostat=info) idim
    if (info /= 0) then
      idim = 40
      info = psb_success_
    end if
  end if
  call get_command_argument(2,arg)
  if (len_trim(arg) > 0) then
    read(arg,*,iostat=info) nrep
    if ((info /= 0).or.(nrep <= 0)) then
      nrep = 7
      info = psb_success_
    end if
  end if
  call get_command_argument(3,arg)
  if (len_trim(arg) > 0) then
    read(arg,*,iostat=info) nwarm
    if ((info /= 0).or.(nwarm < 0)) then
      nwarm = 1
      info = psb_success_
    end if
  end if
  call get_command_argument(4,arg)
  if (len_trim(arg) > 0) then
    read(arg,*,iostat=info) itmax
    if ((info /= 0).or.(itmax <= 0)) then
      itmax = 1000
      info = psb_success_
    end if
  end if
  ! call psb_set_debug_level(psb_debug_ext_)


  ! call probe_ieee('before psb_init')
  call psb_init(ctxt)
  ! call probe_ieee('after psb_init')
  ! call clear_ieee_flags()
  ! call probe_ieee('after clear_ieee_flags')
  call psb_info(ctxt, iam, np)

  allocate(setup_time(n_precs,n_schemes,nrep), solve_time(n_precs,n_schemes,nrep), &
       & total_time(n_precs,n_schemes,nrep), final_error(n_precs,n_schemes,nrep), &
      & krylov_it_time(n_precs,n_schemes,nrep), total_it_time(n_precs,n_schemes,nrep), &
       & iter_count(n_precs,n_schemes,nrep), solve_info(n_precs,n_schemes,nrep), &
       & prec_init_time(n_precs,n_schemes,nrep), prec_bld_time(n_precs,n_schemes,nrep), &
       & comm_set_time(n_precs,n_schemes,nrep), krylov_time(n_precs,n_schemes,nrep), stat=info)
  if (info /= psb_success_) stop 1

  if (iam == psb_root_) then
    write(psb_out_unit,*) 'Welcome to PSBLAS version: ', psb_version_string_
    write(psb_out_unit,*) 'This is the comm/cg test program'
    write(psb_out_unit,'("Grid dimensions      : ",i4," x ",i4," x ",i4)') idim,idim,idim
    write(psb_out_unit,'("Number of processors : ",i0)') np
    write(psb_out_unit,'("Iterative method     : CG")')
    write(psb_out_unit,'("Preconditioners      : NONE, DIAG")')
    write(psb_out_unit,'("Max iterations (CG)  : ",i0)') itmax
    write(psb_out_unit,'("Repetitions          : ",i0)') nrep
    write(psb_out_unit,'("Warmup solves        : ",i0)') nwarm
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Usage: ./psb_comm_cg_test [idim] [nrep] [nwarm] [itmax]")')
    write(psb_out_unit,'(" ")')
  end if

  call psb_barrier(ctxt)
  t_start = psb_wtime()
  ! call probe_ieee('before psb_d_gen_pde3d')
  call psb_d_gen_pde3d(ctxt,idim,a,b,x,desc_a,afmt,info)
  ! call probe_ieee('after psb_d_gen_pde3d')
  if (info /= psb_success_) goto 9999

  do prec_idx = 1, n_precs
    do scheme_idx = 1, n_schemes
      do rep = 1, nrep
        call psb_geaxpby(dzero,b,dzero,x,desc_a,info)
        if (info /= psb_success_) goto 9999

        call psb_barrier(ctxt)

        t_start = psb_wtime()
        call prec%init(ctxt,trim(prec_type(prec_idx)),info)
        if (info /= psb_success_) goto 9999

        prec_init_time(prec_idx,scheme_idx,rep) = psb_wtime() - t_start
        call psb_amx(ctxt,prec_init_time(prec_idx,scheme_idx,rep))

        t_start = psb_wtime()
        call prec%build(a,desc_a,info)
        if (info /= psb_success_) goto 9999
        prec_bld_time(prec_idx,scheme_idx,rep) = psb_wtime() - t_start
        call psb_amx(ctxt,prec_bld_time(prec_idx,scheme_idx,rep))

        if (.not.allocated(prec%prec)) then
          info = psb_err_internal_error_
          write(psb_err_unit,*) 'Preconditioner object not allocated after build'
          goto 9999
        end if

        t_start = psb_wtime()
        call psb_comm_set(scheme_type(scheme_idx),x%v%comm_handle,info)
        comm_set_time(prec_idx,scheme_idx,rep) = psb_wtime() - t_start
        call psb_amx(ctxt,comm_set_time(prec_idx,scheme_idx,rep))

        if (info /= psb_success_) goto 9999
        if (.not.allocated(prec%prec)) then
          info = psb_err_internal_error_
          write(psb_err_unit,*) 'Preconditioner object lost after psb_comm_set'
          goto 9999
        end if

        setup_time(prec_idx,scheme_idx,rep) = prec_init_time(prec_idx,scheme_idx,rep) + &
             & prec_bld_time(prec_idx,scheme_idx,rep) + comm_set_time(prec_idx,scheme_idx,rep)

        do iter = 1, nwarm
          call psb_geaxpby(dzero,b,dzero,x,desc_a,info)
          if (info /= psb_success_) goto 9999
          call psb_krylov('CG',a,prec,b,x,eps,desc_a,info,&
               & itmax=itmax,itrace=itrace,istop=istop)
          if (info /= psb_success_) goto 9999
        end do

        call psb_geaxpby(dzero,b,dzero,x,desc_a,info)
        if (info /= psb_success_) goto 9999

        call psb_barrier(ctxt)
        t_start = psb_wtime()
        call psb_krylov('CG',a,prec,b,x,eps,desc_a,info,&
             & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istop)
        krylov_time(prec_idx,scheme_idx,rep) = psb_wtime() - t_start
        call psb_amx(ctxt,krylov_time(prec_idx,scheme_idx,rep))

        if (info /= psb_success_) goto 9999

        call psb_geaxpby(dzero,b,dzero,x,desc_a,info)
        if (info /= psb_success_) goto 9999

        call psb_barrier(ctxt)
        if (.not.allocated(prec%prec)) then
          info = psb_err_internal_error_
          write(psb_err_unit,*) 'Preconditioner object lost before psb_krylov'
          goto 9999
        end if

        call prec%free(info)
        if (info /= psb_success_) goto 9999

        solve_time(prec_idx,scheme_idx,rep) = krylov_time(prec_idx,scheme_idx,rep)
        total_time(prec_idx,scheme_idx,rep) = setup_time(prec_idx,scheme_idx,rep) + &
             & solve_time(prec_idx,scheme_idx,rep)
        iter_count(prec_idx,scheme_idx,rep) = iter
        iter_denom = real(max(iter,1_psb_ipk_),psb_dpk_)
        krylov_it_time(prec_idx,scheme_idx,rep) = solve_time(prec_idx,scheme_idx,rep)/iter_denom
        total_it_time(prec_idx,scheme_idx,rep) = total_time(prec_idx,scheme_idx,rep)/iter_denom
        final_error(prec_idx,scheme_idx,rep) = err
        solve_info(prec_idx,scheme_idx,rep) = info

        if (iam == psb_root_) then
          select type(ch => x%v%comm_handle)
          type is(psb_comm_neighbor_handle)
            write(psb_out_unit,'("DIAG_COMM scheme=",a,", prec=",a,", rep=",i0)') &
                 & trim(scheme_name(scheme_idx)), trim(prec_name(prec_idx)), rep
            write(psb_out_unit,'("DIAG_COMM counters: init=",i0,", start=",i0,", wait=",i0,", realloc=",i0)') &
              & ch%diag_init_calls, ch%diag_start_calls, ch%diag_wait_calls, &
              & ch%diag_buffer_reallocs
            write(psb_out_unit,'("DIAG_COMM state: ready=",l1,", bsz=",i0)') &
              & ch%persistent_request_ready, ch%persistent_buffer_size
          class default
            continue
          end select
        end if

        if (info /= psb_success_) goto 9999
      end do
    end do
  end do

  if (iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'(100("="))')
    write(psb_out_unit,'("CG TIMING STATISTICS - FINAL RESULTS TABLE")')
    write(psb_out_unit,'(100("="))')
    
    do prec_idx = 1, n_precs
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Preconditioner: ",a)') trim(prec_name(prec_idx))
      write(psb_out_unit,'(100("-"))')
      
      ! Print header
      write(psb_out_unit,'(a25,a18,a18,a18,a18,a18)') &
        & 'Scheme', 'Phase', 'Minimum [s]', 'Average [s]', 'Maximum [s]', 'Std Dev [s]'
      write(psb_out_unit,'(100("-"))')
      
      do scheme_idx = 1, n_schemes
        ! Preconditioner init phase
        call compute_stats(prec_init_time(prec_idx,scheme_idx,:),avg_t,std_t,med_t,p10_t,p90_t,min_t,max_t)
        write(psb_out_unit,'(a25,a18,es18.10,es18.10,es18.10,es18.10)') &
             & trim(scheme_name(scheme_idx)), 'Prec Init', min_t, avg_t, max_t, std_t
        
        ! Preconditioner build phase
        call compute_stats(prec_bld_time(prec_idx,scheme_idx,:),avg_t,std_t,med_t,p10_t,p90_t,min_t,max_t)
        write(psb_out_unit,'(a25,a18,es18.10,es18.10,es18.10,es18.10)') &
             & ' ', 'Prec Build', min_t, avg_t, max_t, std_t
        
        ! Communication setup phase
        call compute_stats(comm_set_time(prec_idx,scheme_idx,:),avg_t,std_t,med_t,p10_t,p90_t,min_t,max_t)
        write(psb_out_unit,'(a25,a18,es18.10,es18.10,es18.10,es18.10)') &
             & ' ', 'Comm Setup', min_t, avg_t, max_t, std_t
        
        ! Total setup phase
        call compute_stats(setup_time(prec_idx,scheme_idx,:),avg_t,std_t,med_t,p10_t,p90_t,min_t,max_t)
        write(psb_out_unit,'(a25,a18,es18.10,es18.10,es18.10,es18.10)') &
             & ' ', 'SetupTotal', min_t, avg_t, max_t, std_t

        ! Krylov solve phase (this is the actual solver)
        call compute_stats(krylov_time(prec_idx,scheme_idx,:),avg_t,std_t,med_t,p10_t,p90_t,min_t,max_t)
        write(psb_out_unit,'(a25,a18,es18.10,es18.10,es18.10,es18.10)') &
             & ' ', 'KrylovSolve', min_t, avg_t, max_t, std_t

           ! Krylov solve normalized per actual CG iteration
           call compute_stats(krylov_it_time(prec_idx,scheme_idx,:),avg_t,std_t,med_t,p10_t,p90_t,min_t,max_t)
           write(psb_out_unit,'(a25,a18,es18.10,es18.10,es18.10,es18.10)') &
             & ' ', 'KrylovPerIter', min_t, avg_t, max_t, std_t
        
        ! Total phase statistics
        call compute_stats(total_time(prec_idx,scheme_idx,:),avg_t,std_t,med_t,p10_t,p90_t,min_t,max_t)
        write(psb_out_unit,'(a25,a18,es18.10,es18.10,es18.10,es18.10)') &
             & ' ', 'TotalTime', min_t, avg_t, max_t, std_t

           ! Total (setup+solve) normalized per actual CG iteration
           call compute_stats(total_it_time(prec_idx,scheme_idx,:),avg_t,std_t,med_t,p10_t,p90_t,min_t,max_t)
           write(psb_out_unit,'(a25,a18,es18.10,es18.10,es18.10,es18.10)') &
             & ' ', 'TotalPerIter', min_t, avg_t, max_t, std_t
        
        ! Final error and convergence info
        call compute_stats(final_error(prec_idx,scheme_idx,:),avg_t,std_t,med_t,p10_t,p90_t,min_t,max_t)
        write(psb_out_unit,'(a25,a18,es18.10,es18.10,es18.10,es18.10)') &
             & ' ', 'Final Error', min_t, avg_t, max_t, std_t
        write(psb_out_unit,'(a25,a18,"Iterations: ",i8,"  Info code: ",i6)') &
             & ' ', ' ', iter_count(prec_idx,scheme_idx,nrep), solve_info(prec_idx,scheme_idx,nrep)
        
        write(psb_out_unit,'(100("-"))')
      end do
    end do
    
    write(psb_out_unit,'(100("="))')
    write(psb_out_unit,'(" ")')
  end if

  call psb_gefree(b,desc_a,info)
  call psb_gefree(x,desc_a,info)
  call psb_spfree(a,desc_a,info)
  call psb_precfree(prec,info)
  call psb_cdfree(desc_a,info)
  deallocate(setup_time,solve_time,total_time,final_error,iter_count,solve_info, &
      & prec_init_time,prec_bld_time,comm_set_time,krylov_time, &
      & krylov_it_time,total_it_time)

  call psb_exit(ctxt)
  stop

9999 call psb_error(ctxt)
  stop 1

contains

  subroutine sort_real_inplace(v)
    real(psb_dpk_), intent(inout) :: v(:)
    integer(psb_ipk_) :: i, j
    real(psb_dpk_) :: key

    do i = 2, size(v)
      key = v(i)
      j = i - 1
      do while (j >= 1)
        if (v(j) <= key) exit
        v(j+1) = v(j)
        j = j - 1
      end do
      v(j+1) = key
    end do
  end subroutine sort_real_inplace

  subroutine compute_stats(vals,mean_v,std_v,median_v,p10_v,p90_v,min_v,max_v)
    real(psb_dpk_), intent(in)  :: vals(:)
    real(psb_dpk_), intent(out) :: mean_v,std_v,median_v,p10_v,p90_v,min_v,max_v
    real(psb_dpk_), allocatable :: tmp(:)
    integer(psb_ipk_) :: n, idx10, idx90

    n = size(vals)
    if (n <= 0) then
      mean_v = dzero; std_v = dzero; median_v = dzero
      p10_v  = dzero; p90_v = dzero; min_v = dzero; max_v = dzero
      return
    end if

    mean_v = sum(vals)/real(n,psb_dpk_)
    if (n > 1) then
      std_v = sqrt(sum((vals-mean_v)**2)/real(n-1,psb_dpk_))
    else
      std_v = dzero
    end if

    allocate(tmp(n))
    tmp = vals
    call sort_real_inplace(tmp)

    if (mod(n,2) == 0) then
      median_v = (tmp(n/2)+tmp(n/2+1))/2.0_psb_dpk_
    else
      median_v = tmp((n+1)/2)
    end if

    idx10 = int(ceiling(0.10_psb_dpk_*real(n,psb_dpk_)),kind=psb_ipk_)
    idx90 = int(ceiling(0.90_psb_dpk_*real(n,psb_dpk_)),kind=psb_ipk_)
    idx10 = max(1_psb_ipk_,min(n,idx10))
    idx90 = max(1_psb_ipk_,min(n,idx90))
    p10_v = tmp(idx10)
    p90_v = tmp(idx90)
    min_v = tmp(1)
    max_v = tmp(n)

    deallocate(tmp)
  end subroutine compute_stats

  function b1(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = dzero
  end function b1

  function b2(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = dzero
  end function b2

  function b3(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = dzero
  end function b3

  function cfun(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = dzero
  end function cfun

  function a1(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = done/80
  end function a1

  function a2(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = done/80
  end function a2

  function a3(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = done/80
  end function a3

  function gfun(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = dzero
    if (x == done) then
      val = done
    else if (x == dzero) then
      val = exp(y**2-z**2)
    end if
  end function gfun

  subroutine probe_ieee(where)
    character(len=*), intent(in) :: where
    logical :: invalid_flag, divzero_flag, overflow_flag, underflow_flag

    call ieee_get_flag(ieee_invalid, invalid_flag)
    call ieee_get_flag(ieee_divide_by_zero, divzero_flag)
    call ieee_get_flag(ieee_overflow, overflow_flag)
    call ieee_get_flag(ieee_underflow, underflow_flag)

    if (invalid_flag .or. divzero_flag .or. overflow_flag .or. underflow_flag) then
      write(psb_out_unit,'("IEEE probe [",a,"] invalid=",l1,", div0=",l1,", overflow=",l1,", underflow=",l1)') &
           trim(where), invalid_flag, divzero_flag, overflow_flag, underflow_flag
    end if
  end subroutine probe_ieee

  subroutine clear_ieee_flags()
    call ieee_set_flag(ieee_invalid, .false.)
    call ieee_set_flag(ieee_divide_by_zero, .false.)
    call ieee_set_flag(ieee_overflow, .false.)
    call ieee_set_flag(ieee_underflow, .false.)
  end subroutine clear_ieee_flags

  subroutine psb_d_gen_pde3d(ctxt,idim,a,bv,xv,desc_a,afmt,info)
    implicit none
    integer(psb_ipk_), intent(in)     :: idim
    type(psb_dspmat_type), intent(out) :: a
    type(psb_d_vect_type), intent(out) :: xv,bv
    type(psb_desc_type), intent(out)   :: desc_a
    type(psb_ctxt_type), intent(in)    :: ctxt
    integer(psb_ipk_), intent(out)     :: info
    character(len=*), intent(in)       :: afmt

    integer(psb_ipk_), parameter :: nb=20
    real(psb_dpk_) :: zt(nb),x,y,z
    integer(psb_lpk_) :: m,n,glob_row
    integer(psb_ipk_) :: nnz,nlr,i,ii,ib,k
    integer(psb_ipk_) :: ix,iy,iz
    integer(psb_ipk_) :: np, iam, nr, nt
    integer(psb_ipk_) :: icoeff
    integer(psb_lpk_), allocatable :: irow(:),icol(:),myidx(:)
    real(psb_dpk_), allocatable :: val(:)
    real(psb_dpk_) :: deltah, sqdeltah, deltah2
    real(psb_dpk_) :: t0, t1, t2, t3, tasb, talc, ttot, tgen, tcdasb
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name, ch_err, tmpfmt

    info = psb_success_
    name = 'create_matrix'
    call psb_erractionsave(err_act)

    call psb_info(ctxt, iam, np)

    if (idim <= 0) then
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='idim must be > 0')
      goto 9999
    end if
    if (np <= 0) then
      info = psb_err_context_error_
      call psb_errpush(info,name,a_err='invalid context: np <= 0')
      goto 9999
    end if
    if (iam < 0) then
      info = psb_err_context_error_
      call psb_errpush(info,name,a_err='invalid context: iam < 0')
      goto 9999
    end if

    deltah   = done/(idim+2)
    sqdeltah = deltah*deltah
    deltah2  = 2.d0*deltah

    if (abs(deltah) <= tiny(deltah)) then
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='invalid mesh spacing: deltah ~ 0')
      goto 9999
    end if
    if (abs(sqdeltah) <= tiny(sqdeltah)) then
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='invalid mesh spacing: sqdeltah ~ 0')
      goto 9999
    end if
    if (abs(deltah2) <= tiny(deltah2)) then
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='invalid mesh spacing: deltah2 ~ 0')
      goto 9999
    end if

    m   = idim*idim*idim
    n   = m
    if (n <= 0) then
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='invalid global size: n <= 0')
      goto 9999
    end if
    nnz = ((n*9)/(np))
    if (nnz <= 0) then
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='invalid local nnz estimate: nnz <= 0')
      goto 9999
    end if
    if(iam == psb_root_) write(psb_out_unit,'("Generating Matrix (size=",i0,")...")')n

    nt = (m+np-1)/np
    nr = max(0,min(nt,m-(iam*nt)))
    nt = nr
    call psb_sum(ctxt,nt)
    if (nt /= m) then
      write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,m
      info = -1
      call psb_barrier(ctxt)
      call psb_abort(ctxt)
      return
    end if

    call psb_barrier(ctxt)
    t0 = psb_wtime()
    ! call probe_ieee('enter psb_cdall')
    call psb_cdall(ctxt,desc_a,info,nl=nr)
    ! call probe_ieee('after psb_cdall')
    if (info == psb_success_) call psb_spall(a,desc_a,info,nnz=nnz)
    ! call probe_ieee('after psb_spall')
    if (info == psb_success_) call psb_geall(xv,desc_a,info)
    if (info == psb_success_) call psb_geall(bv,desc_a,info)
    call psb_barrier(ctxt)
    talc = psb_wtime()-t0

    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='allocation rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    allocate(val(20*nb),irow(20*nb),icol(20*nb),stat=info)
    if (info /= psb_success_) then
      info=psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    endif

    myidx = desc_a%get_global_indices()
    nlr = size(myidx)

    call psb_barrier(ctxt)
    t1 = psb_wtime()
    do ii=1, nlr,nb
      ib = min(nb,nlr-ii+1)
      icoeff = 1
      do k=1,ib
        i=ii+k-1
        glob_row=myidx(i)

        if (mod(glob_row,(idim*idim)) == 0) then
          ix = glob_row/(idim*idim)
        else
          ix = glob_row/(idim*idim)+1
        endif
        if (mod((glob_row-(ix-1)*idim*idim),idim) == 0) then
          iy = (glob_row-(ix-1)*idim*idim)/idim
        else
          iy = (glob_row-(ix-1)*idim*idim)/idim+1
        endif
        iz = glob_row-(ix-1)*idim*idim-(iy-1)*idim

        x = (ix-1)*deltah
        y = (iy-1)*deltah
        z = (iz-1)*deltah
        zt(k) = dzero

        val(icoeff) = -a1(x,y,z)/sqdeltah-b1(x,y,z)/deltah2
        if (ix == 1) then
          zt(k) = gfun(dzero,y,z)*(-val(icoeff)) + zt(k)
        else
          icol(icoeff) = (ix-2)*idim*idim+(iy-1)*idim+iz
          irow(icoeff) = glob_row
          icoeff = icoeff+1
        endif

        val(icoeff)  = -a2(x,y,z)/sqdeltah-b2(x,y,z)/deltah2
        if (iy == 1) then
          zt(k) = gfun(x,dzero,z)*(-val(icoeff)) + zt(k)
        else
          icol(icoeff) = (ix-1)*idim*idim+(iy-2)*idim+iz
          irow(icoeff) = glob_row
          icoeff = icoeff+1
        endif

        val(icoeff) = -a3(x,y,z)/sqdeltah-b3(x,y,z)/deltah2
        if (iz == 1) then
          zt(k) = gfun(x,y,dzero)*(-val(icoeff)) + zt(k)
        else
          icol(icoeff) = (ix-1)*idim*idim+(iy-1)*idim+(iz-1)
          irow(icoeff) = glob_row
          icoeff = icoeff+1
        endif

        val(icoeff)=2.d0*(a1(x,y,z)+a2(x,y,z)+a3(x,y,z))/sqdeltah + cfun(x,y,z)
        icol(icoeff) = (ix-1)*idim*idim+(iy-1)*idim+iz
        irow(icoeff) = glob_row
        icoeff = icoeff+1

        val(icoeff) = -a3(x,y,z)/sqdeltah+b3(x,y,z)/deltah2
        if (iz == idim) then
          zt(k) = gfun(x,y,done)*(-val(icoeff)) + zt(k)
        else
          icol(icoeff) = (ix-1)*idim*idim+(iy-1)*idim+(iz+1)
          irow(icoeff) = glob_row
          icoeff = icoeff+1
        endif

        val(icoeff) = -a2(x,y,z)/sqdeltah+b2(x,y,z)/deltah2
        if (iy == idim) then
          zt(k) = gfun(x,done,z)*(-val(icoeff)) + zt(k)
        else
          icol(icoeff) = (ix-1)*idim*idim+iy*idim+iz
          irow(icoeff) = glob_row
          icoeff = icoeff+1
        endif

        val(icoeff) = -a1(x,y,z)/sqdeltah+b1(x,y,z)/deltah2
        if (ix == idim) then
          zt(k) = gfun(done,y,z)*(-val(icoeff)) + zt(k)
        else
          icol(icoeff) = ix*idim*idim+(iy-1)*idim+iz
          irow(icoeff) = glob_row
          icoeff = icoeff+1
        endif
      end do

      call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
      ! call probe_ieee('after psb_spins')
      if(info /= psb_success_) then
        write(psb_err_unit,'("INSERT FAIL rank=",i0,", call=psb_spins, ii=",i0,", ib=",i0,", icoeff=",i0)') &
             iam, ii, ib, icoeff
        write(psb_err_unit,'("  glob_row=",i0,", ix=",i0,", iy=",i0,", iz=",i0)') glob_row, ix, iy, iz
        exit
      end if
      call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),bv,desc_a,info)
      ! call probe_ieee('after psb_geins bv')
      if(info /= psb_success_) then
        write(psb_err_unit,'("INSERT FAIL rank=",i0,", call=psb_geins bv, ii=",i0,", ib=",i0,", icoeff=",i0)') &
             iam, ii, ib, icoeff
        write(psb_err_unit,'("  glob_row=",i0,", ix=",i0,", iy=",i0,", iz=",i0)') glob_row, ix, iy, iz
        exit
      end if
      zt(:)=dzero
      call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),xv,desc_a,info)
      ! call probe_ieee('after psb_geins xv')
      if(info /= psb_success_) then
        write(psb_err_unit,'("INSERT FAIL rank=",i0,", call=psb_geins xv, ii=",i0,", ib=",i0,", icoeff=",i0)') &
             iam, ii, ib, icoeff
        write(psb_err_unit,'("  glob_row=",i0,", ix=",i0,", iy=",i0,", iz=",i0)') glob_row, ix, iy, iz
        exit
      end if
    end do

    tgen = psb_wtime()-t1
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='insert rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    deallocate(val,irow,icol)

    call psb_barrier(ctxt)
    t1 = psb_wtime()
    ! call probe_ieee('before psb_cdasb')
    call psb_cdasb(desc_a,info)
    ! call probe_ieee('after psb_cdasb')
    tcdasb = psb_wtime()-t1

    call psb_barrier(ctxt)
    t1 = psb_wtime()
    ! call probe_ieee('before psb_spasb')
    if (info == psb_success_) call psb_spasb(a,desc_a,info,afmt=afmt)
    ! call probe_ieee('after psb_spasb')
    call psb_barrier(ctxt)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='asb rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    if (info == psb_success_) call psb_geasb(xv,desc_a,info)
    if (info == psb_success_) call psb_geasb(bv,desc_a,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='asb rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    tasb = psb_wtime()-t1
    call psb_barrier(ctxt)
    ttot = psb_wtime() - t0

    call psb_amx(ctxt,talc)
    call psb_amx(ctxt,tgen)
    call psb_amx(ctxt,tasb)
    call psb_amx(ctxt,ttot)
    if(iam == psb_root_) then
      tmpfmt = a%get_fmt()
      write(psb_out_unit,'("The matrix has been generated and assembled in ",a3," format.")') tmpfmt
      write(psb_out_unit,'("-allocation  time : ",es12.5)') talc
      write(psb_out_unit,'("-coeff. gen. time : ",es12.5)') tgen
      write(psb_out_unit,'("-desc asbly  time : ",es12.5)') tcdasb
      write(psb_out_unit,'("- mat asbly  time : ",es12.5)') tasb
      write(psb_out_unit,'("-total       time : ",es12.5)') ttot
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)
    return
  end subroutine psb_d_gen_pde3d

end program psb_comm_cg_test
