!
! Test program for z-type halo exchange: baseline vs neighbor topology.
!
! This test exercises the lower-level psi_swapdata interface directly
! to compare the two communication paths implemented in psi_dswapdata.F90:
!
!   1. Baseline  (Isend/Irecv)         : flag = IOR(psb_swap_send_, psb_swap_recv_)
!   2. Neighbor topology (Ineighbor_alltoallv) : flag = psb_swap_start_ then psb_swap_wait_
!
! It builds a 3D block-partitioned descriptor with a 7-point stencil,
! fills owned entries with their global index, performs halo exchange
! via both paths, then checks:
!   (a) The two paths produce identical results (cross-check)
!   (b) Every halo entry equals the global index of its source (absolute check)
!
! Run with:  mpirun -np <P> ./psb_z_comm_test
!
program psb_z_comm_test
  use psb_base_mod
  use psb_util_mod
  use psb_error_mod, only: psb_set_debug_level, psb_debug_ext_
  use psi_mod
  use psb_comm_factory_mod, only: psb_comm_set, psb_comm_free
  use psb_comm_schemes_mod, only: psb_comm_ineighbor_alltoallv_, psb_comm_persistent_ineighbor_alltoallv_, &
    & psb_comm_isend_irecv_, psb_comm_rma_pull_, psb_comm_rma_push_
  use psb_comm_schemes_mod, only: psb_comm_status_start_, psb_comm_status_wait_, psb_comm_status_unknown_
  implicit none

  ! ---- parameters ----
  integer(psb_ipk_)                     :: idim
  integer(psb_ipk_)                     :: argc
  integer(psb_ipk_)                     :: iters
  character(len=256)                    :: arg
  character(len=16)                     :: mode
  character(len=256)                    :: matrix_file
  character(len=2)                      :: matrix_fmt
  logical                               :: debug_swapdata
  logical                               :: use_external_matrix

  ! ---- descriptor / context ----
  type(psb_ctxt_type)                   :: ctxt
  type(psb_desc_type)                   :: desc_a
  integer(psb_ipk_)                     :: my_rank, np, info, i, nr, number_of_local_rows
  integer(psb_lpk_)                     :: m, nt
  integer(psb_lpk_), allocatable        :: myidx(:)
  type(psb_zspmat_type)               :: a_mat
  type(psb_lzspmat_type)              :: aux_a

  ! ---- vectors ----
  type(psb_z_vect_type)               :: v_baseline, v_neighbor, v_neighbor_persistent, v_rma_get, v_rma_put

  ! ---- temporary / comparison arrays ----
  complex(psb_dpk_), allocatable          :: vals(:)
  complex(psb_dpk_), allocatable          :: result_baseline(:), result_neighbor(:), result_persistent(:), &
    & result_rma_get(:), result_rma_put(:)
  complex(psb_dpk_), allocatable          :: expected(:)

  ! ---- halo index bookkeeping ----
  integer(psb_ipk_)                     :: nrow, ncol, num_neighbors, send_indexes, receive_indexes
  class(psb_i_base_vect_type), pointer  :: halo_indexes

  ! ---- error / reporting ----
  integer(psb_ipk_)                     :: n_pass, n_total, imode
  logical                               :: run_baseline, run_neighbor, run_persistent, run_rma_get, run_rma_put
  logical                               :: mat_allocated
  logical                               :: comm_ok
  real(psb_dpk_)                      :: err, tol
  real(psb_dpk_)                      :: first_swap_baseline, first_swap_neighbor, first_swap_persistent, &
    & first_swap_rma_get, first_swap_rma_put
  real(psb_dpk_)                      :: comm_setup_time_baseline, comm_setup_time_neighbor, comm_setup_time_persistent, &
    & comm_setup_time_rma_get, comm_setup_time_rma_put
  real(psb_dpk_)                      :: t0, t1, dt, tsum_baseline, tsum_neighbor, tsum_neighbor_persistent, &
    & tsum_rma_get, tsum_rma_put
  integer(psb_lpk_), allocatable        :: glob_col(:)
  character(len=40)                     :: name
  real(psb_dpk_)                      :: huge_d

  name    = 'psb_z_comm_test'
  tol     = 1.0e-12_psb_dpk_
  huge_d  = huge(done)
  n_pass  = 0
  n_total = 0
  iters = 5
  mode = 'both'
  debug_swapdata = .false.
  matrix_file = ''
  matrix_fmt = 'MM'
  use_external_matrix = .false.
  mat_allocated = .false.

  ! ---- parse command-line argument for idim ----
  idim = 10
  argc = command_argument_count()
  do i = 1, argc
    call get_command_argument(i, arg)
    if (trim(arg) == '--dim') then
      if (i < argc) then
        call get_command_argument(i+1, arg)
        read(arg, *) idim
      end if
    else if (trim(arg) == '--iters') then
      if (i < argc) then
        call get_command_argument(i+1, arg)
        read(arg, *) iters
      end if
    else if (index(psb_toupper(trim(arg)),'--MATRIX=') == 1) then
      matrix_file = adjustl(arg(10:len_trim(arg)))
    else if (trim(psb_toupper(arg)) == '--MATRIX') then
      if (i < argc) then
        call get_command_argument(i+1, matrix_file)
      end if
    else if (index(psb_toupper(trim(arg)),'--FMT=') == 1) then
      arg = psb_toupper(adjustl(arg(7:len_trim(arg))))
      if ((trim(arg) == 'MM') .or. (trim(arg) == 'HB')) matrix_fmt = trim(arg)
    else if (trim(psb_toupper(arg)) == '--FMT') then
      if (i < argc) then
        call get_command_argument(i+1, arg)
        arg = psb_toupper(trim(arg))
        if ((trim(arg) == 'MM') .or. (trim(arg) == 'HB')) matrix_fmt = trim(arg)
      end if
    end if
  end do
  use_external_matrix = (len_trim(matrix_file) > 0)

  ! parse optional mode flag
  do i = 1, argc
    call get_command_argument(i, arg)
    if (trim(arg) == '--mode') then
      if (i < argc) then
        call get_command_argument(i+1, arg)
        read(arg, *) mode
      end if
    else if (trim(arg) == '--debug') then
      debug_swapdata = .true.
    end if
  end do

  if (debug_swapdata) then
    call psb_set_debug_level(psb_debug_ext_)
  end if

  run_baseline    = .false.
  run_neighbor    = .false.
  run_persistent  = .false.
  run_rma_get     = .false.
  run_rma_put     = .false.


  select case (trim(adjustl(mode)))
  case ('both','all')
    run_baseline    = .true.
    run_neighbor    = .true.
    run_persistent  = .true.
    run_rma_get     = .true.
    run_rma_put     = .true.
  case ('baseline')
    run_baseline   = .true.
  case ('neighbor')
    run_neighbor   = .true.
  case ('persistent','persistent_neighbor','persistent-neighbor')
    run_persistent = .true.
  case ('rma_get')
    run_rma_get = .true.
  case ('rma_put')
    run_rma_put = .true.
  case default
    run_baseline   = .true.
    run_neighbor   = .true.
    run_persistent = .true.
    run_rma_get    = .true.
    run_rma_put    = .true.
  end select

  if ((.not.use_external_matrix) .and. (idim <= 0)) then
      write(*,*) 'Invalid dimension specified. Usage: --dim <positive integer>'
      call psb_abort(ctxt)
  end if
  ! ==================================================================
  !  1. Initialise MPI / PSBLAS context
  ! ==================================================================
  call psb_init(ctxt)
  call psb_info(ctxt, my_rank, np)

  if (my_rank == 0) then
    write(psb_out_unit,'("================================================")')
    write(psb_out_unit,'("  Test: z-type halo  baseline vs neighbor topo")')
    write(psb_out_unit,'("  Processes : ",i0)') np
    if (use_external_matrix) then
      write(psb_out_unit,'("  Matrix    : ",a)') trim(matrix_file)
      write(psb_out_unit,'("  Format    : ",a)') trim(matrix_fmt)
    else
      write(psb_out_unit,'("  Grid      : ",i0," x ",i0," x ",i0)') idim,idim,idim
    end if
        write(psb_out_unit,'("  Usage     : ./psb_z_comm_test [--dim N] [--iters N] [--mode ...] ",&
          &"[--matrix <path>] [--fmt MM|HB]")')
    write(psb_out_unit,'("================================================")')
  end if

  ! ==================================================================
  !  2. Build descriptor with 7-point stencil connectivity
  ! ==================================================================
  if (use_external_matrix) then
    select case(psb_toupper(trim(matrix_fmt)))
    case('MM')
      call mm_mat_read(aux_a,info,filename=trim(matrix_file))
    case('HB')
      call hb_read(aux_a,info,filename=trim(matrix_file))
    case default
      info = psb_err_internal_error_
    end select
    if (info /= psb_success_) then
      write(psb_err_unit,*) my_rank, 'matrix read error:', info
      call psb_abort(ctxt)
    end if
    if (aux_a%get_nrows() /= aux_a%get_ncols()) then
      write(psb_err_unit,*) my_rank, 'matrix must be square for this test'
      call psb_abort(ctxt)
    end if
    m = aux_a%get_nrows()
    call psb_matdist(aux_a, a_mat, ctxt, desc_a, info, fmt='CSR', parts=part_block)
    if (info /= psb_success_) then
      write(psb_err_unit,*) my_rank, 'matdist error:', info
      call psb_abort(ctxt)
    end if
    mat_allocated = .true.
    myidx = desc_a%get_global_indices()
    number_of_local_rows = size(myidx)
  else
    m  = (1_psb_lpk_ * idim) * idim * idim
    nt = (m + np - 1) / np
    nr = max(0, min(int(nt,psb_ipk_), int(m - (my_rank * nt),psb_ipk_)))

    call psb_cdall(ctxt, desc_a, info, nl=nr)
    if (info /= psb_success_) then
      write(psb_err_unit,*) my_rank, 'cdall error:', info
      call psb_abort(ctxt)
    end if

    myidx = desc_a%get_global_indices()
    number_of_local_rows   = size(myidx)

    do i = 1, number_of_local_rows
      call psb_cdins(1_psb_ipk_, (/myidx(i)/), (/myidx(i)/), desc_a, info)
      if (myidx(i) > 1) &
           & call psb_cdins(1_psb_ipk_, (/myidx(i)/), (/myidx(i)-1/), desc_a, info)
      if (myidx(i) < m) &
           & call psb_cdins(1_psb_ipk_, (/myidx(i)/), (/myidx(i)+1/), desc_a, info)
      if (myidx(i) > idim) &
           & call psb_cdins(1_psb_ipk_, (/myidx(i)/), (/myidx(i)-idim/), desc_a, info)
      if (myidx(i) + idim <= m) &
           & call psb_cdins(1_psb_ipk_, (/myidx(i)/), (/myidx(i)+idim/), desc_a, info)
      if (myidx(i) > int(idim,psb_lpk_)*idim) &
           & call psb_cdins(1_psb_ipk_, (/myidx(i)/), &
           &   (/myidx(i) - int(idim,psb_lpk_)*idim/), desc_a, info)
      if (myidx(i) + int(idim,psb_lpk_)*idim <= m) &
           & call psb_cdins(1_psb_ipk_, (/myidx(i)/), &
           &   (/myidx(i) + int(idim,psb_lpk_)*idim/), desc_a, info)
    end do

    call psb_cdasb(desc_a, info)
    if (info /= psb_success_) then
      write(psb_err_unit,*) my_rank, 'cdasb error:', info
      call psb_abort(ctxt)
    end if
  end if

  nrow = desc_a%get_local_rows()   ! owned
  ncol = desc_a%get_local_cols()   ! owned + halo

  ! ==================================================================
  !  3. Allocate two D vectors (scratch) and fill owned entries
  ! ==================================================================
  call psb_geall(v_baseline, desc_a, info)
  if (info /= psb_success_) then
    write(psb_err_unit,*) my_rank, 'geall baseline error:', info
    call psb_abort(ctxt)
  end if
  call psb_geall(v_neighbor, desc_a, info)
  if (info /= psb_success_) then
    write(psb_err_unit,*) my_rank, 'geall neighbor error:', info
    call psb_abort(ctxt)
  end if
  call psb_geall(v_neighbor_persistent, desc_a, info)
  if (info /= psb_success_) then
    write(psb_err_unit,*) my_rank, 'geall persistent-neighbor error:', info
    call psb_abort(ctxt)
  end if
  call psb_geall(v_rma_get, desc_a, info)
  if (info /= psb_success_) then
    write(psb_err_unit,*) my_rank, 'geall rma-get error:', info
    call psb_abort(ctxt)
  end if
  call psb_geall(v_rma_put, desc_a, info)
  if (info /= psb_success_) then
    write(psb_err_unit,*) my_rank, 'geall rma-put error:', info
    call psb_abort(ctxt)
  end if
  
  
  call psb_geasb(v_baseline, desc_a, info, scratch=.true.)
  if (info /= psb_success_) then
    write(psb_err_unit,*) my_rank, 'geasb baseline error:', info
    call psb_abort(ctxt)
  end if
  call psb_geasb(v_neighbor, desc_a, info, scratch=.true.)
  if (info /= psb_success_) then
    write(psb_err_unit,*) my_rank, 'geasb neighbor error:', info
    call psb_abort(ctxt)
  end if
  call psb_geasb(v_neighbor_persistent, desc_a, info, scratch=.true.)
  if (info /= psb_success_) then
    write(psb_err_unit,*) my_rank, 'geasb persistent-neighbor error:', info
    call psb_abort(ctxt)
  end if
  call psb_geasb(v_rma_get, desc_a, info, scratch=.true.)
  if (info /= psb_success_) then
    write(psb_err_unit,*) my_rank, 'geasb rma-get error:', info
    call psb_abort(ctxt)
  end if
  call psb_geasb(v_rma_put, desc_a, info, scratch=.true.)
  if (info /= psb_success_) then
    write(psb_err_unit,*) my_rank, 'geasb rma-put error:', info
    call psb_abort(ctxt)
  end if



  ! Fill owned entries with the global index value
  allocate(vals(ncol))
  vals = zzero
  do i = 1, number_of_local_rows
    vals(i) = real(myidx(i),psb_dpk_)
  end do
  call v_baseline%set_vect(vals)
  call v_neighbor%set_vect(vals)
  call v_neighbor_persistent%set_vect(vals)
  call v_rma_get%set_vect(vals)
  call v_rma_put%set_vect(vals)
  deallocate(vals)

  ! ==================================================================
  !  4. Build the expected result for halo positions
  !     glob_col(j) = global index of local column j
  !     After halo exchange every position j should hold glob_col(j).
  ! ==================================================================
  allocate(glob_col(ncol), expected(ncol))
  glob_col = desc_a%get_global_indices(owned=.false.)
  do i = 1, ncol
    expected(i) = real(glob_col(i),psb_dpk_)
  end do
  allocate(result_baseline(ncol), result_neighbor(ncol), result_persistent(ncol), &
    & result_rma_get(ncol), result_rma_put(ncol))
  result_baseline   = huge_d
  result_neighbor   = huge_d
  result_persistent = huge_d
  result_rma_get    = huge_d
  result_rma_put    = huge_d

  first_swap_baseline = dzero
  first_swap_neighbor = dzero
  first_swap_persistent = dzero
  first_swap_rma_get = dzero
  first_swap_rma_put = dzero

  comm_setup_time_baseline = dzero
  comm_setup_time_neighbor = dzero
  comm_setup_time_persistent = dzero
  comm_setup_time_rma_get = dzero
  comm_setup_time_rma_put = dzero

  ! ==================================================================
  !  6. Baseline halo exchange  (Isend/Irecv in one call)
  ! ==================================================================
  if (run_baseline) then
    comm_setup_time_baseline = psb_wtime()
    call psb_comm_set(psb_comm_isend_irecv_, v_baseline%v%comm_handle, info)
    if (info /= 0) then
      write(psb_err_unit,*) my_rank, 'psb_comm_set baseline error:', info
      call psb_abort(ctxt)
    end if
    comm_setup_time_baseline = psb_wtime() - comm_setup_time_baseline

    first_swap_baseline = psb_wtime()
    call psi_swapdata( &
      swap_status=psb_comm_status_start_, &
      beta=zzero, &
      y=v_baseline%v, &
      desc_a=desc_a, &
      info=info, &
      data=psb_comm_halo_)
    if (info /= psb_success_) then
      write(psb_err_unit,*) my_rank, 'baseline swap error:', info
      call psb_abort(ctxt)
    end if

    call psi_swapdata( &
        swap_status=psb_comm_status_wait_, &
        beta=zzero, &
        y=v_baseline%v, &
        desc_a=desc_a, &
        info=info, &
        data=psb_comm_halo_)
    if (info /= psb_success_) then
      write(psb_err_unit,*) my_rank, 'baseline swap error:', info
      call psb_abort(ctxt)
    end if
    first_swap_baseline = psb_wtime() - first_swap_baseline
  end if


  ! ==================================================================
  !  7. Neighbor topology halo exchange  (start + wait)
  ! ==================================================================
  if (run_neighbor) then
    comm_setup_time_neighbor = psb_wtime()
    call psb_comm_set(psb_comm_ineighbor_alltoallv_, v_neighbor%v%comm_handle, info)
    if (info /= 0) then
      write(psb_err_unit,*) my_rank, 'psb_comm_set neighbor error:', info
      call psb_abort(ctxt)
    end if
    comm_setup_time_neighbor = psb_wtime() - comm_setup_time_neighbor


    first_swap_neighbor = psb_wtime()
    call psi_swapdata(psb_comm_status_start_, zzero, v_neighbor%v, desc_a, info, data=psb_comm_halo_)
    if (info /= psb_success_) then
      write(psb_err_unit,*) my_rank, 'neighbor start error:', info
      call psb_abort(ctxt)
    end if

    call psi_swapdata(psb_comm_status_wait_, zzero, v_neighbor%v, desc_a, info, data=psb_comm_halo_)
    if (info /= psb_success_) then
      write(psb_err_unit,*) my_rank, 'neighbor wait error:', info
      call psb_abort(ctxt)
    end if
    first_swap_neighbor = psb_wtime() - first_swap_neighbor
  end if

  ! ==================================================================
  !  7b. Persistent-neighbor halo exchange  (start + wait)
  ! ==================================================================
  if (run_persistent) then
    comm_setup_time_persistent = psb_wtime()
    call psb_comm_set(psb_comm_persistent_ineighbor_alltoallv_, v_neighbor_persistent%v%comm_handle, info)
    if (info /= 0) then
      write(psb_err_unit,*) my_rank, 'psb_comm_set persistent-neighbor error:', info
      call psb_abort(ctxt)
    end if
    comm_setup_time_persistent = psb_wtime() - comm_setup_time_persistent


    first_swap_persistent = psb_wtime()
    call psi_swapdata(psb_comm_status_start_, zzero, v_neighbor_persistent%v, desc_a, info, data=psb_comm_halo_)
    if (info /= psb_success_) then
      write(psb_err_unit,*) my_rank, 'persistent-neighbor start error:', info
      call psb_abort(ctxt)
    end if
    call psi_swapdata(psb_comm_status_wait_, zzero, v_neighbor_persistent%v, desc_a, info, data=psb_comm_halo_)
    if (info /= psb_success_) then
      write(psb_err_unit,*) my_rank, 'persistent-neighbor wait error:', info
      call psb_abort(ctxt)
    end if
    first_swap_persistent = psb_wtime() - first_swap_persistent
  end if



  if(run_rma_get) then 
    comm_setup_time_rma_get = psb_wtime()
    call psb_comm_set(psb_comm_rma_pull_, v_rma_get%v%comm_handle, info)
    if (info /= 0) then
      write(psb_err_unit,*) my_rank, 'psb_comm_set RMA get error:', info
      call psb_abort(ctxt)
    end if
    comm_setup_time_rma_get = psb_wtime() - comm_setup_time_rma_get

    first_swap_rma_get = psb_wtime()
    call psi_swapdata(psb_comm_status_start_, zzero, v_rma_get%v, desc_a, info, data=psb_comm_halo_)
    if (info /= psb_success_) then
      write(psb_err_unit,*) my_rank, 'RMA get start error:', info
      call psb_abort(ctxt)
    end if
    call psi_swapdata(psb_comm_status_wait_, zzero, v_rma_get%v, desc_a, info, data=psb_comm_halo_)
    if (info /= psb_success_) then
      write(psb_err_unit,*) my_rank, 'RMA get wait error:', info
      call psb_abort(ctxt)
    end if
    first_swap_rma_get = psb_wtime() - first_swap_rma_get
  end if




    if(run_rma_put) then 
    comm_setup_time_rma_put = psb_wtime()
    call psb_comm_set(psb_comm_rma_push_, v_rma_put%v%comm_handle, info)
    if (info /= 0) then
      write(psb_err_unit,*) my_rank, 'psb_comm_set RMA put error:', info
      call psb_abort(ctxt)
    end if
    comm_setup_time_rma_put = psb_wtime() - comm_setup_time_rma_put

    first_swap_rma_put = psb_wtime()
    call psi_swapdata(psb_comm_status_start_, zzero, v_rma_put%v, desc_a, info, data=psb_comm_halo_)
    if (info /= psb_success_) then
      write(psb_err_unit,*) my_rank, 'RMA put start error:', info
      call psb_abort(ctxt)
    end if
    call psi_swapdata(psb_comm_status_wait_, zzero, v_rma_put%v, desc_a, info, data=psb_comm_halo_)
    if (info /= psb_success_) then
      write(psb_err_unit,*) my_rank, 'RMA put wait error:', info
      call psb_abort(ctxt)
    end if
    first_swap_rma_put = psb_wtime() - first_swap_rma_put
  end if



  ! ==================================================================
  !  8. Performance: repeat exchanges and measure timings
  ! ==================================================================
  if (my_rank == 0) then
    write(psb_out_unit,'("Timing: running ",i0," iterations for selected exchange mode(s)")') iters
  end if

  tsum_baseline = dzero
  tsum_neighbor = dzero
  tsum_neighbor_persistent = dzero
  tsum_rma_get = dzero
  tsum_rma_put = dzero

  do i = 1, iters
    if (run_baseline) then
      t0 = psb_wtime()
      call psi_swapdata(psb_comm_status_start_, zzero, v_baseline%v, desc_a, info, data=psb_comm_halo_)
      call psi_swapdata(psb_comm_status_wait_, zzero, v_baseline%v, desc_a, info, data=psb_comm_halo_)
      t1 = psb_wtime()
      dt = t1 - t0
      call psb_amx(ctxt, dt)
      tsum_baseline = tsum_baseline + dt
    end if

    if (run_neighbor) then
      t0 = psb_wtime()
      call psi_swapdata(psb_comm_status_start_, zzero, v_neighbor%v, desc_a, info, data=psb_comm_halo_)
      call psi_swapdata(psb_comm_status_wait_, zzero, v_neighbor%v, desc_a, info, data=psb_comm_halo_)
      t1 = psb_wtime()
      dt = t1 - t0
      call psb_amx(ctxt, dt)
      tsum_neighbor = tsum_neighbor + dt
    end if

    if (run_persistent) then
      t0 = psb_wtime()
      call psi_swapdata(psb_comm_status_start_, zzero, v_neighbor_persistent%v, desc_a, info, data=psb_comm_halo_)
      call psi_swapdata(psb_comm_status_wait_, zzero, v_neighbor_persistent%v, desc_a, info, data=psb_comm_halo_)
      t1 = psb_wtime()
      dt = t1 - t0
      call psb_amx(ctxt, dt)
      tsum_neighbor_persistent = tsum_neighbor_persistent + dt
    end if


    if (run_rma_get) then
      t0 = psb_wtime()
      call psi_swapdata(psb_comm_status_start_, zzero, v_rma_get%v, desc_a, info, data=psb_comm_halo_)
      call psi_swapdata(psb_comm_status_wait_, zzero, v_rma_get%v, desc_a, info, data=psb_comm_halo_)
      t1 = psb_wtime()
      dt = t1 - t0
      call psb_amx(ctxt, dt)
      tsum_rma_get = tsum_rma_get + dt
    end if


    if (run_rma_put) then
      t0 = psb_wtime()
      call psi_swapdata(psb_comm_status_start_, zzero, v_rma_put%v, desc_a, info, data=psb_comm_halo_)
      call psi_swapdata(psb_comm_status_wait_, zzero, v_rma_put%v, desc_a, info, data=psb_comm_halo_)
      t1 = psb_wtime()
      dt = t1 - t0
      call psb_amx(ctxt, dt)
      tsum_rma_put = tsum_rma_put + dt
    end if


  end do



  call psb_amx(ctxt, tsum_baseline)
  call psb_amx(ctxt, tsum_neighbor)
  call psb_amx(ctxt, tsum_neighbor_persistent)
  call psb_amx(ctxt, first_swap_baseline)
  call psb_amx(ctxt, first_swap_neighbor)
  call psb_amx(ctxt, first_swap_persistent)
  call psb_amx(ctxt, comm_setup_time_baseline)
  call psb_amx(ctxt, comm_setup_time_neighbor)
  call psb_amx(ctxt, comm_setup_time_persistent)
  call psb_amx(ctxt, comm_setup_time_rma_get)
  call psb_amx(ctxt, comm_setup_time_rma_put)



  
  if (my_rank == 0) then
    if (run_baseline) then
      write(psb_out_unit,'("  Avg baseline time  : ",es12.5)') (tsum_baseline / real(iters,psb_dpk_))
      write(psb_out_unit,'("  Tot baseline time  : ",es12.5)') tsum_baseline
      write(psb_out_unit,'("  First baseline time: ",es12.5)') first_swap_baseline
      write(psb_out_unit,'("  Baseline comm setup: ",es12.5)') comm_setup_time_baseline
    end if
    if (run_neighbor) then
      write(psb_out_unit,'("  Avg neighbor time  : ",es12.5)') (tsum_neighbor / real(iters,psb_dpk_))
      write(psb_out_unit,'("  Tot neighbor time  : ",es12.5)') tsum_neighbor
      write(psb_out_unit,'("  First neighbor time: ",es12.5)') first_swap_neighbor
      write(psb_out_unit,'("  Neighbor comm setup: ",es12.5)') comm_setup_time_neighbor
    end if
    if (run_persistent) then
      write(psb_out_unit,'("  Avg pers-neigh time: ",es12.5)') (tsum_neighbor_persistent / real(iters,psb_dpk_))
      write(psb_out_unit,'("  Tot pers-neigh time: ",es12.5)') tsum_neighbor_persistent
      write(psb_out_unit,'("  First pers-neigh time: ",es12.5)') first_swap_persistent
      write(psb_out_unit,'("  Persistent comm setup: ",es12.5)') comm_setup_time_persistent
    end if
    if (run_rma_get) then
      write(psb_out_unit,'("  Avg RMA get time   : ",es12.5)') (tsum_rma_get / real(iters,psb_dpk_))
      write(psb_out_unit,'("  Tot RMA get time   : ",es12.5)') tsum_rma_get
      write(psb_out_unit,'("  First RMA get time : ",es12.5)') first_swap_rma_get
      write(psb_out_unit,'("  RMA get comm setup : ",es12.5)') comm_setup_time_rma_get
    end if
    if (run_rma_put) then
      write(psb_out_unit,'("  Avg RMA put time   : ",es12.5)') (tsum_rma_put / real(iters,psb_dpk_))
      write(psb_out_unit,'("  Tot RMA put time   : ",es12.5)') tsum_rma_put
      write(psb_out_unit,'("  First RMA put time : ",es12.5)') first_swap_rma_put
      write(psb_out_unit,'("  RMA put comm setup : ",es12.5)') comm_setup_time_rma_put
    end if
  end if

  ! ==================================================================
  !  8. Extract results and compare
  ! ==================================================================
  result_baseline   = v_baseline%v%v
  result_neighbor   = v_neighbor%v%v
  result_persistent = v_neighbor_persistent%v%v
  result_rma_get    = v_rma_get%v%v
  result_rma_put    = v_rma_put%v%v

  ! --- Cross Checks ---
  if (run_baseline .and. run_neighbor) &
    call check_result("cross-check baseline vs neighbor", result_baseline, result_neighbor, ncol, tol)

  if (run_baseline .and. run_persistent) &
    call check_result("cross-check baseline vs pers-nei", result_baseline, result_persistent, ncol, tol)

  if (run_baseline .and. run_rma_get) &
    call check_result("cross-check baseline vs rma-get", result_baseline, result_rma_get, ncol, tol)
  
  if (run_baseline .and. run_rma_put) &
    call check_result("cross-check baseline vs rma-put", result_baseline, result_rma_put, ncol, tol)

  ! --- Absolute Correctness Checks against Expected ---
  if (run_baseline) &
    call check_result("baseline absolute correctness", result_baseline, expected, ncol, tol)
    
  if (run_neighbor) &
    call check_result("neighbor absolute correctness", result_neighbor, expected, ncol, tol)
    
  if (run_persistent) &
    call check_result("pers-neigh absolute correctness", result_persistent, expected, ncol, tol)

  if (run_rma_get) &
    call check_result("rma_get absolute correctness", result_rma_get, expected, ncol, tol)

  if (run_rma_put) &
    call check_result("rma_put absolute correctness", result_rma_put, expected, ncol, tol)

  ! ==================================================================
  !  9. Summary
  ! ==================================================================
  call psb_barrier(ctxt)
  if (my_rank == 0) then
    write(psb_out_unit,'("================================================")')
    write(psb_out_unit,'("  Results: ",i0," / ",i0," tests passed")') n_pass, n_total
    if (n_pass == n_total) then
      write(psb_out_unit,'("  STATUS: ALL PASSED")')
    else
      write(psb_out_unit,'("  STATUS: SOME FAILURES")')
    end if
    write(psb_out_unit,'("================================================")')
  end if
  call psb_barrier(ctxt)

  ! ==================================================================
  !  10. Cleanup
  ! ==================================================================
  deallocate(result_baseline, result_neighbor, result_persistent, expected, glob_col)
  
9999 call psb_gefree(v_baseline, desc_a, info)
  call psb_gefree(v_neighbor, desc_a, info)
  call psb_gefree(v_neighbor_persistent, desc_a, info)
  if (mat_allocated) call psb_spfree(a_mat, desc_a, info)
  call psb_cdfree(desc_a, info)
  call psb_exit(ctxt)


contains 

  ! Helper routine to compare two arrays, reduce the error across MPI ranks, and print the result
  subroutine check_result(test_name, arr1, arr2, n, tolerance)
    character(len=*), intent(in) :: test_name
    complex(psb_dpk_), intent(in)  :: arr1(:), arr2(:)
    integer(psb_ipk_), intent(in):: n
    real(psb_dpk_), intent(in) :: tolerance
    real(psb_dpk_)             :: err_val

    n_total = n_total + 1
    
    if (n > 0) then
      err_val = maxval(abs(arr1(1:n) - arr2(1:n)))
    else
      err_val = dzero
    end if
    
    ! Get max error across all processes
    call psb_amx(ctxt, err_val)

    if (my_rank == 0) then
      if ((err_val >= dzero) .and. (err_val < tolerance)) then
        write(psb_out_unit,'("  [PASS] ", a, t45, ": err = ",es12.5)') test_name, err_val
        n_pass = n_pass + 1
      else
        write(psb_out_unit,'("  [FAIL] ", a, t45, ": err = ",es12.5)') test_name, err_val
      end if
    end if
  end subroutine check_result

end program psb_z_comm_test
