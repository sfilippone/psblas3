!
! Test program for D-type halo exchange: baseline vs neighbor topology.
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
! Run with:  mpirun -np <P> ./test_halo_new
!
program psb_comm_test
  use psb_base_mod
  use psi_mod
  use psb_comm_factory_mod, only: psb_comm_set, psb_comm_free
  use psb_comm_schemes_mod, only: psb_comm_ineighbor_alltoallv_, psb_comm_persistent_ineighbor_alltoallv_, &
    & psb_comm_isend_irecv_
  use psb_comm_schemes_mod, only: psb_comm_status_start_, psb_comm_status_wait_, psb_comm_status_unknown_
  implicit none

  ! ---- parameters ----
  integer(psb_ipk_)                     :: idim
  integer(psb_ipk_)                     :: argc
  integer(psb_ipk_)                     :: iters
  character(len=32)                     :: arg
  character(len=16)                     :: mode

  ! ---- descriptor / context ----
  type(psb_ctxt_type)                   :: ctxt
  type(psb_desc_type)                   :: desc_a
  integer(psb_ipk_)                     :: my_rank, np, info, i, nr, number_of_local_rows
  integer(psb_lpk_)                     :: m, nt
  integer(psb_lpk_), allocatable        :: myidx(:)

  ! ---- vectors ----
  type(psb_d_vect_type)                 :: v_baseline, v_neighbor, v_neighbor_persistent

  ! ---- temporary / comparison arrays ----
  real(psb_dpk_), allocatable           :: vals(:)
  real(psb_dpk_), allocatable           :: result_baseline(:), result_neighbor(:), result_persistent(:)
  real(psb_dpk_), allocatable           :: expected(:)

  ! ---- halo index bookkeeping ----
  integer(psb_ipk_)                     :: nrow, ncol, num_neighbors, send_indexes, receive_indexes
  class(psb_i_base_vect_type), pointer  :: halo_indexes

  ! ---- error / reporting ----
  integer(psb_ipk_)                     :: n_pass, n_total, imode
  logical                               :: comm_ok
  real(psb_dpk_)                        :: err, tol
  real(psb_dpk_)                        :: t0, t1, dt, tsum_baseline, tsum_neighbor, tsum_neighbor_persistent
  integer(psb_lpk_), allocatable        :: glob_col(:)
  character(len=40)                     :: name

  name    = 'test_halo_new'
  tol     = 1.0d-12
  n_pass  = 0
  n_total = 0
  iters = 5
  mode = 'both'

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
    end if
  end do

  ! parse optional mode flag
  do i = 1, argc
    call get_command_argument(i, arg)
    if (trim(arg) == '--mode') then
      if (i < argc) then
        call get_command_argument(i+1, arg)
        read(arg, *) mode
      end if
    end if
  end do

  if (idim <= 0) then
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
    write(psb_out_unit,'("  Test: D-type halo  baseline vs neighbor topo")')
    write(psb_out_unit,'("  Processes : ",i0)') np
    write(psb_out_unit,'("  Grid      : ",i0," x ",i0," x ",i0)') idim,idim,idim
    write(psb_out_unit,'("================================================")')
  end if

  ! ==================================================================
  !  2. Build descriptor with 7-point stencil connectivity
  ! ==================================================================
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

  nrow = desc_a%get_local_rows()   ! owned
  ncol = desc_a%get_local_cols()   ! owned + halo

  ! ==================================================================
  !  3. Allocate two D vectors (scratch) and fill owned entries
  ! ==================================================================
  call psb_geall(v_baseline, desc_a, info)
  call psb_geall(v_neighbor, desc_a, info)
  call psb_geall(v_neighbor_persistent, desc_a, info)
  call psb_geasb(v_baseline, desc_a, info, scratch=.true.)
  call psb_geasb(v_neighbor, desc_a, info, scratch=.true.)
  call psb_geasb(v_neighbor_persistent, desc_a, info, scratch=.true.)

  ! Fill owned entries with the global index value
  allocate(vals(ncol))
  vals = dzero
  do i = 1, number_of_local_rows
    vals(i) = real(myidx(i), psb_dpk_)
  end do
  call v_baseline%set_vect(vals)
  call v_neighbor%set_vect(vals)
  call v_neighbor_persistent%set_vect(vals)
  deallocate(vals)

  ! ==================================================================
  !  4. Build the expected result for halo positions
  !     glob_col(j) = global index of local column j
  !     After halo exchange every position j should hold glob_col(j).
  ! ==================================================================
  allocate(glob_col(ncol), expected(ncol))
  glob_col = desc_a%get_global_indices(owned=.false.)
  do i = 1, ncol
    expected(i) = real(glob_col(i), psb_dpk_)
  end do

  ! ==================================================================
  !  6. Baseline halo exchange  (Isend/Irecv in one call)
  ! ==================================================================
  ! v_baseline%v is a psb_d_base_vect_type
  call psi_swapdata( &
    swap_status=psb_comm_status_start_, &
    beta=dzero, &
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
      beta=dzero, &
      y=v_baseline%v, &
      desc_a=desc_a, &
      info=info, &
      data=psb_comm_halo_)
  if (info /= psb_success_) then
    write(psb_err_unit,*) my_rank, 'baseline swap error:', info
    call psb_abort(ctxt)
  end if


  ! ==================================================================
  !  7. Neighbor topology halo exchange  (start + wait)
  ! ==================================================================
  call psb_comm_set(psb_comm_ineighbor_alltoallv_, v_neighbor%v%comm_handle, info)
  if (info /= 0) then
    write(psb_err_unit,*) my_rank, 'psb_comm_set neighbor error:', info
    call psb_abort(ctxt)
  end if
  call psi_swapdata(psb_comm_status_start_, dzero, v_neighbor%v, desc_a, info, data=psb_comm_halo_)
  if (info /= psb_success_) then
    write(psb_err_unit,*) my_rank, 'neighbor start error:', info
    call psb_abort(ctxt)
  end if

  call psi_swapdata(psb_comm_status_wait_, dzero, v_neighbor%v, desc_a, info, data=psb_comm_halo_)
  if (info /= psb_success_) then
    write(psb_err_unit,*) my_rank, 'neighbor wait error:', info
    call psb_abort(ctxt)
  end if

  ! ==================================================================
  !  7b. Persistent-neighbor halo exchange  (start + wait)
  ! ==================================================================
  call psb_comm_set(psb_comm_persistent_ineighbor_alltoallv_, v_neighbor_persistent%v%comm_handle, info)
  if (info /= 0) then
    write(psb_err_unit,*) my_rank, 'psb_comm_set persistent-neighbor error:', info
    call psb_abort(ctxt)
  end if
  call psi_swapdata(psb_comm_status_start_, dzero, v_neighbor_persistent%v, desc_a, info, data=psb_comm_halo_)
  if (info /= psb_success_) then
    write(psb_err_unit,*) my_rank, 'persistent-neighbor start error:', info
    call psb_abort(ctxt)
  end if
  call psi_swapdata(psb_comm_status_wait_, dzero, v_neighbor_persistent%v, desc_a, info, data=psb_comm_halo_)
  if (info /= psb_success_) then
    write(psb_err_unit,*) my_rank, 'persistent-neighbor wait error:', info
    call psb_abort(ctxt)
  end if

  ! ==================================================================
  !  8. Performance: repeat exchanges and measure timings
  ! ==================================================================
  if (my_rank == 0) then
    write(psb_out_unit,'("Timing: running ",i0," iterations for baseline, neighbor and persistent-neighbor")') iters
  end if

  tsum_baseline = 0.0_psb_dpk_
  tsum_neighbor = 0.0_psb_dpk_
  tsum_neighbor_persistent = 0.0_psb_dpk_

  call psb_comm_set(psb_comm_isend_irecv_, v_baseline%v%comm_handle, info)
  call psb_comm_set(psb_comm_ineighbor_alltoallv_, v_neighbor%v%comm_handle, info)
  call psb_comm_set(psb_comm_persistent_ineighbor_alltoallv_, v_neighbor_persistent%v%comm_handle, info)

  ! ---- Comm check: verify selected communication schemes ----
  n_total = n_total + 1
  comm_ok = allocated(v_baseline%v%comm_handle) .and. allocated(v_neighbor%v%comm_handle) .and. &
    & allocated(v_neighbor_persistent%v%comm_handle)

  if (comm_ok) then
    comm_ok = (v_baseline%v%comm_handle%comm_type == psb_comm_isend_irecv_) .and. &
      & (v_neighbor%v%comm_handle%comm_type == psb_comm_ineighbor_alltoallv_) .and. &
      & (v_neighbor_persistent%v%comm_handle%comm_type == psb_comm_persistent_ineighbor_alltoallv_)
  end if

  if (my_rank == 0) then
    if (comm_ok) then
      write(psb_out_unit,'("  [PASS] comm scheme selection            : baseline/neighbor/persistent OK")')
      n_pass = n_pass + 1
    else
      write(psb_out_unit,'("  [FAIL] comm scheme selection            : unexpected comm_type mapping")')
    end if
  end if

  do i = 1, iters
    ! baseline timing
    t0 = psb_wtime()
    call psi_swapdata( &
      swap_status=psb_comm_status_start_, &
      beta=dzero, &
      y=v_baseline%v, &
      desc_a=desc_a, &
      info=info, &
      data=psb_comm_halo_)
    call psi_swapdata( &
      swap_status=psb_comm_status_wait_, &
      beta=dzero, &
      y=v_baseline%v, &
      desc_a=desc_a, &
      info=info, &
      data=psb_comm_halo_)
    t1 = psb_wtime()
    dt = t1 - t0
    call psb_amx(ctxt, dt)
    tsum_baseline = tsum_baseline + dt

    ! neighbor timing (start + wait)
    t0 = psb_wtime()
    call psi_swapdata(psb_comm_status_start_, dzero, v_neighbor%v, desc_a, info, data=psb_comm_halo_)
    call psi_swapdata(psb_comm_status_wait_, dzero, v_neighbor%v, desc_a, info, data=psb_comm_halo_)
    t1 = psb_wtime()
    dt = t1 - t0
    call psb_amx(ctxt, dt)
    tsum_neighbor = tsum_neighbor + dt

    ! persistent-neighbor timing (start + wait)
    t0 = psb_wtime()
    call psi_swapdata(psb_comm_status_start_, dzero, v_neighbor_persistent%v, desc_a, info, data=psb_comm_halo_)
    call psi_swapdata(psb_comm_status_wait_, dzero, v_neighbor_persistent%v, desc_a, info, data=psb_comm_halo_)
    t1 = psb_wtime()
    dt = t1 - t0
    call psb_amx(ctxt, dt)
    tsum_neighbor_persistent = tsum_neighbor_persistent + dt
  end do

  if (my_rank == 0) then
    write(psb_out_unit,'("  Avg baseline time  : ",es12.5)') (tsum_baseline / real(iters,psb_dpk_))
    write(psb_out_unit,'("  Tot baseline time  : ",es12.5)') tsum_baseline
    write(psb_out_unit,'("  Avg neighbor time  : ",es12.5)') (tsum_neighbor / real(iters,psb_dpk_))
    write(psb_out_unit,'("  Tot neighbor time  : ",es12.5)') tsum_neighbor
    write(psb_out_unit,'("  Avg pers-neigh time: ",es12.5)') (tsum_neighbor_persistent / real(iters,psb_dpk_))
    write(psb_out_unit,'("  Tot pers-neigh time: ",es12.5)') tsum_neighbor_persistent
  end if

  ! ==================================================================
  !  8. Extract results and compare
  ! ==================================================================
  result_baseline = v_baseline%get_vect()
  result_neighbor = v_neighbor%get_vect()
  result_persistent = v_neighbor_persistent%get_vect()

  ! ---- Test 1: cross-check baseline vs neighbor (all entries) ----
  n_total = n_total + 1
  err = maxval(abs(result_baseline(1:ncol) - result_neighbor(1:ncol)))
  call psb_amx(ctxt, err)
  if (my_rank == 0) then
    if (err < tol) then
      write(psb_out_unit,'("  [PASS] cross-check baseline vs neighbor : err = ",es12.5)') err
      n_pass = n_pass + 1
    else
      write(psb_out_unit,'("  [FAIL] cross-check baseline vs neighbor : err = ",es12.5)') err
    end if
  end if

  ! ---- Test 2: baseline absolute correctness (halo = global index) ----
  n_total = n_total + 1
  err = maxval(abs(result_baseline(1:ncol) - expected(1:ncol)))
  call psb_amx(ctxt, err)
  if (my_rank == 0) then
    if (err < tol) then
      write(psb_out_unit,'("  [PASS] baseline absolute correctness    : err = ",es12.5)') err
      n_pass = n_pass + 1
    else
      write(psb_out_unit,'("  [FAIL] baseline absolute correctness    : err = ",es12.5)') err
    end if
  end if

  ! ---- Test 3: neighbor absolute correctness (halo = global index) ----
  n_total = n_total + 1
  err = maxval(abs(result_neighbor(1:ncol) - expected(1:ncol)))
  call psb_amx(ctxt, err)
  if (my_rank == 0) then
    if (err < tol) then
      write(psb_out_unit,'("  [PASS] neighbor absolute correctness    : err = ",es12.5)') err
      n_pass = n_pass + 1
    else
      write(psb_out_unit,'("  [FAIL] neighbor absolute correctness    : err = ",es12.5)') err
    end if
  end if

  ! ---- Test 4: cross-check baseline vs persistent-neighbor (all entries) ----
  n_total = n_total + 1
  err = maxval(abs(result_baseline(1:ncol) - result_persistent(1:ncol)))
  call psb_amx(ctxt, err)
  if (my_rank == 0) then
    if (err < tol) then
      write(psb_out_unit,'("  [PASS] cross-check baseline vs pers-nei : err = ",es12.5)') err
      n_pass = n_pass + 1
    else
      write(psb_out_unit,'("  [FAIL] cross-check baseline vs pers-nei : err = ",es12.5)') err
    end if
  end if

  ! ---- Test 5: persistent-neighbor absolute correctness ----
  n_total = n_total + 1
  err = maxval(abs(result_persistent(1:ncol) - expected(1:ncol)))
  call psb_amx(ctxt, err)
  if (my_rank == 0) then
    if (err < tol) then
      write(psb_out_unit,'("  [PASS] pers-neigh absolute correctness  : err = ",es12.5)') err
      n_pass = n_pass + 1
    else
      write(psb_out_unit,'("  [FAIL] pers-neigh absolute correctness  : err = ",es12.5)') err
    end if
  end if

  ! ---- Test 6: repeat neighbor exchange (topology reuse) ----
  ! Reset halo entries to zero, run again, and check
  do i = nrow+1, ncol
    result_neighbor(i) = dzero
  end do
  call v_neighbor%set_vect(result_neighbor)

  call psi_swapdata(psb_comm_status_start_, dzero, v_neighbor%v, desc_a, info, data=psb_comm_halo_)
  call psi_swapdata(psb_comm_status_wait_, dzero, v_neighbor%v, desc_a, info, data=psb_comm_halo_)

  result_neighbor = v_neighbor%get_vect()
  n_total = n_total + 1
  err = maxval(abs(result_neighbor(1:ncol) - expected(1:ncol)))
  call psb_amx(ctxt, err)
  if (my_rank == 0) then
    if (err < tol) then
      write(psb_out_unit,'("  [PASS] neighbor topology reuse          : err = ",es12.5)') err
      n_pass = n_pass + 1
    else
      write(psb_out_unit,'("  [FAIL] neighbor topology reuse          : err = ",es12.5)') err
    end if
  end if

  ! ---- Test 7: repeat persistent-neighbor exchange (buffer reuse) ----
  do i = nrow+1, ncol
    result_persistent(i) = dzero
  end do
  call v_neighbor_persistent%set_vect(result_persistent)

  call psi_swapdata(psb_comm_status_start_, dzero, v_neighbor_persistent%v, desc_a, info, data=psb_comm_halo_)
  call psi_swapdata(psb_comm_status_wait_, dzero, v_neighbor_persistent%v, desc_a, info, data=psb_comm_halo_)

  result_persistent = v_neighbor_persistent%get_vect()
  n_total = n_total + 1
  err = maxval(abs(result_persistent(1:ncol) - expected(1:ncol)))
  call psb_amx(ctxt, err)
  if (my_rank == 0) then
    if (err < tol) then
      write(psb_out_unit,'("  [PASS] pers-neigh buffer reuse          : err = ",es12.5)') err
      n_pass = n_pass + 1
    else
      write(psb_out_unit,'("  [FAIL] pers-neigh buffer reuse          : err = ",es12.5)') err
    end if
  end if

  ! ==================================================================
  !  9. Summary
  ! ==================================================================
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

  ! ==================================================================
  !  10. Cleanup
  ! ==================================================================
  deallocate(result_baseline, result_neighbor, result_persistent, expected, glob_col)
  
9999 call psb_gefree(v_baseline, desc_a, info)
  call psb_gefree(v_neighbor, desc_a, info)
  call psb_gefree(v_neighbor_persistent, desc_a, info)
  call psb_cdfree(desc_a, info)
  call psb_exit(ctxt)

end program psb_comm_test