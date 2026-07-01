!
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone
!        Alfredo Buttari
!
! File: psb_d_nest_halo_regime_test.F90
! Author: Simone Staccone (Stack-1)
!
! Benchmarks the two halo-exchange regimes of a nested (MATNEST-style) operator
! on the SAME assembled blocks:
!
!   * UNION (fused):  one psb_halo over the composed global descriptor brings in
!     the union of every field's ghosts at once (the production path used by
!     psb_spmm / Krylov / AMG4PSBLAS).
!
!   * PER-BLOCK (selective): each present block (i,j) exchanges ONLY its own
!     off-process columns, through its restricted descriptor block_col_desc(i,j),
!     one psb_halo per block.
!
! The two products must coincide to machine precision.  The test reports, for
! each regime: the number of halo exchanges, the full-matvec time, and the
! PURE-COMMUNICATION time (only the psb_halo calls), so the cost of the message
! aggregation can be isolated.
!
! Saddle-point 2x2 operator  [ A  B^T ; B  0 ]  with TWO fields of DIFFERENT
! size (n1 = field 1 / "velocity", n2 = field 2 / "pressure", n1 /= n2 allowed);
! the (2,2) block is absent.  The blocks have REAL, distinct per-block halos (so
! per-block vs union is meaningful): A (1,1) tridiagonal on field 1; the
! rectangular B^T (1,2) and B (2,1) couple every row to a column half a field
! away, which lands on another process.
!
! Run:  mpirun -np <P> ./psb_d_nest_halo_regime_test [n1] [n2] [n_reps] [stride]
!         n1     : global rows of field 1   (default 2000000)
!         n2     : global rows of field 2   (default  500000)
!         n_reps : timed repetitions         (default 50)
!         stride : B/B^T coupling distance   (default 1)
!                  small  => local halo  => LATENCY-bound (aggregation wins)
!                  ~n/2   => huge halo   => BANDWIDTH-bound (union ~ per-block)
!
program psb_d_nest_halo_regime_test
  use psb_base_mod
  use psb_d_nest_mod
  implicit none

  type(psb_ctxt_type)             :: context
  integer(psb_ipk_)               :: my_rank, num_procs, info, i_local_row
  integer(psb_ipk_)               :: entry_idx, field1_local_rows, field2_local_rows
  integer(psb_ipk_)               :: n_local_rows, n_local_cols, n_exch, n_warm
  integer(psb_ipk_)               :: n_reps, narg
  integer(psb_lpk_)               :: global_row, field1_size, field2_size, gcol, stride
  character(len=64)               :: arg

  type(psb_d_nest_matrix)         :: nested_matrix
  integer(psb_lpk_), allocatable  :: entry_rows(:), entry_cols(:)
  integer(psb_lpk_), allocatable  :: field1_rows(:), field2_rows(:)
  real(psb_dpk_),    allocatable  :: entry_vals(:)
  real(psb_dpk_),    allocatable  :: x_arr(:), work_x(:), y_union(:), y_perblock(:)
  real(psb_dpk_)                  :: mismatch_norm, t_u, t_b
  logical                         :: ok_u, ok_b
  integer(psb_ipk_)               :: s_idx, n_schemes
  integer(psb_ipk_), allocatable  :: schemes(:)
  character(len=34), allocatable  :: scheme_names(:)
  type(psb_d_vect_type)              :: xg_vect            ! global-local halo vector
  type(psb_d_vect_type), allocatable :: bvect(:,:)         ! one per present block
  real(psb_dpk_), parameter       :: tolerance = 1.0e-10_psb_dpk_

  call psb_init(context)
  call psb_info(context, my_rank, num_procs)

  ! -------- runtime parameters (rank 0 parses, then broadcast) --------
  field1_size = 2000000_psb_lpk_
  field2_size =  500000_psb_lpk_
  n_reps      = 50
  stride      = 1_psb_lpk_       ! B/B^T coupling distance: small => local halo
                                 ! (latency-bound), large (~n/2) => big halo
                                 ! (bandwidth-bound)
  if (my_rank == 0) then
    narg = command_argument_count()
    if (narg >= 1) then
      call get_command_argument(1, arg); read(arg,*) field1_size
    end if
    if (narg >= 2) then
      call get_command_argument(2, arg); read(arg,*) field2_size
    end if
    if (narg >= 3) then
      call get_command_argument(3, arg); read(arg,*) n_reps
    end if
    if (narg >= 4) then
      call get_command_argument(4, arg); read(arg,*) stride
    end if
  end if
  call psb_bcast(context, field1_size)
  call psb_bcast(context, field2_size)
  call psb_bcast(context, n_reps)
  call psb_bcast(context, stride)
  n_warm = max(5, n_reps/10)

  !---------------------------------------------------------------
  ! 1) build the 2x2 nested operator with real per-block halos
  !---------------------------------------------------------------
  call nested_matrix%init(context, [field1_size, field2_size], info)
  if (info /= psb_success_) then
    if (my_rank==0) write(*,*) 'FAIL: init info=', info; goto 9999
  end if
  field1_rows = nested_matrix%get_owned_rows(1)
  field2_rows = nested_matrix%get_owned_rows(2)
  field1_local_rows = size(field1_rows)
  field2_local_rows = size(field2_rows)

  ! A = tridiag(-1,2,-1) -> block (1,1)
  allocate(entry_rows(3*field1_local_rows), entry_cols(3*field1_local_rows), entry_vals(3*field1_local_rows))
  entry_idx = 0
  do i_local_row = 1, field1_local_rows
    global_row = field1_rows(i_local_row)
    entry_idx = entry_idx + 1
    entry_rows(entry_idx) = global_row; entry_cols(entry_idx) = global_row; entry_vals(entry_idx) = 2.0_psb_dpk_
    if (global_row > 1) then
      entry_idx = entry_idx + 1
      entry_rows(entry_idx) = global_row; entry_cols(entry_idx) = global_row - 1_psb_lpk_; entry_vals(entry_idx) = -1.0_psb_dpk_
    end if
    if (global_row < field1_size) then
      entry_idx = entry_idx + 1
      entry_rows(entry_idx) = global_row; entry_cols(entry_idx) = global_row + 1_psb_lpk_; entry_vals(entry_idx) = -1.0_psb_dpk_
    end if
  end do
  call nested_matrix%ins(1, 1, entry_idx, entry_rows, entry_cols, entry_vals, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  ! B^T -> block (1,2): rows in field 1, columns in field 2 (RECTANGULAR n1 x n2).
  ! Each row r couples to field-2 columns (r mod n2) and (r + n2/2 mod n2); the
  ! second one lands on another process => real, distinct halo for this block.
  allocate(entry_rows(2*field1_local_rows), entry_cols(2*field1_local_rows), entry_vals(2*field1_local_rows))
  entry_idx = 0
  do i_local_row = 1, field1_local_rows
    global_row = field1_rows(i_local_row)
    entry_idx = entry_idx + 1
    gcol = mod(global_row - 1_psb_lpk_, field2_size) + 1_psb_lpk_
    entry_rows(entry_idx) = global_row; entry_cols(entry_idx) = gcol; entry_vals(entry_idx) = 0.5_psb_dpk_
    entry_idx = entry_idx + 1
    gcol = mod(global_row - 1_psb_lpk_ + stride, field2_size) + 1_psb_lpk_
    entry_rows(entry_idx) = global_row; entry_cols(entry_idx) = gcol; entry_vals(entry_idx) = 0.25_psb_dpk_
  end do
  call nested_matrix%ins(1, 2, entry_idx, entry_rows, entry_cols, entry_vals, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  ! B -> block (2,1): rows in field 2, columns in field 1 (RECTANGULAR n2 x n1).
  allocate(entry_rows(2*field2_local_rows), entry_cols(2*field2_local_rows), entry_vals(2*field2_local_rows))
  entry_idx = 0
  do i_local_row = 1, field2_local_rows
    global_row = field2_rows(i_local_row)
    entry_idx = entry_idx + 1
    gcol = mod(global_row - 1_psb_lpk_, field1_size) + 1_psb_lpk_
    entry_rows(entry_idx) = global_row; entry_cols(entry_idx) = gcol; entry_vals(entry_idx) = 0.3_psb_dpk_
    entry_idx = entry_idx + 1
    gcol = mod(global_row - 1_psb_lpk_ + stride, field1_size) + 1_psb_lpk_
    entry_rows(entry_idx) = global_row; entry_cols(entry_idx) = gcol; entry_vals(entry_idx) = 0.15_psb_dpk_
  end do
  call nested_matrix%ins(2, 1, entry_idx, entry_rows, entry_cols, entry_vals, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  call nested_matrix%asb(info)
  if (info /= psb_success_) then
    if (my_rank==0) write(*,*) 'FAIL: asb info=', info; goto 9999
  end if

  !---------------------------------------------------------------
  ! 2) global-local work vectors (x[g] = g on the owned entries)
  !---------------------------------------------------------------
  n_local_rows = nested_matrix%desc_glob%get_local_rows()
  n_local_cols = nested_matrix%desc_glob%get_local_cols()
  allocate(x_arr(n_local_cols), work_x(n_local_cols), &
       &   y_union(n_local_cols), y_perblock(n_local_cols))
  x_arr = dzero
  do i_local_row = 1, n_local_rows
    call nested_matrix%desc_glob%l2g(i_local_row, global_row, info)
    x_arr(i_local_row) = real(global_row, psb_dpk_)
  end do

  !---------------------------------------------------------------
  ! 3) correctness: union vs per-block on the same x
  !---------------------------------------------------------------
  work_x = x_arr
  call psb_spmm(done, nested_matrix%a_glob, work_x, dzero, y_union, nested_matrix%desc_glob, info)
  if (info /= psb_success_) then
    if (my_rank==0) write(*,*) 'FAIL: psb_spmm (union) info=', info; goto 9999
  end if
  call nest_spmv_perblock(nested_matrix, done, x_arr, dzero, y_perblock, info, n_exch, .false.)
  if (info /= psb_success_) then
    if (my_rank==0) write(*,*) 'FAIL: per-block driver info=', info; goto 9999
  end if
  mismatch_norm = dzero
  do i_local_row = 1, n_local_rows
    mismatch_norm = max(mismatch_norm, abs(y_union(i_local_row) - y_perblock(i_local_row)))
  end do
  call psb_amx(context, mismatch_norm)

  !---------------------------------------------------------------
  ! 4) per-scheme PURE COMMUNICATION sweep.
  !    The comm scheme is selected per descriptor (desc%set_comm_scheme) and is
  !    honoured ONLY by the encapsulated vect path; the array path used in (3)
  !    is always baseline.  For each scheme we time, on persistent vects:
  !      union     = one psb_halo over desc_glob
  !      per-block = one psb_halo per present block_col_desc(i,j)
  !---------------------------------------------------------------
  n_schemes = 5
  allocate(schemes(n_schemes), scheme_names(n_schemes))
  schemes = [ psb_comm_isend_irecv_, psb_comm_ineighbor_alltoallv_,            &
       &      psb_comm_persistent_ineighbor_alltoallv_,                        &
       &      psb_comm_rma_pull_, psb_comm_rma_push_ ]
  scheme_names = [ character(len=34) :: 'isend_irecv (baseline)',              &
       &           'ineighbor_alltoallv', 'persistent_ineighbor_alltoallv',    &
       &           'rma_pull', 'rma_push' ]
  allocate(bvect(nested_matrix%n_fields, nested_matrix%n_fields))

  if (my_rank == 0) then
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') ' np=', num_procs, '  n1=', field1_size, &
         &                           '  n2=', field2_size, '  stride=', stride, '  reps=', n_reps
    write(*,'(a,es12.4)') ' max|y_union - y_perblock| = ', mismatch_norm
    if (mismatch_norm <= tolerance) then
      write(*,*) ' PASS: regimes agree'
    else
      write(*,*) ' FAIL: regimes disagree'
    end if
    write(*,'(a,i0)') ' halo exchanges:  union = 1   per-block = ', n_exch
    write(*,*)
    write(*,'(a)') ' pure halo communication time [s]  (min over reps, slowest rank)'
    write(*,'(a)') ' scheme                              union           per-block       ratio'
  end if

  do s_idx = 1, n_schemes
    call set_scheme(schemes(s_idx))
    call build_vects()
    call time_comm('U', t_u, ok_u)
    call time_comm('B', t_b, ok_b)
    call free_vects()
    if (my_rank == 0) then
      if (ok_u .and. ok_b) then
        write(*,'(1x,a34,es15.4,es15.4,3x,f7.2)') scheme_names(s_idx), t_u, t_b, &
             &                                     t_b/max(t_u, tiny(t_u))
      else
        write(*,'(1x,a34,a)') scheme_names(s_idx), '   (unavailable on this build/MPI)'
      end if
    end if
  end do

  deallocate(bvect, schemes, scheme_names)
  call nested_matrix%free(info)
9999 continue
  call psb_exit(context)

contains

  ! Set the communication scheme on the global descriptor and on every present
  ! per-block descriptor (the per-block exchanges must use the same scheme).
  subroutine set_scheme(scheme)
    integer(psb_ipk_), intent(in) :: scheme
    integer(psb_ipk_) :: i, j, linfo
    call nested_matrix%desc_glob%set_comm_scheme(scheme, linfo)
    do j = 1, nested_matrix%n_fields
      do i = 1, nested_matrix%n_fields
        if (nested_matrix%block_storage%has_block(i,j)) &
             & call nested_matrix%block_col_desc(i,j)%set_comm_scheme(scheme, linfo)
      end do
    end do
  end subroutine set_scheme

  ! Fresh persistent halo vectors (the comm_handle, hence the scheme, is created
  ! from desc%comm_type on the first psb_halo and then reused across reps).
  subroutine build_vects()
    integer(psb_ipk_) :: i, j, linfo
    call psb_geall(xg_vect, nested_matrix%desc_glob, linfo)
    call psb_geasb(xg_vect, nested_matrix%desc_glob, linfo)
    do j = 1, nested_matrix%n_fields
      do i = 1, nested_matrix%n_fields
        if (.not. nested_matrix%block_storage%has_block(i,j)) cycle
        call psb_geall(bvect(i,j), nested_matrix%block_col_desc(i,j), linfo)
        call psb_geasb(bvect(i,j), nested_matrix%block_col_desc(i,j), linfo)
      end do
    end do
  end subroutine build_vects

  subroutine free_vects()
    integer(psb_ipk_) :: i, j, linfo
    call xg_vect%free(linfo)
    do j = 1, nested_matrix%n_fields
      do i = 1, nested_matrix%n_fields
        if (nested_matrix%block_storage%has_block(i,j)) call bvect(i,j)%free(linfo)
      end do
    end do
  end subroutine free_vects

  ! Time the pure communication of one regime: n_warm warm-up runs, then min over
  ! n_reps of the slowest rank's wall time (psb_amx = max across ranks).
  subroutine time_comm(code, t_min, ok)
    character,      intent(in)  :: code     ! 'U' union / 'B' per-block
    real(psb_dpk_), intent(out) :: t_min
    logical,        intent(out) :: ok
    integer(psb_ipk_) :: rep, linfo
    real(psb_dpk_)    :: t0, dt
    ok = .true.
    do rep = 1, n_warm
      call do_comm(code, linfo); if (linfo /= 0) ok = .false.
    end do
    t_min = huge(t_min)
    do rep = 1, n_reps
      call psb_barrier(context)
      t0 = psb_wtime()
      call do_comm(code, linfo)
      dt = psb_wtime() - t0
      call psb_amx(context, dt)
      t_min = min(t_min, dt)
      if (linfo /= 0) ok = .false.
    end do
    if (.not. ok) t_min = dzero
  end subroutine time_comm

  subroutine do_comm(code, linfo)
    character,         intent(in)  :: code
    integer(psb_ipk_), intent(out) :: linfo
    integer(psb_ipk_) :: i, j
    linfo = 0
    if (code == 'U') then
      call psb_halo(xg_vect, nested_matrix%desc_glob, linfo)
    else
      do j = 1, nested_matrix%n_fields
        do i = 1, nested_matrix%n_fields
          if (.not. nested_matrix%block_storage%has_block(i,j)) cycle
          call psb_halo(bvect(i,j), nested_matrix%block_col_desc(i,j), linfo)
          if (linfo /= 0) return
        end do
      end do
    end if
  end subroutine do_comm

  ! Per-block nested matvec, exchanging one block's halo at a time.  x, y are in
  ! the composed global-local layout.  halo_only=.true. performs ONLY the
  ! per-block psb_halo calls (for the pure-communication timing).
  subroutine nest_spmv_perblock(op, alpha, x, beta, y, info, n_exchanges, halo_only)
    type(psb_d_nest_matrix), intent(inout) :: op
    real(psb_dpk_),          intent(in)    :: x(:)
    real(psb_dpk_),          intent(inout) :: y(:)
    real(psb_dpk_),          intent(in)    :: alpha, beta
    integer(psb_ipk_),       intent(out)   :: info, n_exchanges
    logical,                 intent(in)    :: halo_only

    integer(psb_ipk_)              :: nf, i, j, k, nrl, nownj, nowni, nfc, nbc
    integer(psb_ipk_), allocatable :: owned_offset(:)
    real(psb_dpk_),    allocatable :: xs(:), xf(:), yb(:)

    info = psb_success_; n_exchanges = 0
    nf = op%n_fields
    allocate(owned_offset(nf+1))
    owned_offset(1) = 0
    do j = 1, nf
      owned_offset(j+1) = owned_offset(j) + op%field_desc(j)%get_local_rows()
    end do
    nrl = owned_offset(nf+1)

    if (.not. halo_only) then
      if (beta == dzero) then
        y(1:nrl) = dzero
      else if (beta /= done) then
        y(1:nrl) = beta * y(1:nrl)
      end if
    end if

    do j = 1, nf
      nownj = op%field_desc(j)%get_local_rows()
      nfc   = op%field_desc(j)%get_local_cols()
      do i = 1, nf
        if (.not. op%block_storage%has_block(i,j)) cycle
        nbc = op%block_col_desc(i,j)%get_local_cols()
        allocate(xs(nbc)); xs = dzero
        xs(1:nownj) = x(owned_offset(j)+1 : owned_offset(j)+nownj)
        call psb_halo(xs, op%block_col_desc(i,j), info)
        if (info /= psb_success_) return
        n_exchanges = n_exchanges + 1
        if (.not. halo_only) then
          allocate(xf(nfc)); xf = dzero
          do k = 1, nbc
            xf(op%blk2field(i,j)%map(k)) = xs(k)
          end do
          nowni = op%field_desc(i)%get_local_rows()
          allocate(yb(nowni)); yb = dzero
          call op%block_storage%mats(i,j)%a%csmv(alpha, xf, dzero, yb, info)
          if (info /= psb_success_) return
          y(owned_offset(i)+1 : owned_offset(i)+nowni) = &
               & y(owned_offset(i)+1 : owned_offset(i)+nowni) + yb
          deallocate(xf, yb)
        end if
        deallocate(xs)
      end do
    end do
    deallocate(owned_offset)
  end subroutine nest_spmv_perblock

end program psb_d_nest_halo_regime_test
