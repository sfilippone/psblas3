!
! File: psb_d_nest_petsc_stokes_test.F90
!
! Read the PETSc binary Stokes files used by PETSc ex87/ex87_stokes,
! assemble the 2x2 operator as a PSBLAS nested matrix, and solve it through
! the standard PSBLAS Krylov/preconditioner interface.
!
! Expected files in PETSC_STOKES_DIR, or in
!   ${DATAFILESPATH}/matrices/hpddm/GENEO
! when PETSC_STOKES_DIR is not set:
!   B00.dat, B10.dat, B11.dat, rhs_B.dat
!
program psb_d_nest_petsc_stokes_test
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use psb_base_mod
  use psb_prec_mod
  use psb_linsolve_mod
  use psb_d_nest_mod
  implicit none

  integer(int32), parameter :: petsc_mat_classid = 1211216_int32
  integer(int32), parameter :: petsc_vec_classid = 1211214_int32

  type petsc_matrix
    integer(psb_lpk_) :: nrow = 0_psb_lpk_
    integer(psb_lpk_) :: ncol = 0_psb_lpk_
    integer(psb_lpk_) :: nnz  = 0_psb_lpk_
    integer(psb_lpk_), allocatable :: row(:), col(:)
    real(psb_dpk_),    allocatable :: val(:)
  end type petsc_matrix

  type petsc_vector
    integer(psb_lpk_) :: n = 0_psb_lpk_
    real(psb_dpk_), allocatable :: val(:)
  end type petsc_vector

  type(psb_ctxt_type)     :: context
  type(psb_d_nest_matrix) :: nested_matrix
  type(psb_dprec_type)    :: preconditioner
  type(psb_d_vect_type)   :: rhs, x_solution, residual

  type(petsc_matrix) :: b00, b10, b11
  type(petsc_vector) :: rhs_full

  integer(psb_ipk_) :: my_rank, num_procs, info
  integer(psb_ipk_) :: max_iter, trace_level, stop_criterion, n_iter, n_insert, schur_maxit, sub_fillin
  integer(psb_lpk_) :: n_u, n_p, n_total
  integer(psb_lpk_), allocatable :: rows(:), cols(:)
  real(psb_dpk_),    allocatable :: vals(:)
  real(psb_dpk_) :: eps, check_tol, final_residual, residual_norm, rhs_norm
  real(psb_dpk_) :: t0, t_read, t_assemble, t_prec, t_solve
  character(len=256) :: input_dir, method, ptype, composition, schur_solve, block_solve, sub_solve
  logical :: empty_a11

  call psb_init(context)
  call psb_info(context, my_rank, num_procs)

  call get_stokes_dir(input_dir)
  call get_string_env('PETSC_STOKES_METHOD', method, 'BICGSTAB')
  call get_string_env('PETSC_STOKES_PREC', ptype, 'NEST')
  call get_string_env('PETSC_STOKES_COMPOSITION', composition, 'SCHUR_FULL')
  call get_string_env('PETSC_STOKES_SCHUR_SOLVE', schur_solve, 'MATRIX_FREE')
  call get_string_env('PETSC_STOKES_BLOCK_SOLVE', block_solve, 'BJAC')
  call get_string_env('PETSC_STOKES_SUB_SOLVE', sub_solve, 'ILU')
  call get_int_env('PETSC_STOKES_ITMAX', max_iter, 2000)
  call get_int_env('PETSC_STOKES_TRACE', trace_level, 0)
  call get_int_env('PETSC_STOKES_SCHUR_MAXIT', schur_maxit,20)
  call get_int_env('PETSC_STOKES_SUB_FILLIN', sub_fillin, 1)
  call get_real_env('PETSC_STOKES_EPS', eps, 1.0e-4_psb_dpk_)
  call get_real_env('PETSC_STOKES_CHECK_TOL', check_tol, 1.0e-4_psb_dpk_)
  stop_criterion = 2

  if (my_rank == 0) then
    write(*,'(a)') 'PETSc Stokes nested binary test'
    write(*,'(a,a)') '  input dir  : ', trim(input_dir)
    write(*,'(a,a)') '  method     : ', trim(method)
    write(*,'(a,a)') '  prec       : ', trim(ptype)
    write(*,'(a,a)') '  composition: ', trim(composition)
    write(*,'(a,a)') '  schur solve: ', trim(schur_solve)
    write(*,'(a,a)') '  block solve: ', trim(block_solve)
    write(*,'(a,a)') '  sub solve  : ', trim(sub_solve)
    write(*,'(a,i0)') '  schur maxit: ', schur_maxit
    write(*,'(a,i0)') '  sub fillin : ', sub_fillin
    write(*,'(a,es12.4)') '  tolerance  : ', eps
    write(*,'(a,es12.4)') '  check tol  : ', check_tol
  end if

  t0 = psb_wtime()
  call read_petsc_matrix(path_join(input_dir, 'B00.dat'), b00)
  call read_petsc_matrix(path_join(input_dir, 'B10.dat'), b10)
  call read_petsc_matrix(path_join(input_dir, 'B11.dat'), b11)
  call read_petsc_vector(path_join(input_dir, 'rhs_B.dat'), rhs_full)
  t_read = psb_wtime() - t0

  n_u = b00%nrow
  n_p = b10%nrow
  n_total = n_u + n_p
  empty_a11 = (b11%nnz == 0_psb_lpk_)

  if ((b00%ncol /= n_u) .or. (b10%ncol /= n_u) .or. &
      (b11%nrow /= n_p) .or. (b11%ncol /= n_p) .or. &
      (rhs_full%n /= n_total)) then
    if (my_rank == 0) write(*,*) 'FAIL: inconsistent PETSc Stokes dimensions'
    call psb_abort(context)
  end if

  t0 = psb_wtime()
  call nested_matrix%init(context, [n_u, n_p], info)
  call check_info(info, 'nested_matrix%init')

  call select_owned_entries(b00, nested_matrix%get_owned_rows(1), rows, cols, vals, n_insert)
  call nested_matrix%ins(1, 1, n_insert, rows, cols, vals, info)
  call check_info(info, 'insert B00')
  call clear_triplets(rows, cols, vals)

  call select_owned_entries_transpose(b10, nested_matrix%get_owned_rows(1), rows, cols, vals, n_insert)
  call nested_matrix%ins(1, 2, n_insert, rows, cols, vals, info)
  call check_info(info, 'insert B10 transpose')
  call clear_triplets(rows, cols, vals)

  call select_owned_entries(b10, nested_matrix%get_owned_rows(2), rows, cols, vals, n_insert)
  call nested_matrix%ins(2, 1, n_insert, rows, cols, vals, info)
  call check_info(info, 'insert B10')
  call clear_triplets(rows, cols, vals)

  if (.not. empty_a11) then
    call select_owned_entries(b11, nested_matrix%get_owned_rows(2), rows, cols, vals, n_insert)
    call nested_matrix%ins(2, 2, n_insert, rows, cols, vals, info)
    call check_info(info, 'insert B11')
    call clear_triplets(rows, cols, vals)
  end if

  call nested_matrix%asb(info)
  call check_info(info, 'nested_matrix%asb')
  t_assemble = psb_wtime() - t0

  call psb_geall(rhs, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_geall(rhs)')
  call insert_owned_rhs(rhs_full, n_u, nested_matrix%get_owned_rows(1), nested_matrix%get_owned_rows(2), &
       & rhs, nested_matrix%desc_glob, info)
  call check_info(info, 'insert rhs')
  call psb_geasb(rhs, nested_matrix%desc_glob, info, dupl=psb_dupl_add_)
  call check_info(info, 'psb_geasb(rhs)')

  call psb_geall(x_solution, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_geall(x)')
  call x_solution%zero()
  call psb_geasb(x_solution, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_geasb(x)')

  t0 = psb_wtime()
  call preconditioner%init(context, trim(ptype), info)
  call check_info(info, 'prec%init')
  if (psb_toupper(trim(ptype)) == 'NEST') then
    call preconditioner%set('COMPOSITION', trim(composition), info)
    call preconditioner%set('SCHUR_SOLVE', trim(schur_solve), info)
    call preconditioner%set('SCHUR_MAXIT', schur_maxit, info)
    call preconditioner%set('BLOCK_SOLVE', trim(block_solve), info)
    call preconditioner%set('SUB_SOLVE', trim(sub_solve), info)
    call preconditioner%set('SUB_FILLIN', sub_fillin, info)
  end if
  call preconditioner%build(nested_matrix%a_glob, nested_matrix%desc_glob, info)
  call check_info(info, 'prec%build')
  t_prec = psb_wtime() - t0

  t0 = psb_wtime()
  call psb_krylov(trim(method), nested_matrix%a_glob, preconditioner, rhs, x_solution, eps, &
       & nested_matrix%desc_glob, info, itmax=max_iter, iter=n_iter, err=final_residual, &
       & itrace=trace_level, istop=stop_criterion)
  call check_info(info, 'psb_krylov')
  t_solve = psb_wtime() - t0

  call psb_geall(residual, nested_matrix%desc_glob, info)
  call psb_geasb(residual, nested_matrix%desc_glob, info)
  call psb_geaxpby(done, rhs, dzero, residual, nested_matrix%desc_glob, info)
  call psb_spmm(-done, nested_matrix%a_glob, x_solution, done, residual, nested_matrix%desc_glob, info)
  call check_info(info, 'residual spmm')
  residual_norm = psb_genrm2(residual, nested_matrix%desc_glob, info)
  rhs_norm      = psb_genrm2(rhs,      nested_matrix%desc_glob, info)

  if (my_rank == 0) then
    write(*,'(a,i0,a,i0,a,i0,a,i0)') '  block sizes: n_u=', n_u, ' n_p=', n_p, ' total=', n_total, ' np=', num_procs
    write(*,'(a,l1)') '  empty B11  : ', empty_a11
    write(*,'(a,i8,a,es12.4)') '  iterations=', n_iter, '  solver err=', final_residual
    write(*,'(a,es12.4)') '  ||b-Ax||_2 / ||b||_2=', residual_norm / max(rhs_norm, tiny(done))
    if (residual_norm / max(rhs_norm, tiny(done)) > check_tol) then
      write(*,'(a)') '[FAIL] psb_d_nest_petsc_stokes_test: residual above PETSC_STOKES_CHECK_TOL'
      call psb_abort(context)
    end if
    write(*,'(a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
         & '  times read=', t_read, ' assemble=', t_assemble, ' prec=', t_prec, ' solve=', t_solve
    write(*,'(a)') '[PASS] psb_d_nest_petsc_stokes_test'
  end if

9999 continue
  call preconditioner%free(info)
  call psb_gefree(residual, nested_matrix%desc_glob, info)
  call psb_gefree(x_solution, nested_matrix%desc_glob, info)
  call psb_gefree(rhs, nested_matrix%desc_glob, info)
  call nested_matrix%free(info)
  call psb_exit(context)

contains

  subroutine check_info(info_value, label)
    integer(psb_ipk_), intent(in) :: info_value
    character(len=*), intent(in) :: label
    if (info_value /= psb_success_) then
      if (my_rank == 0) write(*,*) 'FAIL: ', trim(label), ' info=', info_value
      call psb_abort(context)
    end if
  end subroutine check_info

  function path_join(dir, name) result(path)
    character(len=*), intent(in) :: dir, name
    character(len=512) :: path
    integer :: n
    n = len_trim(dir)
    if (n > 0 .and. dir(n:n) == '/') then
      path = trim(dir) // trim(name)
    else
      path = trim(dir) // '/' // trim(name)
    end if
  end function path_join

  subroutine get_stokes_dir(value)
    character(len=*), intent(out) :: value
    character(len=256) :: datafiles
    integer :: status

    call get_environment_variable('PETSC_STOKES_DIR', value, status=status)
    if (status == 0 .and. len_trim(value) > 0) return

    call get_environment_variable('DATAFILESPATH', datafiles, status=status)
    if (status == 0 .and. len_trim(datafiles) > 0) then
      value = trim(datafiles) // '/matrices/hpddm/GENEO'
    else
      value = '.'
    end if
  end subroutine get_stokes_dir

  subroutine get_string_env(name, value, default_value)
    character(len=*), intent(in) :: name, default_value
    character(len=*), intent(out) :: value
    integer :: status
    call get_environment_variable(name, value, status=status)
    if (status /= 0 .or. len_trim(value) == 0) value = default_value
  end subroutine get_string_env

  subroutine get_int_env(name, value, default_value)
    character(len=*), intent(in) :: name
    integer(psb_ipk_), intent(out) :: value
    integer(psb_ipk_), intent(in) :: default_value
    character(len=64) :: text
    integer :: status
    call get_environment_variable(name, text, status=status)
    if (status == 0 .and. len_trim(text) > 0) then
      read(text,*) value
    else
      value = default_value
    end if
  end subroutine get_int_env

  subroutine get_real_env(name, value, default_value)
    character(len=*), intent(in) :: name
    real(psb_dpk_), intent(out) :: value
    real(psb_dpk_), intent(in) :: default_value
    character(len=64) :: text
    integer :: status
    call get_environment_variable(name, text, status=status)
    if (status == 0 .and. len_trim(text) > 0) then
      read(text,*) value
    else
      value = default_value
    end if
  end subroutine get_real_env

  subroutine read_petsc_matrix(filename, mat)
    character(len=*), intent(in) :: filename
    type(petsc_matrix), intent(out) :: mat
    integer :: unit, ios
    integer(int32) :: classid, nrow32, ncol32, nnz32
    integer(int32), allocatable :: row_nnz(:), cols0(:)
    real(real64), allocatable :: vals64(:)
    integer(psb_lpk_) :: i, j, k

    open(newunit=unit, file=trim(filename), status='old', action='read', access='stream', &
         & form='unformatted', convert='big_endian', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Could not open PETSc matrix file: ', trim(filename)
      stop 2
    end if
    read(unit) classid
    if (classid /= petsc_mat_classid) stop 'Bad PETSc matrix classid'
    read(unit) nrow32
    read(unit) ncol32
    read(unit) nnz32
    mat%nrow = int(nrow32, psb_lpk_)
    mat%ncol = int(ncol32, psb_lpk_)
    mat%nnz  = int(nnz32,  psb_lpk_)

    allocate(row_nnz(mat%nrow))
    if (mat%nrow > 0_psb_lpk_) read(unit) row_nnz
    allocate(cols0(mat%nnz), vals64(mat%nnz))
    if (mat%nnz > 0_psb_lpk_) then
      read(unit) cols0
      read(unit) vals64
    end if
    close(unit)

    allocate(mat%row(mat%nnz), mat%col(mat%nnz), mat%val(mat%nnz))
    k = 0_psb_lpk_
    do i = 1_psb_lpk_, mat%nrow
      do j = 1_psb_lpk_, int(row_nnz(i), psb_lpk_)
        k = k + 1_psb_lpk_
        mat%row(k) = i
        mat%col(k) = int(cols0(k), psb_lpk_) + 1_psb_lpk_
        mat%val(k) = real(vals64(k), psb_dpk_)
      end do
    end do
    if (k /= mat%nnz) stop 'Bad PETSc matrix row lengths'
  end subroutine read_petsc_matrix

  subroutine read_petsc_vector(filename, vec)
    character(len=*), intent(in) :: filename
    type(petsc_vector), intent(out) :: vec
    integer :: unit, ios
    integer(int32) :: classid, n32
    real(real64), allocatable :: vals64(:)

    open(newunit=unit, file=trim(filename), status='old', action='read', access='stream', &
         & form='unformatted', convert='big_endian', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Could not open PETSc vector file: ', trim(filename)
      stop 2
    end if
    read(unit) classid
    if (classid /= petsc_vec_classid) stop 'Bad PETSc vector classid'
    read(unit) n32
    vec%n = int(n32, psb_lpk_)
    allocate(vals64(vec%n), vec%val(vec%n))
    if (vec%n > 0_psb_lpk_) read(unit) vals64
    close(unit)
    vec%val = real(vals64, psb_dpk_)
  end subroutine read_petsc_vector

  subroutine select_owned_entries(mat, owned_rows, out_rows, out_cols, out_vals, count)
    type(petsc_matrix), intent(in) :: mat
    integer(psb_lpk_), intent(in) :: owned_rows(:)
    integer(psb_lpk_), allocatable, intent(out) :: out_rows(:), out_cols(:)
    real(psb_dpk_), allocatable, intent(out) :: out_vals(:)
    integer(psb_ipk_), intent(out) :: count
    integer(psb_lpk_) :: i, first_owned, last_owned
    integer(psb_ipk_) :: k

    count = 0
    if (size(owned_rows) == 0) then
      allocate(out_rows(0), out_cols(0), out_vals(0))
      return
    end if

    first_owned = minval(owned_rows)
    last_owned  = maxval(owned_rows)
    do i = 1_psb_lpk_, mat%nnz
      if (mat%row(i) >= first_owned .and. mat%row(i) <= last_owned) count = count + 1
    end do

    allocate(out_rows(count), out_cols(count), out_vals(count))
    k = 0
    do i = 1_psb_lpk_, mat%nnz
      if (mat%row(i) >= first_owned .and. mat%row(i) <= last_owned) then
        k = k + 1
        out_rows(k) = mat%row(i)
        out_cols(k) = mat%col(i)
        out_vals(k) = mat%val(i)
      end if
    end do
  end subroutine select_owned_entries

  subroutine select_owned_entries_transpose(mat, owned_rows, out_rows, out_cols, out_vals, count)
    type(petsc_matrix), intent(in) :: mat
    integer(psb_lpk_), intent(in) :: owned_rows(:)
    integer(psb_lpk_), allocatable, intent(out) :: out_rows(:), out_cols(:)
    real(psb_dpk_), allocatable, intent(out) :: out_vals(:)
    integer(psb_ipk_), intent(out) :: count
    integer(psb_lpk_) :: i, first_owned, last_owned
    integer(psb_ipk_) :: k

    count = 0
    if (size(owned_rows) == 0) then
      allocate(out_rows(0), out_cols(0), out_vals(0))
      return
    end if

    first_owned = minval(owned_rows)
    last_owned  = maxval(owned_rows)
    do i = 1_psb_lpk_, mat%nnz
      if (mat%col(i) >= first_owned .and. mat%col(i) <= last_owned) count = count + 1
    end do

    allocate(out_rows(count), out_cols(count), out_vals(count))
    k = 0
    do i = 1_psb_lpk_, mat%nnz
      if (mat%col(i) >= first_owned .and. mat%col(i) <= last_owned) then
        k = k + 1
        out_rows(k) = mat%col(i)
        out_cols(k) = mat%row(i)
        out_vals(k) = mat%val(i)
      end if
    end do
  end subroutine select_owned_entries_transpose

  subroutine insert_owned_rhs(vec, velocity_size, owned_u, owned_p, rhs_vec, desc, info)
    type(petsc_vector), intent(in) :: vec
    integer(psb_lpk_), intent(in) :: velocity_size
    integer(psb_lpk_), intent(in) :: owned_u(:), owned_p(:)
    type(psb_d_vect_type), intent(inout) :: rhs_vec
    type(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(out) :: info
    integer(psb_lpk_), allocatable :: idx(:)
    real(psb_dpk_), allocatable :: values(:)
    integer :: i, n

    n = size(owned_u) + size(owned_p)
    allocate(idx(n), values(n))
    do i = 1, size(owned_u)
      idx(i) = owned_u(i)
      values(i) = vec%val(owned_u(i))
    end do
    do i = 1, size(owned_p)
      idx(size(owned_u) + i) = velocity_size + owned_p(i)
      values(size(owned_u) + i) = vec%val(velocity_size + owned_p(i))
    end do
    call psb_geins(n, idx, values, rhs_vec, desc, info)
    deallocate(idx, values)
  end subroutine insert_owned_rhs

  subroutine clear_triplets(in_rows, in_cols, in_vals)
    integer(psb_lpk_), allocatable, intent(inout) :: in_rows(:), in_cols(:)
    real(psb_dpk_), allocatable, intent(inout) :: in_vals(:)
    if (allocated(in_rows)) deallocate(in_rows)
    if (allocated(in_cols)) deallocate(in_cols)
    if (allocated(in_vals)) deallocate(in_vals)
  end subroutine clear_triplets

end program psb_d_nest_petsc_stokes_test