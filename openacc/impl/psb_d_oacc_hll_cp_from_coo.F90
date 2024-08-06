submodule (psb_d_oacc_hll_mat_mod) psb_d_oacc_hll_cp_from_coo_impl
  use psb_base_mod
contains
  module subroutine psb_d_oacc_hll_cp_from_coo(a, b, info)
    implicit none

    class(psb_d_oacc_hll_sparse_mat), intent(inout) :: a
    class(psb_d_coo_sparse_mat), intent(in)         :: b
    integer(psb_ipk_), intent(out)                  :: info

    integer(psb_ipk_) :: i, j, k, row, col, nz_per_row
    real(psb_dpk_) :: value
    integer(psb_ipk_), allocatable :: row_counts(:)
    integer(psb_ipk_) :: hacksize, nza

    info = psb_success_
    hacksize = 32  ! Assuming a default hack size of 32

    call a%set_nrows(b%get_nrows())
    call a%set_ncols(b%get_ncols())
    nz_per_row = a%nzt

    if (.not.allocated(a%val)) then
      allocate(a%val(nz_per_row * a%get_nrows()))
      allocate(a%ja(nz_per_row * a%get_nrows()))
      allocate(a%irn(a%get_nrows()))
      allocate(a%idiag(a%get_nrows()))
      allocate(a%hkoffs((a%get_nrows() + hacksize - 1) / hacksize))
    end if
    a%val = 0.0_psb_dpk_
    a%ja = -1
    a%irn = 0
    a%idiag = 0

    allocate(row_counts(a%get_nrows()))
    row_counts = 0

    nza = b%get_nzeros()

    !$acc parallel loop present(b, a, row_counts)
    do k = 1, nza
      row = b%ia(k)
      col = b%ja(k)
      value = b%val(k)
      if (row_counts(row) < nz_per_row) then
        a%val(row_counts(row) + 1 + (row - 1) * nz_per_row) = value
        a%ja(row_counts(row) + 1 + (row - 1) * nz_per_row) = col
        row_counts(row) = row_counts(row) + 1
      else
        info = psb_err_invalid_mat_state_
        !goto 9999
      end if
    end do

    a%irn = row_counts

    !$acc parallel loop present(a)
    do i = 1, a%get_nrows()
      do j = 1, nz_per_row
        if (a%ja(j + (i - 1) * nz_per_row) == i) then
          a%idiag(i) = j
          exit
        end if
      end do
    end do

    ! Calculate hkoffs for HLL format
    !$acc parallel loop present(a)
    do i = 1, size(a%hkoffs)
      a%hkoffs(i) = (i - 1) * hacksize
    end do

    deallocate(row_counts)

    call a%set_dev()
    if (info /= 0) goto 9999

    return

9999 continue
    info = psb_err_alloc_dealloc_
    return

  end subroutine psb_d_oacc_hll_cp_from_coo
end submodule psb_d_oacc_hll_cp_from_coo_impl
