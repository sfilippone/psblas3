submodule (psb_s_oacc_ell_mat_mod) psb_s_oacc_ell_cp_from_coo_impl
  use psb_base_mod
contains
  module subroutine psb_s_oacc_ell_cp_from_coo(a, b, info)
    implicit none

    class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a
    class(psb_s_coo_sparse_mat), intent(in)         :: b
    integer(psb_ipk_), intent(out)                  :: info

    integer(psb_ipk_) :: i, j, k, row, col, nz_per_row
    real(psb_spk_) :: value
    integer(psb_ipk_), allocatable :: row_counts(:)
    integer(psb_ipk_) :: hacksize, nza

    info = psb_success_
    hacksize = 1 

    call a%set_nrows(b%get_nrows())
    call a%set_ncols(b%get_ncols())
    nz_per_row = a%nzt

    if (.not.allocated(a%val)) then
      allocate(a%val(a%get_nrows(), nz_per_row))
      allocate(a%ja(a%get_nrows(), nz_per_row))
      allocate(a%irn(a%get_nrows()))
      allocate(a%idiag(a%get_nrows()))
    end if
    a%val = szero
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
        a%val(row, row_counts(row) + 1) = value
        a%ja(row, row_counts(row) + 1) = col
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
        if (a%ja(i, j) == i) then
          a%idiag(i) = j
          exit
        end if
      end do
    end do

    deallocate(row_counts)

    call a%set_dev()
    if (info /= 0) goto 9999

    return

9999 continue
    info = psb_err_alloc_dealloc_
    return

  end subroutine psb_s_oacc_ell_cp_from_coo
end submodule psb_s_oacc_ell_cp_from_coo_impl
