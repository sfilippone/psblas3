submodule (psb_s_oacc_ell_mat_mod) psb_s_oacc_ell_allocate_mnnz_impl
  use psb_base_mod
contains
  module subroutine psb_s_oacc_ell_allocate_mnnz(m, n, a, nz)
    implicit none 
    integer(psb_ipk_), intent(in) :: m, n
    class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a
    integer(psb_ipk_), intent(in), optional :: nz
    integer(psb_ipk_) :: info
    integer(psb_ipk_) :: err_act, nz_
    character(len=20)  :: name='allocate_mnnz'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = psb_success_

    if (present(nz)) then
      nz_ = nz
    else
      nz_ = 10  
    end if

    call a%psb_s_ell_sparse_mat%allocate(m, n, nz_)

    if (.not.allocated(a%val)) then
      allocate(a%val(m, nz_))
      allocate(a%ja(m, nz_))
      allocate(a%irn(m))
      allocate(a%idiag(m))
    end if

    a%val = szero
    a%ja = -1
    a%irn = 0
    a%idiag = 0

    call a%set_dev()
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_s_oacc_ell_allocate_mnnz
end submodule psb_s_oacc_ell_allocate_mnnz_impl
