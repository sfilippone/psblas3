submodule (psb_d_oacc_csr_mat_mod) psb_d_oacc_csr_allocate_mnnz_impl
  use psb_base_mod
contains
  module subroutine psb_d_oacc_csr_allocate_mnnz(m, n, a, nz)
    implicit none 
    integer(psb_ipk_), intent(in) :: m, n
    class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
    integer(psb_ipk_), intent(in), optional :: nz
    integer(psb_ipk_) :: info
    integer(psb_ipk_) :: err_act, nz_
    character(len=20)  :: name='allocate_mnz'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = psb_success_

    call a%psb_d_csr_sparse_mat%allocate(m, n, nz)
    call a%set_host()
    call a%sync_dev_space()
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_d_oacc_csr_allocate_mnnz
end submodule psb_d_oacc_csr_allocate_mnnz_impl
