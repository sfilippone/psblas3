submodule (psb_d_oacc_ell_mat_mod) psb_d_oacc_ell_reallocate_nz_impl
  use psb_base_mod
contains
  module subroutine psb_d_oacc_ell_reallocate_nz(nz, a)
    implicit none 
    integer(psb_ipk_), intent(in) :: nz
    class(psb_d_oacc_ell_sparse_mat), intent(inout) :: a
    integer(psb_ipk_) :: info
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='d_oacc_ell_reallocate_nz'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = psb_success_

    call a%psb_d_ell_sparse_mat%reallocate(nz)

    call a%set_dev()
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_d_oacc_ell_reallocate_nz
end submodule psb_d_oacc_ell_reallocate_nz_impl
