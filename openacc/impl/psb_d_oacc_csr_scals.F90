submodule (psb_d_oacc_csr_mat_mod) psb_d_oacc_csr_scals_impl
  use psb_base_mod
contains
  module subroutine psb_d_oacc_csr_scals(d, a, info)
    implicit none 
    class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
    real(psb_dpk_), intent(in)      :: d
    integer(psb_ipk_), intent(out)  :: info

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='scal'
    logical, parameter :: debug=.false.
    integer(psb_ipk_) :: i

    info  = psb_success_
    call psb_erractionsave(err_act)

    if (a%is_host()) call a%sync()

    !$acc parallel loop present(a)
    do i = 1, size(a%val)
      a%val(i) = a%val(i) * d
    end do

    call a%set_dev()

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_d_oacc_csr_scals
end submodule psb_d_oacc_csr_scals_impl
