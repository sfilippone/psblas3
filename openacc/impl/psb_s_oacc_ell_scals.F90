submodule (psb_s_oacc_ell_mat_mod) psb_s_oacc_ell_scals_impl
  use psb_base_mod
contains
  module subroutine psb_s_oacc_ell_scals(d, a, info)
    implicit none 
    class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a
    real(psb_spk_), intent(in)      :: d
    integer(psb_ipk_), intent(out)  :: info

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='scal'
    logical, parameter :: debug=.false.
    integer(psb_ipk_) :: i, j, nzt, m

    info  = psb_success_
    call psb_erractionsave(err_act)

    if (a%is_host()) call a%sync()

    m = a%get_nrows()
    nzt = a%nzt

    !$acc parallel loop collapse(2) present(a)
    do i = 1, m
      do j = 1, nzt
        a%val(i, j) = a%val(i, j) * d
      end do
    end do

    call a%set_dev()

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_s_oacc_ell_scals
end submodule psb_s_oacc_ell_scals_impl
