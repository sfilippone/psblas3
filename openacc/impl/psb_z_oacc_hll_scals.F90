submodule (psb_z_oacc_hll_mat_mod) psb_z_oacc_hll_scals_impl
  use psb_base_mod
contains
  module subroutine psb_z_oacc_hll_scals(d, a, info)
    implicit none 
    class(psb_z_oacc_hll_sparse_mat), intent(inout) :: a
    complex(psb_dpk_), intent(in)      :: d
    integer(psb_ipk_), intent(out)  :: info

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name = 'scal'
    integer(psb_ipk_) :: i, j, k, hksz, nzt, nhacks

    info  = psb_success_
    call psb_erractionsave(err_act)

    if (a%is_host()) call a%sync()

    hksz = a%hksz
    nhacks = (a%get_nrows() + hksz - 1) / hksz
    nzt = a%nzt

    ! $ a parallel loop collapse(2) present(a)
    !$acc parallel loop  present(a)
    do i = 1, nhacks
      do j = a%hkoffs(i), a%hkoffs(i + 1) - 1
        a%val(j) = a%val(j) * d
      end do
    end do

    call a%set_dev()

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_z_oacc_hll_scals
end submodule psb_z_oacc_hll_scals_impl
