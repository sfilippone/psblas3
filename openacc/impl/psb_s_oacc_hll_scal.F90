submodule (psb_s_oacc_hll_mat_mod) psb_s_oacc_hll_scal_impl
  use psb_base_mod
contains
  module subroutine psb_s_oacc_hll_scal(d, a, info, side)
    implicit none 
    class(psb_s_oacc_hll_sparse_mat), intent(inout) :: a
    real(psb_spk_), intent(in)      :: d(:)
    integer(psb_ipk_), intent(out)  :: info
    character, intent(in), optional :: side

    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name = 'scal'
    integer(psb_ipk_)  :: i, j, k, hksz, nzt, nhacks

    info  = psb_success_
    call psb_erractionsave(err_act)

    if (a%is_host()) call a%sync()

    hksz = a%hksz
    nhacks = (a%get_nrows() + hksz - 1) / hksz
    nzt = a%nzt

    if (present(side)) then
      if (side == 'L') then
        ! $ a    parallel loop collapse(2) present(a, d)
        !$acc parallel loop  present(a, d)
        do i = 1, nhacks
          do j = a%hkoffs(i), a%hkoffs(i + 1) - 1
            k = (j - a%hkoffs(i)) / nzt + (i - 1) * hksz + 1
            a%val(j) = a%val(j) * d(k)
          end do
        end do
      else if (side == 'R') then
        ! $ a  parallel loop collapse(2) present(a, d)
        !$acc parallel loop present(a, d)
        do i = 1, nhacks
          do j = a%hkoffs(i), a%hkoffs(i + 1) - 1
            a%val(j) = a%val(j) * d(a%ja(j))
          end do
        end do
      end if
    else
      ! $ a parallel loop collapse(2) present(a, d)
      !$acc parallel loop  present(a, d)
      do i = 1, nhacks
        do j = a%hkoffs(i), a%hkoffs(i + 1) - 1
          a%val(j) = a%val(j) * d(j - a%hkoffs(i) + 1)
        end do
      end do
    end if

    call a%set_dev()

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_s_oacc_hll_scal
end submodule psb_s_oacc_hll_scal_impl
