submodule (psb_c_oacc_ell_mat_mod) psb_c_oacc_ell_scal_impl
  use psb_base_mod
contains
  module subroutine psb_c_oacc_ell_scal(d, a, info, side)
    implicit none 
    class(psb_c_oacc_ell_sparse_mat), intent(inout) :: a
    complex(psb_spk_), intent(in)      :: d(:)
    integer(psb_ipk_), intent(out)  :: info
    character, intent(in), optional :: side

    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='scal'
    logical, parameter :: debug=.false.
    integer(psb_ipk_)  :: i, j, m, nzt

    info  = psb_success_
    call psb_erractionsave(err_act)

    if (a%is_host()) call a%sync()

    m = a%get_nrows()
    nzt = a%nzt

    if (present(side)) then
      if (side == 'L') then
        !$acc parallel loop collapse(2) present(a, d)
        do i = 1, m
          do j = 1, nzt
            a%val(i, j) = a%val(i, j) * d(i)
          end do
        end do
      else if (side == 'R') then
        !$acc parallel loop collapse(2) present(a, d)
        do i = 1, m
          do j = 1, nzt
            a%val(i, j) = a%val(i, j) * d(a%ja(i, j))
          end do
        end do
      end if
    else
      !$acc parallel loop collapse(2) present(a, d)
      do i = 1, m
        do j = 1, nzt
          a%val(i, j) = a%val(i, j) * d(j)
        end do
      end do
    end if

    call a%set_dev()

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_c_oacc_ell_scal
end submodule psb_c_oacc_ell_scal_impl
