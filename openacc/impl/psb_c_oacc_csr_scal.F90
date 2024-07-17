submodule (psb_c_oacc_csr_mat_mod) psb_c_oacc_csr_scal_impl
  use psb_base_mod
contains
  module subroutine psb_c_oacc_csr_scal(d, a, info, side)
    implicit none 
    class(psb_c_oacc_csr_sparse_mat), intent(inout) :: a
    complex(psb_spk_), intent(in)      :: d(:)
    integer(psb_ipk_), intent(out)  :: info
    character, intent(in), optional :: side

    integer(psb_ipk_)  :: err_act
    character(len=20)  :: name='scal'
    logical, parameter :: debug=.false.
    integer(psb_ipk_)  :: i, j

    info  = psb_success_
    call psb_erractionsave(err_act)

    if (a%is_host()) call a%sync()

    if (present(side)) then
      if (side == 'L') then
        !$acc parallel loop present(a, d)
        do i = 1, a%get_nrows()
          do j = a%irp(i), a%irp(i+1) - 1
            a%val(j) = a%val(j) * d(i)
          end do
        end do
      else if (side == 'R') then
        !$acc parallel loop present(a, d)
        do i = 1, a%get_ncols()
          do j = a%irp(i), a%irp(i+1) - 1
            a%val(j) = a%val(j) * d(a%ja(j))
          end do
        end do
      end if
    else
      !$acc parallel loop present(a, d)
      do i = 1, size(a%val)
        a%val(i) = a%val(i) * d(i)
      end do
    end if

    call a%set_dev()

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_c_oacc_csr_scal
end submodule psb_c_oacc_csr_scal_impl
