submodule (psb_c_oacc_ell_mat_mod) psb_c_oacc_ell_csmv_impl
  use psb_base_mod
contains
  module subroutine psb_c_oacc_ell_csmv(alpha, a, x, beta, y, info, trans)
    implicit none 
    class(psb_c_oacc_ell_sparse_mat), intent(in) :: a
    complex(psb_spk_), intent(in)          :: alpha, beta
    complex(psb_spk_), intent(in)          :: x(:)
    complex(psb_spk_), intent(inout)       :: y(:)
    integer(psb_ipk_), intent(out)       :: info
    character, optional, intent(in)      :: trans

    character :: trans_
    integer(psb_ipk_) :: i, j, m, n, nzt
    logical :: tra
    integer(psb_ipk_) :: err_act
    character(len=20) :: name = 'c_oacc_ell_csmv'
    logical, parameter :: debug = .false.

    call psb_erractionsave(err_act)
    info = psb_success_

    if (present(trans)) then
      trans_ = trans
    else
      trans_ = 'N'
    end if

    if (.not.a%is_asb()) then 
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, name)
      goto 9999
    endif

    tra = (psb_toupper(trans_) == 'T') .or. (psb_toupper(trans_) == 'C')

    if (tra) then 
      m = a%get_ncols()
      n = a%get_nrows()
    else
      n = a%get_ncols()
      m = a%get_nrows()
    end if

    if (size(x,1) < n) then 
      info = 36
      call psb_errpush(info, name, i_err = (/3 * ione, n, izero, izero, izero/))
      goto 9999
    end if

    if (size(y,1) < m) then 
      info = 36
      call psb_errpush(info, name, i_err = (/5 * ione, m, izero, izero, izero/))
      goto 9999
    end if

    if (tra) then 
      call a%psb_c_ell_sparse_mat%spmm(alpha, x, beta, y, info, trans) 
    else
      nzt = a%nzt

      !$acc parallel loop present(a, x, y)
      do i = 1, m
        y(i) = beta * y(i)
      end do

      !$acc parallel loop present(a, x, y)
      do i = 1, m
        do j = 1, nzt
          y(i) = y(i) + alpha * a%val(i, j) * x(a%ja(i, j))
        end do
      end do
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_c_oacc_ell_csmv
end submodule psb_c_oacc_ell_csmv_impl
