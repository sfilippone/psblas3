subroutine psb_d_oacc_csr_inner_vect_sv(alpha, a, x, beta, y, info, trans)
    use psb_base_mod
    !    use elldev_mod
    !   use psb_vectordev_mod
    use psb_d_oacc_csr_mat_mod, psb_protect_name => psb_d_oacc_csr_inner_vect_sv
    use psb_d_oacc_vect_mod
    implicit none
    class(psb_d_oacc_csr_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in) :: alpha, beta
    class(psb_d_base_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(out) :: info
    character, optional, intent(in) :: trans
  
    real(psb_dpk_), allocatable :: rx(:), ry(:)
    logical :: tra
    character :: trans_
    integer(psb_ipk_) :: err_act
    character(len=20) :: name = 'd_oacc_csr_inner_vect_sv'
    logical, parameter :: debug = .false.
    integer(psb_ipk_) :: i
  
    call psb_get_erraction(err_act)
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
  
    if (tra .or. (beta /= dzero)) then
      call x%sync()
      call y%sync()
      call a%psb_d_csr_sparse_mat%inner_spsm(alpha, x, beta, y, info, trans)
      call y%set_host()
    else
      select type (xx => x)
      type is (psb_d_vect_oacc)
        select type(yy => y)
        type is (psb_d_vect_oacc)
          if (xx%is_host()) call xx%sync()
          if (beta /= dzero) then
            if (yy%is_host()) call yy%sync()
          end if
          !$acc parallel loop present(a, xx, yy)
          do i = 1, size(a%val)
            yy%v(i) = alpha * a%val(i) * xx%v(a%ja(i)) + beta * yy%v(i)
          end do
          call yy%set_dev()
        class default
          rx = xx%get_vect()
          ry = y%get_vect()
          call a%psb_d_csr_sparse_mat%inner_spsm(alpha, rx, beta, ry, info)
          call y%bld(ry)
        end select
      class default
        rx = x%get_vect()
        ry = y%get_vect()
        call a%psb_d_csr_sparse_mat%inner_spsm(alpha, rx, beta, ry, info)
        call y%bld(ry)
      end select
    endif
  
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info, name, a_err = 'csrg_vect_sv')
      goto 9999
    endif
  
    call psb_erractionrestore(err_act)
    return
  
  9999 call psb_error_handler(err_act)
    return
  end subroutine psb_d_oacc_csr_inner_vect_sv
  
