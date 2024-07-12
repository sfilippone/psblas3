subroutine psb_d_oacc_csr_csmm(alpha, a, x, beta, y, info, trans)
    use psb_base_mod
    !use elldev_mod
    !use psb_vectordev_mod
    use psb_d_oacc_csr_mat_mod, psb_protect_name => psb_d_oacc_csr_csmm
    implicit none
    class(psb_d_oacc_csr_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)          :: alpha, beta
    real(psb_dpk_), intent(in)          :: x(:,:)
    real(psb_dpk_), intent(inout)       :: y(:,:)
    integer(psb_ipk_), intent(out)       :: info
    character, optional, intent(in)      :: trans
  
    character :: trans_
    integer(psb_ipk_) :: i, j, m, n,k, nxy
    logical :: tra
    integer(psb_ipk_) :: err_act
    character(len=20) :: name = 'd_oacc_csmm'
    logical, parameter :: debug = .false.
  
    info = psb_success_
    call psb_erractionsave(err_act)
  
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
      call a%psb_d_csr_sparse_mat%spmm(alpha, x, beta, y, info, trans)
    else
      nxy = min(size(x,2), size(y,2))
      
      !$acc parallel loop collapse(2) present(a, x, y)
      do j = 1, nxy
        do i = 1, m
          y(i,j) = beta * y(i,j)
        end do
      end do
      
      !$acc parallel loop collapse(2) present(a, x, y)
      do j = 1, nxy
        do i = 1, n
          do k = a%irp(i), a%irp(i+1) - 1
            y(a%ja(k), j) = y(a%ja(k), j) + alpha * a%val(k) * x(i, j)
          end do
        end do
      end do
    endif
  
    call psb_erractionrestore(err_act)
    return
  
  9999 call psb_error_handler(err_act)
    return
  
  end subroutine psb_d_oacc_csr_csmm
  
