subroutine psb_s_oacc_mlt_v_2(alpha, x, y, beta, z, info, conjgx, conjgy)
  use psb_s_oacc_vect_mod, psb_protect_name  => psb_s_oacc_mlt_v_2
  use psb_string_mod
  implicit none 
  real(psb_spk_), intent(in)                 :: alpha, beta
  class(psb_s_base_vect_type), intent(inout) :: x
  class(psb_s_base_vect_type), intent(inout) :: y
  class(psb_s_vect_oacc), intent(inout)      :: z
  integer(psb_ipk_), intent(out)             :: info
  character(len=1), intent(in), optional     :: conjgx, conjgy
  integer(psb_ipk_) :: i, n
  logical :: conjgx_, conjgy_, device_done

  conjgx_ = .false.
  conjgy_ = .false.
  device_done = .false.
  if (present(conjgx)) conjgx_ = (psb_toupper(conjgx) == 'C')
  if (present(conjgy)) conjgy_ = (psb_toupper(conjgy) == 'C')

  n = min(x%get_nrows(), y%get_nrows(), z%get_nrows())
  info = 0    
  select type(xx  => x)
  class is (psb_s_vect_oacc)
    select type (yy  => y)
    class is (psb_s_vect_oacc)
      if (xx%is_host()) call xx%sync()
      if (yy%is_host()) call yy%sync()
      if ((beta /= szero) .and. (z%is_host())) call z%sync()
      call s_inner_oacc_mlt_v_2(n,alpha, xx%v, yy%v, beta, z%v, info, conjgx_, conjgy_)
      call z%set_dev()
      device_done = .true.
    end select
  end select
  if (.not.device_done) then 
    if (x%is_dev()) call x%sync()
    if (y%is_dev()) call y%sync()
    if ((beta /= szero) .and. (z%is_dev())) call z%sync()
    if (conjgx_.and.conjgy_) then 
      do i = 1, n
        z%v(i) = alpha * (x%v(i)) * (y%v(i)) + beta * z%v(i)
      end do
    else if (conjgx_.and.(.not.conjgy_)) then 
      do i = 1, n
        z%v(i) = alpha * (x%v(i)) * (y%v(i)) + beta * z%v(i)
      end do
    else if ((.not.conjgx_).and.(conjgy_)) then 
      do i = 1, n
        z%v(i) = alpha * (x%v(i)) * (y%v(i)) + beta * z%v(i)
      end do
    else
      do i = 1, n
        z%v(i) = alpha * (x%v(i)) * (y%v(i)) + beta * z%v(i)
      end do
    end if
    call z%set_host()
  end if

contains
  subroutine s_inner_oacc_mlt_v_2(n,alpha, x, y, beta, z, info, conjgx, conjgy)
    implicit none
    integer(psb_ipk_), intent(in) :: n
    real(psb_spk_), intent(in)                 :: alpha, beta
    real(psb_spk_), intent(inout) :: x(:), y(:), z(:)
    integer(psb_ipk_), intent(out)             :: info
    logical, intent(in)         :: conjgx, conjgy

    integer(psb_ipk_) :: i
    if (conjgx.and.conjgy) then 
      !$acc parallel loop present(x,y,z)
      do i = 1, n
        z(i) = alpha * (x(i)) * (y(i)) + beta * z(i)
      end do
    else if (conjgx.and.(.not.conjgy)) then 
      !$acc parallel loop present(x,y,z)
      do i = 1, n
        z(i) = alpha * (x(i)) * (y(i)) + beta * z(i)
      end do
    else if ((.not.conjgx).and.(conjgy)) then 
      !$acc parallel loop present(x,y,z)
      do i = 1, n
        z(i) = alpha * (x(i)) * (y(i)) + beta * z(i)
      end do
    else
      !$acc parallel loop present(x,y,z)
      do i = 1, n
        z(i) = alpha * (x(i)) * (y(i)) + beta * z(i)
      end do
    end if
  end subroutine s_inner_oacc_mlt_v_2
end subroutine psb_s_oacc_mlt_v_2

