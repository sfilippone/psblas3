
subroutine s_oacc_mlt_v(x, y, info)
  use psb_s_oacc_vect_mod, psb_protect_name => s_oacc_mlt_v

  implicit none 
  class(psb_s_base_vect_type), intent(inout) :: x
  class(psb_s_vect_oacc), intent(inout)       :: y
  integer(psb_ipk_), intent(out)             :: info

  integer(psb_ipk_) :: i, n

  info = 0    
  n = min(x%get_nrows(), y%get_nrows())
  select type(xx => x)
  class is (psb_s_vect_oacc)
    if (y%is_host()) call y%sync()
    if (xx%is_host()) call xx%sync()
    !$acc parallel loop
    do i = 1, n
      y%v(i) = y%v(i) * xx%v(i)
    end do
    call y%set_dev()
  class default
    if (xx%is_dev()) call xx%sync()
    if (y%is_dev()) call y%sync()
    do i = 1, n
      y%v(i) = y%v(i) * xx%v(i)
    end do
    call y%set_host()
  end select
end subroutine s_oacc_mlt_v
