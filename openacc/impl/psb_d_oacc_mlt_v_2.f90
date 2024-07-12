submodule (psb_d_oacc_vect_mod) d_oacc_mlt_v_2_impl
  use psb_string_mod
contains
  module subroutine d_oacc_mlt_v_2(alpha, x, y, beta, z, info, conjgx, conjgy)
    implicit none 
    real(psb_dpk_), intent(in)                 :: alpha, beta
    class(psb_d_base_vect_type), intent(inout) :: x
    class(psb_d_base_vect_type), intent(inout) :: y
    class(psb_d_vect_oacc), intent(inout)      :: z
    integer(psb_ipk_), intent(out)             :: info
    character(len=1), intent(in), optional     :: conjgx, conjgy
    integer(psb_ipk_) :: i, n
    logical :: conjgx_, conjgy_

    conjgx_ = .false.
    conjgy_ = .false.
    if (present(conjgx)) conjgx_ = (psb_toupper(conjgx) == 'C')
    if (present(conjgy)) conjgy_ = (psb_toupper(conjgy) == 'C')

    n = min(x%get_nrows(), y%get_nrows(), z%get_nrows())

    info = 0    
    select type(xx => x)
    class is (psb_d_vect_oacc)
      select type (yy => y)
      class is (psb_d_vect_oacc)
        if (xx%is_host()) call xx%sync()
        if (yy%is_host()) call yy%sync()
        if ((beta /= dzero) .and. (z%is_host())) call z%sync()
        !$acc parallel loop
        do i = 1, n
          z%v(i) = alpha * xx%v(i) * yy%v(i) + beta * z%v(i)
        end do
        call z%set_dev()
      class default
        if (xx%is_dev()) call xx%sync()
        if (yy%is_dev()) call yy%sync()
        if ((beta /= dzero) .and. (z%is_dev())) call z%sync()
        do i = 1, n
          z%v(i) = alpha * xx%v(i) * yy%v(i) + beta * z%v(i)
        end do
        call z%set_host()
      end select
    class default
      if (x%is_dev()) call x%sync()
      if (y%is_dev()) call y%sync()
      if ((beta /= dzero) .and. (z%is_dev())) call z%sync()
      do i = 1, n
        z%v(i) = alpha * x%v(i) * y%v(i) + beta * z%v(i)
      end do
      call z%set_host()
    end select
  end subroutine d_oacc_mlt_v_2
end submodule d_oacc_mlt_v_2_impl

