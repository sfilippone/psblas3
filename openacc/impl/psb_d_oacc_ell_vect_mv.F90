submodule (psb_d_oacc_ell_mat_mod) psb_d_oacc_ell_vect_mv_impl
  use psb_base_mod
contains
  module subroutine psb_d_oacc_ell_vect_mv(alpha, a, x, beta, y, info, trans)
    implicit none

    real(psb_dpk_), intent(in)        :: alpha, beta
    class(psb_d_oacc_ell_sparse_mat), intent(in) :: a
    class(psb_d_base_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(out)    :: info
    character, optional, intent(in)   :: trans

    integer(psb_ipk_) :: m, n, nzt, nc

    info = psb_success_
    m = a%get_nrows()
    n = a%get_ncols()
    nzt = a%nzt
    nc = size(a%ja,2)
    if ((n /= size(x%v)) .or. (m /= size(y%v))) then
      write(0,*) 'Size error ', m, n, size(x%v), size(y%v)
      info = psb_err_invalid_mat_state_
      return
    end if

    if (a%is_host()) call a%sync()
    if (x%is_host()) call x%sync()
    if (y%is_host()) call y%sync()

    call inner_spmv(m, n, nc, alpha, a%val, a%ja, x%v, beta, y%v, info)

    call y%set_dev()

  contains

    subroutine inner_spmv(m, n, nc, alpha, val, ja, x, beta, y, info)
      implicit none
      integer(psb_ipk_) :: m, n, nc
      real(psb_dpk_), intent(in) :: alpha, beta
      real(psb_dpk_) :: val(:,:), x(:), y(:)
      integer(psb_ipk_) :: ja(:,:)
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_) :: i, j, ii, isz
      real(psb_dpk_) :: tmp
      integer(psb_ipk_), parameter :: vsz = 256

      info = 0

      !$acc parallel loop vector_length(vsz) private(isz)
      do ii = 1, m, vsz
        isz = min(vsz, m - ii + 1)
        !$acc loop independent private(tmp)
        do i = ii, ii + isz - 1
          tmp = 0.0_psb_dpk_
          !$acc loop seq
          do j = 1, nc 
            if (ja(i,j) > 0) then
              tmp = tmp + val(i,j) * x(ja(i,j))
            end if
          end do
          y(i) = alpha * tmp + beta * y(i)
        end do
      end do
    end subroutine inner_spmv

  end subroutine psb_d_oacc_ell_vect_mv
end submodule psb_d_oacc_ell_vect_mv_impl
