submodule (psb_z_oacc_hll_mat_mod) psb_z_oacc_hll_vect_mv_impl
  use psb_base_mod
contains
  module subroutine psb_z_oacc_hll_vect_mv(alpha, a, x, beta, y, info, trans)
    implicit none

    complex(psb_dpk_), intent(in)        :: alpha, beta
    class(psb_z_oacc_hll_sparse_mat), intent(in) :: a
    class(psb_z_base_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(out)    :: info
    character, optional, intent(in)   :: trans

    integer(psb_ipk_) :: m, n, nhacks, hksz

    info = psb_success_
    m = a%get_nrows()
    n = a%get_ncols()
    nhacks = size(a%hkoffs) - 1
    hksz = a%hksz

    if ((n /= size(x%v)) .or. (m /= size(y%v))) then
      write(0,*) 'Size error ', m, n, size(x%v), size(y%v)
      info = psb_err_invalid_mat_state_
      return
    end if

    if (a%is_host()) call a%sync()
    if (x%is_host()) call x%sync()
    if (y%is_host()) call y%sync()

    call inner_spmv(m, nhacks, hksz, alpha, a%val, a%ja, a%hkoffs, x%v, beta, y%v, info)
    call y%set_dev()

  contains

    subroutine inner_spmv(m, nhacks, hksz, alpha, val, ja, hkoffs, x, beta, y, info)
      implicit none
      integer(psb_ipk_) :: m, nhacks, hksz
      complex(psb_dpk_), intent(in) :: alpha, beta
      complex(psb_dpk_) :: val(:), x(:), y(:)
      integer(psb_ipk_) :: ja(:), hkoffs(:)
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_) :: i, j, idx, k
      complex(psb_dpk_) :: tmp

      info = 0

      !$acc parallel loop present(val, ja, hkoffs, x, y)
      do i = 1, nhacks
        do k = 0, hksz - 1
          idx = hkoffs(i) + k
          if (idx <= hkoffs(i + 1) - 1) then
            tmp = 0.0_psb_dpk_
            !$acc loop seq
            do j = hkoffs(i) + k, hkoffs(i + 1) - 1, hksz
              if (ja(j) > 0) then
                tmp = tmp + val(j) * x(ja(j))
              end if
            end do
            y(k + 1 + (i - 1) * hksz) = alpha * tmp + beta * y(k + 1 + (i - 1) * hksz)
          end if
        end do
      end do
    end subroutine inner_spmv

  end subroutine psb_z_oacc_hll_vect_mv
end submodule psb_z_oacc_hll_vect_mv_impl
