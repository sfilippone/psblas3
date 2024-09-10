submodule (psb_c_oacc_hll_mat_mod) psb_c_oacc_hll_vect_mv_impl
  use psb_base_mod
contains
  module subroutine psb_c_oacc_hll_vect_mv(alpha, a, x, beta, y, info, trans)
    implicit none

    complex(psb_spk_), intent(in)        :: alpha, beta
    class(psb_c_oacc_hll_sparse_mat), intent(in) :: a
    class(psb_c_base_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(out)    :: info
    character, optional, intent(in)   :: trans

    integer(psb_ipk_) :: m, n, nhacks, hksz
    character :: trans_
    logical :: device_done, tra

    info = psb_success_
    m = a%get_nrows()
    n = a%get_ncols()
    nhacks = size(a%hkoffs) - 1
    hksz = a%hksz

    if ((n > size(x%v)) .or. (m > size(y%v))) then
      write(0,*) 'Size error ', m, n, size(x%v), size(y%v)
      info = psb_err_invalid_mat_state_
      return
    end if
    device_done = .false.
    if (present(trans)) then
      trans_ = trans
    else
      trans_ = 'N'
    end if
    tra = (psb_toupper(trans_) == 'T') .or. (psb_toupper(trans_) == 'C')

    if (.not.tra) then 
      select type(xx  => x)
      class is (psb_c_vect_oacc)
        select type (yy  => y)
        class is (psb_c_vect_oacc)
          if (a%is_host()) call a%sync()
          if (xx%is_host()) call xx%sync()
          if (yy%is_host()) call yy%sync()
          call inner_spmv(m, nhacks, hksz, alpha, a%val, a%ja, a%hkoffs, x%v, beta, y%v, info)
          call y%set_dev()
          device_done = .true.
        end select
      end select
    end if
    
    if (.not.device_done) then
      if (x%is_dev()) call x%sync()
      if (y%is_dev()) call y%sync()
      call a%psb_c_hll_sparse_mat%spmm(alpha, x%v, beta, y%v, info, trans)
      call y%set_host()
    end if
  contains

    subroutine inner_spmv(m, nhacks, hksz, alpha, val, ja, hkoffs, x, beta, y, info)
      implicit none
      integer(psb_ipk_) :: m, nhacks, hksz
      complex(psb_spk_), intent(in) :: alpha, beta
      complex(psb_spk_) :: val(:), x(:), y(:)
      integer(psb_ipk_) :: ja(:), hkoffs(:)
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_) :: i, j, idx, k, ipnt,ir,nr,nlc,isz,ii
      complex(psb_spk_) :: tmp

      info = 0
      !$acc parallel loop  private(nlc, isz,ir,nr)
      do i = 1, nhacks
        isz = hkoffs(i + 1) - hkoffs(i) 
        nlc = isz/hksz
        ir  = (i-1)*hksz
        nr  = min(hksz,m-ir)
        !$acc loop independent private(tmp,ii,ipnt)
        do ii = 1, nr
          ipnt = hkoffs(i) + ii 
          tmp = czero
          !$acc loop seq
          do j = 1, nlc
            tmp = tmp + val(ipnt) * x(ja(ipnt))
            ipnt = ipnt + hksz
          end do
          y(ii+ir) = alpha * tmp + beta * y(ii+ir)
        end do
      end do
    end subroutine inner_spmv
  end subroutine psb_c_oacc_hll_vect_mv
end submodule psb_c_oacc_hll_vect_mv_impl
