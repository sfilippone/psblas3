subroutine psb_d_oacc_csr_vect_mv(alpha, a, x, beta, y, info, trans)
    use psb_base_mod
    use psb_d_oacc_csr_mat_mod, psb_protect_name => psb_d_oacc_csr_vect_mv
    implicit none

    real(psb_dpk_), intent(in)        :: alpha, beta
    class(psb_d_oacc_csr_sparse_mat), intent(in) :: a
    class(psb_d_base_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(out)    :: info
    character, optional, intent(in)   :: trans

    integer(psb_ipk_) :: m, n

    info = psb_success_
    m = a%get_nrows()
    n = a%get_ncols()

    if ((n /= size(x%v)) .or. (n /= size(y%v))) then
        write(0,*) 'Size error ', m, n, size(x%v), size(y%v)
        info = psb_err_invalid_mat_state_
        return
    end if

    if (a%is_host()) call a%sync()
    if (x%is_host()) call x%sync()
    if (y%is_host()) call y%sync()

    call inner_spmv(m, n, alpha, a%val, a%ja, a%irp, x%v, beta, y%v, info)
    call y%set_dev()

contains

    subroutine inner_spmv(m, n, alpha, val, ja, irp, x, beta, y, info)
        implicit none
        integer(psb_ipk_) :: m, n
        real(psb_dpk_), intent(in) :: alpha, beta
        real(psb_dpk_) :: val(:), x(:), y(:)
        integer(psb_ipk_) :: ja(:), irp(:)
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
                do j = irp(i), irp(i + 1) - 1
                    tmp = tmp + val(j) * x(ja(j))
                end do
                y(i) = alpha * tmp + beta * y(i)
            end do
        end do
    end subroutine inner_spmv

end subroutine psb_d_oacc_csr_vect_mv
