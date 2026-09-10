subroutine psb_z_powermethod(a, prec, lambda, desc, info, x, flag, itmax, iter, tol)
    use psb_base_mod
    use psb_prec_mod
    implicit none
    type(psb_zspmat_type), intent(in)     :: a
    class(psb_zprec_type), intent(inout)  :: prec 
    complex(psb_dpk_), intent(out)           :: lambda
    type(psb_desc_type), intent(in)       :: desc
    integer(psb_ipk_), intent(out)        :: info
    type(psb_z_vect_type), intent(inout), optional  :: x
    logical, intent(in), optional                   :: flag
    integer(psb_ipk_), intent(in), optional         :: itmax
    integer(psb_ipk_), intent(out), optional        :: iter
    real(psb_dpk_), intent(in), optional            :: tol ! def 10^-3

    type(psb_z_vect_type)   :: z, q 
    integer(psb_ipk_)       :: i, itmax_
    real(psb_dpk_)          :: tol_
    complex(psb_dpk_)          :: lambda_old, norm_factor
    
    if(present(itmax)) then
        itmax_ = itmax
    else
        itmax_ = 20_psb_ipk_
    end if

    if(present(tol)) then
        tol_ = tol
    else
        tol_ = real(1.0e-3, psb_dpk_)
    end if

    call psb_geall(z, desc, info)
    call psb_geall(q, desc, info)
    call psb_geasb(z, desc, info)
    call psb_geasb(q, desc, info)
    call q%set(zone)

    if(present(x) .and. present(flag)) then     !TO DO: can we avoid the allocation of one vector in this case?
        if(flag) call psb_geaxpby(zone, x, zzero, q, desc, info)
    end if

    lambda_old = zzero
    do i = 1, itmax_
        norm_factor = zone / psb_genrm2(q, desc, info)
        call psb_geaxpby(norm_factor, q, zzero, z, desc, info)      ! z_k = q_k / |q_k|_2
        call psb_spmm(zone, a, z, zzero, q, desc, info)             ! q_k = A z_k 
        call prec%apply(q, desc, info)                              ! q_k = B q_k
        lambda = psb_gedot(z, q, desc, info)                        ! lambda = <z_k, q_k> = z_k^T BA z_k
        if(abs(lambda - lambda_old) < tol_ * abs(lambda)) exit
        lambda_old = lambda
    end do

    if(present(iter)) iter = i
    if(present(x)) call psb_geaxpby(zone, z, zzero, x, desc, info)

    call psb_gefree(z, desc, info)
    call psb_gefree(q, desc, info)    
end subroutine