subroutine psb_c_powermethod(a, prec, lambda, desc, info, x, flag, itmax, iter, tol)
    use psb_base_mod
    use psb_prec_mod
    implicit none
    type(psb_cspmat_type), intent(in)     :: a
    class(psb_cprec_type), intent(inout)  :: prec 
    complex(psb_spk_), intent(out)           :: lambda
    type(psb_desc_type), intent(in)       :: desc
    integer(psb_ipk_), intent(out)        :: info
    type(psb_c_vect_type), intent(inout), optional  :: x
    logical, intent(in), optional                   :: flag
    integer(psb_ipk_), intent(in), optional         :: itmax
    integer(psb_ipk_), intent(out), optional        :: iter
    complex(psb_spk_), intent(in), optional            :: tol ! def 10^-3

    type(psb_c_vect_type)   :: z, q 
    integer(psb_ipk_)       :: i, itmax_
    complex(psb_spk_)          :: tol_
    complex(psb_spk_)          :: lambda_old, norm_factor
    
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
    call q%set(cone)

    if(present(x) .and. present(flag)) then     !TO DO: can we avoid the allocation of one vector in this case?
        if(flag) call psb_geaxpby(cone, x, czero, q, desc, info)
    end if

    lambda_old = czero
    do i = 1, itmax_
        norm_factor = cone / psb_genrm2(q, desc, info)
        call psb_geaxpby(norm_factor, q, czero, z, desc, info)      ! z_k = q_k / |q_k|_2
        call psb_spmm(cone, a, z, czero, q, desc, info)             ! q_k = A z_k 
        call prec%apply(q, desc, info)                              ! q_k = B q_k
        lambda = psb_gedot(z, q, desc, info)                        ! lambda = <z_k, q_k> = z_k^T BA z_k
        if(abs(lambda - lambda_old) < tol_ * abs(lambda)) exit
        lambda_old = lambda
    end do

    if(present(iter)) iter = i
    if(present(x)) call psb_geaxpby(cone, z, czero, x, desc, info)

    call psb_gefree(z, desc, info)
    call psb_gefree(q, desc, info)    
end subroutine