subroutine psb_d_powermethod(a, prec, lambda, desc, info, x, flag, itmax, iter, tol)
    use psb_base_mod, only : psb_ipk_, psb_dpk_, psb_desc_type, &
                            & psb_dspmat_type, psb_d_vect_type
    use psb_prec_mod, only : psb_dprec_type
    type(psb_dspmat_type), intent(in)     :: a
    class(psb_dprec_type), intent(inout)  :: prec 
    real(psb_dpk_), intent(out)           :: lambda
    type(psb_desc_type), intent(in)       :: desc
    integer(psb_ipk_), intent(out)        :: info
    type(psb_d_vect_type), intent(inout), optional  :: x
    logical, intent(in), optional                   :: flag
    integer(psb_ipk_), intent(in), optional         :: itmax
    integer(psb_ipk_), intent(out), optional        :: iter
    real(psb_dpk_), intent(in), optional            :: tol ! def 10^-3

    ! call psb_geasb(tz,desc_a,info,mold=vmold,scratch=.true.)
    ! call psb_geasb(tt,desc_a,info,mold=vmold,scratch=.true.)
    ! call psb_geasb(wv(1),desc_a,info,mold=vmold,scratch=.true.)
    ! call psb_geasb(wv(2),desc_a,info,mold=vmold,scratch=.true.)
    ! call psb_geall(tq,desc_a,info)
    ! call tq%set(done)
    ! call psb_geasb(tq,desc_a,info,mold=vmold)
    ! call psb_spmm(done,a,tq,dzero,tt,desc_a,info) !
    ! call sm%sv%apply_v(done,tt,dzero,tz,desc_a,'NoTrans',work,wv,info) ! z_{k+1} = BA q_k
    ! do i=1,sm%rho_estimate_iterations
    !     znrm = psb_genrm2(tz,desc_a,info)               ! znrm = |z_k|_2
    !     call psb_geaxpby((done/znrm),tz,dzero,tq,desc_a,info)  ! q_k = z_k/znrm
    !     call psb_spmm(done,a,tq,dzero,tt,desc_a,info) ! t_{k+1} = BA q_k
    !     call sm%sv%apply_v(done,tt,dzero,tz,desc_a,'NoTrans',work,wv,info) ! z_{k+1} = B t_{k+1}
    !     lambda = psb_gedot(tq,tz,desc_a,info)      ! lambda = q_k^T z_{k+1} = q_k^T BA q_k
    !     !write(0,*) 'BLD: lambda estimate ',i,lambda
    ! end do
    ! sm%rho_ba = lambda
end subroutine