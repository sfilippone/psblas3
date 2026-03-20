module psb_d_pMPK_mod
    use psb_base_mod
    use psb_d_prec_mod

    interface psb_pMPK
        subroutine psb_d_pMPK_packd(spmat, prec, vec_in, mvec_out, s, desc, info, base_type, alpha, beta, gamma)
            import psb_dspmat_type, psb_dprec_type, psb_d_vect_type, psb_d_multivect_type, &
                        & psb_ipk_, psb_dpk_, psb_desc_type 
            implicit none
            type(psb_dspmat_type), intent(in)           :: spmat
            class(psb_dprec_type), intent(inout)        :: prec
            type(psb_d_vect_type), intent(inout)        :: vec_in
            type(psb_d_multivect_type), intent(inout)   :: mvec_out
            integer(psb_ipk_), intent(in)               :: s
            type(psb_desc_type), intent(in)             :: desc
            integer(psb_ipk_), intent(out)              :: info
            character, optional, intent(in)             :: base_type
            real(psb_dpk_), optional, intent(in)        :: alpha, beta, gamma
        end subroutine psb_d_pMPK_packd

        subroutine psb_d_pMPK_split(spmat, prec, vec_in, Z, Q, s, desc, info, base_type, alpha, beta, gamma)
            import psb_dspmat_type, psb_dprec_type, psb_d_vect_type, psb_d_multivect_type, &
                        & psb_ipk_, psb_dpk_, psb_desc_type 
            implicit none
            type(psb_dspmat_type), intent(in)           :: spmat
            class(psb_dprec_type), intent(inout)        :: prec
            type(psb_d_vect_type), intent(inout)        :: vec_in
            type(psb_d_multivect_type), intent(inout)   :: Z, Q
            integer(psb_ipk_), intent(in)               :: s
            type(psb_desc_type), intent(in)             :: desc
            integer(psb_ipk_), intent(out)              :: info
            character, optional, intent(in)             :: base_type
            real(psb_dpk_), optional, intent(in)        :: alpha, beta, gamma
        end subroutine psb_d_pMPK_split
    end interface

end module