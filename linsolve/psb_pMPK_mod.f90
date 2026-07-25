module psb_pMPK_mod
  use psb_base_mod
  use psb_prec_mod

  interface psb_pMPK
    subroutine psb_s_pMPK_packd(spmat, prec, vec_in, mvec_out, s, desc, info, & 
                                    base_type, alpha, beta, gamma, mvec_temp, farr_temp)
      import psb_sspmat_type, psb_sprec_type, psb_s_vect_type, psb_s_multivect_type, &
            & psb_ipk_, psb_spk_, psb_desc_type 
      implicit none
      type(psb_sspmat_type), intent(in)           :: spmat
      class(psb_sprec_type), intent(inout)        :: prec
      type(psb_s_vect_type), intent(inout)        :: vec_in
      type(psb_s_multivect_type), intent(inout)   :: mvec_out
      integer(psb_ipk_), intent(in)               :: s
      type(psb_desc_type), intent(in)             :: desc
      integer(psb_ipk_), intent(out)              :: info
      character, optional, intent(in)                             :: base_type
      real(psb_spk_), optional, intent(in)                        :: alpha, beta, gamma
      type(psb_s_multivect_type), optional, target, intent(inout) :: mvec_temp
      real(psb_spk_), optional, target, intent(inout)             :: farr_temp(:)
    end subroutine psb_s_pMPK_packd

    subroutine psb_s_pMPK_split(spmat, prec, vec_in, Z, Q, s, desc, info, & 
                                    base_type, alpha, beta, gamma, mvec_temp, farr_temp)
      import psb_sspmat_type, psb_sprec_type, psb_s_vect_type, psb_s_multivect_type, &
            & psb_ipk_, psb_spk_, psb_desc_type 
      implicit none
      type(psb_sspmat_type), intent(in)           :: spmat
      class(psb_sprec_type), intent(inout)        :: prec
      type(psb_s_vect_type), intent(inout)        :: vec_in
      type(psb_s_multivect_type), intent(inout)   :: Z, Q
      integer(psb_ipk_), intent(in)               :: s
      type(psb_desc_type), intent(in)             :: desc
      integer(psb_ipk_), intent(out)              :: info
      character, optional, intent(in)                             :: base_type
      real(psb_spk_), optional, intent(in)                        :: alpha, beta, gamma
      type(psb_s_multivect_type), optional, target, intent(inout) :: mvec_temp
      real(psb_spk_), optional, target, intent(inout)             :: farr_temp(:)
    end subroutine psb_s_pMPK_split

    subroutine psb_c_pMPK_packd(spmat, prec, vec_in, mvec_out, s, desc, info, & 
                                    base_type, alpha, beta, gamma, mvec_temp, farr_temp)
      import psb_cspmat_type, psb_cprec_type, psb_c_vect_type, psb_c_multivect_type, &
            & psb_ipk_, psb_spk_, psb_desc_type 
      implicit none
      type(psb_cspmat_type), intent(in)           :: spmat
      class(psb_cprec_type), intent(inout)        :: prec
      type(psb_c_vect_type), intent(inout)        :: vec_in
      type(psb_c_multivect_type), intent(inout)   :: mvec_out
      integer(psb_ipk_), intent(in)               :: s
      type(psb_desc_type), intent(in)             :: desc
      integer(psb_ipk_), intent(out)              :: info
      character, optional, intent(in)                             :: base_type
      real(psb_spk_), optional, intent(in)                        :: alpha, beta, gamma
      type(psb_c_multivect_type), optional, target, intent(inout) :: mvec_temp
      real(psb_spk_), optional, target, intent(inout)             :: farr_temp(:)
    end subroutine psb_c_pMPK_packd

    subroutine psb_c_pMPK_split(spmat, prec, vec_in, Z, Q, s, desc, info, & 
                                    base_type, alpha, beta, gamma, mvec_temp, farr_temp)
      import psb_cspmat_type, psb_cprec_type, psb_c_vect_type, psb_c_multivect_type, &
            & psb_ipk_, psb_spk_, psb_desc_type 
      implicit none
      type(psb_cspmat_type), intent(in)           :: spmat
      class(psb_cprec_type), intent(inout)        :: prec
      type(psb_c_vect_type), intent(inout)        :: vec_in
      type(psb_c_multivect_type), intent(inout)   :: Z, Q
      integer(psb_ipk_), intent(in)               :: s
      type(psb_desc_type), intent(in)             :: desc
      integer(psb_ipk_), intent(out)              :: info
      character, optional, intent(in)                             :: base_type
      real(psb_spk_), optional, intent(in)                        :: alpha, beta, gamma
      type(psb_c_multivect_type), optional, target, intent(inout) :: mvec_temp
      real(psb_spk_), optional, target, intent(inout)             :: farr_temp(:)
    end subroutine psb_d_pMPK_split

    subroutine psb_d_pMPK_packd(spmat, prec, vec_in, mvec_out, s, desc, info, & 
                                    base_type, alpha, beta, gamma, mvec_temp, farr_temp)
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
      character, optional, intent(in)                             :: base_type
      real(psb_dpk_), optional, intent(in)                        :: alpha, beta, gamma
      type(psb_d_multivect_type), optional, target, intent(inout) :: mvec_temp
      real(psb_dpk_), optional, target, intent(inout)             :: farr_temp(:)
    end subroutine psb_d_pMPK_packd

    subroutine psb_d_pMPK_split(spmat, prec, vec_in, Z, Q, s, desc, info, & 
                                    base_type, alpha, beta, gamma, mvec_temp, farr_temp)
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
      character, optional, intent(in)                             :: base_type
      real(psb_dpk_), optional, intent(in)                        :: alpha, beta, gamma
      type(psb_d_multivect_type), optional, target, intent(inout) :: mvec_temp
      real(psb_dpk_), optional, target, intent(inout)             :: farr_temp(:)
    end subroutine psb_d_pMPK_split

    subroutine psb_z_pMPK_packd(spmat, prec, vec_in, mvec_out, s, desc, info, & 
                                    base_type, alpha, beta, gamma, mvec_temp, farr_temp)
      import psb_zspmat_type, psb_dprec_type, psb_z_vect_type, psb_z_multivect_type, &
            & psb_ipk_, psb_dpk_, psb_desc_type 
      implicit none
      type(psb_zspmat_type), intent(in)           :: spmat
      class(psb_zprec_type), intent(inout)        :: prec
      type(psb_z_vect_type), intent(inout)        :: vec_in
      type(psb_z_multivect_type), intent(inout)   :: mvec_out
      integer(psb_ipk_), intent(in)               :: s
      type(psb_desc_type), intent(in)             :: desc
      integer(psb_ipk_), intent(out)              :: info
      character, optional, intent(in)                             :: base_type
      real(psb_dpk_), optional, intent(in)                        :: alpha, beta, gamma
      type(psb_z_multivect_type), optional, target, intent(inout) :: mvec_temp
      real(psb_dpk_), optional, target, intent(inout)             :: farr_temp(:)
    end subroutine psb_z_pMPK_packd

    subroutine psb_z_pMPK_split(spmat, prec, vec_in, Z, Q, s, desc, info, & 
                                    base_type, alpha, beta, gamma, mvec_temp, farr_temp)
      import psb_zspmat_type, psb_zprec_type, psb_z_vect_type, psb_z_multivect_type, &
            & psb_ipk_, psb_dpk_, psb_desc_type 
      implicit none
      type(psb_zspmat_type), intent(in)           :: spmat
      class(psb_zprec_type), intent(inout)        :: prec
      type(psb_z_vect_type), intent(inout)        :: vec_in
      type(psb_z_multivect_type), intent(inout)   :: Z, Q
      integer(psb_ipk_), intent(in)               :: s
      type(psb_desc_type), intent(in)             :: desc
      integer(psb_ipk_), intent(out)              :: info
      character, optional, intent(in)                             :: base_type
      real(psb_dpk_), optional, intent(in)                        :: alpha, beta, gamma
      type(psb_z_multivect_type), optional, target, intent(inout) :: mvec_temp
      real(psb_dpk_), optional, target, intent(inout)             :: farr_temp(:)
    end subroutine psb_z_pMPK_split
  end interface
end module