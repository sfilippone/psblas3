module psb_eigsolve_mod
  use psb_base_mod
  use psb_prec_mod
  
  interface psb_powermethod
    subroutine psb_s_powermethod(a, prec, lambda, desc, info, x, flag, itmax, iter, tol)
      use psb_base_mod, only : psb_ipk_, psb_desc_type, psb_sspmat_type, psb_spk_, psb_s_vect_type
      use psb_prec_mod, only : psb_sprec_type
      type(psb_sspmat_type), intent(in)     :: a
      class(psb_sprec_type), intent(inout)  :: prec 
      real(psb_spk_), intent(out)           :: lambda
      type(psb_desc_type), intent(in)       :: desc
      integer(psb_ipk_), intent(out)        :: info
      type(psb_s_vect_type), intent(inout), optional  :: x
      logical, intent(in), optional                   :: flag
      integer(psb_ipk_), intent(in), optional         :: itmax
      integer(psb_ipk_), intent(out), optional        :: iter
      real(psb_spk_), intent(in), optional            :: tol ! def 10^-3
    end subroutine psb_s_powermethod

    subroutine psb_c_powermethod(a, prec, lambda, desc, info, x, flag, itmax, iter, tol)
      use psb_base_mod, only : psb_ipk_, psb_desc_type, psb_cspmat_type, psb_spk_, psb_c_vect_type
      use psb_prec_mod, only : psb_cprec_type
      type(psb_cspmat_type), intent(in)     :: a
      class(psb_cprec_type), intent(inout)  :: prec 
      complex(psb_spk_), intent(out)        :: lambda
      type(psb_desc_type), intent(in)       :: desc
      integer(psb_ipk_), intent(out)        :: info
      type(psb_c_vect_type), intent(inout), optional  :: x
      logical, intent(in), optional                   :: flag
      integer(psb_ipk_), intent(in), optional         :: itmax
      integer(psb_ipk_), intent(out), optional        :: iter
      real(psb_spk_), intent(in), optional            :: tol ! def 10^-3
    end subroutine psb_c_powermethod

    subroutine psb_d_powermethod(a, prec, lambda, desc, info, x, flag, itmax, iter, tol)
      use psb_base_mod, only : psb_ipk_, psb_desc_type, psb_dspmat_type, psb_dpk_, psb_d_vect_type
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
    end subroutine psb_d_powermethod

    subroutine psb_z_powermethod(a, prec, lambda, desc, info, x, flag, itmax, iter, tol)
      use psb_base_mod, only : psb_ipk_, psb_desc_type, psb_zspmat_type, psb_dpk_, psb_z_vect_type
      use psb_prec_mod, only : psb_zprec_type
      type(psb_zspmat_type), intent(in)     :: a
      class(psb_zprec_type), intent(inout)  :: prec 
      complex(psb_dpk_), intent(out)        :: lambda
      type(psb_desc_type), intent(in)       :: desc
      integer(psb_ipk_), intent(out)        :: info
      type(psb_z_vect_type), intent(inout), optional  :: x
      logical, intent(in), optional                   :: flag
      integer(psb_ipk_), intent(in), optional         :: itmax
      integer(psb_ipk_), intent(out), optional        :: iter
      real(psb_dpk_), intent(in), optional            :: tol ! def 10^-3
    end subroutine psb_z_powermethod
  end interface
end module psb_eigsolve_mod
 