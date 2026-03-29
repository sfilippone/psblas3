module psb_eigsolve_mod
  use psb_base_mod
  use psb_prec_mod
  interface psb_powermethod
    subroutine psb_d_powermethod(a, prec, lambda, desc, info, x, flag, itmax, iter)
      import :: psb_ipk_, psb_dpk_, psb_desc_type, &
            psb_dspmat_type, psb_d_vect_type, psb_dprec_type
      type(psb_dspmat_type), intent(in)     :: a
      class(psb_dprec_type), intent(inout)  :: prec 
      real(psb_dpk_), intent(out)           :: lambda
      type(psb_desc_type), intent(in)       :: desc
      integer(psb_ipk_), intent(out)        :: info
      type(psb_d_vect_type), intent(inout), optional  :: x
      logical, intent(in), optional                   :: flag
      integer(psb_ipk_), intent(in), optional         :: itmax
      integer(psb_ipk_), intent(out), optional        :: iter
    end subroutine

    !TO DO: add s, c, z versions

  end interface
end module psb_eigsolve_mod
 