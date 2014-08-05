module psb_eigen_mod
  use psb_base_mod
  use psb_krylov_mod
  use psb_prec_mod

  interface 
    Subroutine psb_d_power_vect(method,a,prec,b,x,eps,desc_a,info,&
         & itmax,iter,err,itrace,irst,istop,cond)
      
      use psb_base_mod, only  : psb_ipk_, psb_desc_type, psb_dspmat_type, &
           & psb_dpk_, psb_d_vect_type
      use psb_prec_mod, only : psb_dprec_type
      
      character(len=*)                      :: method
      Type(psb_dspmat_type), Intent(in)     :: a
      Type(psb_desc_type), Intent(in)       :: desc_a
      class(psb_dprec_type), intent(inout)  :: prec 
      type(psb_d_vect_type), Intent(inout)  :: b
      type(psb_d_vect_type), Intent(inout)  :: x
      Real(psb_dpk_), Intent(in)            :: eps
      integer(psb_ipk_), intent(out)                  :: info
      integer(psb_ipk_), Optional, Intent(in)         :: itmax, itrace, irst,istop
      integer(psb_ipk_), Optional, Intent(out)        :: iter
      Real(psb_dpk_), Optional, Intent(out) :: err,cond

    end Subroutine psb_d_power_vect
  end interface
  
end module psb_eigen_mod
