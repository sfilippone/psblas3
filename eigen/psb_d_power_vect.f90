Subroutine psb_d_power_vect(method,a,prec,b,x,eps,desc_a,info,&
     & itmax,iter,err,itrace,irst,istop,cond)

  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod
  use psb_eigen_mod, psb_protect_name =>psb_d_power_vect

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

