subroutine psb_dscg_vect(a, prec, b, x, s, eps, desc_a, info, itmax, iter, err, itrace, istop, cond)
  use psb_base_mod
  use psb_prec_mod
  use psb_d_linsolve_conv_mod
  use psb_linsolve_mod

  implicit none
  type(psb_dspmat_type), intent(in)     :: a
  class(psb_dprec_type), intent(inout)  :: prec
  type(psb_d_vect_type), intent(inout)  :: b, x
  integer(psb_ipk_), intent(in)         :: s
  real(psb_dpk_), intent(in)            :: eps
  type(psb_desc_type), intent(in)       :: desc_a
  integer(psb_ipk_), intent(out)        :: info
  integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, istop
  integer(psb_ipk_), optional, intent(out)  :: iter
  real(psb_dpk_), optional, intent(out)     :: err, cond


  integer(psb_ipk_) :: err_act
  character(len=20)           :: name = 'psb_dscg'
  character(len=*), parameter :: methdname = 'sCG'

  info = psb_success_
  call psb_erractionsave(err_act)
  info = psb_err_from_subroutine_
  call psb_errpush(info, name // ' (not implemented - empty placeholder)')
  goto 9999


9999 call psb_error_handler(err_act)
  return
end subroutine psb_dscg_vect