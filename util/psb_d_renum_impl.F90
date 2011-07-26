subroutine psb_d_mat_renum(alg,mat,info)
  use psb_base_mod
  use psb_renum_mod, psb_protect_name => psb_d_mat_renum
  use psb_gps_mod
  integer, intent(in) :: alg
  type(psb_dspmat_type), intent(inout) :: mat
  integer, intent(out) :: info
  
  integer :: err_act
  character(len=20)           :: name, ch_err

  info = psb_success_
  name = 'mat_distf'
  call psb_erractionsave(err_act)

  info = psb_success_


  


  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return


end subroutine psb_d_mat_renum
