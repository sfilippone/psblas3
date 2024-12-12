module psb_oacc_env_mod
  use psb_penv_mod
  use psb_const_mod
  use psb_error_mod
  logical, private :: oacc_do_maybe_free_buffer = .false.

contains  
  function psb_oacc_get_maybe_free_buffer() result(res)
    logical :: res
    res = oacc_do_maybe_free_buffer
  end function psb_oacc_get_maybe_free_buffer

  subroutine psb_oacc_set_maybe_free_buffer(val)
    logical, intent(in) :: val
    oacc_do_maybe_free_buffer = val
  end subroutine psb_oacc_set_maybe_free_buffer

  subroutine psb_oacc_init(ctxt, dev)
    type(psb_ctxt_type), intent(in) :: ctxt
    integer, intent(in), optional :: dev
    oacc_do_maybe_free_buffer = .false.    
  end subroutine psb_oacc_init

  subroutine psb_oacc_exit()
    integer :: res
    
  end subroutine psb_oacc_exit

end module psb_oacc_env_mod
