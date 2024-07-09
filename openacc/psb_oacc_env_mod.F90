module psb_oacc_env_mod
contains

    subroutine psb_oacc_init(ctxt, dev)
        use psb_penv_mod
        use psb_const_mod
        use psb_error_mod
        type(psb_ctxt_type), intent(in) :: ctxt
        integer, intent(in), optional :: dev

    end subroutine psb_oacc_init

    subroutine psb_oacc_exit()
        integer :: res

    end subroutine psb_oacc_exit

end module psb_oacc_env_mod
