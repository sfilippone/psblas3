subroutine dspmm_row_by_row_ub(a,b,c,info)
    use psb_error_mod
    use psb_base_mat_mod
    use psb_d_mat_mod, only : psb_dspmat_type
    use psb_objhandle_mod, only: spmat_t, config_t
    implicit none 
    type(psb_dspmat_type), intent(in)   :: a,b
    type(psb_dspmat_type), intent(out)  :: c
    integer(psb_ipk_), intent(out)      :: info

    ! TODO : implement the C interface
end subroutine dspmm_row_by_row_ub