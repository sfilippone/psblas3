module sp3mm_mod
    use iso_c_binding
    ! use psb_const_mod
    ! use psb_error_mod

    interface psb_dspmm
        subroutine dspmm(a,b,c,info)
            use psb_mat_mod
            import :: psb_ipk_
            type(psb_dspmat_type), intent(in) :: a,b
            type(psb_dspmat_type), intent(out):: c
            integer(psb_ipk_), intent(out)    :: info
        end subroutine dspmm
    end interface psb_dspmm

    ! interface spmm_row_by_row
    !     subroutine dspmm_row_by_row_ub(a,b,c,info)
    !         use psb_d_mat_mod, only : psb_dspmat_type
    !         import :: psb_ipk_
    !         implicit none 
    !         type(psb_dspmat_type), intent(in)   :: a,b
    !         type(psb_dspmat_type), intent(out)  :: c
    !         integer(psb_ipk_), intent(out)      :: info
    !     end subroutine dspmm_row_by_row_ub

    !     subroutine dspmm_row_by_row_symb_num(a,b,c,info)
    !         use psb_d_mat_mod, only : psb_dspmat_type
    !         import :: psb_ipk_
    !         implicit none 
    !         type(psb_dspmat_type), intent(in)   :: a,b
    !         type(psb_dspmat_type), intent(out)  :: c
    !         integer(psb_ipk_), intent(out)      :: info
    !     end subroutine dspmm_row_by_row_symb_num

    !     subroutine dspmm_row_by_row_1d_blocks_symb_num(a,b,c,info)
    !         use psb_d_mat_mod, only : psb_dspmat_type
    !         import :: psb_ipk_
    !         implicit none 
    !         type(psb_dspmat_type), intent(in)   :: a,b
    !         type(psb_dspmat_type), intent(out)  :: c
    !         integer(psb_ipk_), intent(out)      :: info
    !     end subroutine dspmm_row_by_row_1d_blocks_symb_num

    !     subroutine dspmm_row_by_row_2d_blocks_symb_num(a,b,c,info)
    !         use psb_d_mat_mod, only : psb_dspmat_type
    !         import :: psb_ipk_
    !         implicit none 
    !         type(psb_dspmat_type), intent(in)   :: a,b
    !         type(psb_dspmat_type), intent(out)  :: c
    !         integer(psb_ipk_), intent(out)      :: info
    !     end subroutine dspmm_row_by_row_2d_blocks_symb_num
    ! end interface spmm_row_by_row

end module sp3mm_mod