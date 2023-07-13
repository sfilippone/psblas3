subroutine dspmm(a,b,c,info, impl_choice)
    ! TODO :
    ! * Split the C function into two subroutines:
    !    - Memory estimation
    !    - Actual computation
    ! * Call estimation function
    ! * Allocate c matrix arrays
    ! * pass the arrays as c pointers to the computation routine
    use psb_d_mat_mod
    use iso_c_binding
    implicit none
    type(psb_d_csr_sparse_mat), intent(in), target :: a,b
    type(psb_d_csr_sparse_mat), intent(inout), target :: c
    integer(psb_ipk_), intent(out)      :: info
    integer(psb_ipk_), intent(in)       :: impl_choice

    ! Internal variables
    integer(c_size_t)   :: a_m,a_n,a_nz
    type(c_ptr)         :: a_as,a_ja,a_irp
    integer(c_size_t)   :: a_max_row_nz
    integer(c_size_t)   :: b_m,b_n,b_nz
    type(c_ptr)         :: b_as,b_ja,b_irp
    integer(c_size_t)   :: b_max_row_nz
    integer(c_int)      :: impl_choice_
    type(c_ptr)         :: accumul, rows_sizes, tmp_matrix
    integer(c_size_t)   :: nnz
    type(c_ptr)         :: c_as,c_ja,c_irp

    interface spmm_build_spacc
        subroutine psb_f_spmm_build_spacc(c_a_m,c_a_n,c_a_nz,&
                            c_a_as,c_a_ja,c_a_irp,&
                            c_a_max_row_nz,&
                            c_b_m,c_b_n,c_b_nz,&
                            c_b_as,c_b_ja,c_b_irp,&
                            c_b_max_row_nz,&
                            c_impl_choice,&
                            c_accumul,&
                            c_rows_sizes,&
                            c_tmp_matrix,&
                            c_info,&
                            c_nnz) bind(C)
            use iso_c_binding
            use psb_base_mod
            integer(c_size_t), intent(in), value :: c_a_m,c_a_n,c_a_nz
            type(c_ptr), intent(in) :: c_a_as,c_a_ja,c_a_irp
            integer(c_size_t), intent(in), value :: c_a_max_row_nz
            integer(c_size_t), intent(in), value :: c_b_m,c_b_n,c_b_nz
            type(c_ptr), intent(in) :: c_b_as,c_b_ja,c_b_irp
            integer(c_size_t), intent(in), value :: c_b_max_row_nz
            integer(c_int), intent(in), value :: c_impl_choice
            type(c_ptr), intent(out) :: c_accumul,c_rows_sizes,c_tmp_matrix
            integer(psb_ipk_), intent(out) :: c_info
            integer(c_size_t), intent(out) :: c_nnz
        end subroutine psb_f_spmm_build_spacc
    end interface spmm_build_spacc

    interface spmm_row_by_row_populate
        subroutine psb_f_spmm_merge_spacc(c_accumul,&
                                        c_rows_sizes,&
                                        c_tmp_matrix,&
                                        c_impl_choice,&
                                        as,ja,irp,&
                                        c_info) bind(C)
            use iso_c_binding
            use psb_base_mod
            type(c_ptr), intent(in) :: c_accumul,c_rows_sizes,c_tmp_matrix
            integer(c_int), intent(in), value :: c_impl_choice
            integer(psb_ipk_), intent(out) :: c_info
            type(c_ptr), intent(out) :: as,ja,irp
        end subroutine psb_f_spmm_merge_spacc
    end interface spmm_row_by_row_populate

    ! Initializing internal variables
    a_m = a%get_nrows()
    a_n = a%get_ncols()
    a_nz = a%get_nzeros()
    write(*,*) 'IRP(1:5) ',a%irp(1:5)
    a_as = c_loc(a%val)
    a_ja = c_loc(a%ja)
    a_irp = c_loc(a%irp)
    ! a_max_row_nz

    b_m = b%get_nrows()
    b_n = b%get_ncols()
    b_nz = b%get_nzeros()
    b_as = c_loc(b%val)
    b_ja = c_loc(b%ja)
    b_irp = c_loc(b%irp)

    ! call calculateSize
    call psb_f_spmm_build_spacc(a_m,a_n,a_nz,a_as,&
                                a_ja,a_irp,a_max_row_nz,&
                                b_m,b_n,b_nz,b_as,b_ja,&
                                b_irp,b_max_row_nz,&
                                impl_choice_,accumul,&
                                rows_sizes,tmp_matrix,&
                                info,nnz)

    ! allocate c%val, c%ja and c%irp
    allocate(c%val(nnz))
    allocate(c%ja(nnz))
    allocate(c%irp(a_m + 1))
    
    ! c_as = c_loc(c%val)
    ! c_ja = c_loc(c%ja)
    ! c_irp = c_loc(c%irp)

    ! c%set_nrows(a_m)
    ! c%set_ncols(b_n)

    ! call spmmRowByRowPopulate
    call psb_f_spmm_merge_spacc(accumul,&
                                rows_sizes,&
                                tmp_matrix,&
                                impl_choice_,&
                                c_as,&
                                c_ja,&
                                c_irp,&
                                info)

end subroutine dspmm