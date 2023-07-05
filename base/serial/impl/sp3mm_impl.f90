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
    type(psb_d_csr_sparse_mat), intent(in) :: a,b
    type(psb_d_csr_sparse_mat), intent(inout):: c
    integer(psb_ipk_), intent(out)      :: info
    integer(psb_ipk_), intent(in)       :: impl_choice

    ! Internal variables
    integer(c_size_t):: a_m,a_n,a_nz
    real(c_double), pointer :: a_as(:)
    integer(c_size_t), pointer :: a_ja(:),a_irp(:)
    type(c_ptr) :: a_as_ptr,a_ja_ptr,a_irp_ptr
    integer(c_size_t) :: a_max_row_nz
    integer(c_size_t) :: b_m,b_n,b_nz
    real(c_double), pointer :: b_as(:)
    integer(c_size_t), pointer :: b_ja(:),b_irp(:)
    type(c_ptr) :: b_as_ptr,b_ja_ptr,b_irp_ptr
    integer(c_size_t) :: b_max_row_nz
    integer(c_int) :: impl_choice_
    type(c_ptr) :: accumul, rows_sizes, tmp_matrix
    integer(c_size_t) :: nnz
    real(c_double), pointer :: c_as(:)
    integer(c_size_t), pointer :: c_ja(:),c_irp(:)
    type(c_ptr) :: c_as_ptr,c_ja_ptr,c_irp_ptr

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
                                        c_as,c_ja,c_irp,&
                                        c_info) bind(C)
            use iso_c_binding
            use psb_base_mod
            type(c_ptr), intent(in) :: c_accumul,c_rows_sizes,c_tmp_matrix
            integer(c_int), intent(in), value :: c_impl_choice
            integer(psb_ipk_), intent(out) :: c_info
            type(c_ptr), intent(out) :: c_as,c_ja,c_irp
        end subroutine psb_f_spmm_merge_spacc
    end interface spmm_row_by_row_populate

    ! Initializing internal variables
    a_m = a%get_nrows()
    a_n = a%get_ncols()
    a_nz = a%get_nzeros()
    a_as = a%val
    a_as_ptr = c_loc(a_as)
    a_ja = a%ja
    a_ja_ptr = c_loc(a_ja)
    a_irp = a%irp
    a_irp_ptr = c_loc(a_irp)
    ! ! a_max_row_nz
    b_m = b%get_nrows()
    b_n = b%get_ncols()
    b_nz = b%get_nzeros()
    b_as = b%val
    b_as_ptr = c_loc(b_as)
    b_ja = b%ja
    b_ja_ptr = c_loc(b_ja)
    b_irp = b%irp
    b_irp_ptr = c_loc(b_irp)

    ! call calculateSize
    call psb_f_spmm_build_spacc(a_m,a_n,a_nz,a_as_ptr,&
                                a_ja_ptr,a_irp_ptr,a_max_row_nz,&
                                b_m,b_n,b_nz,b_as_ptr,b_ja_ptr,&
                                b_irp_ptr,b_max_row_nz,&
                                impl_choice_,accumul,&
                                rows_sizes,tmp_matrix,&
                                info,nnz)

    ! allocate c%val, c%ja and c%irp
    allocate(c%val(nnz))
    allocate(c%ja(nnz))
    allocate(c%irp(a_m + 1))
    
    c_as = c%val
    c_as_ptr = c_loc(c_as)
    c_ja = c%ja
    c_ja_ptr = c_loc(c_ja)
    c_irp = c%irp
    c_irp_ptr = c_loc(c_irp)

    ! c%set_nrows(a_m)
    ! c%set_ncols(b_n)

    ! call spmmRowByRowPopulate
    call psb_f_spmm_merge_spacc(accumul,&
                                rows_sizes,&
                                tmp_matrix,&
                                impl_choice_,&
                                c_as_ptr,&
                                c_ja_ptr,&
                                c_irp_ptr,&
                                info)

end subroutine dspmm