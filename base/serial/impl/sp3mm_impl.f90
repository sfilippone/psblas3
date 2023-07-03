subroutine dspmm(a,b,c,info, impl_choice)
    use psb_d_mat_mod
    implicit none
    type(psb_d_csr_sparse_mat), intent(in) :: a,b
    type(psb_d_csr_sparse_mat), intent(out):: c
    integer(psb_ipk_), intent(out)    :: info
    integer(psb_ipk_), intent(in), optional :: impl_choice

    ! Internal variables
    integer(c_size_t), value :: a_m,a_n,a_nz
    real(c_double) :: a_as
    integer(c_size_t) :: a_ja,a_irp

    integer(c_size_t), value :: a_max_row_nz
    integer(c_size_t), value :: b_m,b_n,b_nz
    real(c_double) :: b_as
    integer(c_size_t) :: b_ja,b_irp
    integer(c_size_t), value :: b_max_row_nz
    integer(c_int), value :: impl_choice_
    integer(c_size_t), value :: c_m,c_n,c_nz
    real(c_double) :: c_as
    integer(c_size_t) :: c_ja,c_irp
    integer(c_size_t), value :: c_max_row_nz

    ! Initializing internal variables
    a_m = a%get_nrows()
    a_n = a%get_ncols()
    a_nz = a%get_nzeros()
    a_as = a%val
    a_ja = a%ja
    a_irp = a%irp
    ! a_max_row_nz
    b_m = b%get_nrows()
    b_n = b%get_ncols()
    b_nz = b%get_nzeros()
    b_as = b%val
    b_ja = b%ja
    b_irp = b%irp

    if (present(impl_choice)) then
        impl_choice_ = impl_choice
    else
        impl_choice_ = 0
    end if

    ! Calling psb_f_spmm
    psb_f_spmm(a_m,a_n,a_nz,&
            & a_as,a_ja,a_irp,&
            & a_max_row_nz,&
            & b_m,b_n,b_nz,&
            & b_as,b_ja,b_irp,&
            & impl_choice_,&
            & c_m,c_n,c_nz,&
            & c_as,c_ja,c_irp,&
            & c_max_row_nz,&
            & info)

    ! Putting results in c
    c%m = c_m
    c%n = c_n
    c%val = c_as
    c%ja = c_ja
    c%irp = c_irp

    interface
        subroutine psb_f_spmm(c_a_m,c_a_n,c_a_nz,&
                            & c_a_as,c_a_ja,c_a_irp,&
                            & c_a_max_row_nz,
                            & c_b_m,c_b_n,c_b_nz,&
                            & c_b_as,c_b_ja,c_b_irp,&
                            & c_b_max_row_nz,
                            & c_impl_choice,&
                            & c_c_m,c_c_n,c_c_nz,&
                            & c_c_as,c_c_ja,c_c_irp,
                            & c_c_max_row_nz,&
                            & c_info) bind(c)
            use iso_c_binding
            use psb_base_mod
            integer(c_size_t), intent(in), value :: c_a_m,c_a_n,c_a_nz
            real(c_double), intent(in) :: c_a_as
            integer(c_size_t), intent(in) :: c_a_ja,c_a_irp
            integer(c_size_t), intent(in), value :: c_a_max_row_nz
            integer(c_size_t), intent(in), value :: c_b_m,c_b_n,c_b_nz
            real(c_double), intent(in) :: c_b_as
            integer(c_size_t), intent(in) :: c_b_ja,c_b_irp
            integer(c_size_t), intent(in), value :: c_b_max_row_nz
            integer(c_int), intent(in), value :: c_impl_choice
            integer(c_size_t), intent(out), value :: c_c_m,c_c_n,c_c_nz
            real(c_double), intent(out) :: c_c_as
            integer(c_size_t), intent(out) :: c_c_ja,c_c_irp
            integer(c_size_t), intent(out), value :: c_c_max_row_nz
            integer(psb_ipk_), intent(out) :: c_info
        end subroutine psb_f_spmm
    end interface
end subroutine dspmm