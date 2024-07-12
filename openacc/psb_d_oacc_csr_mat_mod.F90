module psb_d_oacc_csr_mat_mod

    use iso_c_binding
    use psb_d_mat_mod
    use psb_d_oacc_vect_mod
    !use oaccsparse_mod

    integer(psb_ipk_), parameter, private :: is_host = -1
    integer(psb_ipk_), parameter, private :: is_sync = 0 
    integer(psb_ipk_), parameter, private :: is_dev  = 1 

    type, extends(psb_d_csr_sparse_mat) :: psb_d_oacc_csr_sparse_mat
        integer(psb_ipk_) :: devstate = is_host
    contains
        procedure, nopass  :: get_fmt       => d_oacc_csr_get_fmt
        procedure, pass(a) :: sizeof        => d_oacc_csr_sizeof
        procedure, pass(a) :: vect_mv       => psb_d_oacc_csr_vect_mv
        procedure, pass(a) :: in_vect_sv    => psb_d_oacc_csr_inner_vect_sv
        procedure, pass(a) :: csmm          => psb_d_oacc_csr_csmm
        procedure, pass(a) :: csmv          => psb_d_oacc_csr_csmv
        procedure, pass(a) :: scals         => psb_d_oacc_csr_scals
        procedure, pass(a) :: scalv         => psb_d_oacc_csr_scal
        procedure, pass(a) :: reallocate_nz => psb_d_oacc_csr_reallocate_nz
        procedure, pass(a) :: allocate_mnnz => psb_d_oacc_csr_allocate_mnnz
        procedure, pass(a) :: cp_from_coo   => psb_d_oacc_csr_cp_from_coo
        procedure, pass(a) :: cp_from_fmt   => psb_d_oacc_csr_cp_from_fmt
        procedure, pass(a) :: mv_from_coo   => psb_d_oacc_csr_mv_from_coo
        procedure, pass(a) :: mv_from_fmt   => psb_d_oacc_csr_mv_from_fmt
        procedure, pass(a) :: free          => d_oacc_csr_free
        procedure, pass(a) :: mold          => psb_d_oacc_csr_mold
        procedure, pass(a) :: all           => d_oacc_csr_all
        procedure, pass(a) :: is_host       => d_oacc_csr_is_host
        procedure, pass(a) :: is_sync       => d_oacc_csr_is_sync
        procedure, pass(a) :: is_dev        => d_oacc_csr_is_dev
        procedure, pass(a) :: set_host      => d_oacc_csr_set_host
        procedure, pass(a) :: set_sync      => d_oacc_csr_set_sync
        procedure, pass(a) :: set_dev       => d_oacc_csr_set_dev
        procedure, pass(a) :: sync_space    => d_oacc_csr_sync_space
        procedure, pass(a) :: sync          => d_oacc_csr_sync
    end type psb_d_oacc_csr_sparse_mat

    interface 
        subroutine psb_d_oacc_csr_mold(a,b,info)
            import :: psb_d_oacc_csr_sparse_mat, psb_d_base_sparse_mat, psb_ipk_
            class(psb_d_oacc_csr_sparse_mat), intent(in) :: a
            class(psb_d_base_sparse_mat), intent(inout), allocatable :: b
            integer(psb_ipk_), intent(out) :: info
        end subroutine psb_d_oacc_csr_mold
    end interface

    interface
        subroutine psb_d_oacc_csr_cp_from_fmt(a,b,info)
            import :: psb_d_oacc_csr_sparse_mat, psb_d_base_sparse_mat, psb_ipk_
            class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
            class(psb_d_base_sparse_mat), intent(in) :: b
            integer(psb_ipk_), intent(out) :: info
        end subroutine psb_d_oacc_csr_cp_from_fmt
    end interface

    interface 
        subroutine psb_d_oacc_csr_mv_from_coo(a,b,info)
            import :: psb_d_oacc_csr_sparse_mat, psb_d_coo_sparse_mat, psb_ipk_
            class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
            class(psb_d_coo_sparse_mat), intent(inout) :: b
            integer(psb_ipk_), intent(out) :: info
        end subroutine psb_d_oacc_csr_mv_from_coo
    end interface

    interface
        subroutine psb_d_oacc_csr_mv_from_fmt(a,b,info)
            import :: psb_d_oacc_csr_sparse_mat, psb_d_base_sparse_mat, psb_ipk_
            class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
            class(psb_d_base_sparse_mat), intent(inout) :: b
            integer(psb_ipk_), intent(out) :: info
        end subroutine psb_d_oacc_csr_mv_from_fmt
    end interface

    interface 
        subroutine psb_d_oacc_csr_vect_mv(alpha, a, x, beta, y, info, trans)
            import :: psb_d_oacc_csr_sparse_mat, psb_dpk_, psb_d_base_vect_type, psb_ipk_
            class(psb_d_oacc_csr_sparse_mat), intent(in) :: a
            real(psb_dpk_), intent(in) :: alpha, beta
            class(psb_d_base_vect_type), intent(inout) :: x, y
            integer(psb_ipk_), intent(out) :: info
            character, optional, intent(in) :: trans
        end subroutine psb_d_oacc_csr_vect_mv
    end interface

    interface
        subroutine psb_d_oacc_csr_inner_vect_sv(alpha, a, x, beta, y, info, trans)
            import :: psb_d_oacc_csr_sparse_mat, psb_dpk_, psb_d_base_vect_type, psb_ipk_
            class(psb_d_oacc_csr_sparse_mat), intent(in) :: a
            real(psb_dpk_), intent(in) :: alpha, beta
            class(psb_d_base_vect_type), intent(inout) :: x,y
            integer(psb_ipk_), intent(out) :: info
            character, optional, intent(in) :: trans
        end subroutine psb_d_oacc_csr_inner_vect_sv
    end interface

    interface
        subroutine psb_d_oacc_csr_csmm(alpha, a, x, beta, y, info, trans)
            import :: psb_d_oacc_csr_sparse_mat, psb_dpk_, psb_d_base_vect_type, psb_ipk_
            class(psb_d_oacc_csr_sparse_mat), intent(in) :: a
            real(psb_dpk_), intent(in) :: alpha, beta, x(:,:)
            real(psb_dpk_), intent(inout) :: y(:,:)
            integer(psb_ipk_), intent(out) :: info
            character, optional, intent(in) :: trans
        end subroutine psb_d_oacc_csr_csmm
    end interface

    interface
        subroutine psb_d_oacc_csr_csmv(alpha, a, x, beta, y, info, trans)
            import :: psb_d_oacc_csr_sparse_mat, psb_dpk_, psb_d_base_vect_type, psb_ipk_
            class(psb_d_oacc_csr_sparse_mat), intent(in) :: a
            real(psb_dpk_), intent(in) :: alpha, beta, x(:)
            real(psb_dpk_), intent(inout) :: y(:)
            integer(psb_ipk_), intent(out) :: info
            character, optional, intent(in) :: trans
        end subroutine psb_d_oacc_csr_csmv
    end interface
            
    interface
        subroutine psb_d_oacc_csr_scals(d, a, info)
            import :: psb_d_oacc_csr_sparse_mat, psb_dpk_, psb_ipk_
            class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
            real(psb_dpk_), intent(in) :: d
            integer(psb_ipk_), intent(out) :: info
        end subroutine psb_d_oacc_csr_scals
    end interface

    interface 
        subroutine psb_d_oacc_csr_scal(d,a,info,side)
            import :: psb_d_oacc_csr_sparse_mat, psb_dpk_, psb_ipk_
            class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
            real(psb_dpk_), intent(in) :: d(:)
            integer(psb_ipk_), intent(out) :: info
            character, optional, intent(in) :: side
        end subroutine psb_d_oacc_csr_scal
    end interface

    interface
        subroutine psb_d_oacc_csr_reallocate_nz(nz,a)
            import :: psb_d_oacc_csr_sparse_mat, psb_ipk_
            class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
            integer(psb_ipk_), intent(in) :: nz
        end subroutine psb_d_oacc_csr_reallocate_nz
    end interface

    interface
        subroutine psb_d_oacc_csr_allocate_mnnz(m,n,a,nz)
            import :: psb_d_oacc_csr_sparse_mat, psb_ipk_
            class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
            integer(psb_ipk_), intent(in) :: m,n
            integer(psb_ipk_), intent(in), optional :: nz
        end subroutine psb_d_oacc_csr_allocate_mnnz
    end interface

    interface 
        subroutine psb_d_oacc_csr_cp_from_coo(a,b,info)
            import :: psb_d_oacc_csr_sparse_mat, psb_d_coo_sparse_mat, psb_ipk_
            class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
            class(psb_d_coo_sparse_mat), intent(in) :: b
            integer(psb_ipk_), intent(out) :: info
        end subroutine psb_d_oacc_csr_cp_from_coo
    end interface

contains
    

    subroutine d_oacc_csr_free(a)
        use psb_base_mod
        implicit none 
        class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
        integer(psb_ipk_) :: info
    
        if (allocated(a%val)) then
        !$acc exit data delete(a%val)
        end if
        if (allocated(a%ja)) then
        !$acc exit data delete(a%ja)
        end if
        if (allocated(a%irp)) then
        !$acc exit data delete(a%irp)
        end if
    
        call a%psb_d_csr_sparse_mat%free()
        
        return
    end subroutine d_oacc_csr_free
    
    

    function d_oacc_csr_sizeof(a) result(res)
        implicit none
        class(psb_d_oacc_csr_sparse_mat), intent(in) :: a
        integer(psb_epk_) :: res

        if (a%is_dev()) call a%sync()
        
        res = 8
        res = res + psb_sizeof_dp * size(a%val)
        res = res + psb_sizeof_ip * size(a%irp)
        res = res + psb_sizeof_ip * size(a%ja)

    end function d_oacc_csr_sizeof


    function d_oacc_csr_get_fmt() result(res)
        implicit none
        character(len=5) :: res
        res = 'CSR_oacc'
    end function d_oacc_csr_get_fmt

    subroutine d_oacc_csr_all(m, n, nz, a, info)
        implicit none 
        integer(psb_ipk_), intent(in)      :: m, n, nz
        class(psb_d_oacc_csr_sparse_mat), intent(out) :: a
        integer(psb_ipk_), intent(out)     :: info

        info = 0 
        if (allocated(a%val)) then
        !$acc exit data delete(a%val) finalize
        deallocate(a%val, stat=info)
        end if
        if (allocated(a%ja)) then
        !$acc exit data delete(a%ja) finalize
        deallocate(a%ja, stat=info)
        end if
        if (allocated(a%irp)) then
        !$acc exit data delete(a%irp) finalize
        deallocate(a%irp, stat=info)
        end if

        call a%set_nrows(m)
        call a%set_ncols(n)
        
        allocate(a%val(nz),stat=info)
        allocate(a%ja(nz),stat=info)
        allocate(a%irp(m+1),stat=info)
        if (info == 0) call a%set_host()
        if (info == 0) call a%sync_space()
    end subroutine d_oacc_csr_all

    function d_oacc_csr_is_host(a) result(res)
        implicit none
        class(psb_d_oacc_csr_sparse_mat), intent(in) :: a
        logical :: res

        res = (a%devstate == is_host)
    end function d_oacc_csr_is_host

    function d_oacc_csr_is_sync(a) result(res)
        implicit none
        class(psb_d_oacc_csr_sparse_mat), intent(in) :: a
        logical :: res

        res = (a%devstate == is_sync)
    end function d_oacc_csr_is_sync

    function d_oacc_csr_is_dev(a) result(res)
        implicit none
        class(psb_d_oacc_csr_sparse_mat), intent(in) :: a
        logical :: res

        res = (a%devstate == is_dev)
    end function d_oacc_csr_is_dev

    subroutine d_oacc_csr_set_host(a)
        implicit none
        class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a

        a%devstate = is_host
    end subroutine d_oacc_csr_set_host

    subroutine d_oacc_csr_set_sync(a)
        implicit none
        class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a

        a%devstate = is_sync
    end subroutine d_oacc_csr_set_sync

    subroutine d_oacc_csr_set_dev(a)
        implicit none
        class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a

        a%devstate = is_dev
    end subroutine d_oacc_csr_set_dev

    subroutine d_oacc_csr_sync_space(a)
        implicit none
        class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
        if (allocated(a%val)) then
            call d_oacc_create_dev(a%val)
        end if
        if (allocated(a%ja)) then
            call i_oacc_create_dev(a%ja)
        end if
        if (allocated(a%irp)) then
            call i_oacc_create_dev(a%irp)
        end if
    contains
        subroutine d_oacc_create_dev(v)
            implicit none
            real(psb_dpk_), intent(in) :: v(:)
            !$acc enter data copyin(v)          
        end subroutine d_oacc_create_dev
        subroutine i_oacc_create_dev(v)
            implicit none
            integer(psb_ipk_), intent(in) :: v(:)
            !$acc enter data copyin(v)          
        end subroutine i_oacc_create_dev
    end subroutine d_oacc_csr_sync_space

    subroutine d_oacc_csr_sync(a)
        implicit none
        class(psb_d_oacc_csr_sparse_mat), target, intent(in) :: a
        class(psb_d_oacc_csr_sparse_mat), pointer :: tmpa
        integer(psb_ipk_) :: info
    
        tmpa => a
        if (a%is_dev()) then
            call d_oacc_csr_to_host(a%val)
            call i_oacc_csr_to_host(a%ja)
            call i_oacc_csr_to_host(a%irp)
        else if (a%is_host()) then
            call d_oacc_csr_to_dev(a%val)
            call i_oacc_csr_to_dev(a%ja)
            call i_oacc_csr_to_dev(a%irp)
        end if
        call tmpa%set_sync()
    end subroutine d_oacc_csr_sync
    
    subroutine d_oacc_csr_to_dev(v)
        implicit none
        real(psb_dpk_), intent(in) :: v(:)
        !$acc update device(v)          
    end subroutine d_oacc_csr_to_dev
    
    subroutine d_oacc_csr_to_host(v)
        implicit none
        real(psb_dpk_), intent(in) :: v(:)
        !$acc update self(v)          
    end subroutine d_oacc_csr_to_host
    
    subroutine i_oacc_csr_to_dev(v)
        implicit none
        integer(psb_ipk_), intent(in) :: v(:)
        !$acc update device(v)          
    end subroutine i_oacc_csr_to_dev
    
    subroutine i_oacc_csr_to_host(v)
        implicit none
        integer(psb_ipk_), intent(in) :: v(:)
        !$acc update self(v)          
    end subroutine i_oacc_csr_to_host
        
    
end module psb_d_oacc_csr_mat_mod
        
