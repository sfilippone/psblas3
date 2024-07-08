module psb_d_oacc_csr_mat_mod

    use iso_c_binding
    use psb_d_mat_mod
    use psb_d_oacc_vect_mod
    use oaccsparse_mod

    integer(psb_ipk_), parameter, private :: is_host = -1
    integer(psb_ipk_), parameter, private :: is_sync = 0 
    integer(psb_ipk_), parameter, private :: is_dev  = 1 

    type, extends(psb_d_csr_sparse_mat) :: psb_d_oacc_csr_sparse_mat
        integer(psb_ipk_) :: devstate = is_host
    contains
        procedure, pass(a) :: all         => d_oacc_csr_all
        procedure, pass(a) :: is_host     => d_oacc_csr_is_host
        procedure, pass(a) :: is_sync     => d_oacc_csr_is_sync
        procedure, pass(a) :: is_dev      => d_oacc_csr_is_dev
        procedure, pass(a) :: set_host    => d_oacc_csr_set_host
        procedure, pass(a) :: set_sync    => d_oacc_csr_set_sync
        procedure, pass(a) :: set_dev     => d_oacc_csr_set_dev
        procedure, pass(a) :: sync_space  => d_oacc_csr_sync_space
        procedure, pass(a) :: sync        => d_oacc_csr_sync
        procedure, pass(a) :: vect_mv     => psb_d_oacc_csr_vect_mv
    end type psb_d_oacc_csr_sparse_mat

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

contains

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
        
