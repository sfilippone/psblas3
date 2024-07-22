module psb_d_oacc_ell_mat_mod
    use iso_c_binding
    use psb_d_mat_mod
    use psb_d_ell_mat_mod
    use psb_d_oacc_vect_mod

    integer(psb_ipk_), parameter, private :: is_host = -1
    integer(psb_ipk_), parameter, private :: is_sync = 0 
    integer(psb_ipk_), parameter, private :: is_dev  = 1 

    type, extends(psb_d_ell_sparse_mat) :: psb_d_oacc_ell_sparse_mat
        integer(psb_ipk_) :: devstate = is_host
    contains
        procedure, nopass  :: get_fmt        => d_oacc_ell_get_fmt
        procedure, pass(a) :: sizeof         => d_oacc_ell_sizeof
        procedure, pass(a) :: is_host        => d_oacc_ell_is_host
        procedure, pass(a) :: is_sync        => d_oacc_ell_is_sync
        procedure, pass(a) :: is_dev         => d_oacc_ell_is_dev
        procedure, pass(a) :: set_host       => d_oacc_ell_set_host
        procedure, pass(a) :: set_sync       => d_oacc_ell_set_sync
        procedure, pass(a) :: set_dev        => d_oacc_ell_set_dev
        procedure, pass(a) :: sync_space     => d_oacc_ell_sync_space
        procedure, pass(a) :: sync           => d_oacc_ell_sync
        procedure, pass(a) :: free           => d_oacc_ell_free

    end type psb_d_oacc_ell_sparse_mat

    contains

        subroutine d_oacc_ell_free(a)
            use psb_base_mod
            implicit none 
            class(psb_d_oacc_ell_sparse_mat), intent(inout) :: a
            integer(psb_ipk_) :: info
        
            if (allocated(a%val)) then
                !$acc exit data delete(a%val)
            end if
            if (allocated(a%ja)) then
                !$acc exit data delete(a%ja)
            end if
            if (allocated(a%irn)) then
                !$acc exit data delete(a%irn)
            end if
            if (allocated(a%idiag)) then
                !$acc exit data delete(a%idiag)
            end if
        
            call a%psb_d_ell_sparse_mat%free()
        
            return
        end subroutine d_oacc_ell_free
        

        function d_oacc_ell_sizeof(a) result(res)
            implicit none
            class(psb_d_oacc_ell_sparse_mat), intent(in) :: a
            integer(psb_epk_) :: res
        
            if (a%is_dev()) call a%sync()
        
            res = 8
            res = res + psb_sizeof_dp * size(a%val)
            res = res + psb_sizeof_ip * size(a%ja)
            res = res + psb_sizeof_ip * size(a%irn)
            res = res + psb_sizeof_ip * size(a%idiag)
        
        end function d_oacc_ell_sizeof    
            
        subroutine d_oacc_ell_sync_space(a)
            implicit none
            class(psb_d_oacc_ell_sparse_mat), intent(inout) :: a
        
            if (allocated(a%val)) then
                call d_oacc_create_dev(a%val)
            end if
            if (allocated(a%ja)) then
                call i_oacc_create_dev(a%ja)
            end if
            if (allocated(a%irn)) then
                call i_oacc_create_dev_scalar(a%irn)
            end if
            if (allocated(a%idiag)) then
                call i_oacc_create_dev_scalar(a%idiag)
            end if
        
        contains
            subroutine d_oacc_create_dev(v)
                implicit none
                real(psb_dpk_), intent(in) :: v(:,:)
                !$acc enter data copyin(v)
            end subroutine d_oacc_create_dev
        
            subroutine i_oacc_create_dev(v)
                implicit none
                integer(psb_ipk_), intent(in) :: v(:,:)
                !$acc enter data copyin(v)
            end subroutine i_oacc_create_dev
        
            subroutine i_oacc_create_dev_scalar(v)
                implicit none
                integer(psb_ipk_), intent(in) :: v(:)
                !$acc enter data copyin(v)
            end subroutine i_oacc_create_dev_scalar
        end subroutine d_oacc_ell_sync_space
        
        

        function d_oacc_ell_is_host(a) result(res)
            implicit none
            class(psb_d_oacc_ell_sparse_mat), intent(in) :: a
            logical :: res

            res = (a%devstate == is_host)
        end function d_oacc_ell_is_host

        function d_oacc_ell_is_sync(a) result(res)
            implicit none
            class(psb_d_oacc_ell_sparse_mat), intent(in) :: a
            logical :: res

            res = (a%devstate == is_sync)
        end function d_oacc_ell_is_sync

        function d_oacc_ell_is_dev(a) result(res)
            implicit none
            class(psb_d_oacc_ell_sparse_mat), intent(in) :: a
            logical :: res

            res = (a%devstate == is_dev)
        end function d_oacc_ell_is_dev

        subroutine d_oacc_ell_set_host(a)
            implicit none
            class(psb_d_oacc_ell_sparse_mat), intent(inout) :: a

            a%devstate = is_host
        end subroutine d_oacc_ell_set_host

        subroutine d_oacc_ell_set_sync(a)
            implicit none
            class(psb_d_oacc_ell_sparse_mat), intent(inout) :: a

            a%devstate = is_sync
        end subroutine d_oacc_ell_set_sync

        subroutine d_oacc_ell_set_dev(a)
            implicit none
            class(psb_d_oacc_ell_sparse_mat), intent(inout) :: a

            a%devstate = is_dev
        end subroutine d_oacc_ell_set_dev

        function d_oacc_ell_get_fmt() result(res)
            implicit none
            character(len=5) :: res
            res = 'ELL_oacc'
        end function d_oacc_ell_get_fmt

        subroutine d_oacc_ell_sync(a)
            implicit none
            class(psb_d_oacc_ell_sparse_mat), target, intent(in) :: a
            class(psb_d_oacc_ell_sparse_mat), pointer :: tmpa
            integer(psb_ipk_) :: info
        
            tmpa => a
            if (a%is_dev()) then
                call d_oacc_ell_to_host(a%val)
                call i_oacc_ell_to_host(a%ja)
                call i_oacc_ell_to_host_scalar(a%irn)
                call i_oacc_ell_to_host_scalar(a%idiag)
            else if (a%is_host()) then
                call d_oacc_ell_to_dev(a%val)
                call i_oacc_ell_to_dev(a%ja)
                call i_oacc_ell_to_dev_scalar(a%irn)
                call i_oacc_ell_to_dev_scalar(a%idiag)
            end if
            call tmpa%set_sync()
        end subroutine d_oacc_ell_sync

        subroutine d_oacc_ell_to_host(v)
            implicit none
            real(psb_dpk_), intent(in) :: v(:,:)
            !$acc update self(v)
        end subroutine d_oacc_ell_to_host
        
        subroutine d_oacc_ell_to_host_scalar(v)
            implicit none
            real(psb_dpk_), intent(in) :: v(:)
            !$acc update self(v)
        end subroutine d_oacc_ell_to_host_scalar
        
        subroutine i_oacc_ell_to_dev(v)
            implicit none
            integer(psb_ipk_), intent(in) :: v(:,:)
            !$acc update device(v)
        end subroutine i_oacc_ell_to_dev
        
        subroutine i_oacc_ell_to_dev_scalar(v)
            implicit none
            integer(psb_ipk_), intent(in) :: v(:)
            !$acc update device(v)
        end subroutine i_oacc_ell_to_dev_scalar
        
        subroutine i_oacc_ell_to_host(v)
            implicit none
            integer(psb_ipk_), intent(in) :: v(:,:)
            !$acc update self(v)
        end subroutine i_oacc_ell_to_host
        
        subroutine i_oacc_ell_to_host_scalar(v)
            implicit none
            integer(psb_ipk_), intent(in) :: v(:)
            !$acc update self(v)
        end subroutine i_oacc_ell_to_host_scalar        

       
end module psb_d_oacc_ell_mat_mod