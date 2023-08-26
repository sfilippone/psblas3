! Red black tree implementation ordered by index
!
! Each node contains and index and a double precision value
!
! The tree should always be well balanced
!
! inserting a node with an existing index will 
! add up the new value to the old one
module psb_d_rb_idx_tree_mod
    use psb_const_mod
    implicit none
    type :: psb_d_rb_idx_node
        integer(psb_ipk_) :: idx
        real(psb_dpk_) :: val
        type(psb_d_rb_idx_node), pointer :: left, right, parent
        logical :: is_red
    end type psb_d_rb_idx_node

    type :: psb_d_rb_idx_tree
        type(psb_d_rb_idx_node), pointer :: root
        integer(psb_ipk_) :: nnz

        contains

        procedure :: insert => psb_d_rb_idx_tree_insert
    end type psb_d_rb_idx_tree

    interface psb_rb_idx_tree_insert
    subroutine psb_d_rb_idx_tree_insert(this, idx, val)
        import :: psb_ipk_, psb_dpk_, psb_d_rb_idx_tree
        implicit none
        class(psb_d_rb_idx_tree), intent(inout) :: this
        integer(psb_ipk_), intent(in) :: idx
        real(psb_dpk_), intent(in) :: val
    end subroutine psb_d_rb_idx_tree_insert
    end interface psb_rb_idx_tree_insert
    
    interface psb_rb_idx_tree_scalar_sparse_row_mul
    subroutine psb_d_rb_idx_tree_scalar_sparse_row_mul(tree, scalar, mat, row_num)
        use psb_d_csr_mat_mod, only : psb_d_csr_sparse_mat
        import :: psb_ipk_, psb_dpk_, psb_d_rb_idx_tree
        implicit none
        type(psb_d_rb_idx_tree), intent(inout) :: tree
        real(psb_dpk_), intent(in) :: scalar
        type(psb_d_csr_sparse_mat), intent(in) :: mat
        integer(psb_ipk_), intent(in) :: row_num
    end subroutine psb_d_rb_idx_tree_scalar_sparse_row_mul
    end interface psb_rb_idx_tree_scalar_sparse_row_mul
 
    interface psb_rb_idx_tree_merge
    subroutine psb_d_rb_idx_tree_merge(trees, mat)
        use psb_d_csr_mat_mod, only : psb_d_csr_sparse_mat
        import :: psb_d_rb_idx_tree
        type(psb_d_rb_idx_tree), allocatable, intent(inout) :: trees(:)
        type(psb_d_csr_sparse_mat), intent(inout) :: mat
    end subroutine psb_d_rb_idx_tree_merge
    end interface psb_rb_idx_tree_merge
    
    interface psb_rb_idx_tree_fix_insertion
    subroutine psb_d_rb_idx_tree_fix_insertion(this, node)
        import :: psb_d_rb_idx_tree, psb_d_rb_idx_node
        implicit none
        class(psb_d_rb_idx_tree), intent(inout) :: this
        type(psb_d_rb_idx_node), pointer, intent(inout) :: node
    end subroutine psb_d_rb_idx_tree_fix_insertion
    end interface psb_rb_idx_tree_fix_insertion

    interface psb_rb_idx_tree_swap_colors
    subroutine psb_d_rb_idx_tree_swap_colors(n1, n2)
        import :: psb_d_rb_idx_node
        implicit none
        type(psb_d_rb_idx_node), pointer, intent(inout) :: n1, n2
    end subroutine psb_d_rb_idx_tree_swap_colors
    end interface psb_rb_idx_tree_swap_colors

    interface psb_rb_idx_tree_rotate_right
    subroutine psb_d_rb_idx_tree_rotate_right(node)
        import :: psb_d_rb_idx_node
        implicit none
        type(psb_d_rb_idx_node), pointer, intent(inout) :: node
    end subroutine psb_d_rb_idx_tree_rotate_right
    end interface psb_rb_idx_tree_rotate_right

    interface psb_rb_idx_tree_rotate_left
    subroutine psb_d_rb_idx_tree_rotate_left(node)
        import :: psb_d_rb_idx_node
        implicit none
        type(psb_d_rb_idx_node), pointer, intent(inout) :: node
    end subroutine psb_d_rb_idx_tree_rotate_left
    end interface psb_rb_idx_tree_rotate_left
end module psb_d_rb_idx_tree_mod