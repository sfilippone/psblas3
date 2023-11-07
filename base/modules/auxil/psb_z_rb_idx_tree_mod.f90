!
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone
!        Alfredo Buttari
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!
!
!
! package: psb_z_rb_idx_tree_mod
! 
! Red black tree implementation ordered by index
!
! Each node contains and index and a double precision value
!
! The tree should always be well balanced
!
! inserting a node with an existing index will 
! add up the new value to the old one
! Contributed by Dimitri Walther
! 
module psb_z_rb_idx_tree_mod
  use psb_const_mod
  implicit none

  type :: psb_z_rb_idx_node
    integer(psb_ipk_) :: idx
    complex(psb_dpk_) :: val
    type(psb_z_rb_idx_node), pointer :: left, right, parent
    logical :: is_red
  end type psb_z_rb_idx_node

  type :: psb_z_rb_idx_tree
    type(psb_z_rb_idx_node), pointer :: root
    integer(psb_ipk_) :: nnz

  contains

    procedure :: insert => psb_z_rb_idx_tree_insert
  end type psb_z_rb_idx_tree

  interface psb_rb_idx_tree_insert
    subroutine psb_z_rb_idx_tree_insert(this, idx, val)
      import :: psb_ipk_, psb_dpk_, psb_z_rb_idx_tree
      implicit none
      class(psb_z_rb_idx_tree), intent(inout) :: this
      integer(psb_ipk_), intent(in) :: idx
      complex(psb_dpk_), intent(in) :: val
    end subroutine psb_z_rb_idx_tree_insert
  end interface psb_rb_idx_tree_insert

  interface psb_rb_idx_tree_scalar_sparse_row_mul
    subroutine psb_z_rb_idx_tree_scalar_sparse_row_mul(tree, scalar, mat, row_num)
      use psb_z_csr_mat_mod, only : psb_z_csr_sparse_mat
      import :: psb_ipk_, psb_dpk_, psb_z_rb_idx_tree
      implicit none
      type(psb_z_rb_idx_tree), intent(inout) :: tree
      complex(psb_dpk_), intent(in) :: scalar
      type(psb_z_csr_sparse_mat), intent(in) :: mat
      integer(psb_ipk_), intent(in) :: row_num
    end subroutine psb_z_rb_idx_tree_scalar_sparse_row_mul
  end interface psb_rb_idx_tree_scalar_sparse_row_mul

  interface psb_rb_idx_tree_merge
    subroutine psb_z_rb_idx_tree_merge(trees, mat)
      use psb_z_csr_mat_mod, only : psb_z_csr_sparse_mat
      import :: psb_z_rb_idx_tree
      type(psb_z_rb_idx_tree), allocatable, intent(inout) :: trees(:)
      type(psb_z_csr_sparse_mat), intent(inout) :: mat
    end subroutine psb_z_rb_idx_tree_merge
  end interface psb_rb_idx_tree_merge

  interface psb_rb_idx_tree_fix_insertion
    subroutine psb_z_rb_idx_tree_fix_insertion(this, node)
      import :: psb_z_rb_idx_tree, psb_z_rb_idx_node
      implicit none
      class(psb_z_rb_idx_tree), intent(inout) :: this
      type(psb_z_rb_idx_node), pointer, intent(inout) :: node
    end subroutine psb_z_rb_idx_tree_fix_insertion
  end interface psb_rb_idx_tree_fix_insertion

  interface psb_rb_idx_tree_swap_colors
    subroutine psb_z_rb_idx_tree_swap_colors(n1, n2)
      import :: psb_z_rb_idx_node
      implicit none
      type(psb_z_rb_idx_node), pointer, intent(inout) :: n1, n2
    end subroutine psb_z_rb_idx_tree_swap_colors
  end interface psb_rb_idx_tree_swap_colors

  interface psb_rb_idx_tree_rotate_right
    subroutine psb_z_rb_idx_tree_rotate_right(node)
      import :: psb_z_rb_idx_node
      implicit none
      type(psb_z_rb_idx_node), pointer, intent(inout) :: node
    end subroutine psb_z_rb_idx_tree_rotate_right
  end interface psb_rb_idx_tree_rotate_right

  interface psb_rb_idx_tree_rotate_left
    subroutine psb_z_rb_idx_tree_rotate_left(node)
      import :: psb_z_rb_idx_node
      implicit none
      type(psb_z_rb_idx_node), pointer, intent(inout) :: node
    end subroutine psb_z_rb_idx_tree_rotate_left
  end interface psb_rb_idx_tree_rotate_left
end module psb_z_rb_idx_tree_mod
