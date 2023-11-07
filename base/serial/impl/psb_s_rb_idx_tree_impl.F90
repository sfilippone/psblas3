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
! package: psb_s_rb_idx_tree_impl
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
subroutine psb_s_rb_idx_tree_insert(this, idx, val)
  use psb_s_rb_idx_tree_mod, psb_protect_name => psb_s_rb_idx_tree_insert
  implicit none
  class(psb_s_rb_idx_tree), intent(inout) :: this
  integer(psb_ipk_), intent(in) :: idx
  real(psb_spk_), intent(in) :: val

  character(len=22) :: name
  type(psb_s_rb_idx_node), pointer :: new_node
  type(psb_s_rb_idx_node), pointer :: current, previous
  name='psb_rb_idx_tree_insert'

  allocate(new_node)
  new_node%idx = idx
  new_node%val = val
  nullify(new_node%left)
  nullify(new_node%right)
  nullify(new_node%parent)
  new_node%is_red = .true.


  if (.not. associated(this%root)) then
    this%root => new_node
    this%nnz = 1
    new_node%is_red = .false.
    return
  end if

  current => this%root

  do while (associated(current))
    previous => current

    if (idx == current%idx) then
      current%val = current%val + val
      deallocate(new_node)
      return
    else if (idx < current%idx) then
      current => current%left
    else

      current => current%right
    end if
  end do

  if (idx < previous%idx) then
    new_node%parent => previous
    previous%left => new_node
  else
    new_node%parent => previous
    previous%right => new_node
  end if

  call psb_s_rb_idx_tree_fix_insertion(this, new_node)

  this%nnz = this%nnz + 1
end subroutine psb_s_rb_idx_tree_insert

subroutine psb_s_rb_idx_tree_fix_insertion(this, node)
  use psb_s_rb_idx_tree_mod, psb_protect_name => psb_s_rb_idx_tree_fix_insertion
  implicit none
  class(psb_s_rb_idx_tree), intent(inout) :: this
  type(psb_s_rb_idx_node), pointer, intent(inout) :: node

  character(len=29) :: name
  type(psb_s_rb_idx_node), pointer :: current, parent, grand_parent, uncle
  name = 'psb_rb_idx_tree_fix_insertion'

  current => node
  parent => current%parent
  do while(associated(parent) .and. parent%is_red)
    ! grand parent exist because root can't be red
    grand_parent => parent%parent
    if (parent%idx < grand_parent%idx) then
      uncle => grand_parent%right
    else
      uncle => grand_parent%left
    end if

    if (associated(uncle) .and. uncle%is_red) then
      parent%is_red = .false.
      uncle%is_red = .false.
      grand_parent%is_red = .true.
      current => grand_parent
      parent => current%parent

      ! Left-Left case
    else if (current%idx < parent%idx .and. &
         parent%idx < grand_parent%idx) then
      call psb_s_rb_idx_tree_rotate_right(grand_parent)
      call psb_s_rb_idx_tree_swap_colors(parent, grand_parent)

      if (this%root%idx == grand_parent%idx) this%root => parent

      return
      ! Left-Right case
    else if (current%idx > parent%idx .and. &
         parent%idx < grand_parent%idx) then
      call psb_s_rb_idx_tree_rotate_left(parent)
      call psb_s_rb_idx_tree_rotate_right(grand_parent)
      call psb_s_rb_idx_tree_swap_colors(current, grand_parent)

      if (this%root%idx == grand_parent%idx) this%root => current

      return
      ! Right-Right case
    else if (current%idx > parent%idx .and. &
         parent%idx > grand_parent%idx) then
      call psb_s_rb_idx_tree_rotate_left(grand_parent)
      call psb_s_rb_idx_tree_swap_colors(parent, grand_parent)

      if (this%root%idx == grand_parent%idx) this%root => parent

      return
      ! Right-Left case
    else
      call psb_s_rb_idx_tree_rotate_right(parent)
      call psb_s_rb_idx_tree_rotate_left(grand_parent)
      call psb_s_rb_idx_tree_swap_colors(current, grand_parent)

      if (this%root%idx == grand_parent%idx) this%root => current

      return
    end if
  end do

  this%root%is_red = .false.
end subroutine psb_s_rb_idx_tree_fix_insertion

subroutine psb_s_rb_idx_tree_swap_colors(n1, n2)
  use psb_s_rb_idx_tree_mod, psb_protect_name => psb_s_rb_idx_tree_swap_colors
  implicit none
  type(psb_s_rb_idx_node), pointer, intent(inout) :: n1, n2

  character(len=27) :: name
  logical :: tmp
  name='psb_rb_idx_tree_swap_colors'

  tmp = n1%is_red
  n1%is_red = n2%is_red
  n2%is_red = tmp
end subroutine psb_s_rb_idx_tree_swap_colors

subroutine psb_s_rb_idx_tree_rotate_right(node)
  use psb_s_rb_idx_tree_mod, psb_protect_name => psb_s_rb_idx_tree_rotate_right
  implicit none
  type(psb_s_rb_idx_node), pointer, intent(inout) :: node

  character(len=28) :: name
  type(psb_s_rb_idx_node), pointer :: l, lr
  name='psb_rb_idx_tree_rotate_right'

  if (.not. associated(node%left)) return

  l => node%left
  lr => l%right
  node%left => lr

  if (associated(lr)) lr%parent => node

  if (associated(node%parent)) then
    if (node%idx < node%parent%idx) then
      node%parent%left => l
    else
      node%parent%right => l
    end if
  end if

  l%parent => node%parent
  node%parent => l

  l%right => node
end subroutine psb_s_rb_idx_tree_rotate_right

subroutine psb_s_rb_idx_tree_rotate_left(node)
  use psb_s_rb_idx_tree_mod, psb_protect_name => psb_s_rb_idx_tree_rotate_left
  implicit none
  type(psb_s_rb_idx_node), pointer, intent(inout) :: node

  character(len=27) :: name
  type(psb_s_rb_idx_node), pointer :: r, rl
  name='psb_rb_idx_tree_rotate_left'

  if (.not. associated(node%right)) return

  r => node%right
  rl => r%left
  node%right => rl

  if (associated(rl)) rl%parent => node

  if (associated(node%parent)) then
    if (node%idx < node%parent%idx) then
      node%parent%left => r
    else
      node%parent%right => r
    end if
  end if

  r%parent => node%parent
  node%parent => r

  r%left => node
end subroutine psb_s_rb_idx_tree_rotate_left

subroutine psb_s_rb_idx_tree_scalar_sparse_row_mul(tree, scalar, mat, row_num)
  use psb_s_rb_idx_tree_mod, psb_protect_name => psb_s_rb_idx_tree_scalar_sparse_row_mul
  use psb_s_csr_mat_mod, only : psb_s_csr_sparse_mat
  implicit none
  type(psb_s_rb_idx_tree), intent(inout) :: tree
  real(psb_spk_), intent(in) :: scalar
  type(psb_s_csr_sparse_mat), intent(in) :: mat
  integer(psb_ipk_), intent(in) :: row_num

  character(len=37) :: name
  integer(psb_ipk_) :: i
  name='psb_rb_idx_tree_scalar_sparse_row_mul'

  do i = mat%irp(row_num), mat%irp(row_num + 1) - 1
    call tree%insert(mat%ja(i),scalar * mat%val(i))
  end do

end subroutine psb_s_rb_idx_tree_scalar_sparse_row_mul

subroutine psb_s_rb_idx_tree_merge(trees, mat)
#if defined(OPENMP)
  use omp_lib
#endif
  use psb_realloc_mod
  use psb_s_rb_idx_tree_mod, psb_protect_name => psb_s_rb_idx_tree_merge
  use psb_s_csr_mat_mod, only : psb_s_csr_sparse_mat
  implicit none
  type(psb_s_rb_idx_tree), allocatable, intent(inout) :: trees(:)
  type(psb_s_csr_sparse_mat), intent(inout) :: mat

  character(len=21) :: name
  integer(psb_ipk_) :: i, j, rows, info, nnz
  type(psb_s_rb_idx_node), pointer :: current, previous
  name='psb_rb_idx_tree_merge'

  rows = size(trees)

  mat%irp(1) = 1

  do i=1, rows
    mat%irp(i + 1) = mat%irp(i) + trees(i)%nnz
  end do

  nnz = mat%irp(rows + 1)
  call psb_realloc(nnz, mat%val, info)
  call psb_realloc(nnz, mat%ja, info)

#if defined(OPENMP)
  !$omp parallel do schedule(static), private(current, previous, j)
#endif
  do i = 1, size(trees)
    j = 0
    current => trees(i)%root
    do while(associated(current))
      ! go to the left-most node
      do while(associated(current%left))
        current => current%left
      end do
      mat%val(j + mat%irp(i)) = current%val
      mat%ja(j + mat%irp(i)) = current%idx
      j = j + 1

      previous => current
      if (associated(current%right)) then
        if (associated(current%parent)) then
          current%parent%left => current%right
        end if
        current%right%parent => current%parent
        current => current%right
      else
        current => current%parent
        if (associated(current)) nullify(current%left)
      end if
      deallocate(previous)
    end do
  end do
#if defined(OPENMP)
  !$omp end parallel do
#endif
end subroutine psb_s_rb_idx_tree_merge
