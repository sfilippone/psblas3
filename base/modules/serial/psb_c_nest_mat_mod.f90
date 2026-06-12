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
!         software without specific without permission.
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
! module: psb_c_nest_mat_mod
! Author: Simone Staccone (Stack-1)
!
! Defines psb_c_nest_sparse_mat: a block-structured distributed sparse
! matrix for double precision real arithmetic.
!
! The matrix is stored as a 2-D array of psb_cspmat_type sub-matrices.
! Block presence is determined directly from the sub-matrix storage: a block
! (i,j) is present iff mats(i,j)%a is allocated (absent blocks contribute zero
! to any product). There is no separate presence flag array.
!
! Descriptor convention (current nested design)
! ---------------------------------------------
! Each matrix block (i,j) is associated with descs(i,j) from the
! corresponding psb_desc_nest_type. Nested tools (psb_spall_nest,
! psb_spins_nest, psb_spasb_nest, psb_spmm) consistently pass
! descs(i,j) together with mats(i,j).
!
! A block may be structurally absent (NULL/zero): this is represented by
! mats(i,j) left unbuilt (mats(i,j)%a not allocated). In that case the
! block contributes zero and is skipped by nested kernels.
!
! Descriptor storage is distinct from matrix presence: descriptors are
! typically defined for all block positions in descs(:,:), while actual
! matrix blocks may be present only on a subset.
!
! Reference examples in test/pdegen:
!   * psb_c_pde_nest.full.F90  (A(2,2) left NULL, mats(2,2)%a not allocated)
!   * psb_c_nest_tools.F90 and psb_c_pde_nest_full_tools.F90
!     (2-D desc_nest%descs(i,j) used in nested allocation/assembly).
!
module psb_c_nest_mat_mod
  use psb_c_mat_mod
  implicit none

  type :: psb_c_nest_sparse_mat
    integer(psb_ipk_)                         :: nrblocks = 0
    integer(psb_ipk_)                         :: ncblocks = 0
    type(psb_cspmat_type), allocatable        :: mats(:,:)
  contains
    procedure :: get_nrblocks  => psb_c_nest_mat_get_nrb
    procedure :: get_ncblocks  => psb_c_nest_mat_get_ncb
    procedure :: has_block     => psb_c_nest_mat_has_block
    procedure :: sizeof        => psb_c_nest_mat_sizeof
    procedure :: free          => psb_c_nest_mat_free
  end type psb_c_nest_sparse_mat

contains

  !  get_nrblocks / get_ncblocks
  function psb_c_nest_mat_get_nrb(a) result(n)
    class(psb_c_nest_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: n
    n = a%nrblocks
  end function psb_c_nest_mat_get_nrb

  function psb_c_nest_mat_get_ncb(a) result(n)
    class(psb_c_nest_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: n
    n = a%ncblocks
  end function psb_c_nest_mat_get_ncb

  !  has_block: return .true. if block (i,j) is non-null
  function psb_c_nest_mat_has_block(a, i_block_row, j_block_col) result(has)
    class(psb_c_nest_sparse_mat), intent(in) :: a
    integer(psb_ipk_),            intent(in) :: i_block_row, j_block_col
    logical :: has

    has = .false.
    if (i_block_row < 1 .or. i_block_row > a%nrblocks) return
    if (j_block_col < 1 .or. j_block_col > a%ncblocks) return
    if (.not. allocated(a%mats)) return
    ! P3: presence is determined solely by whether the sub-matrix has been
    ! built (its polymorphic storage %a is allocated). No parallel flag array.
    has = allocated(a%mats(i_block_row, j_block_col)%a)
  end function psb_c_nest_mat_has_block

  !  sizeof: total storage across all allocated sub-matrices
  function psb_c_nest_mat_sizeof(a) result(total_bytes)
    class(psb_c_nest_sparse_mat), intent(in) :: a
    integer(psb_epk_) :: total_bytes
    integer(psb_ipk_) :: i_block_row, j_block_col

    total_bytes = 0_psb_epk_
    if (allocated(a%mats)) then
      do j_block_col = 1, a%ncblocks
        do i_block_row = 1, a%nrblocks
          if (allocated(a%mats(i_block_row, j_block_col)%a)) &
               & total_bytes = total_bytes + a%mats(i_block_row, j_block_col)%sizeof()
        end do
      end do
    end if
  end function psb_c_nest_mat_sizeof

  !  free: release all sub-matrices
  subroutine psb_c_nest_mat_free(a, info)
    class(psb_c_nest_sparse_mat), intent(inout) :: a
    integer(psb_ipk_),            intent(out)   :: info

    integer(psb_ipk_) :: i_block_row, j_block_col, local_info

    info = 0
    if (allocated(a%mats)) then
      do j_block_col = 1, a%ncblocks
        do i_block_row = 1, a%nrblocks
          if (allocated(a%mats(i_block_row, j_block_col)%a)) then
            call a%mats(i_block_row, j_block_col)%free()
          end if
        end do
      end do
      deallocate(a%mats, stat=local_info)
      if (local_info /= 0 .and. info == 0) info = local_info
    end if
    a%nrblocks = 0
    a%ncblocks = 0
  end subroutine psb_c_nest_mat_free

end module psb_c_nest_mat_mod
