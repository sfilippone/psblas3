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
! module: psb_d_nest_mat_mod
!
! Defines psb_d_nest_sparse_mat: a block-structured distributed sparse
! matrix for double precision real arithmetic.
!
! The matrix is stored as a 2-D array of psb_dspmat_type sub-matrices.
! A companion logical array blk_present(i,j) flags which blocks are
! non-null (absent blocks contribute zero to any product).
!
! Descriptor convention (current nested design)
! ---------------------------------------------
! Each matrix block (i,j) is associated with descs(i,j) from the
! corresponding psb_desc_nest_type. Nested tools (psb_spall_nest,
! psb_spins_nest, psb_spasb_nest, psb_spmm) consistently pass
! descs(i,j) together with mats(i,j).
!
! A block may be structurally absent (NULL/zero): this is represented by
! blk_present(i,j)=.false. and mats(i,j) left unbuilt. In that case the
! block contributes zero and is skipped by nested kernels.
!
! Descriptor storage is distinct from matrix presence: descriptors are
! typically defined for all block positions in descs(:,:), while actual
! matrix blocks may be present only on a subset.
!
! Reference examples in test/pdegen:
!   * psb_d_pde_nest.full.F90  (A(2,2) left NULL, blk_present(2,2)=.false.)
!   * psb_d_nest_tools.F90 and psb_d_pde_nest_full_tools.F90
!     (2-D desc_nest%descs(i,j) used in nested allocation/assembly).
!
module psb_d_nest_mat_mod
  use psb_d_mat_mod
  implicit none

  type :: psb_d_nest_sparse_mat
    integer(psb_ipk_)                         :: nrblocks = 0
    integer(psb_ipk_)                         :: ncblocks = 0
    type(psb_dspmat_type), allocatable        :: mats(:,:)
    logical,               allocatable        :: blk_present(:,:)
  contains
    procedure :: get_nrblocks  => psb_d_nest_mat_get_nrb
    procedure :: get_ncblocks  => psb_d_nest_mat_get_ncb
    procedure :: has_block     => psb_d_nest_mat_has_block
    procedure :: sizeof        => psb_d_nest_mat_sizeof
    procedure :: free          => psb_d_nest_mat_free
  end type psb_d_nest_sparse_mat

contains

  !  get_nrblocks / get_ncblocks
  function psb_d_nest_mat_get_nrb(a) result(n)
    class(psb_d_nest_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: n
    n = a%nrblocks
  end function psb_d_nest_mat_get_nrb

  function psb_d_nest_mat_get_ncb(a) result(n)
    class(psb_d_nest_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: n
    n = a%ncblocks
  end function psb_d_nest_mat_get_ncb

  !  has_block: return .true. if block (i,j) is non-null
  function psb_d_nest_mat_has_block(a, i, j) result(hp)
    class(psb_d_nest_sparse_mat), intent(in) :: a
    integer(psb_ipk_),            intent(in) :: i, j
    logical :: hp

    hp = .false.
    if (i < 1 .or. i > a%nrblocks) return
    if (j < 1 .or. j > a%ncblocks) return
    if (.not. allocated(a%blk_present)) return
    hp = a%blk_present(i, j)
  end function psb_d_nest_mat_has_block

  !  sizeof: total storage across all allocated sub-matrices
  function psb_d_nest_mat_sizeof(a) result(s)
    class(psb_d_nest_sparse_mat), intent(in) :: a
    integer(psb_epk_) :: s
    integer(psb_ipk_) :: i, j

    s = 0_psb_epk_
    if (allocated(a%mats)) then
      do j = 1, a%ncblocks
        do i = 1, a%nrblocks
          if (a%blk_present(i, j)) s = s + a%mats(i, j)%sizeof()
        end do
      end do
    end if
  end function psb_d_nest_mat_sizeof

  !  free: release all sub-matrices
  subroutine psb_d_nest_mat_free(a, info)
    class(psb_d_nest_sparse_mat), intent(inout) :: a
    integer(psb_ipk_),            intent(out)   :: info

    integer(psb_ipk_) :: i, j, linfo

    info = 0
    if (allocated(a%mats)) then
      do j = 1, a%ncblocks
        do i = 1, a%nrblocks
          if (a%blk_present(i, j)) then
            call a%mats(i, j)%free()
          end if
        end do
      end do
      deallocate(a%mats, stat=linfo)
      if (linfo /= 0 .and. info == 0) info = linfo
    end if
    if (allocated(a%blk_present)) then
      deallocate(a%blk_present, stat=linfo)
      if (linfo /= 0 .and. info == 0) info = linfo
    end if
    a%nrblocks = 0
    a%ncblocks = 0
  end subroutine psb_d_nest_mat_free

end module psb_d_nest_mat_mod
