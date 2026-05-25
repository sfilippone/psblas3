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
! module: psb_d_nest_psblas_mod
!
! Parallel BLAS operations for the nested (block-structured) double
! precision real types.
!
! psb_d_nest_spmm
!   Computes  y = alpha * A_nest * x + beta * y  (block SpMV).
!   Three-phase algorithm:
!     Phase 1 — scale y upfront:
!       if beta == 0: zero all y(i)
!       elif beta /= 1: y(i) = beta * y(i)  for each block i
!     Phase 2 — single halo exchange per column block:
!       call psb_d_nest_halo(xnest, descs)
!       Populates ghost entries of x(j) using descs(1,j) for each j.
!       All descs(i,j) for fixed j share the same column space, so
!       one exchange covers all block-rows.
!     Phase 3 — local SpMM accumulation (no further communication):
!       For each present block (i,j):
!         y(i) += alpha * A(i,j) * x(j)
!         (psb_spmm called with doswap=.false. to skip internal halo)
!
! psb_d_nest_geaxpby
!   Computes  y(i) = alpha * x(i) + beta * y(i)  for each block i.
!
! psb_d_nest_genrm2
!   Computes  ||x||_2 = sqrt( sum_i ||x(i)||_2^2 )  with a single
!   global reduction.
!
! psb_d_nest_genrm2s
!   Subroutine form of psb_d_nest_genrm2 (result via intent(out) argument).
!
! psb_d_nest_gedot
!   Computes  dot(x,y) = sum_i dot(x(i), y(i))  with a single global
!   reduction.
!
! psb_d_nest_geamax
!   Computes  ||x||_inf = max_i ||x(i)||_inf  with a single global reduction.
!
! psb_d_nest_geasum
!   Computes  ||x||_1 = sum_i ||x(i)||_1  with a single global reduction.
!
! psb_d_nest_gemin
!   Computes  min(x) = min_i min(x(i))  with a single global reduction.
!
! psb_d_nest_minquotient
!   Computes  min(x/y) = min_i min(x(i)/y(i))  with a single global reduction.
!
! psb_d_nest_gemlt
!   Computes  y(i) = x(i) .* y(i)  element-wise for each block i.
!
! psb_d_nest_gediv
!   Computes  y(i) = x(i) ./ y(i)  element-wise for each block i.
!
! psb_d_nest_geinv
!   Computes  y(i) = 1 / x(i)  element-wise for each block i.
!
! psb_d_nest_geabs
!   Computes  y(i) = |x(i)|  element-wise for each block i.
!
! psb_d_nest_geaddconst
!   Computes  z(i) = x(i) + b  for each block i (b is a scalar).
!
! psb_d_nest_gecmp
!   Computes  z(i) = cmp(x(i), c)  for each block i (c is a scalar).
!
! psb_d_nest_mask
!   Applies mask operation to each block i; t is .true. iff all blocks
!   return .true.
!
! psb_d_nest_upd_xyz
!   Applies psb_upd_xyz(alpha,beta,gamma,delta, x(i),y(i),z(i)) for each block i.
!
! psb_d_nest_spsm
!   Block-diagonal triangular solve: applies psb_spsm to each diagonal
!   block (i,i) of tnest independently.
!
module psb_d_nest_psblas_mod
  use psb_desc_nest_mod
  use psb_d_nest_vect_mod
  use psb_d_nest_mat_mod
  use psb_d_mat_mod,    only : psb_dspmat_type, psb_csmm
  use psb_d_psblas_mod, only : psb_spmm, psb_geaxpby, psb_genrm2, psb_gedot, &
       & psb_geamax, psb_geasum, psb_gemin, psb_minquotient,                 &
       & psb_gemlt, psb_gediv, psb_geinv, psb_geabs, psb_geaddconst,         &
       & psb_gecmp, psb_mask, psb_upd_xyz, psb_spsm
  use psb_d_nest_comm_mod, only : psb_d_nest_halo, psb_d_nest_ovrl
  use psb_penv_mod,     only : psb_sum, psb_max, psb_min, psb_info
  use psb_const_mod,    only : psb_dpk_, psb_ipk_, psb_epk_, psb_ctxt_type, dzero, done
  implicit none

  private
  public :: psb_d_nest_spmm, psb_d_nest_geaxpby,                          &
            psb_d_nest_genrm2, psb_d_nest_genrm2s, psb_d_nest_gedot,     &
            psb_d_nest_geamax, psb_d_nest_geasum, psb_d_nest_gemin,      &
            psb_d_nest_minquotient,                                       &
            psb_d_nest_gemlt, psb_d_nest_gediv, psb_d_nest_geinv,        &
            psb_d_nest_halo, psb_d_nest_ovrl,                            &
            psb_d_nest_geabs, psb_d_nest_geaddconst, psb_d_nest_gecmp,   &
            psb_d_nest_mask, psb_d_nest_upd_xyz, psb_d_nest_spsm

contains

  !  y = alpha * A_nest * x + beta * y
  subroutine psb_d_nest_spmm(alpha, anest, xnest, beta, ynest, descs, info, trans)

    real(psb_dpk_),               intent(in)    :: alpha, beta
    type(psb_d_nest_sparse_mat),  intent(in)    :: anest
    type(psb_d_nest_vect_type),   intent(inout) :: xnest
    type(psb_d_nest_vect_type),   intent(inout) :: ynest
    type(psb_desc_nest_type),     intent(in)    :: descs
    integer(psb_ipk_),            intent(out)   :: info
    character, optional,          intent(in)    :: trans

    integer(psb_ipk_) :: i, j
    character         :: trans_

    info = 0
    if (present(trans)) then
      trans_ = trans
    else
      trans_ = 'N'
    end if

    if (beta == dzero) then
      do i = 1, anest%nrblocks
        call ynest%vects(i)%zero()
      end do
    else if (beta /= done) then
      do i = 1, anest%nrblocks
        call ynest%vects(i)%scal(beta)
      end do
    end if

    call psb_d_nest_halo(xnest, descs, info)
    if (info /= 0) return

    do i = 1, anest%nrblocks
      do j = 1, anest%ncblocks
        if (anest%has_block(i, j)) then
          ! y(i) += alpha * A(i,j) * x(j) (doswap=.false. skips internal halo, already done in Phase 2)
          call psb_spmm(alpha, anest%mats(i,j), xnest%vects(j), &
              &         done, ynest%vects(i), descs%descs(i,j), info, trans=trans_, doswap=.false.)
          if (info /= 0) return
        end if
      end do
    end do
  end subroutine psb_d_nest_spmm

  !  y(i) = alpha * x(i) + beta * y(i)  for each block i
  subroutine psb_d_nest_geaxpby(alpha, xnest, beta, ynest, descs, info)
    real(psb_dpk_),             intent(in)    :: alpha, beta
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_d_nest_vect_type), intent(inout) :: ynest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info

    integer(psb_ipk_) :: i

    info = 0
    do i = 1, xnest%nblocks
      call psb_geaxpby(alpha, xnest%vects(i), beta, ynest%vects(i), &
           &           descs%descs(i,i), info)
      if (info /= 0) return
    end do
  end subroutine psb_d_nest_geaxpby

  !  ||x||_2 = sqrt( sum_i ||x(i)||_2^2 ) 
  !  Uses a single global MPI_Allreduce across all blocks.
  !  global (optional, default .true.): if .false., skips MPI_Allreduce and returns
  !  the process-local partial norm; use when the caller manages the reduction itself.

  function psb_d_nest_genrm2(xnest, descs, info, global) result(res)
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info
    logical, optional,          intent(in)    :: global
    real(psb_dpk_) :: res

    integer(psb_ipk_) :: i
    real(psb_dpk_)    :: loc_sum, blk_nrm
    logical           :: global_
    type(psb_ctxt_type) :: ctxt

    global_ = .true.
    if (present(global)) global_ = global

    info = 0
    loc_sum = dzero
    do i = 1, xnest%nblocks
      ! global=.false. returns local partial norm (sqrt of local partial sum)
      blk_nrm = psb_genrm2(xnest%vects(i), descs%descs(i,i), info, global=.false.)
      if (info /= 0) then
        res = dzero
        return
      end if
      loc_sum = loc_sum + blk_nrm * blk_nrm
    end do

    if (global_) then
      ctxt = descs%descs(1,1)%get_context()
      call psb_sum(ctxt, loc_sum)
    end if

    res = sqrt(loc_sum)
  end function psb_d_nest_genrm2

  !  psb_d_nest_gedot 
  !  dot(x, y) = sum_i dot(x(i), y(i)) 
  !  Uses a single global MPI_Allreduce across all blocks. 
  function psb_d_nest_gedot(xnest, ynest, descs, info, global) result(res)
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_d_nest_vect_type), intent(inout) :: ynest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info
    logical, optional,          intent(in)    :: global
    real(psb_dpk_) :: res

    integer(psb_ipk_) :: i
    real(psb_dpk_)    :: loc_sum, blk_dot
    logical           :: global_
    type(psb_ctxt_type) :: ctxt

    global_ = .true.
    if (present(global)) global_ = global

    info = 0
    loc_sum = dzero
    do i = 1, xnest%nblocks
      blk_dot = psb_gedot(xnest%vects(i), ynest%vects(i), descs%descs(i,i), &
           &              info, global=.false.)
      if (info /= 0) then
        res = dzero
        return
      end if
      loc_sum = loc_sum + blk_dot
    end do

    if (global_) then
      ctxt = descs%descs(1,1)%get_context()
      call psb_sum(ctxt, loc_sum)
    end if

    res = loc_sum
  end function psb_d_nest_gedot

  !  Subroutine form: res = ||x||_2  (single global reduction)
  subroutine psb_d_nest_genrm2s(res, xnest, descs, info, global)
    real(psb_dpk_),             intent(out)   :: res
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info
    logical, optional,          intent(in)    :: global

    res = psb_d_nest_genrm2(xnest, descs, info, global)
  end subroutine psb_d_nest_genrm2s

  !  ||x||_inf = max_i ||x(i)||_inf  (single global reduction)
  function psb_d_nest_geamax(xnest, descs, info, global) result(res)
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info
    logical, optional,          intent(in)    :: global
    real(psb_dpk_) :: res

    integer(psb_ipk_) :: i
    real(psb_dpk_)    :: blk_val
    logical           :: global_
    type(psb_ctxt_type) :: ctxt

    global_ = .true.
    if (present(global)) global_ = global

    info = 0
    res  = dzero
    do i = 1, xnest%nblocks
      blk_val = psb_geamax(xnest%vects(i), descs%descs(i,i), info, global=.false.)
      if (info /= 0) return
      if (blk_val > res) res = blk_val
    end do

    if (global_) then
      ctxt = descs%descs(1,1)%get_context()
      call psb_max(ctxt, res)
    end if
  end function psb_d_nest_geamax
  
  !  ||x||_1 = sum_i ||x(i)||_1  (single global reduction)
  ! function psb_d_nest_geasum(xnest, descs, info, global) result(res)
  !   type(psb_d_nest_vect_type), intent(inout) :: xnest
  !   type(psb_desc_nest_type),   intent(in)    :: descs
  !   integer(psb_ipk_),          intent(out)   :: info
  !   logical, optional,          intent(in)    :: global
  !   real(psb_dpk_) :: res

  !   integer(psb_ipk_) :: i
  !   real(psb_dpk_)    :: blk_val
  !   logical           :: global_
  !   type(psb_ctxt_type) :: ctxt

  !   global_ = .true.
  !   if (present(global)) global_ = global

  !   info = 0
  !   res  = dzero
  !   do i = 1, xnest%nblocks
  !     blk_val = psb_geasum(xnest%vects(i), descs%descs(i,i), info, global=.false.)
  !     if (info /= 0) return
  !     res = res + blk_val
  !   end do

  !   if (global_) then
  !     ctxt = descs%descs(1,1)%get_context()
  !     call psb_sum(ctxt, res)
  !   end if
  ! end function psb_d_nest_geasum

  !  ||x||_1 = sum_i ||x(i)||_1 
  function psb_d_nest_geasum(xnest, descs, info, global) result(res)
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info
    logical, optional,          intent(in)    :: global
    real(psb_dpk_) :: res

    integer(psb_ipk_) :: i
    integer(psb_ipk_) :: nloc
    real(psb_dpk_)    :: blk_val
    real(psb_dpk_), allocatable :: blk_vals(:)
    logical           :: global_
    type(psb_ctxt_type) :: ctxt

    global_ = .true.
    if (present(global)) global_ = global

    info = 0
    res  = dzero
    do i = 1, xnest%nblocks
      nloc = descs%descs(i,i)%get_local_rows()
      blk_vals = xnest%vects(i)%get_vect(nloc)
      if (size(blk_vals) > 0) then
        blk_val = sum(abs(blk_vals))
      else
        blk_val = dzero
      end if
      res = res + blk_val
    end do

    if (global_) then
      ctxt = descs%descs(1,1)%get_context()
      call psb_sum(ctxt, res)
    end if
  end function psb_d_nest_geasum

  !  min(x) = min_i min(x(i))
  function psb_d_nest_gemin(xnest, descs, info, global) result(res)
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info
    logical, optional,          intent(in)    :: global
    real(psb_dpk_) :: res

    integer(psb_ipk_) :: i
    real(psb_dpk_)    :: blk_val
    logical           :: global_
    type(psb_ctxt_type) :: ctxt

    global_ = .true.
    if (present(global)) global_ = global

    info = 0
    res  = huge(dzero)
    do i = 1, xnest%nblocks
      blk_val = psb_gemin(xnest%vects(i), descs%descs(i,i), info, global=.false.)
      if (info /= 0) return
      if (blk_val < res) res = blk_val
    end do

    if (global_) then
      ctxt = descs%descs(1,1)%get_context()
      call psb_min(ctxt, res)
    end if
  end function psb_d_nest_gemin

  !  min(x/y) = min_i min(x(i)/y(i))  (single global reduction)
  function psb_d_nest_minquotient(xnest, ynest, descs, info, global) result(res)
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_d_nest_vect_type), intent(inout) :: ynest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info
    logical, optional,          intent(in)    :: global
    real(psb_dpk_) :: res

    integer(psb_ipk_) :: i
    real(psb_dpk_)    :: blk_val
    logical           :: global_
    type(psb_ctxt_type) :: ctxt

    global_ = .true.
    if (present(global)) global_ = global

    info = 0
    res  = huge(dzero)
    do i = 1, xnest%nblocks
      blk_val = psb_minquotient(xnest%vects(i), ynest%vects(i), &
           &                    descs%descs(i,i), info, global=.false.)
      if (info /= 0) return
      if (blk_val < res) res = blk_val
    end do

    if (global_) then
      ctxt = descs%descs(1,1)%get_context()
      call psb_min(ctxt, res)
    end if
  end function psb_d_nest_minquotient

  !  y(i) = x(i) .* y(i)  element-wise for each block i
  subroutine psb_d_nest_gemlt(xnest, ynest, descs, info)
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_d_nest_vect_type), intent(inout) :: ynest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info

    integer(psb_ipk_) :: i

    info = 0
    do i = 1, xnest%nblocks
      call psb_gemlt(xnest%vects(i), ynest%vects(i), descs%descs(i,i), info)
      if (info /= 0) return
    end do
  end subroutine psb_d_nest_gemlt

  !  y(i) = x(i) ./ y(i)  element-wise for each block i
  subroutine psb_d_nest_gediv(xnest, ynest, descs, info)
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_d_nest_vect_type), intent(inout) :: ynest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info

    integer(psb_ipk_) :: i

    info = 0
    do i = 1, xnest%nblocks
      call psb_gediv(xnest%vects(i), ynest%vects(i), descs%descs(i,i), info)
      if (info /= 0) return
    end do
  end subroutine psb_d_nest_gediv

  !  y(i) = 1/x(i)  element-wise for each block i
  subroutine psb_d_nest_geinv(xnest, ynest, descs, info)
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_d_nest_vect_type), intent(inout) :: ynest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info

    integer(psb_ipk_) :: i

    info = 0
    do i = 1, xnest%nblocks
      call psb_geinv(xnest%vects(i), ynest%vects(i), descs%descs(i,i), info)
      if (info /= 0) return
    end do
  end subroutine psb_d_nest_geinv

  !  y(i) = |x(i)|  element-wise for each block i
  subroutine psb_d_nest_geabs(xnest, ynest, descs, info)
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_d_nest_vect_type), intent(inout) :: ynest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info

    integer(psb_ipk_) :: i

    info = 0
    do i = 1, xnest%nblocks
      call psb_geabs(xnest%vects(i), ynest%vects(i), descs%descs(i,i), info)
      if (info /= 0) return
    end do
  end subroutine psb_d_nest_geabs

  !  z(i) = x(i) + b  for each block i  (b is a scalar)
  subroutine psb_d_nest_geaddconst(xnest, b, znest, descs, info)
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    real(psb_dpk_),             intent(in)    :: b
    type(psb_d_nest_vect_type), intent(inout) :: znest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info

    integer(psb_ipk_) :: i

    info = 0
    do i = 1, xnest%nblocks
      call psb_geaddconst(xnest%vects(i), b, znest%vects(i), descs%descs(i,i), info)
      if (info /= 0) return
    end do
  end subroutine psb_d_nest_geaddconst

  !  z(i) = cmp(x(i), c)  for each block i  (c is a scalar)
  subroutine psb_d_nest_gecmp(xnest, c, znest, descs, info)
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    real(psb_dpk_),             intent(in)    :: c
    type(psb_d_nest_vect_type), intent(inout) :: znest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info

    integer(psb_ipk_) :: i

    info = 0
    do i = 1, xnest%nblocks
      call psb_gecmp(xnest%vects(i), c, znest%vects(i), descs%descs(i,i), info)
      if (info /= 0) return
    end do
  end subroutine psb_d_nest_gecmp

  subroutine psb_d_nest_mask(cnest, xnest, mnest, t, descs, info)
    type(psb_d_nest_vect_type), intent(inout) :: cnest
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_d_nest_vect_type), intent(inout) :: mnest
    logical,                    intent(out)   :: t
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info

    integer(psb_ipk_) :: i
    logical           :: t_dummy
    real(psb_dpk_)    :: mmax

    info = 0
    do i = 1, cnest%nblocks
      ! psb_mask(c, x, m, t, desc) semantics after the double-swap in
      ! d_vect_mask_v → d_base_mask_v → d_base_mask_a:
      !   first arg  → constraint-type selector in d_base_mask_a
      !   second arg → value to test in d_base_mask_a
      ! So pass xnest (constraint types) first, cnest (values) second.
      call psb_mask(xnest%vects(i), cnest%vects(i), mnest%vects(i), &
           &        t_dummy, descs%descs(i,i), info)
      if (info /= 0) return
    end do
    ! t = .true. iff no block has any violated entry (mnest=1).
    mmax = psb_d_nest_geamax(mnest, descs, info)
    t = (mmax < 0.5_psb_dpk_)
  end subroutine psb_d_nest_mask
  

  !  Applies psb_upd_xyz(alpha,beta,gamma,delta,x(i),y(i),z(i))
  !  for each block i.
  subroutine psb_d_nest_upd_xyz(alpha, beta, gamma, delta, xnest, ynest, znest, descs, info)
    real(psb_dpk_),             intent(in)    :: alpha, beta, gamma, delta
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_d_nest_vect_type), intent(inout) :: ynest
    type(psb_d_nest_vect_type), intent(inout) :: znest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info

    integer(psb_ipk_) :: i

    info = 0
    do i = 1, xnest%nblocks
      call psb_upd_xyz(alpha, beta, gamma, delta, &
           &           xnest%vects(i), ynest%vects(i), znest%vects(i), &
           &           descs%descs(i,i), info)
      if (info /= 0) return
    end do
  end subroutine psb_d_nest_upd_xyz

  !  Block-diagonal triangular solve: applies psb_spsm to each
  !  diagonal block (i,i) of tnest independently.
  !  y(i) = alpha * T(i,i)^{-1} x(i) + beta * y(i)
  subroutine psb_d_nest_spsm(alpha, tnest, xnest, beta, ynest, descs, info, &
       &                     trans, scale, choice, work)
    real(psb_dpk_),               intent(in)              :: alpha, beta
    type(psb_d_nest_sparse_mat),  intent(inout)           :: tnest
    type(psb_d_nest_vect_type),   intent(inout)           :: xnest
    type(psb_d_nest_vect_type),   intent(inout)           :: ynest
    type(psb_desc_nest_type),     intent(in)              :: descs
    integer(psb_ipk_),            intent(out)             :: info
    character, optional,          intent(in)              :: trans, scale
    integer(psb_ipk_), optional,  intent(in)              :: choice
    real(psb_dpk_), optional,     intent(inout), target   :: work(:)

    integer(psb_ipk_) :: i

    info = 0
    do i = 1, tnest%nrblocks
      if (.not. tnest%has_block(i, i)) then
        ! No diagonal block: treat as identity => y(i) = alpha*x(i) + beta*y(i)
        call psb_geaxpby(alpha, xnest%vects(i), beta, ynest%vects(i), &
             &           descs%descs(i,i), info)
      else
        call psb_spsm(alpha, tnest%mats(i,i), xnest%vects(i), beta, ynest%vects(i), &
             &        descs%descs(i,i), info, trans=trans, scale=scale, &
             &        choice=choice, work=work)
      end if
      if (info /= 0) return
    end do
  end subroutine psb_d_nest_spsm

end module psb_d_nest_psblas_mod
