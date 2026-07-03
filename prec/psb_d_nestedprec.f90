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
!         notice, this list of conditions and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific prior written permission.
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
module psb_d_nestedprec

  use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_epk_, psb_dpk_, psb_success_, &
       & psb_err_invalid_input_, psb_err_invalid_preca_, psb_err_invalid_mat_state_, &
       & psb_err_alloc_dealloc_, psb_err_transpose_not_n_unsupported_, &
       & psb_root_, psb_toupper, psb_info, psb_errpush, psb_halo, psb_gedot, &
       & psb_spall, psb_spins, psb_spasb, psb_dupl_add_, done, dzero
  use psb_d_base_prec_mod
  use psb_d_nullprec, only : psb_d_null_prec_type
  use psb_d_diagprec, only : psb_d_diag_prec_type
  use psb_d_bjacprec, only : psb_d_bjac_prec_type
  use psb_d_nest_base_mat_mod, only : psb_d_nest_base_mat, psb_d_nest_get_n_fields, &
       & psb_d_nest_get_field_owned, psb_d_nest_get_block, psb_d_nest_get_field_desc, &
       & psb_d_nest_restrict_field, psb_d_nest_restrict_field_local, &
       & psb_d_nest_prolong_field, psb_d_nest_apply_block
  implicit none

  integer(psb_ipk_), parameter :: psb_d_nested_composition_ = 9101
  integer(psb_ipk_), parameter :: psb_d_nested_block_solve_ = 9102
  integer(psb_ipk_), parameter :: psb_d_nested_schur_solve_ = 9103
  integer(psb_ipk_), parameter :: psb_d_nested_schur_maxit_ = 9104
  integer(psb_ipk_), parameter :: psb_d_nested_schur_tol_ = 9105
  integer(psb_ipk_), parameter :: psb_d_nested_inner_solve_ = 9106
  integer(psb_ipk_), parameter :: psb_d_nested_inner_maxit_ = 9107
  integer(psb_ipk_), parameter :: psb_d_nested_inner_tol_ = 9108
  integer(psb_ipk_), parameter :: psb_d_nested_inner_itrace_ = 9109
  integer(psb_ipk_), parameter :: psb_d_nested_inner_istop_ = 9110

  integer(psb_ipk_), parameter, private :: psb_d_nested_diag_ = 1
  integer(psb_ipk_), parameter, private :: psb_d_nested_add_  = 2
  integer(psb_ipk_), parameter, private :: psb_d_nested_mult_ = 3
  integer(psb_ipk_), parameter, private :: psb_d_nested_symm_ = 4
  integer(psb_ipk_), parameter, private :: psb_d_nested_schur_lower_ = 5
  integer(psb_ipk_), parameter, private :: psb_d_nested_schur_upper_ = 6
  integer(psb_ipk_), parameter, private :: psb_d_nested_schur_full_  = 7
  integer(psb_ipk_), parameter, private :: psb_d_nested_schur_pde_control_ = 8
  integer(psb_ipk_), parameter, private :: psb_d_nested_schur_a22_ = 1
  integer(psb_ipk_), parameter, private :: psb_d_nested_schur_mf_  = 2
  integer(psb_ipk_), parameter, private :: psb_d_nested_schur_selfp_ = 3

  type :: psb_d_nested_iopt
    integer(psb_ipk_) :: field
    integer(psb_ipk_) :: what
    integer(psb_ipk_) :: val
  end type psb_d_nested_iopt

  type :: psb_d_nested_ropt
    integer(psb_ipk_) :: field
    integer(psb_ipk_) :: what
    real(psb_dpk_) :: val
  end type psb_d_nested_ropt

  type :: psb_d_nested_copt
    integer(psb_ipk_) :: field
    integer(psb_ipk_) :: what
    character(len=64) :: val
  end type psb_d_nested_copt

  type :: psb_d_nested_block_prec
    character(len=16) :: ptype
    class(psb_d_base_prec_type), allocatable :: pc
  end type psb_d_nested_block_prec

  type(psb_dspmat_type), pointer, save :: selfp_schur_mat => null()
  type(psb_d_nested_block_prec), pointer, save :: selfp_schur_block => null()

  type :: psb_d_nested_krylov_context
    logical :: enabled
    character(len=16) :: method
    integer(psb_ipk_) :: itmax
    integer(psb_ipk_) :: itrace
    integer(psb_ipk_) :: istop
    real(psb_dpk_) :: tol
  end type psb_d_nested_krylov_context

  type, extends(psb_d_base_prec_type) :: psb_d_nested_prec_type
    integer(psb_ipk_) :: composition
    character(len=32) :: composition_name
    character(len=16) :: default_block_ptype
    integer(psb_ipk_) :: schur_solve
    character(len=32) :: schur_solve_name
    integer(psb_ipk_) :: schur_maxit
    real(psb_dpk_) :: schur_tol
    type(psb_d_nested_krylov_context) :: default_krylov
    integer(psb_ipk_) :: nfields
    type(psb_d_nest_base_mat), pointer :: nest_op
    type(psb_d_nested_block_prec), allocatable :: blocks(:)
    logical :: schur_selfp_built
    character(len=16), allocatable :: field_block_ptype(:)
    type(psb_d_nested_krylov_context), allocatable :: field_krylov(:)
    type(psb_d_nested_iopt), allocatable :: field_iopts(:)
    type(psb_d_nested_ropt), allocatable :: field_ropts(:)
    type(psb_d_nested_copt), allocatable :: field_copts(:)
  contains
    procedure, pass(prec) :: d_apply_v => psb_d_nested_apply_vect
    procedure, pass(prec) :: d_apply   => psb_d_nested_apply
    procedure, pass(prec) :: precbld   => psb_d_nested_precbld
    procedure, pass(prec) :: precinit  => psb_d_nested_precinit
    procedure, pass(prec) :: precseti  => psb_d_nested_precseti
    procedure, pass(prec) :: precsetr  => psb_d_nested_precsetr
    procedure, pass(prec) :: precsetc  => psb_d_nested_precsetc
    procedure, pass(prec) :: precdescr => psb_d_nested_precdescr
    procedure, pass(prec) :: dump      => psb_d_nested_dump
    procedure, pass(prec) :: clone     => psb_d_nested_clone
    procedure, pass(prec) :: free      => psb_d_nested_precfree
    procedure, pass(prec) :: sizeof    => psb_d_nested_sizeof
    procedure, pass(prec) :: get_nzeros => psb_d_nested_get_nzeros
  end type psb_d_nested_prec_type

  private :: psb_d_nested_apply_vect, psb_d_nested_apply, psb_d_nested_precbld, &
       & psb_d_nested_precinit, psb_d_nested_precfree, psb_d_nested_precdescr, &
       & psb_d_nested_dump, psb_d_nested_clone, psb_d_nested_sizeof, &
       & psb_d_nested_get_nzeros, psb_d_nested_precseti, &
       & psb_d_nested_precsetr, psb_d_nested_precsetc, &
       & psb_d_nested_clear_built, psb_d_nested_valid_block_solve, &
       & psb_d_nested_get_field_block_ptype, psb_d_nested_field_solve, &
       & psb_d_nested_apply_schur, psb_d_nested_apply_pde_control_schur, &
       & psb_d_nested_replay_field_options, &
       & psb_d_nested_append_iopt, psb_d_nested_append_ropt, &
       & psb_d_nested_append_copt, psb_d_nested_schur_action, &
       & psb_d_nested_schur_solve, psb_d_nested_pde_control_schur_action, &
       & psb_d_nested_pde_control_schur_solve, psb_d_nested_get_field_krylov, &
       & psb_d_nested_inner_solve, psb_d_nested_inner_cg, psb_d_nested_inner_bicgstab, &
       & psb_d_nested_field_matvec, psb_d_nested_ensure_krylov_field, &
       & psb_d_nested_build_selfp, psb_d_nested_replay_schur_options

contains

  ! Reset nested preconditioner options and detach any previous field state.
  subroutine psb_d_nested_precinit(prec, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info

    info = psb_success_
    prec%composition = psb_d_nested_add_
    prec%composition_name = 'ADDITIVE'
    prec%default_block_ptype = 'DIAG'
    prec%schur_solve = psb_d_nested_schur_a22_
    prec%schur_solve_name = 'A22'
    prec%schur_maxit = 4
    prec%schur_tol = 0.0_psb_dpk_
    prec%schur_selfp_built = .false.
    prec%default_krylov%enabled = .false.
    prec%default_krylov%method = 'CG'
    prec%default_krylov%itmax = 20
    prec%default_krylov%itrace = -1
    prec%default_krylov%istop = 2
    prec%default_krylov%tol = 1.0e-6_psb_dpk_
    prec%nfields = 0
    prec%nest_op => null()
    if (allocated(prec%field_block_ptype)) deallocate(prec%field_block_ptype)
    if (allocated(prec%field_krylov)) deallocate(prec%field_krylov)
    if (allocated(prec%field_iopts)) deallocate(prec%field_iopts)
    if (allocated(prec%field_ropts)) deallocate(prec%field_ropts)
    if (allocated(prec%field_copts)) deallocate(prec%field_copts)
  end subroutine psb_d_nested_precinit

  ! Release built field preconditioners while preserving user configuration.
  subroutine psb_d_nested_clear_built(prec, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: i, local_info

    info = psb_success_
    if (allocated(prec%blocks)) then
      do i = 1, size(prec%blocks)
        if (allocated(prec%blocks(i)%pc)) then
          call prec%blocks(i)%pc%free(local_info)
          deallocate(prec%blocks(i)%pc, stat=local_info)
        end if
      end do
      deallocate(prec%blocks, stat=local_info)
      if (local_info /= 0 .and. info == psb_success_) info = local_info
    end if
    if (associated(selfp_schur_block)) then
      if (allocated(selfp_schur_block%pc)) then
        call selfp_schur_block%pc%free(local_info)
        if (local_info /= psb_success_ .and. info == psb_success_) info = local_info
        deallocate(selfp_schur_block%pc, stat=local_info)
        if (local_info /= 0 .and. info == psb_success_) info = local_info
      end if
      deallocate(selfp_schur_block, stat=local_info)
      if (local_info /= 0 .and. info == psb_success_) info = local_info
      nullify(selfp_schur_block)
    end if
    if (associated(selfp_schur_mat)) then
      call selfp_schur_mat%free()
      deallocate(selfp_schur_mat, stat=local_info)
      if (local_info /= 0 .and. info == psb_success_) info = local_info
      nullify(selfp_schur_mat)
    end if
    prec%schur_selfp_built = .false.
    prec%nfields = 0
    prec%nest_op => null()
  end subroutine psb_d_nested_clear_built

  ! Fully release nested preconditioner storage and pending configuration.
  subroutine psb_d_nested_precfree(prec, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: local_info

    call psb_d_nested_clear_built(prec, info)
    if (allocated(prec%field_block_ptype)) then
      deallocate(prec%field_block_ptype, stat=local_info)
      if (local_info /= 0 .and. info == psb_success_) info = local_info
    end if
    if (allocated(prec%field_krylov)) then
      deallocate(prec%field_krylov, stat=local_info)
      if (local_info /= 0 .and. info == psb_success_) info = local_info
    end if
    if (allocated(prec%field_iopts)) deallocate(prec%field_iopts, stat=local_info)
    if (allocated(prec%field_ropts)) deallocate(prec%field_ropts, stat=local_info)
    if (allocated(prec%field_copts)) deallocate(prec%field_copts, stat=local_info)
  end subroutine psb_d_nested_precfree

  ! Set integer-valued nested options or record integer field-block options.
  subroutine psb_d_nested_precseti(prec, what, val, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: what, val
    integer(psb_ipk_), intent(out) :: info

    info = psb_success_
    select case (what)
    case (psb_d_nested_schur_maxit_)
      prec%schur_maxit = max(0, val)
    case (psb_d_nested_inner_maxit_, psb_d_nested_inner_itrace_, &
         & psb_d_nested_inner_istop_)
      call psb_d_nested_set_field_krylov_i(prec, 0, what, val, info)
    case default
      call psb_d_nested_append_iopt(prec, 0, what, val, info)
    end select
  end subroutine psb_d_nested_precseti

  ! Set real-valued nested options or record real field-block options.
  subroutine psb_d_nested_precsetr(prec, what, val, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: what
    real(psb_dpk_), intent(in) :: val
    integer(psb_ipk_), intent(out) :: info

    info = psb_success_
    select case (what)
    case (psb_d_nested_schur_tol_)
      prec%schur_tol = max(dzero, val)
    case (psb_d_nested_inner_tol_)
      call psb_d_nested_set_field_krylov_r(prec, 0, what, val, info)
    case default
      call psb_d_nested_append_ropt(prec, 0, what, val, info)
    end select
  end subroutine psb_d_nested_precsetr

  ! Set character-valued nested options or record character field-block options.
  subroutine psb_d_nested_precsetc(prec, what, val, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: what
    character(len=*), intent(in) :: val
    integer(psb_ipk_), intent(out) :: info

    character(len=64) :: opt

    info = psb_success_
    opt = psb_toupper(trim(val))
    select case (what)
    case (psb_d_nested_composition_)
      select case (trim(opt))
      case ('DIAGONAL','DIAG')
        prec%composition = psb_d_nested_diag_
        prec%composition_name = 'DIAGONAL'
      case ('ADDITIVE','ADD')
        prec%composition = psb_d_nested_add_
        prec%composition_name = 'ADDITIVE'
      case ('MULTIPLICATIVE','MULT')
        prec%composition = psb_d_nested_mult_
        prec%composition_name = 'MULTIPLICATIVE'
      case ('SYMMETRIC_MULTIPLICATIVE','SYMMETRIC','SYM_MULT','SYMM')
        prec%composition = psb_d_nested_symm_
        prec%composition_name = 'SYMMETRIC_MULTIPLICATIVE'
      case ('SCHUR_LOWER','LOWER_SCHUR')
        prec%composition = psb_d_nested_schur_lower_
        prec%composition_name = 'SCHUR_LOWER'
      case ('SCHUR_UPPER','UPPER_SCHUR')
        prec%composition = psb_d_nested_schur_upper_
        prec%composition_name = 'SCHUR_UPPER'
      case ('SCHUR','SCHUR_FULL','FULL_SCHUR')
        prec%composition = psb_d_nested_schur_full_
        prec%composition_name = 'SCHUR_FULL'
      case ('SCHUR_PDE_CONTROL','PDE_CONTROL_SCHUR','SCHUR_2PLUS1', &
           & 'SCHUR3','SCHUR_3FIELD','SCHUR_THREE_FIELD')
        prec%composition = psb_d_nested_schur_pde_control_
        prec%composition_name = 'SCHUR_PDE_CONTROL'
        prec%schur_solve = psb_d_nested_schur_mf_
        prec%schur_solve_name = 'MATRIX_FREE'
      case default
        info = psb_err_invalid_input_
      end select
    case (psb_d_nested_block_solve_)
      if (psb_d_nested_valid_block_solve(opt)) then
        prec%default_block_ptype = trim(opt)
      else
        info = psb_err_invalid_preca_
      end if
    case (psb_d_nested_schur_solve_)
      select case (trim(opt))
      case ('A22','A_22','FIELD','FIELD_BLOCK')
        prec%schur_solve = psb_d_nested_schur_a22_
        prec%schur_solve_name = 'A22'
      case ('A33','A_33','FIELD3','FIELD_3','FIELD_BLOCK3','FIELD_BLOCK_3')
        prec%schur_solve = psb_d_nested_schur_a22_
        prec%schur_solve_name = 'A33'
      case ('MATRIX_FREE','MATFREE','MF')
        prec%schur_solve = psb_d_nested_schur_mf_
        prec%schur_solve_name = 'MATRIX_FREE'
      case ('SELFP','SELF')
        prec%schur_solve = psb_d_nested_schur_selfp_
        prec%schur_solve_name = 'SELFP'
      case default
        info = psb_err_invalid_input_
      end select
    case (psb_d_nested_inner_solve_)
      call psb_d_nested_set_field_krylov_c(prec, 0, what, trim(opt), info)
    case default
      call psb_d_nested_append_copt(prec, 0, what, trim(val), info)
    end select
  end subroutine psb_d_nested_precsetc

  ! Build one field preconditioner for each diagonal block of a nested matrix.
  subroutine psb_d_nested_precbld(a, desc_a, prec, info, amold, vmold, imold)
    type(psb_dspmat_type), intent(in), target :: a
    type(psb_desc_type), intent(inout), target :: desc_a
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info
    class(psb_d_base_sparse_mat), intent(in), optional :: amold
    class(psb_d_base_vect_type), intent(in), optional :: vmold
    class(psb_i_base_vect_type), intent(in), optional :: imold

    integer(psb_ipk_) :: i, nfields
    type(psb_dspmat_type), pointer :: block_ptr
    type(psb_desc_type), pointer :: field_desc
    character(len=24) :: name

    info = psb_success_
    name = 'd_nested_precbld'
    call prec%set_ctxt(desc_a%get_ctxt())

    call psb_d_nested_clear_built(prec, info)
    if (info /= psb_success_) return
    call prec%set_ctxt(desc_a%get_ctxt())

    if (.not. allocated(a%a)) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, name, a_err='matrix storage not allocated')
      return
    end if

    select type (ap => a%a)
    type is (psb_d_nest_base_mat)
      prec%nest_op => ap
    class default
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='NEST preconditioner requires NEST matrix')
      return
    end select

    nfields = psb_d_nest_get_n_fields(prec%nest_op)
    if (nfields <= 0) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, name, a_err='nested matrix has no fields')
      return
    end if

    prec%nfields = nfields
    allocate(prec%blocks(nfields), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name, a_err='blocks')
      return
    end if

    do i = 1, nfields
      field_desc => psb_d_nest_get_field_desc(prec%nest_op, i)
      if (.not. associated(field_desc)) then
        info = psb_err_invalid_mat_state_
        call psb_errpush(info, name, a_err='missing field descriptor')
        return
      end if

      block_ptr => psb_d_nest_get_block(prec%nest_op, i, i)
      if (associated(block_ptr)) then
        prec%blocks(i)%ptype = psb_d_nested_get_field_block_ptype(prec, i)
      else
        prec%blocks(i)%ptype = 'NONE'
      end if
      call psb_d_nested_alloc_block_pc(prec%blocks(i), info)
      if (info /= psb_success_) return
      call prec%blocks(i)%pc%precinit(info)
      if (info /= psb_success_) return
      call psb_d_nested_replay_field_options(prec, i, info)
      if (info /= psb_success_) return

      if (associated(block_ptr)) then
        call prec%blocks(i)%pc%precbld(block_ptr, field_desc, info, &
             & amold=amold, vmold=vmold, imold=imold)
      else
        call prec%blocks(i)%pc%set_ctxt(field_desc%get_ctxt())
      end if
      if (info /= psb_success_) return
    end do

    if ((prec%schur_solve == psb_d_nested_schur_selfp_) .and. &
         & (prec%composition == psb_d_nested_schur_lower_ .or. &
         &  prec%composition == psb_d_nested_schur_upper_ .or. &
         &  prec%composition == psb_d_nested_schur_full_)) then
      call psb_d_nested_build_selfp(prec, info, amold, vmold, imold)
      if (info /= psb_success_) return
    end if
  end subroutine psb_d_nested_precbld

  ! Allocate the concrete PSBLAS preconditioner requested for one field block.
  subroutine psb_d_nested_alloc_block_pc(block, info)
    type(psb_d_nested_block_prec), intent(inout) :: block
    integer(psb_ipk_), intent(out) :: info

    info = psb_success_
    if (allocated(block%pc)) then
      call block%pc%free(info)
      if (info == psb_success_) deallocate(block%pc, stat=info)
      if (info /= psb_success_) return
    end if
    select case (psb_toupper(trim(block%ptype)))
    case ('NONE','NOPREC')
      allocate(psb_d_null_prec_type :: block%pc, stat=info)
    case ('DIAG')
      allocate(psb_d_diag_prec_type :: block%pc, stat=info)
    case ('BJAC')
      allocate(psb_d_bjac_prec_type :: block%pc, stat=info)
    case default
      info = psb_err_invalid_preca_
    end select
  end subroutine psb_d_nested_alloc_block_pc

  logical function psb_d_nested_valid_block_solve(ptype) result(valid)
    character(len=*), intent(in) :: ptype
    character(len=32) :: opt

    opt = psb_toupper(trim(ptype))
    select case (trim(opt))
    case ('NONE','NOPREC','DIAG','BJAC')
      valid = .true.
    case default
      valid = .false.
    end select
  end function psb_d_nested_valid_block_solve

  function psb_d_nested_get_field_block_ptype(prec, field) result(ptype)
    class(psb_d_nested_prec_type), intent(in) :: prec
    integer(psb_ipk_), intent(in) :: field
    character(len=16) :: ptype

    ptype = prec%default_block_ptype
    if (allocated(prec%field_block_ptype)) then
      if (field <= size(prec%field_block_ptype)) then
        if (len_trim(prec%field_block_ptype(field)) > 0) &
             & ptype = prec%field_block_ptype(field)
      end if
    end if
  end function psb_d_nested_get_field_block_ptype

  ! Store a per-field block preconditioner type override.
  subroutine psb_d_nested_set_block_solve_field(prec, field, ptype, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field
    character(len=*), intent(in) :: ptype
    integer(psb_ipk_), intent(out) :: info

    character(len=16), allocatable :: tmp(:)
    character(len=32) :: opt
    integer(psb_ipk_) :: old_size

    info = psb_success_
    opt = psb_toupper(trim(ptype))
    if (field <= 0) then
      info = psb_err_invalid_input_
      call psb_errpush(info, 'd_nested_set_block_solve_field', a_err='field index')
      return
    end if
    if (.not. psb_d_nested_valid_block_solve(opt)) then
      info = psb_err_invalid_preca_
      call psb_errpush(info, 'd_nested_set_block_solve_field', a_err='block solve')
      return
    end if

    old_size = 0
    if (allocated(prec%field_block_ptype)) old_size = size(prec%field_block_ptype)
    if (field > old_size) then
      allocate(tmp(field), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_
        call psb_errpush(info, 'd_nested_set_block_solve_field', a_err='field block options')
        return
      end if
      tmp(:) = ''
      if (old_size > 0) tmp(1:old_size) = prec%field_block_ptype(1:old_size)
      call move_alloc(tmp, prec%field_block_ptype)
    end if
    prec%field_block_ptype(field) = trim(opt)
  end subroutine psb_d_nested_set_block_solve_field

  subroutine psb_d_nested_ensure_krylov_field(prec, field, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field
    integer(psb_ipk_), intent(out) :: info

    type(psb_d_nested_krylov_context), allocatable :: tmp(:)
    integer(psb_ipk_) :: old_size

    info = psb_success_
    if (field <= 0) return

    old_size = 0
    if (allocated(prec%field_krylov)) old_size = size(prec%field_krylov)
    if (field > old_size) then
      allocate(tmp(field), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_
        call psb_errpush(info, 'd_nested_ensure_krylov_field', a_err='field krylov options')
        return
      end if
      tmp(:) = prec%default_krylov
      if (old_size > 0) tmp(1:old_size) = prec%field_krylov(1:old_size)
      call move_alloc(tmp, prec%field_krylov)
    end if
  end subroutine psb_d_nested_ensure_krylov_field

  subroutine psb_d_nested_set_field_krylov_c(prec, field, what, val, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field, what
    character(len=*), intent(in) :: val
    integer(psb_ipk_), intent(out) :: info

    character(len=32) :: opt

    info = psb_success_
    opt = psb_toupper(trim(val))
    if (field < 0) then
      info = psb_err_invalid_input_
      call psb_errpush(info, 'd_nested_set_field_krylov_c', a_err='field index')
      return
    end if
    if (what /= psb_d_nested_inner_solve_) then
      info = psb_err_invalid_input_
      call psb_errpush(info, 'd_nested_set_field_krylov_c', a_err='inner option')
      return
    end if

    select case (trim(opt))
    case ('NONE','NO','OFF','FALSE')
      if (field == 0) then
        prec%default_krylov%enabled = .false.
      else
        call psb_d_nested_ensure_krylov_field(prec, field, info)
        if (info == psb_success_) prec%field_krylov(field)%enabled = .false.
      end if
    case ('CG','BICGSTAB','BICGstab','BiCGSTAB')
      if (field == 0) then
        prec%default_krylov%enabled = .true.
        prec%default_krylov%method = trim(opt)
      else
        call psb_d_nested_ensure_krylov_field(prec, field, info)
        if (info == psb_success_) then
          prec%field_krylov(field)%enabled = .true.
          prec%field_krylov(field)%method = trim(opt)
        end if
      end if
    case default
      info = psb_err_invalid_input_
      call psb_errpush(info, 'd_nested_set_field_krylov_c', a_err='inner method')
    end select
  end subroutine psb_d_nested_set_field_krylov_c

  subroutine psb_d_nested_set_field_krylov_i(prec, field, what, val, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field, what, val
    integer(psb_ipk_), intent(out) :: info

    info = psb_success_
    if (field < 0) then
      info = psb_err_invalid_input_
      call psb_errpush(info, 'd_nested_set_field_krylov_i', a_err='field index')
      return
    end if

    if (field > 0) then
      call psb_d_nested_ensure_krylov_field(prec, field, info)
      if (info /= psb_success_) return
    end if

    select case (what)
    case (psb_d_nested_inner_maxit_)
      if (field == 0) then
        prec%default_krylov%itmax = max(1, val)
      else
        prec%field_krylov(field)%itmax = max(1, val)
      end if
    case (psb_d_nested_inner_itrace_)
      if (field == 0) then
        prec%default_krylov%itrace = val
      else
        prec%field_krylov(field)%itrace = val
      end if
    case (psb_d_nested_inner_istop_)
      if (field == 0) then
        prec%default_krylov%istop = val
      else
        prec%field_krylov(field)%istop = val
      end if
    case default
      info = psb_err_invalid_input_
      call psb_errpush(info, 'd_nested_set_field_krylov_i', a_err='inner option')
    end select
  end subroutine psb_d_nested_set_field_krylov_i

  subroutine psb_d_nested_set_field_krylov_r(prec, field, what, val, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field, what
    real(psb_dpk_), intent(in) :: val
    integer(psb_ipk_), intent(out) :: info

    info = psb_success_
    if (field < 0) then
      info = psb_err_invalid_input_
      call psb_errpush(info, 'd_nested_set_field_krylov_r', a_err='field index')
      return
    end if
    if (what /= psb_d_nested_inner_tol_) then
      info = psb_err_invalid_input_
      call psb_errpush(info, 'd_nested_set_field_krylov_r', a_err='inner option')
      return
    end if

    if (field == 0) then
      prec%default_krylov%tol = max(dzero, val)
    else
      call psb_d_nested_ensure_krylov_field(prec, field, info)
      if (info == psb_success_) prec%field_krylov(field)%tol = max(dzero, val)
    end if
  end subroutine psb_d_nested_set_field_krylov_r

  function psb_d_nested_get_field_krylov(prec, field) result(ctx)
    class(psb_d_nested_prec_type), intent(in) :: prec
    integer(psb_ipk_), intent(in) :: field
    type(psb_d_nested_krylov_context) :: ctx

    ctx = prec%default_krylov
    if (allocated(prec%field_krylov)) then
      if (field <= size(prec%field_krylov)) ctx = prec%field_krylov(field)
    end if
  end function psb_d_nested_get_field_krylov

  ! Append a pending integer option for all fields or one selected field.
  subroutine psb_d_nested_append_iopt(prec, field, what, val, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field, what, val
    integer(psb_ipk_), intent(out) :: info

    type(psb_d_nested_iopt), allocatable :: tmp(:)
    integer(psb_ipk_) :: n

    info = psb_success_
    n = 0
    if (allocated(prec%field_iopts)) n = size(prec%field_iopts)
    allocate(tmp(n+1), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_append_iopt', a_err='field integer options')
      return
    end if
    if (n > 0) tmp(1:n) = prec%field_iopts(1:n)
    tmp(n+1)%field = field
    tmp(n+1)%what = what
    tmp(n+1)%val = val
    call move_alloc(tmp, prec%field_iopts)
  end subroutine psb_d_nested_append_iopt

  ! Append a pending real option for all fields or one selected field.
  subroutine psb_d_nested_append_ropt(prec, field, what, val, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field, what
    real(psb_dpk_), intent(in) :: val
    integer(psb_ipk_), intent(out) :: info

    type(psb_d_nested_ropt), allocatable :: tmp(:)
    integer(psb_ipk_) :: n

    info = psb_success_
    n = 0
    if (allocated(prec%field_ropts)) n = size(prec%field_ropts)
    allocate(tmp(n+1), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_append_ropt', a_err='field real options')
      return
    end if
    if (n > 0) tmp(1:n) = prec%field_ropts(1:n)
    tmp(n+1)%field = field
    tmp(n+1)%what = what
    tmp(n+1)%val = val
    call move_alloc(tmp, prec%field_ropts)
  end subroutine psb_d_nested_append_ropt

  ! Append a pending character option for all fields or one selected field.
  subroutine psb_d_nested_append_copt(prec, field, what, val, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field, what
    character(len=*), intent(in) :: val
    integer(psb_ipk_), intent(out) :: info

    type(psb_d_nested_copt), allocatable :: tmp(:)
    integer(psb_ipk_) :: n

    info = psb_success_
    n = 0
    if (allocated(prec%field_copts)) n = size(prec%field_copts)
    allocate(tmp(n+1), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_append_copt', a_err='field character options')
      return
    end if
    if (n > 0) tmp(1:n) = prec%field_copts(1:n)
    tmp(n+1)%field = field
    tmp(n+1)%what = what
    tmp(n+1)%val = trim(val)
    call move_alloc(tmp, prec%field_copts)
  end subroutine psb_d_nested_append_copt

  ! Record an integer option to replay on a field block preconditioner.
  subroutine psb_d_nested_set_field_precseti(prec, field, what, val, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field, what, val
    integer(psb_ipk_), intent(out) :: info

    info = psb_success_
    if (field < 0) then
      info = psb_err_invalid_input_
      call psb_errpush(info, 'd_nested_set_field_precseti', a_err='field index')
      return
    end if
    call psb_d_nested_append_iopt(prec, field, what, val, info)
  end subroutine psb_d_nested_set_field_precseti

  ! Record a real option to replay on a field block preconditioner.
  subroutine psb_d_nested_set_field_precsetr(prec, field, what, val, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field, what
    real(psb_dpk_), intent(in) :: val
    integer(psb_ipk_), intent(out) :: info

    info = psb_success_
    if (field < 0) then
      info = psb_err_invalid_input_
      call psb_errpush(info, 'd_nested_set_field_precsetr', a_err='field index')
      return
    end if
    call psb_d_nested_append_ropt(prec, field, what, val, info)
  end subroutine psb_d_nested_set_field_precsetr

  ! Record a character option to replay on a field block preconditioner.
  subroutine psb_d_nested_set_field_precsetc(prec, field, what, val, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field, what
    character(len=*), intent(in) :: val
    integer(psb_ipk_), intent(out) :: info

    info = psb_success_
    if (field < 0) then
      info = psb_err_invalid_input_
      call psb_errpush(info, 'd_nested_set_field_precsetc', a_err='field index')
      return
    end if
    call psb_d_nested_append_copt(prec, field, what, val, info)
  end subroutine psb_d_nested_set_field_precsetc

  ! Apply stored field-specific options to a built field block preconditioner.
  subroutine psb_d_nested_replay_field_options(prec, field, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: k

    info = psb_success_
    if (.not. allocated(prec%blocks(field)%pc)) return

    if (allocated(prec%field_iopts)) then
      do k = 1, size(prec%field_iopts)
        if ((prec%field_iopts(k)%field == 0) .or. (prec%field_iopts(k)%field == field)) then
          call prec%blocks(field)%pc%precset(prec%field_iopts(k)%what, &
               & prec%field_iopts(k)%val, info)
          if (info /= psb_success_) return
        end if
      end do
    end if

    if (allocated(prec%field_ropts)) then
      do k = 1, size(prec%field_ropts)
        if ((prec%field_ropts(k)%field == 0) .or. (prec%field_ropts(k)%field == field)) then
          call prec%blocks(field)%pc%precset(prec%field_ropts(k)%what, &
               & prec%field_ropts(k)%val, info)
          if (info /= psb_success_) return
        end if
      end do
    end if

    if (allocated(prec%field_copts)) then
      do k = 1, size(prec%field_copts)
        if ((prec%field_copts(k)%field == 0) .or. (prec%field_copts(k)%field == field)) then
          call prec%blocks(field)%pc%precset(prec%field_copts(k)%what, &
               & trim(prec%field_copts(k)%val), info)
          if (info /= psb_success_) return
        end if
      end do
    end if
  end subroutine psb_d_nested_replay_field_options

  ! Apply the nested preconditioner to PSBLAS vector objects.
  subroutine psb_d_nested_apply_vect(alpha, prec, x, beta, y, desc_data, info, trans)
    type(psb_desc_type), intent(in) :: desc_data
    class(psb_d_nested_prec_type), intent(inout) :: prec
    type(psb_d_vect_type), intent(inout) :: x
    real(psb_dpk_), intent(in) :: alpha, beta
    type(psb_d_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out) :: info
    character(len=1), optional :: trans

    real(psb_dpk_), allocatable :: xbuf(:), ybuf(:)
    integer(psb_ipk_) :: ncol

    info = psb_success_
    ncol = desc_data%get_local_cols()
    xbuf = x%get_vect(ncol)
    allocate(ybuf(max(y%get_nrows(), ncol)), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      return
    end if
    ybuf(:) = dzero
    call psb_d_nested_apply(alpha, prec, xbuf, beta, ybuf, desc_data, info, trans)
    if (info == psb_success_) call y%set(ybuf(1:y%get_nrows()))
    deallocate(ybuf)
  end subroutine psb_d_nested_apply_vect

  ! Apply the selected nested composition to raw array vectors.
  subroutine psb_d_nested_apply(alpha, prec, x, beta, y, desc_data, info, trans, work)
    type(psb_desc_type), intent(in) :: desc_data
    class(psb_d_nested_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(inout) :: x(:)
    real(psb_dpk_), intent(in) :: alpha, beta
    real(psb_dpk_), intent(inout) :: y(:)
    integer(psb_ipk_), intent(out) :: info
    character(len=1), optional :: trans
    real(psb_dpk_), intent(inout), optional, target :: work(:)

    real(psb_dpk_), allocatable :: z(:)
    integer(psb_ipk_) :: nrow, ncol
    character :: trans_
    character(len=24) :: name

    info = psb_success_
    name = 'd_nested_apply'
    trans_ = 'N'
    if (present(trans)) trans_ = psb_toupper(trans)
    if (.not. associated(prec%nest_op) .or. .not. allocated(prec%blocks)) then
      info = 1124
      call psb_errpush(info, name, a_err='nested preconditioner')
      return
    end if

    nrow = desc_data%get_local_rows()
    ncol = desc_data%get_local_cols()
    if (size(x) < nrow .or. size(y) < nrow) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='vector too small')
      return
    end if
    allocate(z(max(nrow,ncol)), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name, a_err='work vector')
      return
    end if
    z(:) = dzero

    if (trans_ /= 'N') then
      info = psb_err_transpose_not_n_unsupported_
      call psb_errpush(info, name)
    else
      select case (prec%composition)
      case (psb_d_nested_diag_, psb_d_nested_add_)
        call psb_d_nested_add_apply(prec, x, z, desc_data, info)
      case (psb_d_nested_mult_)
        call psb_d_nested_sweep(prec, x, z, desc_data, 1, prec%nfields, 1, info)
      case (psb_d_nested_symm_)
        call psb_d_nested_sweep(prec, x, z, desc_data, 1, prec%nfields, 1, info)
        if (info == psb_success_) &
             & call psb_d_nested_sweep(prec, x, z, desc_data, prec%nfields, 1, -1, info)
      case (psb_d_nested_schur_lower_, psb_d_nested_schur_upper_, psb_d_nested_schur_full_)
        call psb_d_nested_apply_schur(prec, x, z, desc_data, info)
      case (psb_d_nested_schur_pde_control_)
        call psb_d_nested_apply_pde_control_schur(prec, x, z, desc_data, info)
      case default
        info = psb_err_invalid_input_
        call psb_errpush(info, name, a_err='unknown composition')
      end select
    end if
    if (info /= psb_success_) then
      deallocate(z)
      return
    end if

    if (beta == dzero) then
      y(1:nrow) = alpha * z(1:nrow)
    else
      y(1:nrow) = alpha * z(1:nrow) + beta * y(1:nrow)
    end if
    deallocate(z)
  end subroutine psb_d_nested_apply

  ! Apply the additive or diagonal field-split composition.
  subroutine psb_d_nested_add_apply(prec, x, z, desc_data, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(in) :: x(:)
    real(psb_dpk_), intent(inout) :: z(:)
    type(psb_desc_type), intent(in) :: desc_data
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: i, n_owned, n_col_field
    real(psb_dpk_), allocatable :: rhs(:), sol(:)
    type(psb_desc_type), pointer :: field_desc

    info = psb_success_
    z(:) = dzero

    do i = 1, prec%nfields
      n_owned = psb_d_nest_get_field_owned(prec%nest_op, i)
      field_desc => psb_d_nest_get_field_desc(prec%nest_op, i)

      if (.not. associated(field_desc)) then
        info = psb_err_invalid_mat_state_
        call psb_errpush(info, 'd_nested_add_apply', a_err='missing field descriptor')
        return
      end if

      n_col_field = field_desc%get_local_cols()

      allocate(rhs(n_col_field), sol(n_col_field), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_
        call psb_errpush(info, 'd_nested_add_apply', a_err='field vectors')
        return
      end if

      rhs(:) = dzero
      sol(:) = dzero

      call psb_d_nest_restrict_field(prec%nest_op, i, x, rhs, info)

      if (info == psb_success_) then
        call psb_d_nested_field_solve(prec, i, rhs, sol, info)
      end if

      if (info == psb_success_) then
        call psb_d_nest_prolong_field(prec%nest_op, i, sol, z, info)
      end if

      deallocate(rhs, sol)

      if (info /= psb_success_) return
    end do

  end subroutine psb_d_nested_add_apply

  ! Apply one forward or backward multiplicative field sweep.
  subroutine psb_d_nested_sweep(prec, x, z, desc_data, first, last, step, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(in) :: x(:)
    real(psb_dpk_), intent(inout) :: z(:)
    type(psb_desc_type), intent(in) :: desc_data
    integer(psb_ipk_), intent(in) :: first, last, step
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: i, j, n_owned, n_col_i, n_col_j
    real(psb_dpk_), allocatable :: rhs(:), sol(:), z_field(:)
    type(psb_desc_type), pointer :: field_desc_i, field_desc_j
    type(psb_dspmat_type), pointer :: block_ptr

    info = psb_success_

    call psb_halo(z, desc_data, info)
    if (info /= psb_success_) return

    i = first
    do
      n_owned = psb_d_nest_get_field_owned(prec%nest_op, i)
      field_desc_i => psb_d_nest_get_field_desc(prec%nest_op, i)

      if (.not. associated(field_desc_i)) then
        info = psb_err_invalid_mat_state_
        call psb_errpush(info, 'd_nested_sweep', a_err='missing row field descriptor')
        return
      end if

      n_col_i = field_desc_i%get_local_cols()

      allocate(rhs(n_col_i), sol(n_col_i), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_
        call psb_errpush(info, 'd_nested_sweep', a_err='field vectors')
        return
      end if

      rhs(:) = dzero
      sol(:) = dzero

      call psb_d_nest_restrict_field(prec%nest_op, i, x, rhs, info)
      if (info /= psb_success_) then
        deallocate(rhs, sol)
        return
      end if

      do j = 1, prec%nfields
        if (j == i) cycle

        block_ptr => psb_d_nest_get_block(prec%nest_op, i, j)
        if (.not. associated(block_ptr)) cycle

        field_desc_j => psb_d_nest_get_field_desc(prec%nest_op, j)
        if (.not. associated(field_desc_j)) then
          info = psb_err_invalid_mat_state_
          call psb_errpush(info, 'd_nested_sweep', a_err='missing column field descriptor')
          deallocate(rhs, sol)
          return
        end if

        n_col_j = field_desc_j%get_local_cols()

        allocate(z_field(n_col_j), stat=info)
        if (info /= 0) then
          info = psb_err_alloc_dealloc_
          call psb_errpush(info, 'd_nested_sweep', a_err='offdiag field vector')
          deallocate(rhs, sol)
          return
        end if

        z_field(:) = dzero

        call psb_d_nest_restrict_field_local(prec%nest_op, j, z, z_field, info)

        if (info == psb_success_) then
          call psb_d_nest_apply_block(prec%nest_op, i, j, -done, z_field, done, rhs, info)
        end if

        deallocate(z_field)

        if (info /= psb_success_) then
          deallocate(rhs, sol)
          return
        end if
      end do

      call psb_d_nested_field_solve(prec, i, rhs, sol, info)

      if (info == psb_success_) then
        call psb_d_nest_prolong_field(prec%nest_op, i, sol, z, info)
      end if

      deallocate(rhs, sol)

      if (info /= psb_success_) return

      call psb_halo(z, desc_data, info)
      if (info /= psb_success_) return

      if (i == last) exit
      i = i + step
    end do
  end subroutine psb_d_nested_sweep

  ! Apply one field block preconditioner to a field-local right-hand side.
  subroutine psb_d_nested_field_solve(prec, field, rhs, sol, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field
    real(psb_dpk_), intent(inout) :: rhs(:)
    real(psb_dpk_), intent(inout) :: sol(:)
    integer(psb_ipk_), intent(out) :: info

    real(psb_dpk_), allocatable :: wrk(:)
    type(psb_desc_type), pointer :: field_desc
    type(psb_d_nested_krylov_context) :: kctx
    integer(psb_ipk_) :: n_col_field

    info = psb_success_
    field_desc => psb_d_nest_get_field_desc(prec%nest_op, field)
    if (.not. associated(field_desc)) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'd_nested_field_solve', a_err='missing field descriptor')
      return
    end if

    n_col_field = field_desc%get_local_cols()
    kctx = prec%default_krylov
    if (allocated(prec%field_krylov)) then
      if (field <= size(prec%field_krylov)) kctx = prec%field_krylov(field)
    end if
    if (kctx%enabled) then
      call psb_d_nested_inner_solve(prec, field, kctx, rhs, sol, info)
      return
    end if

    allocate(wrk(n_col_field), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_field_solve', a_err='work vector')
      return
    end if

    sol(:) = dzero
    call prec%blocks(field)%pc%apply(done, rhs, dzero, sol, field_desc, info, trans='N', work=wrk)
    deallocate(wrk)
  end subroutine psb_d_nested_field_solve

  subroutine psb_d_nested_field_matvec(prec, field, x, y, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field
    real(psb_dpk_), intent(inout) :: x(:)
    real(psb_dpk_), intent(inout) :: y(:)
    integer(psb_ipk_), intent(out) :: info

    real(psb_dpk_), allocatable :: xh(:)
    type(psb_desc_type), pointer :: field_desc
    type(psb_dspmat_type), pointer :: block_ptr
    integer(psb_ipk_) :: n_col_field

    info = psb_success_
    field_desc => psb_d_nest_get_field_desc(prec%nest_op, field)
    block_ptr => psb_d_nest_get_block(prec%nest_op, field, field)
    if (.not. associated(field_desc)) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'd_nested_field_matvec', a_err='missing field descriptor')
      return
    end if
    if (.not. associated(block_ptr)) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'd_nested_field_matvec', a_err='missing diagonal block')
      return
    end if

    n_col_field = field_desc%get_local_cols()
    allocate(xh(n_col_field), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_field_matvec', a_err='work vector')
      return
    end if
    xh(:) = dzero
    xh(1:min(size(x), n_col_field)) = x(1:min(size(x), n_col_field))
    call psb_halo(xh, field_desc, info)
    if (info == psb_success_) then
      y(:) = dzero
      call psb_d_nest_apply_block(prec%nest_op, field, field, done, xh, dzero, y, info)
    end if
    deallocate(xh)
  end subroutine psb_d_nested_field_matvec

  subroutine psb_d_nested_inner_solve(prec, field, kctx, rhs, sol, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field
    type(psb_d_nested_krylov_context), intent(in) :: kctx
    real(psb_dpk_), intent(inout) :: rhs(:)
    real(psb_dpk_), intent(inout) :: sol(:)
    integer(psb_ipk_), intent(out) :: info

    select case (psb_toupper(trim(kctx%method)))
    case ('CG')
      call psb_d_nested_inner_cg(prec, field, kctx, rhs, sol, info)
    case ('BICGSTAB')
      call psb_d_nested_inner_bicgstab(prec, field, kctx, rhs, sol, info)
    case default
      info = psb_err_invalid_input_
      call psb_errpush(info, 'd_nested_inner_solve', a_err='inner Krylov method')
    end select
  end subroutine psb_d_nested_inner_solve

  subroutine psb_d_nested_inner_cg(prec, field, kctx, rhs, sol, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field
    type(psb_d_nested_krylov_context), intent(in) :: kctx
    real(psb_dpk_), intent(inout) :: rhs(:)
    real(psb_dpk_), intent(inout) :: sol(:)
    integer(psb_ipk_), intent(out) :: info

    real(psb_dpk_), allocatable :: r(:), z(:), p(:), q(:), wrk(:)
    type(psb_desc_type), pointer :: field_desc
    integer(psb_ipk_) :: n_col_field, k
    real(psb_dpk_) :: bnorm, rnorm, rz, rz_new, denom, alpha, beta

    info = psb_success_
    field_desc => psb_d_nest_get_field_desc(prec%nest_op, field)
    if (.not. associated(field_desc)) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'd_nested_inner_cg', a_err='missing field descriptor')
      return
    end if
    n_col_field = field_desc%get_local_cols()
    allocate(r(n_col_field), z(n_col_field), p(n_col_field), q(n_col_field), &
         & wrk(n_col_field), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_inner_cg', a_err='inner work vectors')
      return
    end if

    sol(:) = dzero
    r(:) = rhs(:)
    bnorm = sqrt(max(dzero, psb_gedot(rhs, rhs, field_desc, info)))
    if (info /= psb_success_) goto 100
    if (bnorm == dzero) goto 100

    z(:) = dzero
    call prec%blocks(field)%pc%apply(done, r, dzero, z, field_desc, info, trans='N', work=wrk)
    if (info /= psb_success_) goto 100
    p(:) = z(:)
    rz = psb_gedot(r, z, field_desc, info)
    if (info /= psb_success_) goto 100

    do k = 1, max(1, kctx%itmax)
      q(:) = dzero
      call psb_d_nested_field_matvec(prec, field, p, q, info)
      if (info /= psb_success_) exit
      denom = psb_gedot(p, q, field_desc, info)
      if (info /= psb_success_) exit
      if (abs(denom) <= epsilon(done)) exit
      alpha = rz / denom
      sol(:) = sol(:) + alpha * p(:)
      r(:) = r(:) - alpha * q(:)
      rnorm = sqrt(max(dzero, psb_gedot(r, r, field_desc, info)))
      if (info /= psb_success_) exit
      if (rnorm <= kctx%tol * bnorm) exit
      z(:) = dzero
      call prec%blocks(field)%pc%apply(done, r, dzero, z, field_desc, info, trans='N', work=wrk)
      if (info /= psb_success_) exit
      rz_new = psb_gedot(r, z, field_desc, info)
      if (info /= psb_success_) exit
      if (abs(rz) <= epsilon(done)) exit
      beta = rz_new / rz
      p(:) = z(:) + beta * p(:)
      rz = rz_new
    end do

100 continue
    deallocate(r, z, p, q, wrk)
  end subroutine psb_d_nested_inner_cg

  subroutine psb_d_nested_inner_bicgstab(prec, field, kctx, rhs, sol, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: field
    type(psb_d_nested_krylov_context), intent(in) :: kctx
    real(psb_dpk_), intent(inout) :: rhs(:)
    real(psb_dpk_), intent(inout) :: sol(:)
    integer(psb_ipk_), intent(out) :: info

    real(psb_dpk_), allocatable :: r(:), r0(:), p(:), v(:), s(:), t(:), ph(:), sh(:), wrk(:)
    type(psb_desc_type), pointer :: field_desc
    integer(psb_ipk_) :: n_col_field, k
    real(psb_dpk_) :: bnorm, rnorm, rho, rho_old, alpha, beta, omega, denom

    info = psb_success_
    field_desc => psb_d_nest_get_field_desc(prec%nest_op, field)
    if (.not. associated(field_desc)) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'd_nested_inner_bicgstab', a_err='missing field descriptor')
      return
    end if
    n_col_field = field_desc%get_local_cols()
    allocate(r(n_col_field), r0(n_col_field), p(n_col_field), v(n_col_field), &
         & s(n_col_field), t(n_col_field), ph(n_col_field), sh(n_col_field), &
         & wrk(n_col_field), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_inner_bicgstab', a_err='inner work vectors')
      return
    end if

    sol(:) = dzero
    r(:) = rhs(:)
    r0(:) = r(:)
    p(:) = dzero
    v(:) = dzero
    rho_old = done
    alpha = done
    omega = done
    bnorm = sqrt(max(dzero, psb_gedot(rhs, rhs, field_desc, info)))
    if (info /= psb_success_) goto 100
    if (bnorm == dzero) goto 100

    do k = 1, max(1, kctx%itmax)
      rho = psb_gedot(r0, r, field_desc, info)
      if (info /= psb_success_) exit
      if (abs(rho) <= epsilon(done)) exit
      if (k == 1) then
        p(:) = r(:)
      else
        if (abs(omega) <= epsilon(done)) exit
        beta = (rho / rho_old) * (alpha / omega)
        p(:) = r(:) + beta * (p(:) - omega * v(:))
      end if

      ph(:) = dzero
      call prec%blocks(field)%pc%apply(done, p, dzero, ph, field_desc, info, trans='N', work=wrk)
      if (info /= psb_success_) exit
      v(:) = dzero
      call psb_d_nested_field_matvec(prec, field, ph, v, info)
      if (info /= psb_success_) exit
      denom = psb_gedot(r0, v, field_desc, info)
      if (info /= psb_success_) exit
      if (abs(denom) <= epsilon(done)) exit
      alpha = rho / denom
      s(:) = r(:) - alpha * v(:)
      rnorm = sqrt(max(dzero, psb_gedot(s, s, field_desc, info)))
      if (info /= psb_success_) exit
      if (rnorm <= kctx%tol * bnorm) then
        sol(:) = sol(:) + alpha * ph(:)
        exit
      end if

      sh(:) = dzero
      call prec%blocks(field)%pc%apply(done, s, dzero, sh, field_desc, info, trans='N', work=wrk)
      if (info /= psb_success_) exit
      t(:) = dzero
      call psb_d_nested_field_matvec(prec, field, sh, t, info)
      if (info /= psb_success_) exit
      denom = psb_gedot(t, t, field_desc, info)
      if (info /= psb_success_) exit
      if (abs(denom) <= epsilon(done)) exit
      omega = psb_gedot(t, s, field_desc, info) / denom
      if (info /= psb_success_) exit
      sol(:) = sol(:) + alpha * ph(:) + omega * sh(:)
      r(:) = s(:) - omega * t(:)
      rnorm = sqrt(max(dzero, psb_gedot(r, r, field_desc, info)))
      if (info /= psb_success_) exit
      if (rnorm <= kctx%tol * bnorm) exit
      rho_old = rho
    end do

100 continue
    deallocate(r, r0, p, v, s, t, ph, sh, wrk)
  end subroutine psb_d_nested_inner_bicgstab

  ! Apply the matrix-free Schur action Sx on the second field.
  subroutine psb_d_nested_schur_action(prec, x2, y2, gwork, desc_data, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(inout) :: x2(:)
    real(psb_dpk_), intent(inout) :: y2(:)
    real(psb_dpk_), intent(inout) :: gwork(:)
    type(psb_desc_type), intent(in) :: desc_data
    integer(psb_ipk_), intent(out) :: info

    real(psb_dpk_), allocatable :: x2h(:), t1(:), w1(:), w1h(:), t2(:)
    type(psb_desc_type), pointer :: desc1, desc2
    type(psb_dspmat_type), pointer :: block12, block21, block22
    integer(psb_ipk_) :: n_col_1, n_col_2

    info = psb_success_
    desc1 => psb_d_nest_get_field_desc(prec%nest_op, 1)
    desc2 => psb_d_nest_get_field_desc(prec%nest_op, 2)
    if ((.not. associated(desc1)) .or. (.not. associated(desc2))) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'd_nested_schur_action', a_err='missing field descriptor')
      return
    end if

    n_col_1 = desc1%get_local_cols()
    n_col_2 = desc2%get_local_cols()
    allocate(x2h(n_col_2), t1(n_col_1), w1(n_col_1), w1h(n_col_1), &
         & t2(n_col_2), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_schur_action', a_err='field vectors')
      return
    end if

    x2h(:) = dzero
    t1(:) = dzero
    w1(:) = dzero
    w1h(:) = dzero
    t2(:) = dzero
    y2(:) = dzero

    gwork(:) = dzero
    call psb_d_nest_prolong_field(prec%nest_op, 2, x2, gwork, info)
    if (info /= psb_success_) goto 100
    call psb_halo(gwork, desc_data, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field_local(prec%nest_op, 2, gwork, x2h, info)
    if (info /= psb_success_) goto 100

    block22 => psb_d_nest_get_block(prec%nest_op, 2, 2)
    if (associated(block22)) then
      call psb_d_nest_apply_block(prec%nest_op, 2, 2, done, x2h, dzero, y2, info)
      if (info /= psb_success_) goto 100
    end if

    block12 => psb_d_nest_get_block(prec%nest_op, 1, 2)
    block21 => psb_d_nest_get_block(prec%nest_op, 2, 1)
    if (associated(block12) .and. associated(block21)) then
      call psb_d_nest_apply_block(prec%nest_op, 1, 2, done, x2h, dzero, t1, info)
      if (info /= psb_success_) goto 100
      call psb_d_nested_field_solve(prec, 1, t1, w1, info)
      if (info /= psb_success_) goto 100

      gwork(:) = dzero
      call psb_d_nest_prolong_field(prec%nest_op, 1, w1, gwork, info)
      if (info /= psb_success_) goto 100
      call psb_halo(gwork, desc_data, info)
      if (info /= psb_success_) goto 100
      call psb_d_nest_restrict_field_local(prec%nest_op, 1, gwork, w1h, info)
      if (info /= psb_success_) goto 100

      call psb_d_nest_apply_block(prec%nest_op, 2, 1, done, w1h, dzero, t2, info)
      if (info /= psb_success_) goto 100
      y2(:) = y2(:) - t2(:)
    end if

100 continue
    deallocate(x2h, t1, w1, w1h, t2)
  end subroutine psb_d_nested_schur_action


  subroutine psb_d_nested_replay_schur_options(prec, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: k

    info = psb_success_
    if (.not. associated(selfp_schur_block)) return
    if (.not. allocated(selfp_schur_block%pc)) return

    if (allocated(prec%field_iopts)) then
      do k = 1, size(prec%field_iopts)
        if ((prec%field_iopts(k)%field == 0) .or. (prec%field_iopts(k)%field == 2)) then
          call selfp_schur_block%pc%precset(prec%field_iopts(k)%what, &
               & prec%field_iopts(k)%val, info)
          if (info /= psb_success_) return
        end if
      end do
    end if

    if (allocated(prec%field_ropts)) then
      do k = 1, size(prec%field_ropts)
        if ((prec%field_ropts(k)%field == 0) .or. (prec%field_ropts(k)%field == 2)) then
          call selfp_schur_block%pc%precset(prec%field_ropts(k)%what, &
               & prec%field_ropts(k)%val, info)
          if (info /= psb_success_) return
        end if
      end do
    end if

    if (allocated(prec%field_copts)) then
      do k = 1, size(prec%field_copts)
        if ((prec%field_copts(k)%field == 0) .or. (prec%field_copts(k)%field == 2)) then
          call selfp_schur_block%pc%precset(prec%field_copts(k)%what, &
               & trim(prec%field_copts(k)%val), info)
          if (info /= psb_success_) return
        end if
      end do
    end if
  end subroutine psb_d_nested_replay_schur_options
  subroutine psb_d_nested_build_selfp(prec, info, amold, vmold, imold)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info
    class(psb_d_base_sparse_mat), intent(in), optional :: amold
    class(psb_d_base_vect_type), intent(in), optional :: vmold
    class(psb_i_base_vect_type), intent(in), optional :: imold

    type(psb_desc_type), pointer :: desc1, desc2
    type(psb_dspmat_type), pointer :: block11, block12, block21, block22
    real(psb_dpk_), allocatable :: diag1(:), acc(:), vals(:), va21(:), va12(:), va22(:)
    integer(psb_lpk_), allocatable :: owned_rows(:)
    integer(psb_ipk_), allocatable :: ia21(:), ja21(:), ia12(:), ja12(:), ia22(:), ja22(:)
    integer(psb_lpk_), allocatable :: rows(:), cols(:)
    integer(psb_ipk_) :: i, k, m, nz21, nz12, nz22, n_insert
    integer(psb_ipk_) :: n_owned_1, n_owned_2, vloc
    real(psb_dpk_) :: scale, row_norm

    info = psb_success_
    if (prec%nfields < 2) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'd_nested_build_selfp', a_err='SELFP requires two fields')
      return
    end if

    desc1 => psb_d_nest_get_field_desc(prec%nest_op, 1)
    desc2 => psb_d_nest_get_field_desc(prec%nest_op, 2)
    block11 => psb_d_nest_get_block(prec%nest_op, 1, 1)
    block12 => psb_d_nest_get_block(prec%nest_op, 1, 2)
    block21 => psb_d_nest_get_block(prec%nest_op, 2, 1)
    block22 => psb_d_nest_get_block(prec%nest_op, 2, 2)
    if ((.not. associated(desc1)) .or. (.not. associated(desc2)) .or. &
         & (.not. associated(block11)) .or. (.not. associated(block12)) .or. &
         & (.not. associated(block21))) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'd_nested_build_selfp', a_err='missing Schur blocks')
      return
    end if

    diag1 = block11%get_diag(info)
    if (info /= psb_success_) then
      call psb_errpush(info, 'd_nested_build_selfp', a_err='A11 get_diag')
      return
    end if

    n_owned_1 = desc1%get_local_rows()
    n_owned_2 = desc2%get_local_rows()
    owned_rows = desc2%get_global_indices(owned=.true.)
    allocate(acc(n_owned_2), rows(n_owned_2), cols(n_owned_2), vals(n_owned_2), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_build_selfp', a_err='work arrays')
      return
    end if

    if (associated(selfp_schur_mat)) then
      call selfp_schur_mat%free()
      deallocate(selfp_schur_mat, stat=info)
      nullify(selfp_schur_mat)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_
        call psb_errpush(info, 'd_nested_build_selfp', a_err='Schur matrix reset')
        goto 100
      end if
    end if
    allocate(selfp_schur_mat, stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_build_selfp', a_err='Schur matrix')
      goto 100
    end if
    call psb_spall(selfp_schur_mat, desc2, info)
    if (info /= psb_success_) goto 100

    do i = 1, n_owned_2
      acc(:) = dzero

      if (associated(block22)) then
        call block22%csget(i, i, nz22, ia22, ja22, va22, info, jmax=n_owned_2)
        if (info /= psb_success_) exit
        do m = 1, nz22
          if (ja22(m) >= 1 .and. ja22(m) <= n_owned_2) acc(ja22(m)) = acc(ja22(m)) + va22(m)
        end do
      end if

      call block21%csget(i, i, nz21, ia21, ja21, va21, info)
      if (info /= psb_success_) exit
      do k = 1, nz21
        vloc = ja21(k)
        if (vloc < 1 .or. vloc > n_owned_1) cycle
        if (abs(diag1(vloc)) <= tiny(done)) cycle
        scale = -va21(k) / diag1(vloc)
        call block12%csget(vloc, vloc, nz12, ia12, ja12, va12, info, jmax=n_owned_2)
        if (info /= psb_success_) exit
        do m = 1, nz12
          if (ja12(m) >= 1 .and. ja12(m) <= n_owned_2) &
               & acc(ja12(m)) = acc(ja12(m)) + scale * va12(m)
        end do
      end do
      if (info /= psb_success_) exit

      if (abs(acc(i)) <= tiny(done)) then
        row_norm = sum(abs(acc(:)))
        acc(i) = -max(row_norm, done)
      end if

      n_insert = 0
      do m = 1, n_owned_2
        if (abs(acc(m)) > dzero) then
          n_insert = n_insert + 1
          rows(n_insert) = owned_rows(i)
          cols(n_insert) = owned_rows(m)
          vals(n_insert) = acc(m)
        end if
      end do
      if (n_insert > 0) then
        call psb_spins(n_insert, rows(1:n_insert), cols(1:n_insert), vals(1:n_insert), &
             & selfp_schur_mat, desc2, info)
        if (info /= psb_success_) exit
      end if
    end do
    if (info /= psb_success_) goto 100

    call psb_spasb(selfp_schur_mat, desc2, info, dupl=psb_dupl_add_)
    if (info /= psb_success_) goto 100
    call selfp_schur_mat%cscnv(info, type='CSR')
    if (info /= psb_success_) goto 100

    if (associated(selfp_schur_block)) then
      if (allocated(selfp_schur_block%pc)) then
        call selfp_schur_block%pc%free(info)
        if (info /= psb_success_) return
        deallocate(selfp_schur_block%pc, stat=info)
        if (info /= psb_success_) return
      end if
      deallocate(selfp_schur_block, stat=info)
      if (info /= psb_success_) return
      nullify(selfp_schur_block)
    end if
    allocate(selfp_schur_block, stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_build_selfp', a_err='Schur block')
      goto 100
    end if
    selfp_schur_block%ptype = 'BJAC'
    call psb_d_nested_alloc_block_pc(selfp_schur_block, info)
    if (info /= psb_success_) goto 100
    call selfp_schur_block%pc%precinit(info)
    if (info /= psb_success_) goto 100
    if (info /= psb_success_) goto 100
    call selfp_schur_block%pc%precbld(selfp_schur_mat, desc2, info, &
         & amold=amold, vmold=vmold, imold=imold)
    if (info == psb_success_) prec%schur_selfp_built = .true.

100 continue
    if (allocated(diag1)) deallocate(diag1)
    if (allocated(owned_rows)) deallocate(owned_rows)
    if (allocated(acc)) deallocate(acc)
    if (allocated(rows)) deallocate(rows)
    if (allocated(cols)) deallocate(cols)
    if (allocated(vals)) deallocate(vals)
    if (allocated(ia21)) deallocate(ia21)
    if (allocated(ja21)) deallocate(ja21)
    if (allocated(va21)) deallocate(va21)
    if (allocated(ia12)) deallocate(ia12)
    if (allocated(ja12)) deallocate(ja12)
    if (allocated(va12)) deallocate(va12)
    if (allocated(ia22)) deallocate(ia22)
    if (allocated(ja22)) deallocate(ja22)
    if (allocated(va22)) deallocate(va22)
  end subroutine psb_d_nested_build_selfp
  ! Approximately solve the Schur block using A22 or the matrix-free Schur action.
  subroutine psb_d_nested_schur_solve(prec, rhs, sol, gwork, desc_data, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(inout) :: rhs(:)
    real(psb_dpk_), intent(inout) :: sol(:)
    real(psb_dpk_), intent(inout) :: gwork(:)
    type(psb_desc_type), intent(in) :: desc_data
    integer(psb_ipk_), intent(out) :: info

    real(psb_dpk_), allocatable :: sx(:), res(:), dz(:)
    integer(psb_ipk_) :: k, maxit, n_owned
    real(psb_dpk_) :: rnrm
    type(psb_desc_type), pointer :: desc2

    info = psb_success_
    if (prec%schur_solve == psb_d_nested_schur_a22_) then
      call psb_d_nested_field_solve(prec, 2, rhs, sol, info)
      return
    else if (prec%schur_solve == psb_d_nested_schur_selfp_) then
      if (.not. prec%schur_selfp_built .or. (.not. associated(selfp_schur_block)) .or. &
           & (.not. allocated(selfp_schur_block%pc))) then
        info = psb_err_invalid_mat_state_
        call psb_errpush(info, 'd_nested_schur_solve', a_err='SELFP not built')
        return
      end if
      allocate(dz(size(rhs)), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_
        call psb_errpush(info, 'd_nested_schur_solve', a_err='SELFP work vector')
        return
      end if
      desc2 => psb_d_nest_get_field_desc(prec%nest_op, 2)
      if (.not. associated(desc2)) then
        info = psb_err_invalid_mat_state_
        call psb_errpush(info, 'd_nested_schur_solve', a_err='missing field descriptor')
        deallocate(dz)
        return
      end if
      sol(:) = dzero
      call selfp_schur_block%pc%apply(done, rhs, dzero, sol, desc2, info, &
           & trans='N', work=dz)
      deallocate(dz)
      return
    end if

    allocate(sx(size(rhs)), res(size(rhs)), dz(size(rhs)), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_schur_solve', a_err='Schur work vectors')
      return
    end if

    sol(:) = dzero
    maxit = max(1, prec%schur_maxit)
    n_owned = psb_d_nest_get_field_owned(prec%nest_op, 2)
    do k = 1, maxit
      sx(:) = dzero
      call psb_d_nested_schur_action(prec, sol, sx, gwork, desc_data, info)
      if (info /= psb_success_) exit
      res(:) = rhs(:) - sx(:)
      if (prec%schur_tol > dzero) then
        rnrm = sqrt(sum(res(1:n_owned) * res(1:n_owned)))
        if (rnrm <= prec%schur_tol) exit
      end if
      call psb_d_nested_field_solve(prec, 2, res, dz, info)
      if (info /= psb_success_) exit
      sol(:) = sol(:) + dz(:)
    end do

    deallocate(sx, res, dz)
  end subroutine psb_d_nested_schur_solve

  ! Apply the two-field lower, upper, or full Schur-style composition.
  subroutine psb_d_nested_apply_schur(prec, x, z, desc_data, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(in) :: x(:)
    real(psb_dpk_), intent(inout) :: z(:)
    type(psb_desc_type), intent(in) :: desc_data
    integer(psb_ipk_), intent(out) :: info

    real(psb_dpk_), allocatable :: b1(:), b2(:), x1(:), x2(:), t1(:), t2(:), c1(:), gwork(:)
    type(psb_desc_type), pointer :: desc1, desc2
    type(psb_dspmat_type), pointer :: block12, block21
    integer(psb_ipk_) :: n_col_1, n_col_2

    info = psb_success_
    z(:) = dzero

    if (prec%nfields /= 2) then
      info = psb_err_invalid_input_
      call psb_errpush(info, 'd_nested_apply_schur', a_err='Schur composition requires two fields')
      return
    end if

    desc1 => psb_d_nest_get_field_desc(prec%nest_op, 1)
    desc2 => psb_d_nest_get_field_desc(prec%nest_op, 2)
    if ((.not. associated(desc1)) .or. (.not. associated(desc2))) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'd_nested_apply_schur', a_err='missing field descriptor')
      return
    end if

    n_col_1 = desc1%get_local_cols()
    n_col_2 = desc2%get_local_cols()
    allocate(b1(n_col_1), b2(n_col_2), x1(n_col_1), x2(n_col_2), &
         & t1(n_col_1), t2(n_col_2), c1(n_col_1), gwork(size(z)), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_apply_schur', a_err='field vectors')
      return
    end if

    b1(:) = dzero
    b2(:) = dzero
    x1(:) = dzero
    x2(:) = dzero
    t1(:) = dzero
    t2(:) = dzero
    c1(:) = dzero
    gwork(:) = dzero

    call psb_d_nest_restrict_field(prec%nest_op, 1, x, b1, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field(prec%nest_op, 2, x, b2, info)
    if (info /= psb_success_) goto 100

    block12 => psb_d_nest_get_block(prec%nest_op, 1, 2)
    block21 => psb_d_nest_get_block(prec%nest_op, 2, 1)

    select case (prec%composition)
    case (psb_d_nested_schur_lower_, psb_d_nested_schur_full_)
      call psb_d_nested_field_solve(prec, 1, b1, x1, info)
      if (info /= psb_success_) goto 100
      call psb_d_nest_prolong_field(prec%nest_op, 1, x1, z, info)
      if (info /= psb_success_) goto 100
      call psb_halo(z, desc_data, info)
      if (info /= psb_success_) goto 100
      call psb_d_nest_restrict_field_local(prec%nest_op, 1, z, x1, info)
      if (info /= psb_success_) goto 100

      t2(:) = b2(:)
      if (associated(block21)) then
        call psb_d_nest_apply_block(prec%nest_op, 2, 1, -done, x1, done, t2, info)
        if (info /= psb_success_) goto 100
      end if
      call psb_d_nested_schur_solve(prec, t2, x2, gwork, desc_data, info)
      if (info /= psb_success_) goto 100
      call psb_d_nest_prolong_field(prec%nest_op, 2, x2, z, info)
      if (info /= psb_success_) goto 100

      if (prec%composition == psb_d_nested_schur_full_) then
        call psb_halo(z, desc_data, info)
        if (info /= psb_success_) goto 100
        call psb_d_nest_restrict_field_local(prec%nest_op, 2, z, x2, info)
        if (info /= psb_success_) goto 100
        t1(:) = dzero
        if (associated(block12)) then
          call psb_d_nest_apply_block(prec%nest_op, 1, 2, done, x2, dzero, t1, info)
          if (info /= psb_success_) goto 100
        end if
        call psb_d_nested_field_solve(prec, 1, t1, c1, info)
        if (info /= psb_success_) goto 100
        x1(:) = x1(:) - c1(:)
        call psb_d_nest_prolong_field(prec%nest_op, 1, x1, z, info)
        if (info /= psb_success_) goto 100
      end if

    case (psb_d_nested_schur_upper_)
      call psb_d_nested_schur_solve(prec, b2, x2, gwork, desc_data, info)
      if (info /= psb_success_) goto 100
      call psb_d_nest_prolong_field(prec%nest_op, 2, x2, z, info)
      if (info /= psb_success_) goto 100
      call psb_halo(z, desc_data, info)
      if (info /= psb_success_) goto 100
      call psb_d_nest_restrict_field_local(prec%nest_op, 2, z, x2, info)
      if (info /= psb_success_) goto 100

      t1(:) = b1(:)
      if (associated(block12)) then
        call psb_d_nest_apply_block(prec%nest_op, 1, 2, -done, x2, done, t1, info)
        if (info /= psb_success_) goto 100
      end if
      call psb_d_nested_field_solve(prec, 1, t1, x1, info)
      if (info /= psb_success_) goto 100
      call psb_d_nest_prolong_field(prec%nest_op, 1, x1, z, info)
      if (info /= psb_success_) goto 100
    end select

    call psb_halo(z, desc_data, info)

100 continue
    deallocate(b1, b2, x1, x2, t1, t2, c1, gwork)
  end subroutine psb_d_nested_apply_schur

  ! Apply the positive PDE-control Schur action on field 3:
  !   T x3 = [A31 A32] diag(B1,B2) [A13; A23] x3 - A33 x3.
  ! For the KKT matrix [M 0 K; 0 alpha*M -M; K -M 0], T is -S.
  subroutine psb_d_nested_pde_control_schur_action(prec, x3, y3, gwork, desc_data, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(inout) :: x3(:)
    real(psb_dpk_), intent(inout) :: y3(:)
    real(psb_dpk_), intent(inout) :: gwork(:)
    type(psb_desc_type), intent(in) :: desc_data
    integer(psb_ipk_), intent(out) :: info

    real(psb_dpk_), allocatable :: x3h(:), t1(:), t2(:), w1(:), w2(:), &
         & w1h(:), w2h(:), t3(:)
    type(psb_desc_type), pointer :: desc1, desc2, desc3
    type(psb_dspmat_type), pointer :: block13, block23, block31, block32, block33
    integer(psb_ipk_) :: n_col_1, n_col_2, n_col_3

    info = psb_success_
    desc1 => psb_d_nest_get_field_desc(prec%nest_op, 1)
    desc2 => psb_d_nest_get_field_desc(prec%nest_op, 2)
    desc3 => psb_d_nest_get_field_desc(prec%nest_op, 3)
    if ((.not. associated(desc1)) .or. (.not. associated(desc2)) .or. &
         & (.not. associated(desc3))) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'd_nested_pde_schur_action', a_err='missing field descriptor')
      return
    end if

    n_col_1 = desc1%get_local_cols()
    n_col_2 = desc2%get_local_cols()
    n_col_3 = desc3%get_local_cols()
    allocate(x3h(n_col_3), t1(n_col_1), t2(n_col_2), w1(n_col_1), w2(n_col_2), &
         & w1h(n_col_1), w2h(n_col_2), t3(n_col_3), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_pde_schur_action', a_err='field vectors')
      return
    end if

    x3h(:) = dzero
    t1(:) = dzero
    t2(:) = dzero
    w1(:) = dzero
    w2(:) = dzero
    w1h(:) = dzero
    w2h(:) = dzero
    t3(:) = dzero
    y3(:) = dzero

    gwork(:) = dzero
    call psb_d_nest_prolong_field(prec%nest_op, 3, x3, gwork, info)
    if (info /= psb_success_) goto 100
    call psb_halo(gwork, desc_data, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field_local(prec%nest_op, 3, gwork, x3h, info)
    if (info /= psb_success_) goto 100

    block13 => psb_d_nest_get_block(prec%nest_op, 1, 3)
    block23 => psb_d_nest_get_block(prec%nest_op, 2, 3)
    block31 => psb_d_nest_get_block(prec%nest_op, 3, 1)
    block32 => psb_d_nest_get_block(prec%nest_op, 3, 2)
    block33 => psb_d_nest_get_block(prec%nest_op, 3, 3)

    if (associated(block13) .and. associated(block31)) then
      call psb_d_nest_apply_block(prec%nest_op, 1, 3, done, x3h, dzero, t1, info)
      if (info /= psb_success_) goto 100
      call psb_d_nested_field_solve(prec, 1, t1, w1, info)
      if (info /= psb_success_) goto 100

      gwork(:) = dzero
      call psb_d_nest_prolong_field(prec%nest_op, 1, w1, gwork, info)
      if (info /= psb_success_) goto 100
      call psb_halo(gwork, desc_data, info)
      if (info /= psb_success_) goto 100
      call psb_d_nest_restrict_field_local(prec%nest_op, 1, gwork, w1h, info)
      if (info /= psb_success_) goto 100

      call psb_d_nest_apply_block(prec%nest_op, 3, 1, done, w1h, done, y3, info)
      if (info /= psb_success_) goto 100
    end if

    if (associated(block23) .and. associated(block32)) then
      call psb_d_nest_apply_block(prec%nest_op, 2, 3, done, x3h, dzero, t2, info)
      if (info /= psb_success_) goto 100
      call psb_d_nested_field_solve(prec, 2, t2, w2, info)
      if (info /= psb_success_) goto 100

      gwork(:) = dzero
      call psb_d_nest_prolong_field(prec%nest_op, 2, w2, gwork, info)
      if (info /= psb_success_) goto 100
      call psb_halo(gwork, desc_data, info)
      if (info /= psb_success_) goto 100
      call psb_d_nest_restrict_field_local(prec%nest_op, 2, gwork, w2h, info)
      if (info /= psb_success_) goto 100

      call psb_d_nest_apply_block(prec%nest_op, 3, 2, done, w2h, done, y3, info)
      if (info /= psb_success_) goto 100
    end if

    if (associated(block33)) then
      call psb_d_nest_apply_block(prec%nest_op, 3, 3, done, x3h, dzero, t3, info)
      if (info /= psb_success_) goto 100
      y3(:) = y3(:) - t3(:)
    end if

100 continue
    deallocate(x3h, t1, t2, w1, w2, w1h, w2h, t3)
  end subroutine psb_d_nested_pde_control_schur_action

  ! Approximately solve the field-3 PDE-control Schur equation S x = rhs.
  ! Internally this solves T x = -rhs, where T = -S is the positive
  ! PDE-control Schur action assembled matrix-free by
  ! psb_d_nested_pde_control_schur_action.  The old Richardson update used the
  ! field-3 diagonal block as the Schur preconditioner; for KKT systems with
  ! A33 = 0 that is only an identity fallback.  Use a matrix-free CG iteration
  ! instead, with the field-3 block preconditioner as an optional left
  ! preconditioner when it is meaningful.
  subroutine psb_d_nested_pde_control_schur_solve(prec, rhs, sol, gwork, desc_data, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(inout) :: rhs(:)
    real(psb_dpk_), intent(inout) :: sol(:)
    real(psb_dpk_), intent(inout) :: gwork(:)
    type(psb_desc_type), intent(in) :: desc_data
    integer(psb_ipk_), intent(out) :: info

    real(psb_dpk_), allocatable :: q(:), r(:), zc(:), p(:), ap(:)
    type(psb_desc_type), pointer :: desc3
    integer(psb_ipk_) :: k, maxit
    real(psb_dpk_) :: alpha, beta, rz, rz_old, pap, rnorm, rhs_norm, tol

    info = psb_success_
    desc3 => psb_d_nest_get_field_desc(prec%nest_op, 3)
    if (.not. associated(desc3)) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'd_nested_pde_schur_solve', a_err='missing field descriptor')
      return
    end if

    if (prec%schur_solve == psb_d_nested_schur_a22_) then
      call psb_d_nested_field_solve(prec, 3, rhs, sol, info)
      if (info == psb_success_) sol(:) = -sol(:)
      return
    end if

    allocate(q(size(rhs)), r(size(rhs)), zc(size(rhs)), p(size(rhs)), ap(size(rhs)), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_pde_schur_solve', a_err='Schur work vectors')
      return
    end if

    sol(:) = dzero
    q(:) = -rhs(:)
    r(:) = q(:)
    rhs_norm = sqrt(max(dzero, psb_gedot(q, q, desc3, info)))
    if (info /= psb_success_) goto 100
    tol = prec%schur_tol
    if (tol <= dzero) tol = 1.0e-6_psb_dpk_ * max(done, rhs_norm)

    call psb_d_nested_field_solve(prec, 3, r, zc, info)
    if (info /= psb_success_) goto 100
    p(:) = zc(:)
    rz = psb_gedot(r, zc, desc3, info)
    if (info /= psb_success_) goto 100
    if (abs(rz) <= tiny(done)) then
      p(:) = r(:)
      zc(:) = r(:)
      rz = psb_gedot(r, zc, desc3, info)
      if (info /= psb_success_) goto 100
    end if

    maxit = max(1, prec%schur_maxit)
    do k = 1, maxit
      ap(:) = dzero
      call psb_d_nested_pde_control_schur_action(prec, p, ap, gwork, desc_data, info)
      if (info /= psb_success_) exit
      pap = psb_gedot(p, ap, desc3, info)
      if (info /= psb_success_) exit
      if (abs(pap) <= tiny(done)) exit
      alpha = rz / pap
      sol(:) = sol(:) + alpha * p(:)
      r(:) = r(:) - alpha * ap(:)
      rnorm = sqrt(max(dzero, psb_gedot(r, r, desc3, info)))
      if (info /= psb_success_) exit
      if (rnorm <= tol) exit

      call psb_d_nested_field_solve(prec, 3, r, zc, info)
      if (info /= psb_success_) exit
      rz_old = rz
      rz = psb_gedot(r, zc, desc3, info)
      if (info /= psb_success_) exit
      if (abs(rz) <= tiny(done)) then
        zc(:) = r(:)
        rz = psb_gedot(r, zc, desc3, info)
        if (info /= psb_success_) exit
      end if
      beta = rz / rz_old
      p(:) = zc(:) + beta * p(:)
    end do

100 continue
    deallocate(q, r, zc, p, ap)
  end subroutine psb_d_nested_pde_control_schur_solve

  ! Apply a three-field PDE-control Schur preconditioner.
  ! Fields 1 and 2 are grouped into the leading block D=diag(A11,A22);
  ! field 3 is the Schur field.
  subroutine psb_d_nested_apply_pde_control_schur(prec, x, z, desc_data, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(in) :: x(:)
    real(psb_dpk_), intent(inout) :: z(:)
    type(psb_desc_type), intent(in) :: desc_data
    integer(psb_ipk_), intent(out) :: info

    real(psb_dpk_), allocatable :: b1(:), b2(:), b3(:), x1(:), x2(:), x3(:), &
         & t1(:), t2(:), c1(:), c2(:), gwork(:)
    type(psb_desc_type), pointer :: desc1, desc2, desc3
    type(psb_dspmat_type), pointer :: block13, block23, block31, block32
    integer(psb_ipk_) :: n_col_1, n_col_2, n_col_3

    info = psb_success_
    z(:) = dzero

    if (prec%nfields /= 3) then
      info = psb_err_invalid_input_
      call psb_errpush(info, 'd_nested_apply_pde_schur', &
           & a_err='PDE-control Schur composition requires three fields')
      return
    end if

    desc1 => psb_d_nest_get_field_desc(prec%nest_op, 1)
    desc2 => psb_d_nest_get_field_desc(prec%nest_op, 2)
    desc3 => psb_d_nest_get_field_desc(prec%nest_op, 3)
    if ((.not. associated(desc1)) .or. (.not. associated(desc2)) .or. &
         & (.not. associated(desc3))) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'd_nested_apply_pde_schur', a_err='missing field descriptor')
      return
    end if

    n_col_1 = desc1%get_local_cols()
    n_col_2 = desc2%get_local_cols()
    n_col_3 = desc3%get_local_cols()
    allocate(b1(n_col_1), b2(n_col_2), b3(n_col_3), x1(n_col_1), x2(n_col_2), &
         & x3(n_col_3), t1(n_col_1), t2(n_col_2), c1(n_col_1), c2(n_col_2), &
         & gwork(size(z)), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'd_nested_apply_pde_schur', a_err='field vectors')
      return
    end if

    b1(:) = dzero
    b2(:) = dzero
    b3(:) = dzero
    x1(:) = dzero
    x2(:) = dzero
    x3(:) = dzero
    t1(:) = dzero
    t2(:) = dzero
    c1(:) = dzero
    c2(:) = dzero
    gwork(:) = dzero

    call psb_d_nest_restrict_field(prec%nest_op, 1, x, b1, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field(prec%nest_op, 2, x, b2, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field(prec%nest_op, 3, x, b3, info)
    if (info /= psb_success_) goto 100

    block13 => psb_d_nest_get_block(prec%nest_op, 1, 3)
    block23 => psb_d_nest_get_block(prec%nest_op, 2, 3)
    block31 => psb_d_nest_get_block(prec%nest_op, 3, 1)
    block32 => psb_d_nest_get_block(prec%nest_op, 3, 2)

    call psb_d_nested_field_solve(prec, 1, b1, x1, info)
    if (info /= psb_success_) goto 100
    call psb_d_nested_field_solve(prec, 2, b2, x2, info)
    if (info /= psb_success_) goto 100

    call psb_d_nest_prolong_field(prec%nest_op, 1, x1, z, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_prolong_field(prec%nest_op, 2, x2, z, info)
    if (info /= psb_success_) goto 100
    call psb_halo(z, desc_data, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field_local(prec%nest_op, 1, z, x1, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field_local(prec%nest_op, 2, z, x2, info)
    if (info /= psb_success_) goto 100

    if (associated(block31)) then
      call psb_d_nest_apply_block(prec%nest_op, 3, 1, -done, x1, done, b3, info)
      if (info /= psb_success_) goto 100
    end if
    if (associated(block32)) then
      call psb_d_nest_apply_block(prec%nest_op, 3, 2, -done, x2, done, b3, info)
      if (info /= psb_success_) goto 100
    end if

    call psb_d_nested_pde_control_schur_solve(prec, b3, x3, gwork, desc_data, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_prolong_field(prec%nest_op, 3, x3, z, info)
    if (info /= psb_success_) goto 100

    call psb_halo(z, desc_data, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field_local(prec%nest_op, 3, z, x3, info)
    if (info /= psb_success_) goto 100

    if (associated(block13)) then
      call psb_d_nest_apply_block(prec%nest_op, 1, 3, done, x3, dzero, t1, info)
      if (info /= psb_success_) goto 100
      call psb_d_nested_field_solve(prec, 1, t1, c1, info)
      if (info /= psb_success_) goto 100
      x1(:) = x1(:) - c1(:)
      call psb_d_nest_prolong_field(prec%nest_op, 1, x1, z, info)
      if (info /= psb_success_) goto 100
    end if
    if (associated(block23)) then
      call psb_d_nest_apply_block(prec%nest_op, 2, 3, done, x3, dzero, t2, info)
      if (info /= psb_success_) goto 100
      call psb_d_nested_field_solve(prec, 2, t2, c2, info)
      if (info /= psb_success_) goto 100
      x2(:) = x2(:) - c2(:)
      call psb_d_nest_prolong_field(prec%nest_op, 2, x2, z, info)
      if (info /= psb_success_) goto 100
    end if

    call psb_halo(z, desc_data, info)

100 continue
    deallocate(b1, b2, b3, x1, x2, x3, t1, t2, c1, c2, gwork)
  end subroutine psb_d_nested_apply_pde_control_schur

  ! Print a short description of the nested preconditioner configuration.
  subroutine psb_d_nested_precdescr(prec, iout, root, verbosity, prefix)
    class(psb_d_nested_prec_type), intent(in) :: prec
    integer(psb_ipk_), intent(in), optional :: iout
    integer(psb_ipk_), intent(in), optional :: root
    integer(psb_ipk_), intent(in), optional :: verbosity
    character(len=*), intent(in), optional :: prefix

    integer(psb_ipk_) :: iout_, root_, verbosity_, iam, np
    character(len=1024) :: prefix_

    iout_ = 6
    if (present(iout)) iout_ = iout
    root_ = psb_root_
    if (present(root)) root_ = root
    verbosity_ = 0
    if (present(verbosity)) verbosity_ = verbosity
    prefix_ = ''
    if (present(prefix)) prefix_ = prefix
    if (verbosity_ < 0) return

    call psb_info(prec%ctxt, iam, np)
    if (root_ == -1) root_ = iam
    if (iam == root_) then
      write(iout_,*) trim(prefix_), ' ', trim(prec%desc_prefix()), &
           & ' Nested block preconditioner: ', trim(prec%composition_name), &
           & ', block solve: ', trim(prec%default_block_ptype)
      if (prec%default_krylov%enabled) &
           & write(iout_,*) trim(prefix_), ' default inner Krylov: ', &
           & trim(prec%default_krylov%method), ' maxit=', prec%default_krylov%itmax, &
           & ' tol=', prec%default_krylov%tol
    end if
  end subroutine psb_d_nested_precdescr

  ! Provide a placeholder dump hook for the nested preconditioner.
  subroutine psb_d_nested_dump(prec, info, prefix, head)
    class(psb_d_nested_prec_type), intent(in) :: prec
    integer(psb_ipk_), intent(out) :: info
    character(len=*), intent(in), optional :: prefix, head
    info = psb_success_
  end subroutine psb_d_nested_dump

  function psb_d_nested_sizeof(prec) result(val)
    class(psb_d_nested_prec_type), intent(in) :: prec
    integer(psb_epk_) :: val
    integer(psb_ipk_) :: i
    val = 0_psb_epk_
    if (allocated(prec%blocks)) then
      do i = 1, size(prec%blocks)
        if (allocated(prec%blocks(i)%pc)) val = val + prec%blocks(i)%pc%sizeof()
      end do
    end if
  end function psb_d_nested_sizeof

  function psb_d_nested_get_nzeros(prec) result(val)
    class(psb_d_nested_prec_type), intent(in) :: prec
    integer(psb_epk_) :: val
    integer(psb_ipk_) :: i
    val = 0_psb_epk_
    if (allocated(prec%blocks)) then
      do i = 1, size(prec%blocks)
        if (allocated(prec%blocks(i)%pc)) val = val + prec%blocks(i)%pc%get_nzeros()
      end do
    end if
  end function psb_d_nested_get_nzeros

  ! Clone nested preconditioner configuration and any built field blocks.
  subroutine psb_d_nested_clone(prec, precout, info)
    class(psb_d_nested_prec_type), intent(inout) :: prec
    class(psb_d_base_prec_type), allocatable, intent(inout) :: precout
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: i

    info = psb_success_
    if (allocated(precout)) then
      call precout%free(info)
      if (info == psb_success_) deallocate(precout, stat=info)
      if (info /= psb_success_) return
    end if
    allocate(psb_d_nested_prec_type :: precout, stat=info)
    if (info /= 0) return
    select type (pout => precout)
    type is (psb_d_nested_prec_type)
      call pout%set_ctxt(prec%get_ctxt())
      pout%composition = prec%composition
      pout%composition_name = prec%composition_name
      pout%default_block_ptype = prec%default_block_ptype
      pout%schur_solve = prec%schur_solve
      pout%schur_solve_name = prec%schur_solve_name
      pout%schur_maxit = prec%schur_maxit
      pout%schur_tol = prec%schur_tol
      pout%default_krylov = prec%default_krylov
      pout%nfields = prec%nfields
      pout%nest_op => prec%nest_op
      if (allocated(prec%field_block_ptype)) then
        allocate(pout%field_block_ptype(size(prec%field_block_ptype)), stat=info)
        if (info /= 0) return
        pout%field_block_ptype(:) = prec%field_block_ptype(:)
      end if
      if (allocated(prec%field_krylov)) then
        allocate(pout%field_krylov(size(prec%field_krylov)), stat=info)
        if (info /= 0) return
        pout%field_krylov(:) = prec%field_krylov(:)
      end if
      if (allocated(prec%field_iopts)) then
        allocate(pout%field_iopts(size(prec%field_iopts)), stat=info)
        if (info /= 0) return
        pout%field_iopts(:) = prec%field_iopts(:)
      end if
      if (allocated(prec%field_ropts)) then
        allocate(pout%field_ropts(size(prec%field_ropts)), stat=info)
        if (info /= 0) return
        pout%field_ropts(:) = prec%field_ropts(:)
      end if
      if (allocated(prec%field_copts)) then
        allocate(pout%field_copts(size(prec%field_copts)), stat=info)
        if (info /= 0) return
        pout%field_copts(:) = prec%field_copts(:)
      end if
      if (allocated(prec%blocks)) then
        allocate(pout%blocks(size(prec%blocks)), stat=info)
        if (info /= 0) return
        do i = 1, size(prec%blocks)
          pout%blocks(i)%ptype = prec%blocks(i)%ptype
          if (allocated(prec%blocks(i)%pc)) &
               & call prec%blocks(i)%pc%clone(pout%blocks(i)%pc, info)
          if (info /= psb_success_) return
        end do
      end if
    end select
  end subroutine psb_d_nested_clone

end module psb_d_nestedprec
