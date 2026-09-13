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
! package: psb_c_vect_mod
!
! This module contains the definition of the psb_c_vect type which
! is the outer container for dense vectors.
! Therefore all methods simply invoke the corresponding methods of the
! inner component.
!
module psb_c_vect_mod
  use psb_realloc_mod
  use psb_c_base_vect_mod
  use psb_i_vect_mod

  type psb_c_vect_type
    class(psb_c_base_vect_type), allocatable :: v
    integer(psb_ipk_) :: nrmv = 0
    integer(psb_ipk_) :: remote_build = psb_matbld_noremote_
    integer(psb_ipk_) :: dupl = psb_dupl_add_
    complex(psb_spk_), allocatable :: rmtv(:)
    integer(psb_lpk_), allocatable :: rmidx(:)
  contains
    procedure, pass(x) :: get_nrows => c_vect_get_nrows
    procedure, pass(x) :: sizeof   => c_vect_sizeof
    procedure, pass(x) :: get_fmt  => c_vect_get_fmt
    procedure, pass(x) :: is_remote_build => c_vect_is_remote_build
    procedure, pass(x) :: set_remote_build => c_vect_set_remote_build
    procedure, pass(x) :: get_nrmv => c_vect_get_nrmv
    procedure, pass(x) :: set_nrmv => c_vect_set_nrmv
    procedure, pass(x) :: all      => c_vect_all
    procedure, pass(x) :: reall    => c_vect_reall
    procedure, pass(x) :: zero     => c_vect_zero
    procedure, pass(x) :: asb      => c_vect_asb
    procedure, pass(x) :: set_dupl => c_vect_set_dupl 
    procedure, pass(x) :: get_dupl => c_vect_get_dupl
    procedure, pass(x) :: set_ncfs => c_vect_set_ncfs 
    procedure, pass(x) :: get_ncfs => c_vect_get_ncfs
    procedure, pass(x) :: set_state => c_vect_set_state
    procedure, pass(x) :: set_null  => c_vect_set_null
    procedure, pass(x) :: set_bld   => c_vect_set_bld
    procedure, pass(x) :: set_upd   => c_vect_set_upd
    procedure, pass(x) :: set_asb   => c_vect_set_asb
    procedure, pass(x) :: get_state => c_vect_get_state
    procedure, pass(x) :: is_null   => c_vect_is_null
    procedure, pass(x) :: is_bld    => c_vect_is_bld
    procedure, pass(x) :: is_upd    => c_vect_is_upd
    procedure, pass(x) :: is_asb    => c_vect_is_asb
    procedure, pass(x) :: reinit    => c_vect_reinit

    procedure, pass(x) :: gthab    => c_vect_gthab
    procedure, pass(x) :: gthzv    => c_vect_gthzv
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => c_vect_sctb
    generic, public    :: sct      => sctb
    procedure, pass(x) :: free     => c_vect_free
    procedure, pass(x) :: ins_a    => c_vect_ins_a
    procedure, pass(x) :: ins_v    => c_vect_ins_v
    generic, public    :: ins      => ins_v, ins_a
    procedure, pass(x) :: bld_x    => c_vect_bld_x
    procedure, pass(x) :: bld_mn   => c_vect_bld_mn
    procedure, pass(x) :: bld_en   => c_vect_bld_en
    generic, public    :: bld      => bld_x, bld_mn, bld_en
    procedure, pass(x) :: get_vect => c_vect_get_vect
    procedure, pass(x) :: cnv      => c_vect_cnv
    procedure, pass(x) :: set_scal => c_vect_set_scal
    procedure, pass(x) :: set_vect => c_vect_set_vect
    generic, public    :: set      => set_vect, set_scal
    procedure, pass(x) :: clone    => c_vect_clone

    procedure, pass(x) :: sync     => c_vect_sync
    procedure, pass(x) :: is_host  => c_vect_is_host
    procedure, pass(x) :: is_dev   => c_vect_is_dev
    procedure, pass(x) :: is_sync  => c_vect_is_sync
    procedure, pass(x) :: set_host => c_vect_set_host
    procedure, pass(x) :: set_dev  => c_vect_set_dev
    procedure, pass(x) :: set_sync => c_vect_set_sync
    procedure, pass(x) :: check_addr => c_vect_check_addr
    procedure, pass(x) :: get_entry => c_vect_get_entry
    procedure, pass(x) :: set_entry => c_vect_set_entry

    procedure, pass(x) :: dot_v    => c_vect_dot_v
    procedure, pass(x) :: dot_a    => c_vect_dot_a
    generic, public    :: dot      => dot_v, dot_a
    procedure, pass(y) :: axpby_v  => c_vect_axpby_v
    procedure, pass(y) :: axpby_a  => c_vect_axpby_a
    procedure, pass(z) :: axpby_v2 => c_vect_axpby_v2
    procedure, pass(z) :: axpby_a2 => c_vect_axpby_a2
    generic, public    :: axpby    => axpby_v, axpby_a, axpby_v2, axpby_a2
    procedure, pass(z) :: upd_xyz  => c_vect_upd_xyz
    procedure, pass(z) :: xyzw     => c_vect_xyzw
    procedure, pass(y) :: mlt_v    => c_vect_mlt_v
    procedure, pass(y) :: mlt_a    => c_vect_mlt_a
    procedure, pass(z) :: mlt_a_2  => c_vect_mlt_a_2
    procedure, pass(z) :: mlt_v_2  => c_vect_mlt_v_2
    procedure, pass(z) :: mlt_va   => c_vect_mlt_va
    procedure, pass(z) :: mlt_av   => c_vect_mlt_av
    generic, public    :: mlt      => mlt_v, mlt_a, mlt_a_2, &
         & mlt_v_2, mlt_av, mlt_va
    procedure, pass(x) :: div_v    => c_vect_div_v
    procedure, pass(z) :: div_v2    => c_vect_div_v2
    procedure, pass(x) :: div_v_check => c_vect_div_v_check
    procedure, pass(x) :: div_v2_check => c_vect_div_v2_check
    procedure, pass(z) :: div_a2   => c_vect_div_a2
    procedure, pass(z) :: div_a2_check => c_vect_div_a2_check
    generic, public    :: div      => div_v, div_v2, div_v_check, &
                                        div_v2_check, div_a2, div_a2_check
    procedure, pass(y) :: inv_v    => c_vect_inv_v
    procedure, pass(y) :: inv_v_check => c_vect_inv_v_check
    procedure, pass(y) :: inv_a2   => c_vect_inv_a2
    procedure, pass(y) :: inv_a2_check => c_vect_inv_a2_check
    generic, public    :: inv      => inv_v, inv_v_check, inv_a2, inv_a2_check
    procedure, pass(x) :: scal     => c_vect_scal
    procedure, pass(x) :: absval1  => c_vect_absval1
    procedure, pass(x) :: absval2  => c_vect_absval2
    generic, public    :: absval   => absval1, absval2
    procedure, pass(x) :: nrm2std  => c_vect_nrm2
    procedure, pass(x) :: nrm2weight => c_vect_nrm2_weight
    procedure, pass(x) :: nrm2weightmask => c_vect_nrm2_weight_mask
    generic, public    :: nrm2     => nrm2std, nrm2weight, nrm2weightmask
    procedure, pass(x) :: amax     => c_vect_amax
    procedure, pass(x) :: asum     => c_vect_asum
    procedure, pass(z) :: acmp_a2   => c_vect_acmp_a2
    procedure, pass(z) :: acmp_v2   => c_vect_acmp_v2
    generic, public    :: acmp      => acmp_a2, acmp_v2
    procedure, pass(z) :: addconst_a2   => c_vect_addconst_a2
    procedure, pass(z) :: addconst_v2   => c_vect_addconst_v2
    generic, public    :: addconst      => addconst_a2, addconst_v2
  end type psb_c_vect_type

  public  :: psb_c_vect
  private :: constructor, size_const
  interface psb_c_vect
    module procedure constructor, size_const
  end interface psb_c_vect

  private :: c_vect_get_nrows, c_vect_sizeof, c_vect_get_fmt, &
            & c_vect_all, c_vect_reall, c_vect_zero, c_vect_asb, &
            & c_vect_gthab, c_vect_gthzv, c_vect_sctb, &
            & c_vect_free, c_vect_ins_a, c_vect_ins_v, c_vect_blc_x, &
            & c_vect_blc_mn, c_vect_blc_en, c_vect_get_vect, &
            & c_vect_cnv, c_vect_set_scal, &
            & c_vect_set_vect, c_vect_clone, c_vect_sync, c_vect_is_host, &
            & c_vect_is_dev, c_vect_is_sync, c_vect_set_host, &
            & c_vect_set_dev, c_vect_set_sync, &
            & c_vect_set_remote_build, c_is_remote_build, &
            & c_vect_set_dupl, c_get_dupl, c_vect_set_nrmv, c_get_nrmv
  private :: c_vect_dot_v, c_vect_dot_a, c_vect_axpby_v, c_vect_axpby_a, &
            & c_vect_mlt_v, c_vect_mlt_a, c_vect_mlt_a_2, c_vect_mlt_v_2, &
            & c_vect_mlt_va, c_vect_mlt_av, c_vect_scal, c_vect_absval1, &
            & c_vect_absval2, c_vect_nrm2, c_vect_amax, c_vect_asum
  class(psb_c_base_vect_type), allocatable, target, &
       & save, private :: psb_c_base_vect_default

  interface psb_set_vect_default
    module procedure psb_c_set_vect_default
  end interface psb_set_vect_default

  interface psb_get_vect_default
    module procedure psb_c_get_vect_default
  end interface psb_get_vect_default

contains
  function c_vect_get_dupl(x) result(res)
    implicit none
    class(psb_c_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    if(allocated(x%v)) then 
      res = x%v%get_dupl()
    else
      res = psb_dupl_null_
    end if
  end function c_vect_get_dupl

  subroutine c_vect_set_dupl(x, val)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if(allocated(x%v)) then 
      if(present(val)) then
        call x%v%set_dupl(val)
      else
        call x%v%set_dupl(psb_dupl_def_)
      end if
    end if
  end subroutine c_vect_set_dupl

  function c_vect_get_ncfs(x) result(res)
    implicit none
    class(psb_c_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    if(allocated(x%v)) then 
      res = x%v%get_ncfs()
    else
      res = 0
    end if
  end function c_vect_get_ncfs

  subroutine c_vect_set_ncfs(x, val)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if(allocated(x%v)) then 
      if(present(val)) then
        call x%v%set_ncfs(val)
      else
        call x%v%set_ncfs(0)
      end if
    end if
  end subroutine c_vect_set_ncfs

  function c_vect_get_state(x) result(res)
    implicit none
    class(psb_c_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    if(allocated(x%v)) then 
      res = x%v%get_state()
    else
      res = psb_vect_null_
    end if
  end function c_vect_get_state

  function c_vect_is_null(x) result(res)
    implicit none
    class(psb_c_vect_type), intent(in) :: x
    logical :: res

    res = (x%get_state() == psb_vect_null_)
  end function c_vect_is_null

  function c_vect_is_bld(x) result(res)
    implicit none
    class(psb_c_vect_type), intent(in) :: x
    logical :: res

    res = (x%get_state() == psb_vect_bld_)
  end function c_vect_is_bld

  function c_vect_is_upd(x) result(res)
    implicit none
    class(psb_c_vect_type), intent(in) :: x
    logical :: res

    res = (x%get_state() == psb_vect_upd_)
  end function c_vect_is_upd

  function c_vect_is_asb(x) result(res)
    implicit none
    class(psb_c_vect_type), intent(in) :: x
    logical :: res

    res = (x%get_state() == psb_vect_asb_)
  end function c_vect_is_asb

  subroutine c_vect_set_state(n, x)
    implicit none
    integer(psb_ipk_), intent(in)         :: n
    class(psb_c_vect_type), intent(inout) :: x

    if(allocated(x%v)) call x%v%set_state(n)
  end subroutine c_vect_set_state

  subroutine c_vect_set_null(x)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_null_)
  end subroutine c_vect_set_null

  subroutine c_vect_set_bld(x)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_bld_)
  end subroutine c_vect_set_bld

  subroutine c_vect_set_upd(x)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_upd_)
  end subroutine c_vect_set_upd

  subroutine c_vect_set_asb(x)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_asb_)
  end subroutine c_vect_set_asb

  function c_vect_get_nrmv(x) result(res)
    implicit none
    class(psb_c_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = x%nrmv
  end function c_vect_get_nrmv

  subroutine c_vect_set_nrmv(x, val)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: val

    x%nrmv = val
  end subroutine c_vect_set_nrmv

  function c_vect_is_remote_build(x) result(res)
    implicit none
    class(psb_c_vect_type), intent(in) :: x
    logical :: res

    res = (x%remote_build == psb_matbld_remote_)
  end function c_vect_is_remote_build

  subroutine c_vect_set_remote_build(x, val)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if(present(val)) then
      x%remote_build = val
    else
      x%remote_build = psb_matbld_remote_
    end if
  end subroutine c_vect_set_remote_build
        
  subroutine psb_c_set_vect_default(v)
    implicit none
    class(psb_c_base_vect_type), intent(in) :: v

    if(allocated(psb_c_base_vect_default)) deallocate(psb_c_base_vect_default)
    allocate(psb_c_base_vect_default, mold = v)
  end subroutine psb_c_set_vect_default

  function psb_c_get_vect_default(v) result(res)
    implicit none
    class(psb_c_vect_type), intent(in) :: v
    class(psb_c_base_vect_type), pointer :: res

    res => psb_c_get_base_vect_default()
  end function psb_c_get_vect_default

  subroutine psb_c_clear_vect_default()
    implicit none
    if(allocated(psb_c_base_vect_default)) deallocate(psb_c_base_vect_default)
  end subroutine psb_c_clear_vect_default

  function psb_c_get_base_vect_default() result(res)
    implicit none
    class(psb_c_base_vect_type), pointer :: res

    if(.not. allocated(psb_c_base_vect_default)) &
      & allocate(psb_c_base_vect_type :: psb_c_base_vect_default)

    res => psb_c_base_vect_default
  end function psb_c_get_base_vect_default

  subroutine c_vect_clone(x, y, info)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    call y%free(info)
    !
    ! Using sourced allocation here creates
    ! problems with handling of memory allocated
    ! elsewhere (e.g. accelerators), hence delegation
    ! to %bld method
    ! 
    if((info == psb_success_) .and. allocated(x%v)) call y%bld(x%get_vect(), mold = x%v)
  end subroutine c_vect_clone

  subroutine c_vect_bld_x(x, invect, mold, scratch)
    class(psb_c_vect_type), intent(inout) :: x
    complex(psb_spk_), intent(in)            :: invect(:)
    class(psb_c_base_vect_type), intent(in), optional :: mold
    logical, intent(in), optional                     :: scratch

    logical :: scratch_
    integer(psb_ipk_) :: info

    if(present(scratch)) then
      scratch_ = scratch
    else
      scratch_ = .false.
    end if

    info = psb_success_
    if(allocated(x%v)) call x%free(info)

    if(present(mold)) then
      allocate(x%v, stat = info, mold = mold)
    else
      allocate(x%v, stat = info, mold = psb_c_get_base_vect_default())
    endif

    if(info == psb_success_) call x%v%bld(invect, scratch = scratch_)
  end subroutine c_vect_bld_x

  subroutine c_vect_bld_mn(x, n, mold, scratch)
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_mpk_), intent(in)         :: n
    class(psb_c_base_vect_type), intent(in), optional :: mold
    logical, intent(in), optional                     :: scratch

    logical :: scratch_
    integer(psb_ipk_) :: info
    class(psb_c_base_vect_type), pointer :: mld

    if(present(scratch)) then
      scratch_ = scratch
    else
      scratch_ = .false.
    end if

    info = psb_success_
    if(allocated(x%v)) call x%free(info)

    if(present(mold)) then
      allocate(x%v, stat = info, mold = mold)
    else
      allocate(x%v, stat = info, mold = psb_c_get_base_vect_default())
    endif
    if(info == psb_success_) call x%v%bld(n, scratch = scratch_)
  end subroutine c_vect_bld_mn

  subroutine c_vect_bld_en(x, n, mold, scratch)
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_epk_), intent(in)         :: n
    class(psb_c_base_vect_type), intent(in), optional :: mold
    logical, intent(in), optional                     :: scratch

    logical :: scratch_
    integer(psb_ipk_) :: info

    info = psb_success_
    if(present(scratch)) then
      scratch_ = scratch
    else
      scratch_ = .false.
    end if

    if(allocated(x%v)) call x%free(info)

    if(present(mold)) then
      allocate(x%v, stat = info, mold = mold)
    else
      allocate(x%v, stat = info, mold = psb_c_get_base_vect_default())
    endif
    if(info == psb_success_) call x%v%bld(n, scratch = scratch_)
  end subroutine c_vect_bld_en

  function c_vect_get_vect(x, n) result(res)
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), optional :: n
    complex(psb_spk_), allocatable :: res(:)

    integer(psb_ipk_) :: info

    if(allocated(x%v)) res = x%v%get_vect(n)
  end function c_vect_get_vect

  subroutine c_vect_set_scal(x, val, first, last)
    class(psb_c_vect_type), intent(inout) :: x
    complex(psb_spk_), intent(in)            :: val
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info
    if(allocated(x%v)) call x%v%set(val, first, last)
  end subroutine c_vect_set_scal

  subroutine c_vect_set_vect(x, val, first, last)
    class(psb_c_vect_type), intent(inout) :: x
    complex(psb_spk_), intent(in)            :: val(:)
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info
    if(allocated(x%v)) call x%v%set(val, first, last)
  end subroutine c_vect_set_vect

  subroutine c_vect_check_addr(x)
    class(psb_c_vect_type), intent(inout) :: x

    integer(psb_ipk_) :: info
    if(allocated(x%v)) call x%v%check_addr()
  end subroutine c_vect_check_addr

  function constructor(x) result(this)
    complex(psb_spk_) :: x(:)
    type(psb_c_vect_type) :: this

    integer(psb_ipk_) :: info

    call this%bld(x)
    call this%asb(size(x, kind = psb_ipk_), info)
  end function constructor

  function size_const(n) result(this)
    integer(psb_ipk_), intent(in) :: n
    type(psb_c_vect_type) :: this

    integer(psb_ipk_) :: info

    call this%bld(n)
    call this%asb(n, info)
  end function size_const

  function c_vect_get_nrows(x) result(res)
    implicit none
    class(psb_c_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if(allocated(x%v)) res = x%v%get_nrows()
  end function c_vect_get_nrows

  function c_vect_sizeof(x) result(res)
    implicit none
    class(psb_c_vect_type), intent(in) :: x
    integer(psb_epk_) :: res

    res = 0
    if(allocated(x%v)) res = x%v%sizeof()
  end function c_vect_sizeof

  function c_vect_get_fmt(x) result(res)
    implicit none
    class(psb_c_vect_type), intent(in) :: x
    character(len=5) :: res

    res = 'NULL'
    if(allocated(x%v)) res = x%v%get_fmt()
  end function c_vect_get_fmt

  subroutine c_vect_all(n, x, info, mold)
    implicit none
    integer(psb_ipk_), intent(in)         :: n
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info
    class(psb_c_base_vect_type), intent(in), optional :: mold

    if(allocated(x%v)) call x%free(info)

    if(present(mold)) then
      allocate(x%v, stat = info, mold = mold)
    else
      allocate(psb_c_base_vect_type :: x%v, stat = info)
    endif
    if(info == psb_success_) then
      call x%v%all(n, info)
    else
      info = psb_err_alloc_dealloc_
    end if
    call x%set_bld()
  end subroutine c_vect_all

  subroutine c_vect_reinit(x, info, clear)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)      :: info
    logical, intent(in), optional       :: clear

    if(allocated(x%v)) call x%v%reinit(info, clear)
    call x%set_upd()
  end subroutine c_vect_reinit
 
  subroutine c_vect_reall(n, x, info)
    implicit none
    integer(psb_ipk_), intent(in)         :: n
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    if(.not. allocated(x%v)) call x%all(n, info)
    if(info == psb_success_) call x%asb(n, info)
  end subroutine c_vect_reall

  subroutine c_vect_zero(x)
    use psi_serial_mod
    implicit none
    class(psb_c_vect_type), intent(inout) :: x

    if(allocated(x%v)) call x%v%zero()
  end subroutine c_vect_zero

  subroutine c_vect_asb(n, x, info, scratch)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    integer(psb_ipk_), intent(in)         :: n
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info
    logical, intent(in), optional :: scratch

    if(allocated(x%v)) then
      call x%v%asb(n, info, scratch = scratch)
      call x%set_asb()
    end if
  end subroutine c_vect_asb

  subroutine c_vect_gthab(n, idx, alpha, x, beta, y)
    use psi_serial_mod
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    complex(psb_spk_)    :: alpha, beta, y(:)
    class(psb_c_vect_type) :: x

    if(allocated(x%v)) call x%v%gth(n, idx, alpha, beta, y)
  end subroutine c_vect_gthab

  subroutine c_vect_gthzv(n, idx, x, y)
    use psi_serial_mod
    integer(psb_mpk_)       :: n
    integer(psb_ipk_)       :: idx(:)
    class(psb_c_vect_type)  :: x
    complex(psb_spk_)          :: y(:)

    if(allocated(x%v)) call x%v%gth(n, idx, y)
  end subroutine c_vect_gthzv

  subroutine c_vect_sctb(n, idx, x, beta, y)
    use psi_serial_mod
    integer(psb_mpk_)       :: n
    integer(psb_ipk_)       :: idx(:)
    complex(psb_spk_)          :: beta, x(:)
    class(psb_c_vect_type)  :: y

    if(allocated(y%v)) call y%v%sct(n, idx, x, beta)
  end subroutine c_vect_sctb

  subroutine c_vect_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    if(allocated(x%v)) then
      call x%v%free(info)
      if(info == psb_success_) deallocate(x%v, stat = info)
    end if
  end subroutine c_vect_free

  subroutine c_vect_ins_a(n, irl, val, x, maxr, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)         :: n, maxr
    integer(psb_ipk_), intent(in)         :: irl(:)
    complex(psb_spk_), intent(in)            :: val(:)
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_) :: i, dupl

    info = psb_success_
    if(.not. allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if
    dupl = x%get_dupl()
    call x%v%ins(n, irl, val, dupl, maxr, info)
  end subroutine c_vect_ins_a

  subroutine c_vect_ins_v(n, irl, val, x, maxr, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)         :: n, maxr
    class(psb_i_vect_type), intent(inout) :: irl
    class(psb_c_vect_type), intent(inout) :: val
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_) :: i, dupl

    info = psb_success_
    if(.not. (allocated(x%v) .and. allocated(irl%v) .and. allocated(val%v))) then
      info = psb_err_invalid_vect_state_
      return
    end if
    dupl = x%get_dupl()
    call x%v%ins(n, irl%v, val%v, dupl, maxr, info)
  end subroutine c_vect_ins_v

  subroutine c_vect_cnv(x, mold)
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_base_vect_type), intent(in), optional :: mold

    class(psb_c_base_vect_type), allocatable :: tmp
    integer(psb_ipk_) :: info

    info = psb_success_
    if(present(mold)) then
      allocate(tmp, stat = info, mold = mold)
    else
      allocate(tmp, stat = info, mold = psb_c_get_base_vect_default())
    end if
    
    if(allocated(x%v)) then
      if(allocated(x%v%v)) then 
        call x%v%sync()
        if(info == psb_success_) call tmp%bld(x%v%v)
        call x%v%base_cpy(tmp)
        call x%v%free(info)
      endif
    end if
    call move_alloc(tmp, x%v)
  end subroutine c_vect_cnv

  subroutine c_vect_sync(x)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x

    if(allocated(x%v)) call x%v%sync()
  end subroutine c_vect_sync

  subroutine c_vect_set_sync(x)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x

    if(allocated(x%v)) call x%v%set_sync()
  end subroutine c_vect_set_sync

  subroutine c_vect_set_host(x)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x

    if(allocated(x%v)) call x%v%set_host()
  end subroutine c_vect_set_host

  subroutine c_vect_set_dev(x)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x

    if(allocated(x%v)) call x%v%set_dev()
  end subroutine c_vect_set_dev

  function c_vect_is_sync(x) result(res)
    implicit none
    logical :: res
    class(psb_c_vect_type), intent(inout) :: x

    res = .true.
    if(allocated(x%v)) res = x%v%is_sync()
  end function c_vect_is_sync

  function c_vect_is_host(x) result(res)
    implicit none
    logical :: res
    class(psb_c_vect_type), intent(inout) :: x

    res = .true.
    if(allocated(x%v)) res = x%v%is_host()
  end function c_vect_is_host

  function c_vect_is_dev(x) result(res)
    implicit none
    logical :: res
    class(psb_c_vect_type), intent(inout) :: x

    res = .false.
    if(allocated(x%v)) res = x%v%is_dev()
  end function c_vect_is_dev

  function c_vect_get_entry(x, index) result(res)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)         :: index
    complex(psb_spk_) :: res

    res = czero
    if(allocated(x%v)) res = x%v%get_entry(index)
  end function c_vect_get_entry

  subroutine c_vect_set_entry(x, index, val)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)         :: index
    complex(psb_spk_)                        :: val

    if(allocated(x%v)) call x%v%set_entry(index, val)
  end subroutine c_vect_set_entry

  function c_vect_dot_v(n, x, y) result(res)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(in)           :: n
    complex(psb_spk_)  :: res

    res = czero
    if(allocated(x%v) .and. allocated(y%v)) res = x%v%dot(n, y%v)
  end function c_vect_dot_v

  function c_vect_dot_a(n, x, y) result(res)
    implicit none
    integer(psb_ipk_), intent(in)         :: n
    class(psb_c_vect_type), intent(inout) :: x
    complex(psb_spk_), intent(in)            :: y(:)
    complex(psb_spk_) :: res

    res = czero
    if(allocated(x%v)) res = x%v%dot_a(n, y)
  end function c_vect_dot_a

  subroutine c_vect_axpby_v(m, alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)         :: m
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_vect_type), intent(inout) :: y
    complex(psb_spk_), intent(in)            :: alpha, beta
    integer(psb_ipk_), intent(out)        :: info

    if((.not. allocated(x%v)) .or. (.not. allocated(y%v))) then
      info = psb_err_invalid_vect_state_
      return
    endif
    
    call y%v%axpby(m, alpha, x%v, beta, info)
  end subroutine c_vect_axpby_v

  subroutine c_vect_axpby_v2(m, alpha, x, beta, y, z, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)         :: m
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_vect_type), intent(inout) :: y
    class(psb_c_vect_type), intent(inout) :: z
    complex(psb_spk_), intent(in)            :: alpha, beta
    integer(psb_ipk_), intent(out)        :: info

    if((.not. allocated(x%v)) .or. (.not. allocated(y%v))) then
      info = psb_err_invalid_vect_state_
      return
    endif
    
    call z%v%axpby(m, alpha, x%v, beta, y%v, info)
  end subroutine c_vect_axpby_v2

  subroutine c_vect_axpby_a(m, alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)         :: m
    complex(psb_spk_), intent(in)            :: x(:)
    class(psb_c_vect_type), intent(inout) :: y
    complex(psb_spk_), intent(in)            :: alpha, beta
    integer(psb_ipk_), intent(out)        :: info

    if(allocated(y%v)) call y%v%axpby(m, alpha, x, beta, info)
  end subroutine c_vect_axpby_a

  subroutine c_vect_axpby_a2(m, alpha, x, beta, y, z, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)         :: m
    complex(psb_spk_), intent(in)            :: x(:)
    complex(psb_spk_), intent(in)            :: y(:)
    class(psb_c_vect_type), intent(inout) :: z
    complex(psb_spk_), intent(in)            :: alpha, beta
    integer(psb_ipk_), intent(out)        :: info

    if(allocated(z%v)) call z%v%axpby(m, alpha, x, beta, y, info)
  end subroutine c_vect_axpby_a2

  subroutine c_vect_upd_xyz(m, alpha, beta, gamma, delta, x, y, z, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)         :: m
    complex(psb_spk_), intent(in)            :: alpha, beta, gamma, delta
    class(psb_c_vect_type), intent(inout) :: x, y, z
    integer(psb_ipk_), intent(out)        :: info

    if(allocated(z%v)) call z%v%upd_xyz(m, alpha, beta, gamma, delta, x%v, y%v, info)
  end subroutine c_vect_upd_xyz

  subroutine c_vect_xyzw(m, a, b, c, d, e, f, x, y, z, w, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)         :: m
    complex(psb_spk_), intent(in)            :: a, b, c, d, e, f
    class(psb_c_vect_type), intent(inout) :: x, y, z, w
    integer(psb_ipk_), intent(out)        :: info

    if(allocated(w%v)) call w%v%xyzw(m, a, b, c, d, e, f, x%v, y%v, z%v, info)
  end subroutine c_vect_xyzw

  subroutine c_vect_mlt_v(x, y, info)
    use psi_serial_mod
    implicit none
    class(psb_c_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_) :: i, n

    info = psb_success_
    if(allocated(x%v) .and. allocated(y%v)) call y%v%mlt(x%v, info)
  end subroutine c_vect_mlt_v

  subroutine c_vect_mlt_a(x, y, info)
    use psi_serial_mod
    implicit none
    complex(psb_spk_), intent(in)            :: x(:)
    class(psb_c_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)       :: info

    integer(psb_ipk_) :: i, n

    info = psb_success_
    if(allocated(y%v)) call y%v%mlt(x, info)
  end subroutine c_vect_mlt_a

  subroutine c_vect_mlt_a_2(alpha, x, y, beta, z, info)
    use psi_serial_mod
    implicit none
    complex(psb_spk_), intent(in)            :: alpha, beta
    complex(psb_spk_), intent(in)            :: x(:), y(:)
    class(psb_c_vect_type), intent(inout) :: z
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_) :: i, n

    info = psb_success_
    if(allocated(z%v)) call z%v%mlt(alpha, x, y, beta, info)
  end subroutine c_vect_mlt_a_2

  subroutine c_vect_mlt_v_2(alpha, x, y, beta, z, info, conjgx, conjgy)
    use psi_serial_mod
    implicit none
    complex(psb_spk_), intent(in)            :: alpha, beta
    class(psb_c_vect_type), intent(inout) :: x, y, z
    integer(psb_ipk_), intent(out)       :: info
    character(len=1), intent(in), optional  :: conjgx, conjgy

    integer(psb_ipk_) :: i, n

    info = psb_success_
    if(allocated(x%v) .and. allocated(y%v) .and. allocated(z%v)) &
         & call z%v%mlt(alpha, x%v, y%v, beta, info, conjgx, conjgy)
  end subroutine c_vect_mlt_v_2

  subroutine c_vect_mlt_av(alpha, x, y, beta, z, info)
    use psi_serial_mod
    implicit none
    complex(psb_spk_), intent(in)            :: alpha, beta
    complex(psb_spk_), intent(in)            :: x(:)
    class(psb_c_vect_type), intent(inout) :: y, z
    integer(psb_ipk_), intent(out)        :: info
    integer(psb_ipk_) :: i, n

    info = psb_success_
    if(allocated(z%v) .and. allocated(y%v)) call z%v%mlt(alpha, x, y%v, beta, info)
  end subroutine c_vect_mlt_av

  subroutine c_vect_mlt_va(alpha, x, y, beta, z, info)
    use psi_serial_mod
    implicit none
    complex(psb_spk_), intent(in)            :: alpha, beta
    class(psb_c_vect_type), intent(inout) :: x
    complex(psb_spk_), intent(in)            :: y(:)
    class(psb_c_vect_type), intent(inout) :: z
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_) :: i, n

    info = psb_success_
    if(allocated(z%v) .and. allocated(x%v)) call z%v%mlt(alpha, x%v, y, beta, info)
  end subroutine c_vect_mlt_va

  subroutine c_vect_div_v(x, y, info)
    use psi_serial_mod
    implicit none
    class(psb_c_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_) :: i, n

    info = psb_success_
    if(allocated(x%v) .and. allocated(y%v)) call x%v%div(y%v, info)
  end subroutine c_vect_div_v

  subroutine c_vect_div_v2( x, y, z, info)
    use psi_serial_mod
    implicit none
    class(psb_c_vect_type), intent(inout) :: x, y, z
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_) :: i, n

    info = psb_success_
    if(allocated(x%v) .and. allocated(y%v) .and. allocated(z%v)) &
         & call z%v%div(x%v, y%v, info)
  end subroutine c_vect_div_v2

  subroutine c_vect_div_v_check(x, y, info, flag)
    use psi_serial_mod
    implicit none
    class(psb_c_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(out)        :: info
    integer(psb_ipk_) :: i, n
    logical, intent(in) :: flag

    info = psb_success_
    if(allocated(x%v) .and. allocated(y%v)) call x%v%div(y%v, info, flag)
  end subroutine c_vect_div_v_check

  subroutine c_vect_div_v2_check(x, y, z, info, flag)
    use psi_serial_mod
    implicit none
    class(psb_c_vect_type), intent(inout) :: x, y, z
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_)   :: i, n
    logical, intent(in) :: flag

    info = psb_success_
    if(allocated(x%v) .and. allocated(y%v) .and. allocated(z%v)) &
         & call z%v%div(x%v, y%v, info, flag)
  end subroutine c_vect_div_v2_check

  subroutine c_vect_div_a2(x, y, z, info)
    use psi_serial_mod
    implicit none
    complex(psb_spk_), intent(in)            :: x(:), y(:)
    class(psb_c_vect_type), intent(inout) :: z
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_) :: i, n

    info = psb_success_
    if(allocated(z%v)) call z%v%div(x, y, info)
  end subroutine c_vect_div_a2

  subroutine c_vect_div_a2_check(x, y, z, info, flag)
    use psi_serial_mod
    implicit none
    complex(psb_spk_), intent(in)           :: x(:), y(:)
    class(psb_c_vect_type), intent(inout) :: z
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_) :: i, n
    logical, intent(in) :: flag

    info = psb_success_
    if(allocated(z%v)) call z%v%div(x, y, info, flag)
  end subroutine c_vect_div_a2_check

  subroutine c_vect_inv_v(x, y, info)
    use psi_serial_mod
    implicit none
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_) :: i, n

    info = psb_success_
    if(allocated(x%v) .and. allocated(y%v)) call y%v%inv(x%v, info)
  end subroutine c_vect_inv_v

  subroutine c_vect_inv_v_check(x, y, info, flag)
    use psi_serial_mod
    implicit none
    class(psb_c_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_) :: i, n
    logical, intent(in) :: flag

    info = psb_success_
    if(allocated(x%v) .and. allocated(y%v)) call y%v%inv(x%v, info, flag)
  end subroutine c_vect_inv_v_check

  subroutine c_vect_inv_a2(x, y, info)
    use psi_serial_mod
    implicit none
    complex(psb_spk_), intent(inout)         :: x(:)
    class(psb_c_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_) :: i, n

    info = psb_success_
    if(allocated(y%v)) call y%v%inv(x, info)
  end subroutine c_vect_inv_a2

  subroutine c_vect_inv_a2_check(x, y, info, flag)
    use psi_serial_mod
    implicit none
    complex(psb_spk_), intent(inout)         :: x(:)
    class(psb_c_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_) :: i, n
    logical, intent(in) :: flag

    info = psb_success_
    if(allocated(y%v)) call y%v%inv(x, info, flag)
  end subroutine c_vect_inv_a2_check

  subroutine c_vect_acmp_a2(x, c, z, info)
    use psi_serial_mod
    implicit none
    complex(psb_spk_), intent(inout)         :: x(:)
    real(psb_spk_), intent(in)            :: c
    class(psb_c_vect_type), intent(inout) :: z
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    if(allocated(z%v)) call z%acmp(x, c, info)
  end subroutine c_vect_acmp_a2

  subroutine c_vect_acmp_v2(x, c, z, info)
    use psi_serial_mod
    implicit none
    real(psb_spk_), intent(in)            :: c
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_vect_type), intent(inout) :: z
    integer(psb_ipk_), intent(out)         :: info

    info = psb_success_
    if(allocated(x%v) .and. allocated(z%v)) call z%v%acmp(x%v, c, info)
  end subroutine c_vect_acmp_v2

  subroutine c_vect_scal(alpha, x)
    use psi_serial_mod
    implicit none
    complex(psb_spk_), intent(in)            :: alpha
    class(psb_c_vect_type), intent(inout) :: x

    if(allocated(x%v)) call x%v%scal(alpha)
  end subroutine c_vect_scal

  subroutine c_vect_absval1(x)
    class(psb_c_vect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%absval()
  end subroutine c_vect_absval1

  subroutine c_vect_absval2(x, y)
    class(psb_c_vect_type), intent(inout) :: x, y

    if(allocated(x%v)) then
      if(.not. allocated(y%v)) call y%bld(psb_size(x%v%v))
      call x%v%absval(y%v)
    end if
  end subroutine c_vect_absval2

  function c_vect_nrm2(n, x) result(res)
    implicit none
    integer(psb_ipk_), intent(in)         :: n
    class(psb_c_vect_type), intent(inout) :: x
    real(psb_spk_)  :: res

    res = szero
    if(allocated(x%v)) res = x%v%nrm2(n)
  end function c_vect_nrm2

  function c_vect_nrm2_weight(n, x, w, aux) result(res)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)         :: n
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_vect_type), intent(inout) :: w
    class(psb_c_vect_type), intent(inout), optional :: aux
    real(psb_spk_)  :: res

    integer(psb_ipk_) :: info
    ! Temp vectors
    type(psb_c_vect_type) :: wtemp

    info = psb_success_
    if(allocated(w%v)) then
      if(.not. present(aux)) then
        allocate(wtemp%v, mold = w%v)
        call wtemp%v%bld(w%get_vect())
      else
        call psb_geaxpby(n, cone, w%v%v, czero, aux%v%v, info)
      end if
    else
      info = -1
    end if
    if(info /= psb_success_) then
      res = -sone
      return
    end if

    if(allocated(x%v)) then
      if(.not. present(aux)) then
        call wtemp%v%mlt(x%v, info)
        res = wtemp%v%nrm2(n)
      else
        call aux%v%mlt(x%v, info)
        res = aux%v%nrm2(n)
      end if
    else
      res = szero
    end if

    if(.not. present(aux)) call wtemp%free(info)
  end function c_vect_nrm2_weight
 
  function c_vect_nrm2_weight_mask(n, x, w, id, info, aux) result(res)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)         :: n
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_vect_type), intent(inout) :: w, id
    integer(psb_ipk_), intent(out)        :: info
    class(psb_c_vect_type), intent(inout), optional :: aux
    real(psb_spk_)  :: res

    ! Temp vectors
    type(psb_c_vect_type) :: wtemp

    info = psb_success_ 
    if( allocated(w%v) ) then
      if(.not. present(aux)) then
        allocate(wtemp%v, mold = w%v)
        call wtemp%v%bld(w%get_vect())
      else
        call psb_geaxpby(n, cone, w%v%v, czero, aux%v%v, info)
      end if
    else
      info = -1
    end if
    if(info /= psb_success_) then
      res = -sone
      return
    end if

    if(allocated(x%v) .and. allocated(id%v)) then
      if(.not. present(aux)) then
        where(abs(id%v%v) <= szero) wtemp%v%v = szero
        call wtemp%set_host()
        call wtemp%v%mlt(x%v, info)
        res = wtemp%v%nrm2(n)
      else
        where(abs(id%v%v) <= szero) aux%v%v = szero
        call aux%set_host()
        call aux%v%mlt(x%v, info)
        res = aux%v%nrm2(n)
      end if
    else
      res = szero
    end if

    if(.not. present(aux)) call wtemp%free(info)
  end function c_vect_nrm2_weight_mask

  function c_vect_amax(n, x) result(res)
    implicit none
    integer(psb_ipk_), intent(in)         :: n
    class(psb_c_vect_type), intent(inout) :: x
    real(psb_spk_)  :: res

    res = szero
    if(allocated(x%v)) res = x%v%amax(n)
  end function c_vect_amax


  function c_vect_asum(n, x) result(res)
    implicit none
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)  :: res
    
    res = szero
    if(allocated(x%v)) res = x%v%asum(n)
  end function c_vect_asum


  subroutine c_vect_addconst_a2(x, b, z, info)
    use psi_serial_mod
    implicit none
    complex(psb_spk_), intent(inout)         :: x(:)
    real(psb_spk_), intent(in)           :: b
    class(psb_c_vect_type), intent(inout) :: z
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    if(allocated(z%v)) call z%addconst(x, b, info)
  end subroutine c_vect_addconst_a2

  subroutine c_vect_addconst_v2(x, b, z, info)
    use psi_serial_mod
    implicit none
    class(psb_c_vect_type), intent(inout) :: x
    real(psb_spk_), intent(in)           :: b
    class(psb_c_vect_type), intent(inout) :: z
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    if(allocated(x%v) .and. allocated(z%v)) call z%v%addconst(x%v, b, info)
  end subroutine c_vect_addconst_v2
end module psb_c_vect_mod


module psb_c_multivect_mod
  use psb_c_base_multivect_mod
  use psb_const_mod
  use psb_i_vect_mod

  !private
  type psb_c_multivect_type
    class(psb_c_base_multivect_type), allocatable :: v
    integer(psb_ipk_) :: nrmv = 0
    integer(psb_ipk_) :: remote_build = psb_matbld_noremote_
    complex(psb_spk_), allocatable :: rmtv(:, :)
  contains
    procedure, pass(x) :: get_nrows => c_mvect_get_nrows
    procedure, pass(x) :: get_ncols => c_mvect_get_ncols
    procedure, pass(x) :: sizeof   => c_mvect_sizeof
    procedure, pass(x) :: get_fmt  => c_mvect_get_fmt
    procedure, pass(x) :: is_remote_build => c_mvect_is_remote_build
    procedure, pass(x) :: set_remote_build => c_mvect_set_remote_build

    procedure, pass(x) :: all      => c_mvect_all
    procedure, pass(x) :: reall    => c_mvect_reall
    procedure, pass(x) :: zero     => c_mvect_zero
    procedure, pass(x) :: asb      => c_mvect_asb
    procedure, pass(x) :: sync     => c_mvect_sync
    procedure, pass(x) :: free     => c_mvect_free
    procedure, pass(x) :: reinit   => c_mvect_reinit
    procedure, pass(x) :: set_ncfs => c_mvect_set_ncfs 
    procedure, pass(x) :: get_ncfs => c_mvect_get_ncfs
    procedure, pass(x) :: set_dupl => c_mvect_set_dupl 
    procedure, pass(x) :: get_dupl => c_mvect_get_dupl
    procedure, pass(x) :: set_state => c_mvect_set_state
    procedure, pass(x) :: set_null  => c_mvect_set_null
    procedure, pass(x) :: set_bld   => c_mvect_set_bld
    procedure, pass(x) :: set_upd   => c_mvect_set_upd
    procedure, pass(x) :: set_asb   => c_mvect_set_asb
    procedure, pass(x) :: get_state => c_mvect_get_state
    procedure, pass(x) :: is_null   => c_mvect_is_null
    procedure, pass(x) :: is_bld    => c_mvect_is_bld
    procedure, pass(x) :: is_upd    => c_mvect_is_upd
    procedure, pass(x) :: is_asb    => c_mvect_is_asb
    !!$ procedure, pass(x) :: base_cpy  => c_mvect_cpy

    procedure, pass(x) :: ins      => c_mvect_ins
    procedure, pass(x) :: bld_x    => c_mvect_bld_x
    procedure, pass(x) :: bld_n    => c_mvect_bld_n
    generic, public    :: bld      => bld_x, bld_n
    procedure, pass(x) :: get_vect => c_mvect_get_vect
    procedure, pass(x) :: cnv      => c_mvect_cnv

    procedure, pass(x) :: set_scal => c_mvect_set_scal
    procedure, pass(x) :: set_vect => c_mvect_set_vect
    procedure, pass(x) :: set_colm => c_mvect_set_colm
    generic, public    :: set      => set_vect, set_scal, set_colm

    procedure, pass(x) :: clone    => c_mvect_clone
    procedure, pass(x) :: gthab    => c_mvect_gthab
    procedure, pass(x) :: gthzv    => c_mvect_gthzv
    procedure, pass(x) :: gthzv_x  => c_mvect_gthzv_x
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => c_mvect_sctb
    procedure, pass(y) :: sctb_x   => c_mvect_sctb_x
    generic, public    :: sct      => sctb, sctb_x

    ! single column export as vector
    procedure, pass(x) :: extract_col   => c_mvect_extract_col
    ! two term axpy-like operations
    procedure, pass(y) :: axpby_v_i     => c_mvect_axpby_v_idxs
    procedure, pass(y) :: axpby_v_f     => c_mvect_axpby_v_full
    procedure, pass(y) :: axpby_m_i     => c_mvect_axpby_m_idxs
    procedure, pass(y) :: axpby_m_f     => c_mvect_axpby_m_full
    ! two term axpy-like operations with separate output mv
    procedure, pass(z) :: axpby_m_f_o   => c_mvect_axpby_m_full_out
    ! three term axpy-like operations
    procedure, pass(z) :: axpbycz_vv    => c_mvect_axpbycz_vv
    procedure, pass(z) :: axpbycz_mv    => c_mvect_axpbycz_mv
    procedure, pass(z) :: axpbycz_mm_i  => c_mvect_axpbycz_mm_idxs
    procedure, pass(z) :: axpbycz_mm_f  => c_mvect_axpbycz_mm_full
    ! three term axpy-like operations with separate output mv (rename?)
    procedure, pass(w) :: axpbycz_mm_o  => c_mvect_axpbycz_mm_out
    ! linear combinations of columns of the multivector
    procedure, pass(x) :: colspan1D     => c_mvect_colspan1D
    procedure, pass(x) :: colspan2D     => c_mvect_colspan2D
    ! all procedures exported as axpby
    generic, public    :: axpby         => extract_col, &
                                            axpby_v_i, axpby_v_f, & 
                                            axpby_m_i, axpby_m_f, &
                                            axpby_m_f_o, &
                                            axpbycz_vv, axpbycz_mv, axpbycz_mm_i, &
                                            axpbycz_mm_f, axpbycz_mm_o, & 
                                            colspan1D, colspan2D

    ! dot products operations - only full-full version for now
    procedure, pass(x) :: dot_mm  => c_mvect_dot_mm
    procedure, pass(x) :: dot_mv  => c_mvect_dot_mv
    generic, public    :: dot     => dot_mm, dot_mv

    ! Element wise multiplication operations with in place output
    procedure, pass(y) :: mlt_v_f     => c_mvect_mlt_v_full
    procedure, pass(y) :: mlt_v_i     => c_mvect_mlt_v_idxs
    procedure, pass(y) :: mlt_m_f     => c_mvect_mlt_m_full
    procedure, pass(y) :: mlt_m_i     => c_mvect_mlt_m_idxs
    ! Element wise multiplication operations with separate output
    procedure, pass(z) :: mlt_vv_f_o  => c_mvect_mlt_vv_full_out
    procedure, pass(z) :: mlt_vv_i_o  => c_mvect_mlt_vv_idxs_out
    procedure, pass(z) :: mlt_vm_f_o  => c_mvect_mlt_vm_full_out
    procedure, pass(z) :: mlt_vm_i_o  => c_mvect_mlt_vm_idxs_out
    procedure, pass(z) :: mlt_mm_f_o  => c_mvect_mlt_mm_full_out
    procedure, pass(z) :: mlt_mm_i_o  => c_mvect_mlt_mm_idxs_out
    ! Element wise multiplication operations with externel vector output
    procedure, pass(y) :: mlt_vm_e    => c_mvect_mlt_vm_ext
    procedure, pass(y) :: mlt_mm_e    => c_mvect_mlt_mm_ext
    ! All procedures exported as mlt
    generic, public    :: mlt         => mlt_v_f, mlt_v_i, &
                                          mlt_m_f, mlt_m_i, &
                                          mlt_vv_f_o, mlt_vv_i_o, &
                                          mlt_vm_f_o, mlt_vm_i_o, &
                                          mlt_mm_f_o, mlt_mm_i_o, &
                                          mlt_vm_e, mlt_mm_e
    !
    ! Scaling and norms
    !
    procedure, pass(x)  :: nrm2_f => c_mvect_nrm2_full
    procedure, pass(x)  :: nrm2_i => c_mvect_nrm2_idxs
    generic, public     :: nrm2   => nrm2_f, nrm2_i
    
    !!$ procedure, pass(x) :: scal     => c_mvect_scal
    !!$ procedure, pass(x) :: amax     => c_mvect_amax
    !!$ procedure, pass(x) :: asum     => c_mvect_asum
  end type psb_c_multivect_type

  public  :: psb_c_multivect, psb_c_multivect_type, psb_c_base_multivect_type, &
          & psb_set_multivect_default, psb_get_multivect_default

  private
  interface psb_c_multivect
    module procedure constructor, size_const
  end interface psb_c_multivect

  class(psb_c_base_multivect_type), allocatable, target, &
       & save, private :: psb_c_base_multivect_default

  interface psb_set_multivect_default
    module procedure psb_c_set_multivect_default
  end interface psb_set_multivect_default

  interface psb_get_multivect_default
    module procedure psb_c_get_multivect_default
  end interface psb_get_multivect_default

contains
  function c_mvect_is_remote_build(x) result(res)
    implicit none
    class(psb_c_multivect_type), intent(in) :: x
    logical :: res
    
    res = (x%remote_build == psb_matbld_remote_)
  end function c_mvect_is_remote_build

  subroutine c_mvect_set_remote_build(x, val)
    implicit none
    class(psb_c_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if(present(val)) then
      x%remote_build = val
    else
      x%remote_build = psb_matbld_remote_
    end if
  end subroutine c_mvect_set_remote_build

  subroutine psb_c_set_multivect_default(v)
    implicit none
    class(psb_c_base_multivect_type), intent(in) :: v

    if(allocated(psb_c_base_multivect_default)) then
      deallocate(psb_c_base_multivect_default)
    end if
    allocate(psb_c_base_multivect_default, mold = v)
  end subroutine psb_c_set_multivect_default

  function psb_c_get_multivect_default(v) result(res)
    implicit none
    class(psb_c_multivect_type), intent(in) :: v
    class(psb_c_base_multivect_type), pointer :: res

    res => psb_c_get_base_multivect_default()
  end function psb_c_get_multivect_default

  function psb_c_get_base_multivect_default() result(res)
    implicit none
    class(psb_c_base_multivect_type), pointer :: res

    if(.not. allocated(psb_c_base_multivect_default)) then
      allocate(psb_c_base_multivect_type :: psb_c_base_multivect_default)
    end if

    res => psb_c_base_multivect_default
  end function psb_c_get_base_multivect_default

  subroutine c_mvect_clone(x, y, info)
    implicit none
    class(psb_c_multivect_type), intent(inout)  :: x, y
    integer(psb_ipk_), intent(out)                :: info

    info = psb_success_
    call y%free(info)
    if((info == psb_success_) .and. allocated(x%v)) then
      call y%bld_x(x%get_vect(), mold = x%v)
    end if
  end subroutine c_mvect_clone

  subroutine c_mvect_bld_x(x, invect, mold)
    class(psb_c_multivect_type), intent(out)  :: x
    complex(psb_spk_), intent(in)                :: invect(:, :)
    class(psb_c_base_multivect_type), intent(in), optional :: mold

    integer(psb_ipk_) :: info
    info = psb_success_

    if(present(mold)) then
      allocate(x%v, stat = info, mold = mold)
    else
      allocate(x%v, stat = info, mold = psb_c_get_base_multivect_default())
    endif

    if(info == psb_success_) call x%v%bld(invect)
  end subroutine c_mvect_bld_x

  subroutine c_mvect_bld_n(x, m, n, mold, scratch)
    class(psb_c_multivect_type), intent(out)  :: x
    integer(psb_ipk_), intent(in)             :: m, n
    class(psb_c_base_multivect_type), intent(in), optional  :: mold
    logical, intent(in), optional                           :: scratch

    integer(psb_ipk_) :: info
    info = psb_success_

    if(present(mold)) then
      allocate(x%v, stat = info, mold = mold)
    else
      allocate(x%v, stat = info, mold = psb_c_get_base_multivect_default())
    endif

    if(info == psb_success_) call x%v%bld(m, n, scratch = scratch)
  end subroutine c_mvect_bld_n

  function c_mvect_get_vect(x) result(res)
    class(psb_c_multivect_type), intent(inout)  :: x
    complex(psb_spk_), allocatable                 :: res(:, :)
    integer(psb_ipk_) :: info

    if(allocated(x%v)) then
      res = x%v%get_vect()
    end if
  end function c_mvect_get_vect

  subroutine c_mvect_set_scal(x, val, rfirst, rlast)
    class(psb_c_multivect_type), intent(inout)  :: x
    complex(psb_spk_), intent(in)                  :: val
    integer(psb_ipk_), optional :: rfirst, rlast

    integer(psb_ipk_) :: info
    if(allocated(x%v)) call x%v%set(val, rfirst, rlast)
  end subroutine c_mvect_set_scal

  subroutine c_mvect_set_vect(x, val)
    class(psb_c_multivect_type), intent(inout)  :: x
    complex(psb_spk_), intent(in)                  :: val(:, :)

    integer(psb_ipk_) :: info
    if(allocated(x%v)) call x%v%set(val)
  end subroutine c_mvect_set_vect

  subroutine c_mvect_set_colm(x, cidx, val, rfirst, rlast)
    class(psb_c_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: cidx
    complex(psb_spk_), intent(in)                  :: val
    integer(psb_ipk_), optional :: rfirst, rlast

    integer(psb_ipk_) :: info
    if(allocated(x%v)) call x%v%set(cidx, val, rfirst, rlast)
  end subroutine c_mvect_set_colm

  function constructor(x) result(this)
    complex(psb_spk_)  :: x(:, :)
    type(psb_c_multivect_type)  :: this

    integer(psb_ipk_) :: info

    call this%bld_x(x)
    call this%asb(size(x, dim=1, kind=psb_ipk_), size(x, dim=2, kind=psb_ipk_), info)
  end function constructor

  function size_const(m, n) result(this)
    integer(psb_ipk_), intent(in) :: m, n
    type(psb_c_multivect_type)  :: this

    integer(psb_ipk_) :: info

    call this%bld_n(m, n)
    call this%asb(m, n, info)
  end function size_const

  function c_mvect_get_nrows(x) result(res)
    implicit none
    class(psb_c_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if(allocated(x%v)) res = x%v%get_nrows()
  end function c_mvect_get_nrows

  function c_mvect_get_ncols(x) result(res)
    implicit none
    class(psb_c_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if(allocated(x%v)) res = x%v%get_ncols()
  end function c_mvect_get_ncols

  function c_mvect_sizeof(x) result(res)
    implicit none
    class(psb_c_multivect_type), intent(in) :: x
    integer(psb_epk_) :: res

    res = 0
    if(allocated(x%v)) res = x%v%sizeof()
  end function c_mvect_sizeof

  function c_mvect_get_fmt(x) result(res)
    implicit none
    class(psb_c_multivect_type), intent(in) :: x
    character(len=5) :: res

    res = 'NULL'
    if(allocated(x%v)) res = x%v%get_fmt()
  end function c_mvect_get_fmt

  subroutine c_mvect_all(m, n, x, info, mold)
    implicit none
    integer(psb_ipk_), intent(in)             :: m, n
    class(psb_c_multivect_type), intent(out)  :: x
    integer(psb_ipk_), intent(out)            :: info
    class(psb_c_base_multivect_type), intent(in), optional :: mold

    if(present(mold)) then
      allocate(x%v, stat = info, mold = mold)
    else
      allocate(psb_c_base_multivect_type :: x%v, stat = info)
    endif
    if(info == psb_success_) then
      call x%v%all(m, n, info)
    else
      info = psb_err_alloc_dealloc_
    end if
    call x%set_bld()
  end subroutine c_mvect_all

  subroutine c_mvect_reall(m, n, x, info)
    implicit none
    integer(psb_ipk_), intent(in)               :: m, n
    class(psb_c_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = psb_success_
    if(.not. allocated(x%v)) call x%all(m, n, info)
    if(info == psb_success_) call x%asb(m, n, info)
  end subroutine c_mvect_reall

  subroutine c_mvect_reinit(x, info)
    implicit none
    class(psb_c_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = psb_success_
    if(allocated(x%v)) call x%v%reinit(info)
    call x%set_upd()
  end subroutine c_mvect_reinit

  subroutine c_mvect_zero(x)
    use psi_serial_mod
    implicit none
    class(psb_c_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%zero()
  end subroutine c_mvect_zero

  subroutine c_mvect_asb(m, n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m, n
    class(psb_c_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    if(allocated(x%v)) then
      call x%v%asb(m, n, info)
      call x%set_asb()
    end if
  end subroutine c_mvect_asb

  subroutine c_mvect_sync(x)
    implicit none
    class(psb_c_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%sync()
  end subroutine c_mvect_sync

  subroutine c_mvect_gthab(n, idx, alpha, x, beta, y)
    use psi_serial_mod
    integer(psb_mpk_)           :: n
    integer(psb_ipk_)           :: idx(:)
    complex(psb_spk_)              :: alpha, beta, y(:)
    class(psb_c_multivect_type) :: x

    if(allocated(x%v)) call x%v%gth(n, idx, alpha, beta, y)
  end subroutine c_mvect_gthab

  subroutine c_mvect_gthzv(n, idx, x, y)
    use psi_serial_mod
    integer(psb_mpk_)           :: n
    integer(psb_ipk_)           :: idx(:)
    complex(psb_spk_)              :: y(:)
    class(psb_c_multivect_type) :: x

    if(allocated(x%v)) call x%v%gth(n, idx, y)
  end subroutine c_mvect_gthzv

  subroutine c_mvect_gthzv_x(i, n, idx, x, y)
    use psi_serial_mod
    integer(psb_ipk_)           :: i
    integer(psb_mpk_)           :: n
    class(psb_i_base_vect_type) :: idx
    class(psb_c_multivect_type) :: x
    complex(psb_spk_)              :: y(:)

    if(allocated(x%v)) call x%v%gth(i, n, idx, y)
  end subroutine c_mvect_gthzv_x

  subroutine c_mvect_sctb(n, idx, x, beta, y)
    use psi_serial_mod
    integer(psb_mpk_)           :: n
    integer(psb_ipk_)           :: idx(:)
    complex(psb_spk_)              :: beta, x(:)
    class(psb_c_multivect_type) :: y

    if(allocated(y%v)) call y%v%sct(n, idx, x, beta)
  end subroutine c_mvect_sctb

  subroutine c_mvect_sctb_x(i, n, idx, x, beta, y)
    use psi_serial_mod
    integer(psb_ipk_)           :: i
    integer(psb_mpk_)           :: n
    class(psb_i_base_vect_type) :: idx
    complex(psb_spk_)              :: beta, x(:)
    class(psb_c_multivect_type) :: y

    if(allocated(y%v)) call y%v%sct(i, n, idx, x, beta)
  end subroutine c_mvect_sctb_x

  subroutine c_mvect_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    class(psb_c_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = psb_success_
    if(allocated(x%v)) then
      call x%v%free(info)
      if(info == psb_success_) deallocate(x%v, stat = info)
    end if
  end subroutine c_mvect_free

  subroutine c_mvect_set_ncfs(n, x)
    integer(psb_ipk_)                           :: n
    class(psb_c_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%set_ncfs(n)
  end subroutine c_mvect_set_ncfs

  function c_mvect_get_ncfs(n, x) result(res)
    class(psb_c_multivect_type), intent(inout)  :: x
    integer(psb_ipk_) :: res

    if(allocated(x%v)) res = x%v%get_ncfs()
  end function c_mvect_get_ncfs

  subroutine c_mvect_set_dupl(n, x)
    integer(psb_ipk_)                           :: n
    class(psb_c_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%set_dupl(n)
  end subroutine c_mvect_set_dupl

  function c_mvect_get_dupl(x) result(res)
    class(psb_c_multivect_type), intent(inout)  :: x
    integer(psb_ipk_) :: res

    if(allocated(x%v)) res = x%v%get_dupl()
  end function c_mvect_get_dupl

  subroutine c_mvect_set_state(n, x)
    integer(psb_ipk_)                           :: n
    class(psb_c_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%set_state(n)
  end subroutine c_mvect_set_state

  function c_mvect_get_state(n, x) result(res)
    class(psb_c_multivect_type), intent(inout)  :: x
    integer(psb_ipk_) :: res

    if(allocated(x%v)) res = x%v%get_state()
  end function c_mvect_get_state

  subroutine c_mvect_set_null(x)
    class(psb_c_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%set_null()
  end subroutine c_mvect_set_null

  function c_mvect_is_null(x) result(res)
    class(psb_c_multivect_type), intent(inout)  :: x
    logical :: res

    res = .false.
    if(allocated(x%v)) res = x%v%is_null()
  end function c_mvect_is_null

  subroutine c_mvect_set_bld(x)
    class(psb_c_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%set_bld()
  end subroutine c_mvect_set_bld

  function c_mvect_is_bld(x) result(res)
    class(psb_c_multivect_type), intent(inout)  :: x
    logical :: res

    res = .false.
    if(allocated(x%v)) res = x%v%is_bld()
  end function c_mvect_is_bld

  subroutine c_mvect_set_upd(x)
    class(psb_c_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%set_upd()
  end subroutine c_mvect_set_upd

  function c_mvect_is_upd(x) result(res)
    class(psb_c_multivect_type), intent(inout)  :: x
    logical :: res

    res = .false.
    if(allocated(x%v)) res = x%v%is_upd()
  end function c_mvect_is_upd

  subroutine c_mvect_set_asb(x)
    class(psb_c_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%set_asb()
  end subroutine c_mvect_set_asb

  function c_mvect_is_asb(x) result(res)
    class(psb_c_multivect_type), intent(inout)  :: x
    logical :: res

    res = .false.
    if(allocated(x%v)) res = x%v%is_asb()
  end function c_mvect_is_asb
  
  subroutine c_mvect_ins(n, irl, val, x, maxr, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: n, maxr
    integer(psb_ipk_), intent(in)               :: irl(:)
    complex(psb_spk_), intent(in)                  :: val(:, :)
    class(psb_c_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i, dupl

    info = psb_success_
    if(.not. allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if
    dupl = x%get_dupl()
    call x%v%ins(n, irl, val, dupl, maxr, info)
  end subroutine c_mvect_ins

  subroutine c_mvect_cnv(x, mold)
    class(psb_c_multivect_type), intent(inout)  :: x
    class(psb_c_base_multivect_type), intent(in), optional :: mold

    class(psb_c_base_multivect_type), allocatable :: tmp
    integer(psb_ipk_) :: info

    if(present(mold)) then
      allocate(tmp, stat = info, mold = mold)
    else
      allocate(tmp, stat = info, mold = psb_c_get_base_multivect_default())
    endif

    if(allocated(x%v)) then
      call x%v%sync()
      if(info == psb_success_) call tmp%bld(x%v%v)
      call x%v%free(info)
    end if

    call move_alloc(tmp, x%v)
  end subroutine c_mvect_cnv
  subroutine c_mvect_extract_col(m, alpha, x, idx_x, beta, y, info)
    use psi_serial_mod
    use psb_c_vect_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m, idx_x
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_multivect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)       :: y
    integer(psb_ipk_), intent(out)              :: info

    if(.not. allocated(x%v)) then
      info = psb_err_invalid_mvect_state_
      return
    endif

    if(.not. allocated(y%v)) then
      info = psb_err_invalid_vect_state_
      return
    endif

    if(idx_x < 0) then
      info = psb_err_iarg_neg_
      return
    endif

    if(idx_x > x%get_ncols()) then
      info = psb_err_entry_out_of_bounds_
      return
    endif
    
    call x%v%axpby(m, alpha, idx_x, beta, y%v, info)
  end subroutine c_mvect_extract_col

  subroutine c_mvect_axpby_v_idxs(m, alpha, x, beta, y, idx_y, info)
    use psi_serial_mod
    use psb_c_vect_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m, idx_y
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_vect_type), intent(inout)       :: x
    class(psb_c_multivect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info

    if(.not. allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    endif

    if(.not. allocated(y%v)) then
      info = psb_err_invalid_mvect_state_
      return
    endif

    if(idx_y < 0) then
      info = psb_err_iarg_neg_
      return
    endif

    if(idx_y > y%get_ncols()) then
      info = psb_err_entry_out_of_bounds_
      return
    endif
    
    call y%v%axpby(m, alpha, x%v, beta, idx_y, info)
  end subroutine c_mvect_axpby_v_idxs
  
  subroutine c_mvect_axpby_v_full(m, alpha, x, beta, y, info)
    use psi_serial_mod
    use psb_c_vect_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_vect_type), intent(inout)       :: x
    class(psb_c_multivect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info

    if(.not. allocated(x%v)) then
      info = psb_err_invalid_vect_state_ 
      return
    endif

    if(.not. allocated(y%v)) then
      info = psb_err_invalid_mvect_state_ 
      return
    endif
    
    call y%v%axpby(m, alpha, x%v, beta, info)
  end subroutine c_mvect_axpby_v_full
  
  subroutine c_mvect_axpby_m_idxs(m, alpha, x, idx_x, beta, y, idx_y, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m, idx_x, idx_y
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_multivect_type), intent(inout)  :: x, y
    integer(psb_ipk_), intent(out)              :: info

    if((.not. allocated(x%v)) .or. (.not. allocated(y%v))) then
      info = psb_err_invalid_mvect_state_  
      return
    endif

    if(idx_y <= 0 .or. idx_x <= 0) then
      info = psb_err_iarg_neg_
      return
    endif

    if(idx_y > y%get_ncols() .or. idx_x > x%get_ncols()) then
      info = psb_err_entry_out_of_bounds_
      return
    endif
    
    call y%v%axpby(m, alpha, x%v, idx_x, beta, idx_y, info)
  end subroutine c_mvect_axpby_m_idxs
  
  subroutine c_mvect_axpby_m_full(m, alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_multivect_type), intent(inout)  :: x, y
    integer(psb_ipk_), intent(out)              :: info

    if((.not. allocated(x%v)) .or. (.not. allocated(y%v))) then
      info = psb_err_invalid_mvect_state_
      return
    endif

    ! Multivector with different size rise error
    if(x%get_ncols() /= y%get_ncols()) then
      info = psb_err_invalid_mvect_size_
      return
    endif
    
    call y%v%axpby(m, alpha, x%v, beta, info)
  end subroutine c_mvect_axpby_m_full

  subroutine c_mvect_axpby_m_full_out(m, alpha, x, beta, y, z, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_multivect_type), intent(inout)  :: x, y, z
    integer(psb_ipk_), intent(out)              :: info

    if((.not. allocated(x%v)) .or. (.not. allocated(y%v)) .or. (.not. allocated(y%v))) then
      info = psb_err_invalid_mvect_state_
      return
    endif

    ! Multivector with different size rise error
    if((z%get_ncols() /= x%get_ncols()) .or. (z%get_ncols() /= y%get_ncols())) then
      info = psb_err_invalid_mvect_size_
      return
    endif
    
    call z%v%axpby(m, alpha, x%v, beta, y%v, info)
  end subroutine c_mvect_axpby_m_full_out
  
  subroutine c_mvect_axpbycz_vv(m, alpha, x, beta, y, gamma, z, idx_z, info)
    use psi_serial_mod
    use psb_c_vect_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m, idx_z
    complex(psb_spk_), intent(in)                  :: alpha, beta, gamma
    class(psb_c_vect_type), intent(inout)       :: x, y
    class(psb_c_multivect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info

    if((.not. allocated(x%v)) .or. (.not. allocated(y%v))) then
      info = psb_err_invalid_vect_state_
      return
    endif

    if(.not. allocated(z%v)) then
      info = psb_err_invalid_mvect_state_
      return
    endif

    if(idx_z < 0) then
      info = psb_err_iarg_neg_
      return
    endif

    if(idx_z > z%get_ncols()) then
      info = psb_err_entry_out_of_bounds_
      return
    endif
    
    call z%v%axpby(m, alpha, x%v, beta, y%v, gamma, idx_z, info)
  end subroutine c_mvect_axpbycz_vv

  subroutine c_mvect_axpbycz_mv(m, alpha, x, beta, y, idx_y, gamma, z, idx_z, info)
    use psi_serial_mod
    use psb_c_vect_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m, idx_y, idx_z
    complex(psb_spk_), intent(in)                  :: alpha, beta, gamma
    class(psb_c_vect_type), intent(inout)       :: x
    class(psb_c_multivect_type), intent(inout)  :: y, z
    integer(psb_ipk_), intent(out)              :: info

    if(.not. allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    endif

    if((.not. allocated(y%v)) .or. (.not. allocated(z%v))) then
      info = psb_err_invalid_mvect_state_
      return
    endif

    if(idx_y <= 0 .or. idx_z <= 0) then
      info = psb_err_iarg_neg_
      return
    endif

    if(idx_y > y%get_ncols() .or. idx_z > z%get_ncols()) then
      info = psb_err_entry_out_of_bounds_
      return
    endif
    
    call z%v%axpby(m, alpha, x%v, beta, y%v, idx_y, gamma, idx_z, info)
  end subroutine c_mvect_axpbycz_mv

  subroutine c_mvect_axpbycz_mm_idxs(m, alpha, x, idx_x, beta, y, idx_y, gamma, z, idx_z, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m, idx_x, idx_y, idx_z
    complex(psb_spk_), intent(in)                  :: alpha, beta, gamma
    class(psb_c_multivect_type), intent(inout)  :: x, y, z
    integer(psb_ipk_), intent(out)              :: info

    if((.not. allocated(x%v)) .or. (.not. allocated(y%v)) .or. (.not. allocated(z%v))) then
      info = psb_err_invalid_mvect_state_
      return
    endif

    if(idx_x <= 0 .or. idx_y <= 0 .or. idx_z <= 0) then
      info = psb_err_iarg_neg_
      return
    endif

    if(idx_x > x%get_ncols() .or. idx_y > y%get_ncols() .or. idx_z > z%get_ncols()) then
      info = psb_err_entry_out_of_bounds_
      return
    endif
    
    call z%v%axpby(m, alpha, x%v, idx_x, beta, y%v, idx_y, gamma, idx_z, info)
  end subroutine c_mvect_axpbycz_mm_idxs

  subroutine c_mvect_axpbycz_mm_full(m, alpha, x, beta, y, gamma, z, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    class(psb_c_multivect_type), intent(inout)  :: x, y, z
    complex(psb_spk_), intent(in)                  :: alpha, beta, gamma
    integer(psb_ipk_), intent(out)              :: info

    if((.not. allocated(x%v)) .or. (.not. allocated(y%v)) .or. (.not. allocated(z%v))) then
      info = psb_err_invalid_mvect_state_
      return
    endif

    if((x%get_ncols() /= y%get_ncols()) .or. (x%get_ncols() /= z%get_ncols())) then
      info = psb_err_invalid_mvect_size_
      return
    endif
    
    call z%v%axpby(m, alpha, x%v, beta, y%v, gamma, info)
  end subroutine c_mvect_axpbycz_mm_full
  
  subroutine c_mvect_axpbycz_mm_out(m, alpha, x, idx_x, beta, y, idx_y, gamma, z, idx_z, w, idx_w, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m, idx_x, idx_y, idx_z, idx_w
    complex(psb_spk_), intent(in)                  :: alpha, beta, gamma
    class(psb_c_multivect_type), intent(inout)  :: x, y, z, w
    integer(psb_ipk_), intent(out)              :: info

    if((.not. allocated(x%v)) .or. (.not. allocated(y%v)) & 
        .or. (.not. allocated(z%v)).or. (.not. allocated(w%v))) then
      info = psb_err_invalid_mvect_state_
      return
    endif

    if((idx_x <= 0) .or. (idx_y <= 0) .or. (idx_z <= 0) .or. (idx_w <= 0)) then
      info = psb_err_iarg_neg_
      return
    endif

    if((idx_x > x%get_ncols()) .or. (idx_y > y%get_ncols()) & 
          .or. (idx_z > z%get_ncols()) .or. (idx_w > w%get_ncols())) then
      info = psb_err_entry_out_of_bounds_
      return
    endif
    
    call w%v%axpby(m, alpha, x%v, idx_x, beta, y%v, idx_y, gamma, z%v, idx_z, idx_w, info)
  end subroutine c_mvect_axpbycz_mm_out
  
  subroutine c_mvect_colspan1D(m, x, coeff, y, info, upd_flag)
    use psi_serial_mod
    use psb_c_vect_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    class(psb_c_multivect_type), intent(inout)  :: x
    complex(psb_spk_), intent(in)                  :: coeff(:)
    class(psb_c_vect_type), intent(inout)       :: y 
    integer(psb_ipk_), intent(out)              :: info
    logical, intent(in)                         :: upd_flag

    if(.not. allocated(x%v)) then
      info = psb_err_invalid_mvect_state_
      return
    endif

    if(.not. allocated(y%v)) then
      info = psb_err_invalid_vect_state_
      return
    endif

    call x%v%axpby(m, coeff, y%v, info, upd_flag)
  end subroutine c_mvect_colspan1D

  subroutine c_mvect_colspan2D(m, x, coeff, y, info, upd_flag)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    class(psb_c_multivect_type), intent(inout)  :: x, y
    complex(psb_spk_), intent(in)                  :: coeff(:, :)
    integer(psb_ipk_), intent(out)              :: info
    logical, intent(in)                         :: upd_flag

    if((.not. allocated(x%v)) .or. (.not. allocated(y%v))) then
      info = psb_err_invalid_mvect_state_
      return
    endif

    call x%v%axpby(m, coeff, y%v, info, upd_flag)
  end subroutine c_mvect_colspan2D

  subroutine c_mvect_dot_mm(m, x, y, res, info)
    implicit none
    integer(psb_ipk_), intent(in)              :: m
    class(psb_c_multivect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(out)             :: info
    complex(psb_spk_), intent(out)                :: res(:, :)

    if((.not. allocated(x%v)) .or. (.not. allocated(y%v))) then
      info = psb_err_invalid_mvect_state_
      return
    endif

    call x%v%dotsbr(m, y%v, res, info)
  end subroutine c_mvect_dot_mm

  subroutine c_mvect_dot_mv(m, x, y, res, info)
    use psb_c_vect_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    class(psb_c_multivect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)       :: y
    integer(psb_ipk_), intent(out)              :: info
    complex(psb_spk_), intent(out)                 :: res(:)

    if((.not. allocated(x%v)) .or. (.not. allocated(y%v))) then
      info = psb_err_invalid_vect_state_
      return
    endif

    call x%v%dotsbr(m, y%v, res, info)
  end subroutine c_mvect_dot_mv

  subroutine c_mvect_mlt_v_full(m, alpha, x, y, beta, info, conjgx, conjgy)
    use psi_serial_mod
    use psb_c_vect_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_vect_type), intent(inout)       :: x
    class(psb_c_multivect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    character(len=1), intent(in), optional  :: conjgx, conjgy   !TO DO: remove from real cases?

    info = psb_success_
    if(.not. allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if

    if(.not. allocated(y%v)) then
      info = psb_err_invalid_mvect_state_
      return
    end if

    call y%v%mlt(m, alpha, x%v, beta, info, conjgx, conjgy)
  end subroutine c_mvect_mlt_v_full

  subroutine c_mvect_mlt_v_idxs(m, alpha, x, y, idx_y, beta, info, conjgx, conjgy)
    use psi_serial_mod
    use psb_c_vect_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_vect_type), intent(inout)       :: x
    class(psb_c_multivect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(in)               :: idx_y
    integer(psb_ipk_), intent(out)              :: info
    character(len=1), intent(in), optional  :: conjgx, conjgy   !TO DO: remove from real cases?

    info = psb_success_
    if(.not. allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if

    if(.not. allocated(y%v)) then
      info = psb_err_invalid_mvect_state_
      return
    end if
    
    if(idx_y < 0) then
      info = psb_err_iarg_neg_
      return
    endif

    if(idx_y > y%get_ncols()) then
      info = psb_err_entry_out_of_bounds_
      return
    endif

    call y%v%mlt(m, alpha, x%v, idx_y, beta, info, conjgx, conjgy)
  end subroutine c_mvect_mlt_v_idxs

  subroutine c_mvect_mlt_m_full(m, alpha, x, y, beta, info, conjgx, conjgy)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_multivect_type), intent(inout)  :: x, y
    integer(psb_ipk_), intent(out)              :: info
    character(len=1), intent(in), optional  :: conjgx, conjgy   !TO DO: remove from real cases?

    info = psb_success_
    if((.not. allocated(x%v)) .or. (.not. allocated(y%v))) then
      info = psb_err_invalid_mvect_state_
      return
    end if

    ! Multivector with different size rise error
    if(x%get_ncols() /= y%get_ncols()) then
      info = psb_err_invalid_mvect_size_
      return
    endif

    call y%v%mlt(m, alpha, x%v, beta, info, conjgx, conjgy)
  end subroutine c_mvect_mlt_m_full

  subroutine c_mvect_mlt_m_idxs(m, alpha, x, idx_x, y, idx_y, beta, info, conjgx, conjgy)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_multivect_type), intent(inout)  :: x, y
    integer(psb_ipk_), intent(in)               :: idx_x, idx_y
    integer(psb_ipk_), intent(out)              :: info
    character(len=1), intent(in), optional  :: conjgx, conjgy   !TO DO: remove from real cases?

    info = psb_success_
    if((.not. allocated(x%v)) .or. (.not. allocated(y%v))) then
      info = psb_err_invalid_mvect_state_
      return
    end if

    if((idx_y < 0) .or. (idx_x < 0)) then
      info = psb_err_iarg_neg_
      return
    endif

    if((idx_y > y%get_ncols()) .or. (idx_x > x%get_ncols())) then
      info = psb_err_entry_out_of_bounds_
      return
    endif

    call y%v%mlt(m, alpha, x%v, idx_x, idx_y, beta, info, conjgx, conjgy)
  end subroutine c_mvect_mlt_m_idxs

  subroutine c_mvect_mlt_vv_full_out(m, alpha, x, y, beta, z, info, conjgx, conjgy)
    use psi_serial_mod
    use psb_c_vect_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_vect_type), intent(inout)       :: x, y
    class(psb_c_multivect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info
    character(len=1), intent(in), optional  :: conjgx, conjgy   !TO DO: remove from real cases?

    info = psb_success_
    if(.not. allocated(x%v) .or. (.not. allocated(y%v))) then
      info = psb_err_invalid_vect_state_
      return
    end if

    if(.not. allocated(z%v)) then
      info = psb_err_invalid_mvect_state_
      return
    end if

    call z%v%mlt(m, alpha, x%v, y%v, beta, info, conjgx, conjgy)
  end subroutine c_mvect_mlt_vv_full_out

  subroutine c_mvect_mlt_vv_idxs_out(m, alpha, x, y, beta, z, idx_z, info, conjgx, conjgy)
    use psi_serial_mod
    use psb_c_vect_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_vect_type), intent(inout)       :: x, y
    class(psb_c_multivect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(in)               :: idx_z
    integer(psb_ipk_), intent(out)              :: info
    character(len=1), intent(in), optional  :: conjgx, conjgy   !TO DO: remove from real cases?

    info = psb_success_
    if(.not. allocated(x%v) .or. (.not. allocated(y%v))) then
      info = psb_err_invalid_vect_state_
      return
    end if

    if(.not. allocated(z%v)) then
      info = psb_err_invalid_mvect_state_
      return
    end if

    if(idx_z < 0) then
      info = psb_err_iarg_neg_
      return
    endif

    if(idx_z > z%get_ncols()) then
      info = psb_err_entry_out_of_bounds_
      return
    endif

    call z%v%mlt(m, alpha, x%v, y%v, beta, idx_z, info, conjgx, conjgy)
  end subroutine c_mvect_mlt_vv_idxs_out

  subroutine c_mvect_mlt_vm_full_out(m, alpha, x, y, beta, z, info, conjgx, conjgy)
    use psi_serial_mod
    use psb_c_vect_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_vect_type), intent(inout)       :: x
    class(psb_c_multivect_type), intent(inout)  :: y, z
    integer(psb_ipk_), intent(out)              :: info
    character(len=1), intent(in), optional  :: conjgx, conjgy   !TO DO: remove from real cases?

    info = psb_success_
    if(.not. allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if

    if(.not. allocated(y%v) .or. (.not. allocated(z%v))) then
      info = psb_err_invalid_mvect_state_
      return
    end if

    call z%v%mlt(m, alpha, x%v, y%v, beta, info, conjgx, conjgy)
  end subroutine c_mvect_mlt_vm_full_out

  subroutine c_mvect_mlt_vm_idxs_out(m, alpha, x, y, idx_y, beta, z, idx_z, info, conjgx, conjgy)
    use psi_serial_mod
    use psb_c_vect_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_vect_type), intent(inout)       :: x
    class(psb_c_multivect_type), intent(inout)  :: y, z
    integer(psb_ipk_), intent(in)               :: idx_y, idx_z
    integer(psb_ipk_), intent(out)              :: info
    character(len=1), intent(in), optional  :: conjgx, conjgy   !TO DO: remove from real cases?

    info = psb_success_
    if(.not. allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if

    if((.not. allocated(y%v)).or. (.not. allocated(z%v))) then
      info = psb_err_invalid_mvect_state_
      return
    end if
    
    if((idx_y < 0) .or. (idx_z < 0)) then
      info = psb_err_iarg_neg_
      return
    endif

    if((idx_y > y%get_ncols()) .or. (idx_z > z%get_ncols())) then
      info = psb_err_entry_out_of_bounds_
      return
    endif

    call z%v%mlt(m, alpha, x%v, y%v, idx_y, beta, idx_z, info, conjgx, conjgy)
  end subroutine c_mvect_mlt_vm_idxs_out

  subroutine c_mvect_mlt_mm_full_out(m, alpha, x, y, beta, z, info, conjgx, conjgy)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_multivect_type), intent(inout)  :: x, y, z
    integer(psb_ipk_), intent(out)              :: info
    character(len=1), intent(in), optional  :: conjgx, conjgy   !TO DO: remove from real cases?

    info = psb_success_
    if((.not. allocated(x%v)) .or. (.not. allocated(y%v)) &
          .or. (.not. allocated(z%v))) then
      info = psb_err_invalid_mvect_state_
      return
    end if

    ! Multivector with different size rise error
    if((z%get_ncols() /= x%get_ncols()) .or. (z%get_ncols() /= y%get_ncols())) then
      info = psb_err_invalid_mvect_size_
      return
    endif
    
    call z%v%mlt(m, alpha, x%v, y%v, beta, info, conjgx, conjgy)
  end subroutine c_mvect_mlt_mm_full_out

  subroutine c_mvect_mlt_mm_idxs_out(m, alpha, x, idx_x, y, idx_y, beta, z, idx_z, info, conjgx, conjgy)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_multivect_type), intent(inout)  :: x, y, z
    integer(psb_ipk_), intent(in)               :: idx_x, idx_y, idx_z
    integer(psb_ipk_), intent(out)              :: info
    character(len=1), intent(in), optional  :: conjgx, conjgy   !TO DO: remove from real cases?

    info = psb_success_
    if((.not. allocated(x%v)) .or. (.not. allocated(y%v)) &
          .or. (.not. allocated(z%v))) then
      info = psb_err_invalid_mvect_state_
      return
    end if

    if((idx_z < 0) .or. (idx_y < 0) .or. (idx_x < 0)) then
      info = psb_err_iarg_neg_
      return
    endif

    if((idx_z > z%get_ncols()) .or. (idx_y > y%get_ncols()) &
          .or. (idx_x > x%get_ncols())) then
      info = psb_err_entry_out_of_bounds_
      return
    endif

    call z%v%mlt(m, alpha, x%v, idx_x, y%v, idx_y, beta, idx_z, info, conjgx, conjgy)
  end subroutine c_mvect_mlt_mm_idxs_out
  
  subroutine c_mvect_mlt_vm_ext(m, alpha, x, y, idx_y, beta, z, info, conjgx, conjgy)
    use psi_serial_mod
    use psb_c_vect_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_vect_type), intent(inout)       :: x, z
    class(psb_c_multivect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(in)               :: idx_y
    integer(psb_ipk_), intent(out)              :: info
    character(len=1), intent(in), optional  :: conjgx, conjgy   !TO DO: remove from real cases?

    info = psb_success_
    if((.not. allocated(x%v)) .or. (.not. allocated(z%v))) then
      info = psb_err_invalid_vect_state_
      return
    end if

    if(.not. allocated(y%v)) then
      info = psb_err_invalid_mvect_state_
      return
    end if

    if(idx_y < 0) then
      info = psb_err_iarg_neg_
      return
    endif

    if(idx_y > y%get_ncols()) then
      info = psb_err_entry_out_of_bounds_
      return
    endif

    call y%v%mlt(m, alpha, x%v, idx_y, beta, z%v, info, conjgx, conjgy)
  end subroutine c_mvect_mlt_vm_ext
  
  subroutine c_mvect_mlt_mm_ext(m, alpha, x, idx_x, y, idx_y, beta, z, info, conjgx, conjgy)
    use psi_serial_mod
    use psb_c_vect_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                  :: alpha, beta
    class(psb_c_multivect_type), intent(inout)  :: x, y
    integer(psb_ipk_), intent(in)               :: idx_x, idx_y
    class(psb_c_vect_type), intent(inout)       :: z
    integer(psb_ipk_), intent(out)              :: info
    character(len=1), intent(in), optional  :: conjgx, conjgy   !TO DO: remove from real cases?

    info = psb_success_
    if((.not. allocated(x%v)) .or. (.not. allocated(y%v))) then
      info = psb_err_invalid_mvect_state_
      return
    end if

    if(.not. allocated(z%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if

    if((idx_y < 0) .or. (idx_x < 0)) then
      info = psb_err_iarg_neg_
      return
    endif

    if((idx_y > y%get_ncols()) .or. (idx_x > x%get_ncols())) then
      info = psb_err_entry_out_of_bounds_
      return
    endif

    call y%v%mlt(m, alpha, x%v, idx_x, idx_y, beta, z%v, info, conjgx, conjgy)
  end subroutine c_mvect_mlt_mm_ext

  function c_mvect_nrm2_full(m, x) result(res)
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    class(psb_c_multivect_type), intent(inout)  :: x
    real(psb_spk_), allocatable :: res(:)

    if(.not. allocated(x%v)) then
      res = szero
      return
    end if
    
    res = x%v%nrm2(m)
  end function c_mvect_nrm2_full
  
  function c_mvect_nrm2_idxs(m, x, idx) result(res)
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    class(psb_c_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: idx
    real(psb_spk_)  :: res

    if(.not. allocated(x%v)) then
      res = szero
      return
    end if
    
    res = x%v%nrm2(m, idx)
  end function c_mvect_nrm2_idxs

  !!$  subroutine c_mvect_scal(alpha, x)
  !!$    use psi_serial_mod
  !!$    implicit none
  !!$    complex(psb_spk_), intent(in)                 :: alpha
  !!$    class(psb_c_multivect_type), intent(inout) :: x
  !!$
  !!$    if(allocated(x%v)) call x%v%scal(alpha)
  !!$  end subroutine c_mvect_scal
  !!$
  !!$  function c_mvect_amax(n, x) result(res)
  !!$    implicit none
  !!$    integer(psb_ipk_), intent(in)              :: n
  !!$    class(psb_c_multivect_type), intent(inout) :: x
  !!$    real(psb_spk_) :: res
  !!$
  !!$    if(allocated(x%v)) then
  !!$      res = x%v%amax(n)
  !!$    else
  !!$      res = szero
  !!$    end if
  !!$  end function c_mvect_amax
  !!$
  !!$  function c_mvect_asum(n, x) result(res)
  !!$    implicit none
  !!$    integer(psb_ipk_), intent(in)              :: n
  !!$    class(psb_c_multivect_type), intent(inout) :: x
  !!$    real(psb_spk_) :: res
  !!$
  !!$    if(allocated(x%v)) then
  !!$      res = x%v%asum(n)
  !!$    else
  !!$      res = szero
  !!$    end if
  !!$  end function c_mvect_asum
end module psb_c_multivect_mod