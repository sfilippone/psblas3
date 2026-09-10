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
! package: psb_l_vect_mod
!
! This module contains the definition of the psb_l_vect type which
! is the outer container for dense vectors.
! Therefore all methods simply invoke the corresponding methods of the
! inner component.
!
module psb_l_vect_mod
  use psb_realloc_mod
  use psb_l_base_vect_mod
  use psb_i_vect_mod

  type psb_l_vect_type
    class(psb_l_base_vect_type), allocatable :: v
    integer(psb_ipk_) :: nrmv = 0
    integer(psb_ipk_) :: remote_build = psb_matbld_noremote_
    integer(psb_ipk_) :: dupl = psb_dupl_add_
    integer(psb_lpk_), allocatable :: rmtv(:)
    integer(psb_lpk_), allocatable :: rmidx(:)
  contains
    procedure, pass(x) :: get_nrows => l_vect_get_nrows
    procedure, pass(x) :: sizeof   => l_vect_sizeof
    procedure, pass(x) :: get_fmt  => l_vect_get_fmt
    procedure, pass(x) :: is_remote_build => l_vect_is_remote_build
    procedure, pass(x) :: set_remote_build => l_vect_set_remote_build
    procedure, pass(x) :: get_nrmv => l_vect_get_nrmv
    procedure, pass(x) :: set_nrmv => l_vect_set_nrmv
    procedure, pass(x) :: all      => l_vect_all
    procedure, pass(x) :: reall    => l_vect_reall
    procedure, pass(x) :: zero     => l_vect_zero
    procedure, pass(x) :: asb      => l_vect_asb
    procedure, pass(x) :: set_dupl => l_vect_set_dupl 
    procedure, pass(x) :: get_dupl => l_vect_get_dupl
    procedure, pass(x) :: set_ncfs => l_vect_set_ncfs 
    procedure, pass(x) :: get_ncfs => l_vect_get_ncfs
    procedure, pass(x) :: set_state => l_vect_set_state
    procedure, pass(x) :: set_null  => l_vect_set_null
    procedure, pass(x) :: set_bld   => l_vect_set_bld
    procedure, pass(x) :: set_upd   => l_vect_set_upd
    procedure, pass(x) :: set_asb   => l_vect_set_asb
    procedure, pass(x) :: get_state => l_vect_get_state
    procedure, pass(x) :: is_null   => l_vect_is_null
    procedure, pass(x) :: is_bld    => l_vect_is_bld
    procedure, pass(x) :: is_upd    => l_vect_is_upd
    procedure, pass(x) :: is_asb    => l_vect_is_asb
    procedure, pass(x) :: reinit    => l_vect_reinit

    procedure, pass(x) :: gthab    => l_vect_gthab
    procedure, pass(x) :: gthzv    => l_vect_gthzv
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => l_vect_sctb
    generic, public    :: sct      => sctb
    procedure, pass(x) :: free     => l_vect_free
    procedure, pass(x) :: ins_a    => l_vect_ins_a
    procedure, pass(x) :: ins_v    => l_vect_ins_v
    generic, public    :: ins      => ins_v, ins_a
    procedure, pass(x) :: bld_x    => l_vect_bld_x
    procedure, pass(x) :: bld_mn   => l_vect_bld_mn
    procedure, pass(x) :: bld_en   => l_vect_bld_en
    generic, public    :: bld      => bld_x, bld_mn, bld_en
    procedure, pass(x) :: get_vect => l_vect_get_vect
    procedure, pass(x) :: cnv      => l_vect_cnv
    procedure, pass(x) :: set_scal => l_vect_set_scal
    procedure, pass(x) :: set_vect => l_vect_set_vect
    generic, public    :: set      => set_vect, set_scal
    procedure, pass(x) :: clone    => l_vect_clone

    procedure, pass(x) :: sync     => l_vect_sync
    procedure, pass(x) :: is_host  => l_vect_is_host
    procedure, pass(x) :: is_dev   => l_vect_is_dev
    procedure, pass(x) :: is_sync  => l_vect_is_sync
    procedure, pass(x) :: set_host => l_vect_set_host
    procedure, pass(x) :: set_dev  => l_vect_set_dev
    procedure, pass(x) :: set_sync => l_vect_set_sync
    procedure, pass(x) :: check_addr => l_vect_check_addr
  end type psb_l_vect_type

  public  :: psb_l_vect
  private :: constructor, size_const
  interface psb_l_vect
    module procedure constructor, size_const
  end interface psb_l_vect

  private :: l_vect_get_nrows, l_vect_sizeof, l_vect_get_fmt, &
            & l_vect_all, l_vect_reall, l_vect_zero, l_vect_asb, &
            & l_vect_gthab, l_vect_gthzv, l_vect_sctb, &
            & l_vect_free, l_vect_ins_a, l_vect_ins_v, l_vect_bll_x, &
            & l_vect_bll_mn, l_vect_bll_en, l_vect_get_vect, &
            & l_vect_cnv, l_vect_set_scal, &
            & l_vect_set_vect, l_vect_clone, l_vect_sync, l_vect_is_host, &
            & l_vect_is_dev, l_vect_is_sync, l_vect_set_host, &
            & l_vect_set_dev, l_vect_set_sync, &
            & l_vect_set_remote_build, l_is_remote_build, &
            & l_vect_set_dupl, l_get_dupl, l_vect_set_nrmv, l_get_nrmv
  class(psb_l_base_vect_type), allocatable, target, &
       & save, private :: psb_l_base_vect_default

  interface psb_set_vect_default
    module procedure psb_l_set_vect_default
  end interface psb_set_vect_default

  interface psb_get_vect_default
    module procedure psb_l_get_vect_default
  end interface psb_get_vect_default

contains
  function l_vect_get_dupl(x) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    if(allocated(x%v)) then 
      res = x%v%get_dupl()
    else
      res = psb_dupl_null_
    end if
  end function l_vect_get_dupl

  subroutine l_vect_set_dupl(x, val)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if(allocated(x%v)) then 
      if(present(val)) then
        call x%v%set_dupl(val)
      else
        call x%v%set_dupl(psb_dupl_def_)
      end if
    end if
  end subroutine l_vect_set_dupl

  function l_vect_get_ncfs(x) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    if(allocated(x%v)) then 
      res = x%v%get_ncfs()
    else
      res = 0
    end if
  end function l_vect_get_ncfs

  subroutine l_vect_set_ncfs(x, val)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if(allocated(x%v)) then 
      if(present(val)) then
        call x%v%set_ncfs(val)
      else
        call x%v%set_ncfs(0)
      end if
    end if
  end subroutine l_vect_set_ncfs

  function l_vect_get_state(x) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    if(allocated(x%v)) then 
      res = x%v%get_state()
    else
      res = psb_vect_null_
    end if
  end function l_vect_get_state

  function l_vect_is_null(x) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: x
    logical :: res

    res = (x%get_state() == psb_vect_null_)
  end function l_vect_is_null

  function l_vect_is_bld(x) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: x
    logical :: res

    res = (x%get_state() == psb_vect_bld_)
  end function l_vect_is_bld

  function l_vect_is_upd(x) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: x
    logical :: res

    res = (x%get_state() == psb_vect_upd_)
  end function l_vect_is_upd

  function l_vect_is_asb(x) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: x
    logical :: res

    res = (x%get_state() == psb_vect_asb_)
  end function l_vect_is_asb

  subroutine l_vect_set_state(n, x)
    implicit none
    integer(psb_ipk_), intent(in)         :: n
    class(psb_l_vect_type), intent(inout) :: x

    if(allocated(x%v)) call x%v%set_state(n)
  end subroutine l_vect_set_state

  subroutine l_vect_set_null(x)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_null_)
  end subroutine l_vect_set_null

  subroutine l_vect_set_bld(x)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_bld_)
  end subroutine l_vect_set_bld

  subroutine l_vect_set_upd(x)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_upd_)
  end subroutine l_vect_set_upd

  subroutine l_vect_set_asb(x)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_asb_)
  end subroutine l_vect_set_asb

  function l_vect_get_nrmv(x) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = x%nrmv
  end function l_vect_get_nrmv

  subroutine l_vect_set_nrmv(x, val)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: val

    x%nrmv = val
  end subroutine l_vect_set_nrmv

  function l_vect_is_remote_build(x) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: x
    logical :: res

    res = (x%remote_build == psb_matbld_remote_)
  end function l_vect_is_remote_build

  subroutine l_vect_set_remote_build(x, val)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if(present(val)) then
      x%remote_build = val
    else
      x%remote_build = psb_matbld_remote_
    end if
  end subroutine l_vect_set_remote_build
        
  subroutine psb_l_set_vect_default(v)
    implicit none
    class(psb_l_base_vect_type), intent(in) :: v

    if(allocated(psb_l_base_vect_default)) deallocate(psb_l_base_vect_default)
    allocate(psb_l_base_vect_default, mold = v)
  end subroutine psb_l_set_vect_default

  function psb_l_get_vect_default(v) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: v
    class(psb_l_base_vect_type), pointer :: res

    res => psb_l_get_base_vect_default()
  end function psb_l_get_vect_default

  subroutine psb_l_clear_vect_default()
    implicit none
    if(allocated(psb_l_base_vect_default)) deallocate(psb_l_base_vect_default)
  end subroutine psb_l_clear_vect_default

  function psb_l_get_base_vect_default() result(res)
    implicit none
    class(psb_l_base_vect_type), pointer :: res

    if(.not. allocated(psb_l_base_vect_default)) &
      & allocate(psb_l_base_vect_type :: psb_l_base_vect_default)

    res => psb_l_base_vect_default
  end function psb_l_get_base_vect_default

  subroutine l_vect_clone(x, y, info)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x, y
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
  end subroutine l_vect_clone

  subroutine l_vect_bld_x(x, invect, mold, scratch)
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_lpk_), intent(in)            :: invect(:)
    class(psb_l_base_vect_type), intent(in), optional :: mold
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
      allocate(x%v, stat = info, mold = psb_l_get_base_vect_default())
    endif

    if(info == psb_success_) call x%v%bld(invect, scratch = scratch_)
  end subroutine l_vect_bld_x

  subroutine l_vect_bld_mn(x, n, mold, scratch)
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_mpk_), intent(in)         :: n
    class(psb_l_base_vect_type), intent(in), optional :: mold
    logical, intent(in), optional                     :: scratch

    logical :: scratch_
    integer(psb_ipk_) :: info
    class(psb_l_base_vect_type), pointer :: mld

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
      allocate(x%v, stat = info, mold = psb_l_get_base_vect_default())
    endif
    if(info == psb_success_) call x%v%bld(n, scratch = scratch_)
  end subroutine l_vect_bld_mn

  subroutine l_vect_bld_en(x, n, mold, scratch)
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_epk_), intent(in)         :: n
    class(psb_l_base_vect_type), intent(in), optional :: mold
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
      allocate(x%v, stat = info, mold = psb_l_get_base_vect_default())
    endif
    if(info == psb_success_) call x%v%bld(n, scratch = scratch_)
  end subroutine l_vect_bld_en

  function l_vect_get_vect(x, n) result(res)
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), optional :: n
    integer(psb_lpk_), allocatable :: res(:)

    integer(psb_ipk_) :: info

    if(allocated(x%v)) res = x%v%get_vect(n)
  end function l_vect_get_vect

  subroutine l_vect_set_scal(x, val, first, last)
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_lpk_), intent(in)            :: val
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info
    if(allocated(x%v)) call x%v%set(val, first, last)
  end subroutine l_vect_set_scal

  subroutine l_vect_set_vect(x, val, first, last)
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_lpk_), intent(in)            :: val(:)
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info
    if(allocated(x%v)) call x%v%set(val, first, last)
  end subroutine l_vect_set_vect

  subroutine l_vect_check_addr(x)
    class(psb_l_vect_type), intent(inout) :: x

    integer(psb_ipk_) :: info
    if(allocated(x%v)) call x%v%check_addr()
  end subroutine l_vect_check_addr

  function constructor(x) result(this)
    integer(psb_lpk_) :: x(:)
    type(psb_l_vect_type) :: this

    integer(psb_ipk_) :: info

    call this%bld(x)
    call this%asb(size(x, kind = psb_ipk_), info)
  end function constructor

  function size_const(n) result(this)
    integer(psb_ipk_), intent(in) :: n
    type(psb_l_vect_type) :: this

    integer(psb_ipk_) :: info

    call this%bld(n)
    call this%asb(n, info)
  end function size_const

  function l_vect_get_nrows(x) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if(allocated(x%v)) res = x%v%get_nrows()
  end function l_vect_get_nrows

  function l_vect_sizeof(x) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: x
    integer(psb_epk_) :: res

    res = 0
    if(allocated(x%v)) res = x%v%sizeof()
  end function l_vect_sizeof

  function l_vect_get_fmt(x) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: x
    character(len=5) :: res

    res = 'NULL'
    if(allocated(x%v)) res = x%v%get_fmt()
  end function l_vect_get_fmt

  subroutine l_vect_all(n, x, info, mold)
    implicit none
    integer(psb_ipk_), intent(in)         :: n
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info
    class(psb_l_base_vect_type), intent(in), optional :: mold

    if(allocated(x%v)) call x%free(info)

    if(present(mold)) then
      allocate(x%v, stat = info, mold = mold)
    else
      allocate(psb_l_base_vect_type :: x%v, stat = info)
    endif
    if(info == psb_success_) then
      call x%v%all(n, info)
    else
      info = psb_err_alloc_dealloc_
    end if
    call x%set_bld()
  end subroutine l_vect_all

  subroutine l_vect_reinit(x, info, clear)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)      :: info
    logical, intent(in), optional       :: clear

    if(allocated(x%v)) call x%v%reinit(info, clear)
    call x%set_upd()
  end subroutine l_vect_reinit
 
  subroutine l_vect_reall(n, x, info)
    implicit none
    integer(psb_ipk_), intent(in)         :: n
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    if(.not. allocated(x%v)) call x%all(n, info)
    if(info == psb_success_) call x%asb(n, info)
  end subroutine l_vect_reall

  subroutine l_vect_zero(x)
    use psi_serial_mod
    implicit none
    class(psb_l_vect_type), intent(inout) :: x

    if(allocated(x%v)) call x%v%zero()
  end subroutine l_vect_zero

  subroutine l_vect_asb(n, x, info, scratch)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    integer(psb_ipk_), intent(in)         :: n
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info
    logical, intent(in), optional :: scratch

    if(allocated(x%v)) then
      call x%v%asb(n, info, scratch = scratch)
      call x%set_asb()
    end if
  end subroutine l_vect_asb

  subroutine l_vect_gthab(n, idx, alpha, x, beta, y)
    use psi_serial_mod
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_lpk_)    :: alpha, beta, y(:)
    class(psb_l_vect_type) :: x

    if(allocated(x%v)) call x%v%gth(n, idx, alpha, beta, y)
  end subroutine l_vect_gthab

  subroutine l_vect_gthzv(n, idx, x, y)
    use psi_serial_mod
    integer(psb_mpk_)       :: n
    integer(psb_ipk_)       :: idx(:)
    class(psb_l_vect_type)  :: x
    integer(psb_lpk_)          :: y(:)

    if(allocated(x%v)) call x%v%gth(n, idx, y)
  end subroutine l_vect_gthzv

  subroutine l_vect_sctb(n, idx, x, beta, y)
    use psi_serial_mod
    integer(psb_mpk_)       :: n
    integer(psb_ipk_)       :: idx(:)
    integer(psb_lpk_)          :: beta, x(:)
    class(psb_l_vect_type)  :: y

    if(allocated(y%v)) call y%v%sct(n, idx, x, beta)
  end subroutine l_vect_sctb

  subroutine l_vect_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    if(allocated(x%v)) then
      call x%v%free(info)
      if(info == psb_success_) deallocate(x%v, stat = info)
    end if
  end subroutine l_vect_free

  subroutine l_vect_ins_a(n, irl, val, x, maxr, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)         :: n, maxr
    integer(psb_ipk_), intent(in)         :: irl(:)
    integer(psb_lpk_), intent(in)            :: val(:)
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_) :: i, dupl

    info = psb_success_
    if(.not. allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if
    dupl = x%get_dupl()
    call x%v%ins(n, irl, val, dupl, maxr, info)
  end subroutine l_vect_ins_a

  subroutine l_vect_ins_v(n, irl, val, x, maxr, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)         :: n, maxr
    class(psb_i_vect_type), intent(inout) :: irl
    class(psb_l_vect_type), intent(inout) :: val
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    integer(psb_ipk_) :: i, dupl

    info = psb_success_
    if(.not. (allocated(x%v) .and. allocated(irl%v) .and. allocated(val%v))) then
      info = psb_err_invalid_vect_state_
      return
    end if
    dupl = x%get_dupl()
    call x%v%ins(n, irl%v, val%v, dupl, maxr, info)
  end subroutine l_vect_ins_v

  subroutine l_vect_cnv(x, mold)
    class(psb_l_vect_type), intent(inout) :: x
    class(psb_l_base_vect_type), intent(in), optional :: mold

    class(psb_l_base_vect_type), allocatable :: tmp
    integer(psb_ipk_) :: info

    info = psb_success_
    if(present(mold)) then
      allocate(tmp, stat = info, mold = mold)
    else
      allocate(tmp, stat = info, mold = psb_l_get_base_vect_default())
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
  end subroutine l_vect_cnv

  subroutine l_vect_sync(x)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x

    if(allocated(x%v)) call x%v%sync()
  end subroutine l_vect_sync

  subroutine l_vect_set_sync(x)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x

    if(allocated(x%v)) call x%v%set_sync()
  end subroutine l_vect_set_sync

  subroutine l_vect_set_host(x)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x

    if(allocated(x%v)) call x%v%set_host()
  end subroutine l_vect_set_host

  subroutine l_vect_set_dev(x)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x

    if(allocated(x%v)) call x%v%set_dev()
  end subroutine l_vect_set_dev

  function l_vect_is_sync(x) result(res)
    implicit none
    logical :: res
    class(psb_l_vect_type), intent(inout) :: x

    res = .true.
    if(allocated(x%v)) res = x%v%is_sync()
  end function l_vect_is_sync

  function l_vect_is_host(x) result(res)
    implicit none
    logical :: res
    class(psb_l_vect_type), intent(inout) :: x

    res = .true.
    if(allocated(x%v)) res = x%v%is_host()
  end function l_vect_is_host

  function l_vect_is_dev(x) result(res)
    implicit none
    logical :: res
    class(psb_l_vect_type), intent(inout) :: x

    res = .false.
    if(allocated(x%v)) res = x%v%is_dev()
  end function l_vect_is_dev



end module psb_l_vect_mod


module psb_l_multivect_mod
  use psb_l_base_multivect_mod
  use psb_const_mod
  use psb_i_vect_mod

  !private
  type psb_l_multivect_type
    class(psb_l_base_multivect_type), allocatable :: v
    integer(psb_ipk_) :: nrmv = 0
    integer(psb_ipk_) :: remote_build = psb_matbld_noremote_
    integer(psb_lpk_), allocatable :: rmtv(:, :)
  contains
    procedure, pass(x) :: get_nrows => l_mvect_get_nrows
    procedure, pass(x) :: get_ncols => l_mvect_get_ncols
    procedure, pass(x) :: sizeof   => l_mvect_sizeof
    procedure, pass(x) :: get_fmt  => l_mvect_get_fmt
    procedure, pass(x) :: is_remote_build => l_mvect_is_remote_build
    procedure, pass(x) :: set_remote_build => l_mvect_set_remote_build

    procedure, pass(x) :: all      => l_mvect_all
    procedure, pass(x) :: reall    => l_mvect_reall
    procedure, pass(x) :: zero     => l_mvect_zero
    procedure, pass(x) :: asb      => l_mvect_asb
    procedure, pass(x) :: sync     => l_mvect_sync
    procedure, pass(x) :: free     => l_mvect_free
    procedure, pass(x) :: reinit   => l_mvect_reinit
    procedure, pass(x) :: set_ncfs => l_mvect_set_ncfs 
    procedure, pass(x) :: get_ncfs => l_mvect_get_ncfs
    procedure, pass(x) :: set_dupl => l_mvect_set_dupl 
    procedure, pass(x) :: get_dupl => l_mvect_get_dupl
    procedure, pass(x) :: set_state => l_mvect_set_state
    procedure, pass(x) :: set_null  => l_mvect_set_null
    procedure, pass(x) :: set_bld   => l_mvect_set_bld
    procedure, pass(x) :: set_upd   => l_mvect_set_upd
    procedure, pass(x) :: set_asb   => l_mvect_set_asb
    procedure, pass(x) :: get_state => l_mvect_get_state
    procedure, pass(x) :: is_null   => l_mvect_is_null
    procedure, pass(x) :: is_bld    => l_mvect_is_bld
    procedure, pass(x) :: is_upd    => l_mvect_is_upd
    procedure, pass(x) :: is_asb    => l_mvect_is_asb
    !!$ procedure, pass(x) :: base_cpy  => l_mvect_cpy

    procedure, pass(x) :: ins      => l_mvect_ins
    procedure, pass(x) :: bld_x    => l_mvect_bld_x
    procedure, pass(x) :: bld_n    => l_mvect_bld_n
    generic, public    :: bld      => bld_x, bld_n
    procedure, pass(x) :: get_vect => l_mvect_get_vect
    procedure, pass(x) :: cnv      => l_mvect_cnv

    procedure, pass(x) :: set_scal => l_mvect_set_scal
    procedure, pass(x) :: set_vect => l_mvect_set_vect
    generic, public    :: set      => set_vect, set_scal

    procedure, pass(x) :: clone    => l_mvect_clone
    procedure, pass(x) :: gthab    => l_mvect_gthab
    procedure, pass(x) :: gthzv    => l_mvect_gthzv
    procedure, pass(x) :: gthzv_x  => l_mvect_gthzv_x
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => l_mvect_sctb
    procedure, pass(y) :: sctb_x   => l_mvect_sctb_x
    generic, public    :: sct      => sctb, sctb_x

  end type psb_l_multivect_type

  public  :: psb_l_multivect, psb_l_multivect_type, psb_l_base_multivect_type, &
          & psb_set_multivect_default, psb_get_multivect_default

  private
  interface psb_l_multivect
    module procedure constructor, size_const
  end interface psb_l_multivect

  class(psb_l_base_multivect_type), allocatable, target, &
       & save, private :: psb_l_base_multivect_default

  interface psb_set_multivect_default
    module procedure psb_l_set_multivect_default
  end interface psb_set_multivect_default

  interface psb_get_multivect_default
    module procedure psb_l_get_multivect_default
  end interface psb_get_multivect_default

contains
  function l_mvect_is_remote_build(x) result(res)
    implicit none
    class(psb_l_multivect_type), intent(in) :: x
    logical :: res
    
    res = (x%remote_build == psb_matbld_remote_)
  end function l_mvect_is_remote_build

  subroutine l_mvect_set_remote_build(x, val)
    implicit none
    class(psb_l_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if(present(val)) then
      x%remote_build = val
    else
      x%remote_build = psb_matbld_remote_
    end if
  end subroutine l_mvect_set_remote_build

  subroutine psb_l_set_multivect_default(v)
    implicit none
    class(psb_l_base_multivect_type), intent(in) :: v

    if(allocated(psb_l_base_multivect_default)) then
      deallocate(psb_l_base_multivect_default)
    end if
    allocate(psb_l_base_multivect_default, mold = v)
  end subroutine psb_l_set_multivect_default

  function psb_l_get_multivect_default(v) result(res)
    implicit none
    class(psb_l_multivect_type), intent(in) :: v
    class(psb_l_base_multivect_type), pointer :: res

    res => psb_l_get_base_multivect_default()
  end function psb_l_get_multivect_default

  function psb_l_get_base_multivect_default() result(res)
    implicit none
    class(psb_l_base_multivect_type), pointer :: res

    if(.not. allocated(psb_l_base_multivect_default)) then
      allocate(psb_l_base_multivect_type :: psb_l_base_multivect_default)
    end if

    res => psb_l_base_multivect_default
  end function psb_l_get_base_multivect_default

  subroutine l_mvect_clone(x, y, info)
    implicit none
    class(psb_l_multivect_type), intent(inout)  :: x, y
    integer(psb_ipk_), intent(out)                :: info

    info = psb_success_
    call y%free(info)
    if((info == psb_success_) .and. allocated(x%v)) then
      call y%bld_x(x%get_vect(), mold = x%v)
    end if
  end subroutine l_mvect_clone

  subroutine l_mvect_bld_x(x, invect, mold)
    class(psb_l_multivect_type), intent(out)  :: x
    integer(psb_lpk_), intent(in)                :: invect(:, :)
    class(psb_l_base_multivect_type), intent(in), optional :: mold

    integer(psb_ipk_) :: info
    info = psb_success_

    if(present(mold)) then
      allocate(x%v, stat = info, mold = mold)
    else
      allocate(x%v, stat = info, mold = psb_l_get_base_multivect_default())
    endif

    if(info == psb_success_) call x%v%bld(invect)
  end subroutine l_mvect_bld_x

  subroutine l_mvect_bld_n(x, m, n, mold, scratch)
    class(psb_l_multivect_type), intent(out)  :: x
    integer(psb_ipk_), intent(in)             :: m, n
    class(psb_l_base_multivect_type), intent(in), optional  :: mold
    logical, intent(in), optional                           :: scratch

    integer(psb_ipk_) :: info
    info = psb_success_

    if(present(mold)) then
      allocate(x%v, stat = info, mold = mold)
    else
      allocate(x%v, stat = info, mold = psb_l_get_base_multivect_default())
    endif

    if(info == psb_success_) call x%v%bld(m, n, scratch = scratch)
  end subroutine l_mvect_bld_n

  function l_mvect_get_vect(x) result(res)
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_lpk_), allocatable                 :: res(:, :)
    integer(psb_ipk_) :: info

    if(allocated(x%v)) then
      res = x%v%get_vect()
    end if
  end function l_mvect_get_vect

  subroutine l_mvect_set_scal(x, val, rfirst, rlast)
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_lpk_), intent(in)                  :: val
    integer(psb_ipk_), optional :: rfirst, rlast

    integer(psb_ipk_) :: info
    if(allocated(x%v)) call x%v%set(val, rfirst, rlast)
  end subroutine l_mvect_set_scal

  subroutine l_mvect_set_vect(x, val)
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_lpk_), intent(in)                  :: val(:, :)

    integer(psb_ipk_) :: info
    if(allocated(x%v)) call x%v%set(val)
  end subroutine l_mvect_set_vect


  function constructor(x) result(this)
    integer(psb_lpk_)  :: x(:, :)
    type(psb_l_multivect_type)  :: this

    integer(psb_ipk_) :: info

    call this%bld_x(x)
    call this%asb(size(x, dim=1, kind=psb_ipk_), size(x, dim=2, kind=psb_ipk_), info)
  end function constructor

  function size_const(m, n) result(this)
    integer(psb_ipk_), intent(in) :: m, n
    type(psb_l_multivect_type)  :: this

    integer(psb_ipk_) :: info

    call this%bld_n(m, n)
    call this%asb(m, n, info)
  end function size_const

  function l_mvect_get_nrows(x) result(res)
    implicit none
    class(psb_l_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if(allocated(x%v)) res = x%v%get_nrows()
  end function l_mvect_get_nrows

  function l_mvect_get_ncols(x) result(res)
    implicit none
    class(psb_l_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if(allocated(x%v)) res = x%v%get_ncols()
  end function l_mvect_get_ncols

  function l_mvect_sizeof(x) result(res)
    implicit none
    class(psb_l_multivect_type), intent(in) :: x
    integer(psb_epk_) :: res

    res = 0
    if(allocated(x%v)) res = x%v%sizeof()
  end function l_mvect_sizeof

  function l_mvect_get_fmt(x) result(res)
    implicit none
    class(psb_l_multivect_type), intent(in) :: x
    character(len=5) :: res

    res = 'NULL'
    if(allocated(x%v)) res = x%v%get_fmt()
  end function l_mvect_get_fmt

  subroutine l_mvect_all(m, n, x, info, mold)
    implicit none
    integer(psb_ipk_), intent(in)             :: m, n
    class(psb_l_multivect_type), intent(out)  :: x
    integer(psb_ipk_), intent(out)            :: info
    class(psb_l_base_multivect_type), intent(in), optional :: mold

    if(present(mold)) then
      allocate(x%v, stat = info, mold = mold)
    else
      allocate(psb_l_base_multivect_type :: x%v, stat = info)
    endif
    if(info == psb_success_) then
      call x%v%all(m, n, info)
    else
      info = psb_err_alloc_dealloc_
    end if
    call x%set_bld()
  end subroutine l_mvect_all

  subroutine l_mvect_reall(m, n, x, info)
    implicit none
    integer(psb_ipk_), intent(in)               :: m, n
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = psb_success_
    if(.not. allocated(x%v)) call x%all(m, n, info)
    if(info == psb_success_) call x%asb(m, n, info)
  end subroutine l_mvect_reall

  subroutine l_mvect_reinit(x, info)
    implicit none
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = psb_success_
    if(allocated(x%v)) call x%v%reinit(info)
    call x%set_upd()
  end subroutine l_mvect_reinit

  subroutine l_mvect_zero(x)
    use psi_serial_mod
    implicit none
    class(psb_l_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%zero()
  end subroutine l_mvect_zero

  subroutine l_mvect_asb(m, n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m, n
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    if(allocated(x%v)) then
      call x%v%asb(m, n, info)
      call x%set_asb()
    end if
  end subroutine l_mvect_asb

  subroutine l_mvect_sync(x)
    implicit none
    class(psb_l_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%sync()
  end subroutine l_mvect_sync

  subroutine l_mvect_gthab(n, idx, alpha, x, beta, y)
    use psi_serial_mod
    integer(psb_mpk_)           :: n
    integer(psb_ipk_)           :: idx(:)
    integer(psb_lpk_)              :: alpha, beta, y(:)
    class(psb_l_multivect_type) :: x

    if(allocated(x%v)) call x%v%gth(n, idx, alpha, beta, y)
  end subroutine l_mvect_gthab

  subroutine l_mvect_gthzv(n, idx, x, y)
    use psi_serial_mod
    integer(psb_mpk_)           :: n
    integer(psb_ipk_)           :: idx(:)
    integer(psb_lpk_)              :: y(:)
    class(psb_l_multivect_type) :: x

    if(allocated(x%v)) call x%v%gth(n, idx, y)
  end subroutine l_mvect_gthzv

  subroutine l_mvect_gthzv_x(i, n, idx, x, y)
    use psi_serial_mod
    integer(psb_ipk_)           :: i
    integer(psb_mpk_)           :: n
    class(psb_i_base_vect_type) :: idx
    class(psb_l_multivect_type) :: x
    integer(psb_lpk_)              :: y(:)

    if(allocated(x%v)) call x%v%gth(i, n, idx, y)
  end subroutine l_mvect_gthzv_x

  subroutine l_mvect_sctb(n, idx, x, beta, y)
    use psi_serial_mod
    integer(psb_mpk_)           :: n
    integer(psb_ipk_)           :: idx(:)
    integer(psb_lpk_)              :: beta, x(:)
    class(psb_l_multivect_type) :: y

    if(allocated(y%v)) call y%v%sct(n, idx, x, beta)
  end subroutine l_mvect_sctb

  subroutine l_mvect_sctb_x(i, n, idx, x, beta, y)
    use psi_serial_mod
    integer(psb_ipk_)           :: i
    integer(psb_mpk_)           :: n
    class(psb_i_base_vect_type) :: idx
    integer(psb_lpk_)              :: beta, x(:)
    class(psb_l_multivect_type) :: y

    if(allocated(y%v)) call y%v%sct(i, n, idx, x, beta)
  end subroutine l_mvect_sctb_x

  subroutine l_mvect_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = psb_success_
    if(allocated(x%v)) then
      call x%v%free(info)
      if(info == psb_success_) deallocate(x%v, stat = info)
    end if
  end subroutine l_mvect_free

  subroutine l_mvect_set_ncfs(n, x)
    integer(psb_ipk_)                           :: n
    class(psb_l_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%set_ncfs(n)
  end subroutine l_mvect_set_ncfs

  function l_mvect_get_ncfs(n, x) result(res)
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_ipk_) :: res

    if(allocated(x%v)) res = x%v%get_ncfs()
  end function l_mvect_get_ncfs

  subroutine l_mvect_set_dupl(n, x)
    integer(psb_ipk_)                           :: n
    class(psb_l_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%set_dupl(n)
  end subroutine l_mvect_set_dupl

  function l_mvect_get_dupl(x) result(res)
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_ipk_) :: res

    if(allocated(x%v)) res = x%v%get_dupl()
  end function l_mvect_get_dupl

  subroutine l_mvect_set_state(n, x)
    integer(psb_ipk_)                           :: n
    class(psb_l_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%set_state(n)
  end subroutine l_mvect_set_state

  function l_mvect_get_state(n, x) result(res)
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_ipk_) :: res

    if(allocated(x%v)) res = x%v%get_state()
  end function l_mvect_get_state

  subroutine l_mvect_set_null(x)
    class(psb_l_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%set_null()
  end subroutine l_mvect_set_null

  function l_mvect_is_null(x) result(res)
    class(psb_l_multivect_type), intent(inout)  :: x
    logical :: res

    res = .false.
    if(allocated(x%v)) res = x%v%is_null()
  end function l_mvect_is_null

  subroutine l_mvect_set_bld(x)
    class(psb_l_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%set_bld()
  end subroutine l_mvect_set_bld

  function l_mvect_is_bld(x) result(res)
    class(psb_l_multivect_type), intent(inout)  :: x
    logical :: res

    res = .false.
    if(allocated(x%v)) res = x%v%is_bld()
  end function l_mvect_is_bld

  subroutine l_mvect_set_upd(x)
    class(psb_l_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%set_upd()
  end subroutine l_mvect_set_upd

  function l_mvect_is_upd(x) result(res)
    class(psb_l_multivect_type), intent(inout)  :: x
    logical :: res

    res = .false.
    if(allocated(x%v)) res = x%v%is_upd()
  end function l_mvect_is_upd

  subroutine l_mvect_set_asb(x)
    class(psb_l_multivect_type), intent(inout)  :: x

    if(allocated(x%v)) call x%v%set_asb()
  end subroutine l_mvect_set_asb

  function l_mvect_is_asb(x) result(res)
    class(psb_l_multivect_type), intent(inout)  :: x
    logical :: res

    res = .false.
    if(allocated(x%v)) res = x%v%is_asb()
  end function l_mvect_is_asb
  
  subroutine l_mvect_ins(n, irl, val, x, maxr, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: n, maxr
    integer(psb_ipk_), intent(in)               :: irl(:)
    integer(psb_lpk_), intent(in)                  :: val(:, :)
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i, dupl

    info = psb_success_
    if(.not. allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if
    dupl = x%get_dupl()
    call x%v%ins(n, irl, val, dupl, maxr, info)
  end subroutine l_mvect_ins

  subroutine l_mvect_cnv(x, mold)
    class(psb_l_multivect_type), intent(inout)  :: x
    class(psb_l_base_multivect_type), intent(in), optional :: mold

    class(psb_l_base_multivect_type), allocatable :: tmp
    integer(psb_ipk_) :: info

    if(present(mold)) then
      allocate(tmp, stat = info, mold = mold)
    else
      allocate(tmp, stat = info, mold = psb_l_get_base_multivect_default())
    endif

    if(allocated(x%v)) then
      call x%v%sync()
      if(info == psb_success_) call tmp%bld(x%v%v)
      call x%v%free(info)
    end if

    call move_alloc(tmp, x%v)
  end subroutine l_mvect_cnv
end module psb_l_multivect_mod