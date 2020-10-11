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
! package: psb_l_vect_mod
!
! This module contains the definition of the psb_l_vect type which
! is the outer container for dense vectors.
! Therefore all methods simply invoke the corresponding methods of the
! inner component.
!
module psb_l_vect_mod

  use psb_l_base_vect_mod
  use psb_i_vect_mod

  type psb_l_vect_type
    class(psb_l_base_vect_type), allocatable :: v
  contains
    procedure, pass(x) :: get_nrows => l_vect_get_nrows
    procedure, pass(x) :: sizeof   => l_vect_sizeof
    procedure, pass(x) :: get_fmt  => l_vect_get_fmt
    procedure, pass(x) :: all      => l_vect_all
    procedure, pass(x) :: reall    => l_vect_reall
    procedure, pass(x) :: zero     => l_vect_zero
    procedure, pass(x) :: asb      => l_vect_asb
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



  end type psb_l_vect_type

  public  :: psb_l_vect
  private :: constructor, size_const
  interface psb_l_vect
    module procedure constructor, size_const
  end interface psb_l_vect

  private :: l_vect_get_nrows, l_vect_sizeof, l_vect_get_fmt, &
       & l_vect_all, l_vect_reall, l_vect_zero,  l_vect_asb, &
       & l_vect_gthab, l_vect_gthzv, l_vect_sctb, &
       & l_vect_free, l_vect_ins_a, l_vect_ins_v, l_vect_bld_x, &
       & l_vect_bld_mn, l_vect_bld_en, l_vect_get_vect, &
       & l_vect_cnv, l_vect_set_scal, &
       & l_vect_set_vect, l_vect_clone, l_vect_sync, l_vect_is_host, &
       & l_vect_is_dev, l_vect_is_sync, l_vect_set_host, &
       & l_vect_set_dev, l_vect_set_sync


  class(psb_l_base_vect_type), allocatable, target,&
       & save, private :: psb_l_base_vect_default

  interface psb_set_vect_default
    module procedure psb_l_set_vect_default
  end interface psb_set_vect_default

  interface psb_get_vect_default
    module procedure psb_l_get_vect_default
  end interface psb_get_vect_default


contains


  subroutine  psb_l_set_vect_default(v)
    implicit none
    class(psb_l_base_vect_type), intent(in) :: v

    if (allocated(psb_l_base_vect_default)) then
      deallocate(psb_l_base_vect_default)
    end if
    allocate(psb_l_base_vect_default, mold=v)

  end subroutine psb_l_set_vect_default

  function psb_l_get_vect_default(v) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: v
    class(psb_l_base_vect_type), pointer :: res

    res => psb_l_get_base_vect_default()

  end function psb_l_get_vect_default

  subroutine  psb_l_clear_vect_default()
    implicit none

    if (allocated(psb_l_base_vect_default)) then
      deallocate(psb_l_base_vect_default)
    end if

  end subroutine psb_l_clear_vect_default

  function psb_l_get_base_vect_default() result(res)
    implicit none
    class(psb_l_base_vect_type), pointer :: res

    if (.not.allocated(psb_l_base_vect_default)) then
      allocate(psb_l_base_vect_type :: psb_l_base_vect_default)
    end if

    res => psb_l_base_vect_default

  end function psb_l_get_base_vect_default

  subroutine l_vect_clone(x,y,info)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x
    class(psb_l_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    call y%free(info)
    if ((info==0).and.allocated(x%v)) then
      call y%bld(x%get_vect(),mold=x%v)
    end if
  end subroutine l_vect_clone

  subroutine l_vect_bld_x(x,invect,mold)
    integer(psb_lpk_), intent(in)          :: invect(:)
    class(psb_l_vect_type), intent(inout) :: x
    class(psb_l_base_vect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info

    info = psb_success_
    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_l_get_base_vect_default())
    endif

    if (info == psb_success_) call x%v%bld(invect)

  end subroutine l_vect_bld_x


  subroutine l_vect_bld_mn(x,n,mold)
    integer(psb_mpk_), intent(in) :: n
    class(psb_l_vect_type), intent(inout) :: x
    class(psb_l_base_vect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info
    class(psb_l_base_vect_type), pointer :: mld

    info = psb_success_
    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_l_get_base_vect_default())
    endif
    if (info == psb_success_) call x%v%bld(n)

  end subroutine l_vect_bld_mn

  subroutine l_vect_bld_en(x,n,mold)
    integer(psb_epk_), intent(in) :: n
    class(psb_l_vect_type), intent(inout) :: x
    class(psb_l_base_vect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info

    info = psb_success_

    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_l_get_base_vect_default())
    endif
    if (info == psb_success_) call x%v%bld(n)

  end subroutine l_vect_bld_en

  function  l_vect_get_vect(x,n) result(res)
    class(psb_l_vect_type), intent(inout)  :: x
    integer(psb_lpk_), allocatable                 :: res(:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional :: n

    if (allocated(x%v)) then
      res = x%v%get_vect(n)
    end if
  end function l_vect_get_vect

  subroutine l_vect_set_scal(x,val,first,last)
    class(psb_l_vect_type), intent(inout)  :: x
    integer(psb_lpk_), intent(in) :: val
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val,first,last)

  end subroutine l_vect_set_scal

  subroutine l_vect_set_vect(x,val,first,last)
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_lpk_), intent(in)         :: val(:)
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val,first,last)

  end subroutine l_vect_set_vect


  function constructor(x) result(this)
    integer(psb_lpk_)   :: x(:)
    type(psb_l_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(x)
    call this%asb(size(x,kind=psb_ipk_),info)

  end function constructor


  function size_const(n) result(this)
    integer(psb_ipk_), intent(in) :: n
    type(psb_l_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(n)
    call this%asb(n,info)

  end function size_const

  function l_vect_get_nrows(x) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_nrows()
  end function l_vect_get_nrows

  function l_vect_sizeof(x) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: x
    integer(psb_epk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%sizeof()
  end function l_vect_sizeof

  function l_vect_get_fmt(x) result(res)
    implicit none
    class(psb_l_vect_type), intent(in) :: x
    character(len=5) :: res
    res = 'NULL'
    if (allocated(x%v)) res = x%v%get_fmt()
  end function l_vect_get_fmt

  subroutine l_vect_all(n, x, info, mold)

    implicit none
    integer(psb_ipk_), intent(in)           :: n
    class(psb_l_vect_type), intent(inout) :: x
    class(psb_l_base_vect_type), intent(in), optional :: mold
    integer(psb_ipk_), intent(out)      :: info

    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(psb_l_base_vect_type :: x%v,stat=info)
    endif
    if (info == 0) then
      call x%v%all(n,info)
    else
      info = psb_err_alloc_dealloc_
    end if

  end subroutine l_vect_all

  subroutine l_vect_reall(n, x, info)

    implicit none
    integer(psb_ipk_), intent(in)         :: n
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    info = 0
    if (.not.allocated(x%v)) &
         & call x%all(n,info)
    if (info == 0) &
         & call x%asb(n,info)

  end subroutine l_vect_reall

  subroutine l_vect_zero(x)
    use psi_serial_mod
    implicit none
    class(psb_l_vect_type), intent(inout)    :: x

    if (allocated(x%v)) call x%v%zero()

  end subroutine l_vect_zero

  subroutine l_vect_asb(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    integer(psb_ipk_), intent(in)              :: n
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%v)) &
         & call x%v%asb(n,info)

  end subroutine l_vect_asb

  subroutine l_vect_gthab(n,idx,alpha,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    integer(psb_lpk_) :: alpha, beta, y(:)
    class(psb_l_vect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,alpha,beta,y)

  end subroutine l_vect_gthab

  subroutine l_vect_gthzv(n,idx,x,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    integer(psb_lpk_) ::  y(:)
    class(psb_l_vect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,y)

  end subroutine l_vect_gthzv

  subroutine l_vect_sctb(n,idx,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    integer(psb_lpk_) :: beta, x(:)
    class(psb_l_vect_type) :: y

    if (allocated(y%v)) &
         &  call y%v%sct(n,idx,x,beta)

  end subroutine l_vect_sctb

  subroutine l_vect_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    class(psb_l_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (allocated(x%v)) then
      call x%v%free(info)
      if (info == 0) deallocate(x%v,stat=info)
    end if

  end subroutine l_vect_free

  subroutine l_vect_ins_a(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none
    class(psb_l_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    integer(psb_ipk_), intent(in)               :: irl(:)
    integer(psb_lpk_), intent(in)        :: val(:)
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i

    info = 0
    if (.not.allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if

    call  x%v%ins(n,irl,val,dupl,info)

  end subroutine l_vect_ins_a

  subroutine l_vect_ins_v(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none
    class(psb_l_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    class(psb_i_vect_type), intent(inout)       :: irl
    class(psb_l_vect_type), intent(inout)       :: val
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i

    info = 0
    if (.not.(allocated(x%v).and.allocated(irl%v).and.allocated(val%v))) then
      info = psb_err_invalid_vect_state_
      return
    end if

    call  x%v%ins(n,irl%v,val%v,dupl,info)

  end subroutine l_vect_ins_v


  subroutine l_vect_cnv(x,mold)
    class(psb_l_vect_type), intent(inout) :: x
    class(psb_l_base_vect_type), intent(in), optional :: mold
    class(psb_l_base_vect_type), allocatable :: tmp

    integer(psb_ipk_) :: info

    info = psb_success_
    if (present(mold)) then
      allocate(tmp,stat=info,mold=mold)
    else
      allocate(tmp,stat=info,mold=psb_l_get_base_vect_default())
    end if
    if (allocated(x%v)) then
      call x%v%sync()
      if (info == psb_success_) call tmp%bld(x%v%v)
      call x%v%free(info)
    end if
    call move_alloc(tmp,x%v)

  end subroutine l_vect_cnv


  subroutine l_vect_sync(x)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%sync()

  end subroutine l_vect_sync

  subroutine l_vect_set_sync(x)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%set_sync()

  end subroutine l_vect_set_sync

  subroutine l_vect_set_host(x)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%set_host()

  end subroutine l_vect_set_host

  subroutine l_vect_set_dev(x)
    implicit none
    class(psb_l_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%set_dev()

  end subroutine l_vect_set_dev

  function l_vect_is_sync(x) result(res)
    implicit none
    logical :: res
    class(psb_l_vect_type), intent(inout) :: x

    res = .true.
    if (allocated(x%v)) &
         & res = x%v%is_sync()

  end function l_vect_is_sync

  function l_vect_is_host(x) result(res)
    implicit none
    logical :: res
    class(psb_l_vect_type), intent(inout) :: x

    res = .true.
    if (allocated(x%v)) &
         & res = x%v%is_host()

  end function l_vect_is_host

  function l_vect_is_dev(x) result(res)
    implicit none
    logical :: res
    class(psb_l_vect_type), intent(inout) :: x

    res = .false.
    if (allocated(x%v)) &
         & res =  x%v%is_dev()

  end function l_vect_is_dev




end module psb_l_vect_mod



module psb_l_multivect_mod

  use psb_l_base_multivect_mod
  use psb_const_mod
  use psb_i_vect_mod


  !private

  type psb_l_multivect_type
    class(psb_l_base_multivect_type), allocatable :: v
  contains
    procedure, pass(x) :: get_nrows => l_vect_get_nrows
    procedure, pass(x) :: get_ncols => l_vect_get_ncols
    procedure, pass(x) :: sizeof   => l_vect_sizeof
    procedure, pass(x) :: get_fmt  => l_vect_get_fmt

    procedure, pass(x) :: all      => l_vect_all
    procedure, pass(x) :: reall    => l_vect_reall
    procedure, pass(x) :: zero     => l_vect_zero
    procedure, pass(x) :: asb      => l_vect_asb
    procedure, pass(x) :: sync     => l_vect_sync
    procedure, pass(x) :: free     => l_vect_free
    procedure, pass(x) :: ins      => l_vect_ins
    procedure, pass(x) :: bld_x    => l_vect_bld_x
    procedure, pass(x) :: bld_n    => l_vect_bld_n
    generic, public    :: bld      => bld_x, bld_n
    procedure, pass(x) :: get_vect => l_vect_get_vect
    procedure, pass(x) :: cnv      => l_vect_cnv
    procedure, pass(x) :: set_scal => l_vect_set_scal
    procedure, pass(x) :: set_vect => l_vect_set_vect
    generic, public    :: set      => set_vect, set_scal
    procedure, pass(x) :: clone    => l_vect_clone
    procedure, pass(x) :: gthab    => l_vect_gthab
    procedure, pass(x) :: gthzv    => l_vect_gthzv
    procedure, pass(x) :: gthzv_x  => l_vect_gthzv_x
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => l_vect_sctb
    procedure, pass(y) :: sctb_x   => l_vect_sctb_x
    generic, public    :: sct      => sctb, sctb_x
  end type psb_l_multivect_type

  public  :: psb_l_multivect, psb_l_multivect_type,&
       & psb_set_multivect_default, psb_get_multivect_default, &
       & psb_l_base_multivect_type

  private
  interface psb_l_multivect
    module procedure constructor, size_const
  end interface psb_l_multivect

  class(psb_l_base_multivect_type), allocatable, target,&
       & save, private :: psb_l_base_multivect_default

  interface psb_set_multivect_default
    module procedure psb_l_set_multivect_default
  end interface psb_set_multivect_default

  interface psb_get_multivect_default
    module procedure psb_l_get_multivect_default
  end interface psb_get_multivect_default


contains


  subroutine  psb_l_set_multivect_default(v)
    implicit none
    class(psb_l_base_multivect_type), intent(in) :: v

    if (allocated(psb_l_base_multivect_default)) then
      deallocate(psb_l_base_multivect_default)
    end if
    allocate(psb_l_base_multivect_default, mold=v)

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

    if (.not.allocated(psb_l_base_multivect_default)) then
      allocate(psb_l_base_multivect_type :: psb_l_base_multivect_default)
    end if

    res => psb_l_base_multivect_default

  end function psb_l_get_base_multivect_default


  subroutine l_vect_clone(x,y,info)
    implicit none
    class(psb_l_multivect_type), intent(inout) :: x
    class(psb_l_multivect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    call y%free(info)
    if ((info==0).and.allocated(x%v)) then
      call y%bld(x%get_vect(),mold=x%v)
    end if
  end subroutine l_vect_clone

  subroutine l_vect_bld_x(x,invect,mold)
    integer(psb_lpk_), intent(in)          :: invect(:,:)
    class(psb_l_multivect_type), intent(out) :: x
    class(psb_l_base_multivect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info
    class(psb_l_base_multivect_type), pointer :: mld

    info = psb_success_
    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_l_get_base_multivect_default())
    endif

    if (info == psb_success_) call x%v%bld(invect)

  end subroutine l_vect_bld_x


  subroutine l_vect_bld_n(x,m,n,mold)
    integer(psb_ipk_), intent(in) :: m,n
    class(psb_l_multivect_type), intent(out) :: x
    class(psb_l_base_multivect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info

    info = psb_success_
    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_l_get_base_multivect_default())
    endif
    if (info == psb_success_) call x%v%bld(m,n)

  end subroutine l_vect_bld_n

  function  l_vect_get_vect(x) result(res)
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_lpk_), allocatable                 :: res(:,:)
    integer(psb_ipk_) :: info

    if (allocated(x%v)) then
      res = x%v%get_vect()
    end if
  end function l_vect_get_vect

  subroutine l_vect_set_scal(x,val)
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_lpk_), intent(in) :: val

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val)

  end subroutine l_vect_set_scal

  subroutine l_vect_set_vect(x,val)
    class(psb_l_multivect_type), intent(inout) :: x
    integer(psb_lpk_), intent(in)         :: val(:,:)

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val)

  end subroutine l_vect_set_vect


  function constructor(x) result(this)
    integer(psb_lpk_)   :: x(:,:)
    type(psb_l_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(x)
    call this%asb(size(x,dim=1,kind=psb_ipk_),size(x,dim=2,kind=psb_ipk_),info)

  end function constructor


  function size_const(m,n) result(this)
    integer(psb_ipk_), intent(in) :: m,n
    type(psb_l_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(m,n)
    call this%asb(m,n,info)

  end function size_const

  function l_vect_get_nrows(x) result(res)
    implicit none
    class(psb_l_multivect_type), intent(in) :: x
    integer(psb_ipk_)  :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_nrows()
  end function l_vect_get_nrows

  function l_vect_get_ncols(x) result(res)
    implicit none
    class(psb_l_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_ncols()
  end function l_vect_get_ncols

  function l_vect_sizeof(x) result(res)
    implicit none
    class(psb_l_multivect_type), intent(in) :: x
    integer(psb_epk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%sizeof()
  end function l_vect_sizeof

  function l_vect_get_fmt(x) result(res)
    implicit none
    class(psb_l_multivect_type), intent(in) :: x
    character(len=5) :: res
    res = 'NULL'
    if (allocated(x%v)) res = x%v%get_fmt()
  end function l_vect_get_fmt

  subroutine l_vect_all(m,n, x, info, mold)

    implicit none
    integer(psb_ipk_), intent(in)       :: m,n
    class(psb_l_multivect_type), intent(out) :: x
    class(psb_l_base_multivect_type), intent(in), optional :: mold
    integer(psb_ipk_), intent(out)      :: info

    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(psb_l_base_multivect_type :: x%v,stat=info)
    endif
    if (info == 0) then
      call x%v%all(m,n,info)
    else
      info = psb_err_alloc_dealloc_
    end if

  end subroutine l_vect_all

  subroutine l_vect_reall(m,n, x, info)

    implicit none
    integer(psb_ipk_), intent(in)         :: m,n
    class(psb_l_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    info = 0
    if (.not.allocated(x%v)) &
         & call x%all(m,n,info)
    if (info == 0) &
         & call x%asb(m,n,info)

  end subroutine l_vect_reall

  subroutine l_vect_zero(x)
    use psi_serial_mod
    implicit none
    class(psb_l_multivect_type), intent(inout)    :: x

    if (allocated(x%v)) call x%v%zero()

  end subroutine l_vect_zero

  subroutine l_vect_asb(m,n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    integer(psb_ipk_), intent(in)              :: m,n
    class(psb_l_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%v)) &
         & call x%v%asb(m,n,info)

  end subroutine l_vect_asb

  subroutine l_vect_sync(x)
    implicit none
    class(psb_l_multivect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%sync()

  end subroutine l_vect_sync

  subroutine l_vect_gthab(n,idx,alpha,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    integer(psb_lpk_) :: alpha, beta, y(:)
    class(psb_l_multivect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,alpha,beta,y)

  end subroutine l_vect_gthab

  subroutine l_vect_gthzv(n,idx,x,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    integer(psb_lpk_) ::  y(:)
    class(psb_l_multivect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,y)

  end subroutine l_vect_gthzv

  subroutine l_vect_gthzv_x(i,n,idx,x,y)
    use psi_serial_mod
    integer(psb_ipk_) :: i,n
    class(psb_i_base_vect_type) :: idx
    integer(psb_lpk_) ::  y(:)
    class(psb_l_multivect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(i,n,idx,y)

  end subroutine l_vect_gthzv_x

  subroutine l_vect_sctb(n,idx,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    integer(psb_lpk_) :: beta, x(:)
    class(psb_l_multivect_type) :: y

    if (allocated(y%v)) &
         &  call y%v%sct(n,idx,x,beta)

  end subroutine l_vect_sctb

  subroutine l_vect_sctb_x(i,n,idx,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    integer(psb_lpk_) :: beta, x(:)
    class(psb_l_multivect_type) :: y

    if (allocated(y%v)) &
         &  call y%v%sct(i,n,idx,x,beta)

  end subroutine l_vect_sctb_x

  subroutine l_vect_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (allocated(x%v)) then
      call x%v%free(info)
      if (info == 0) deallocate(x%v,stat=info)
    end if

  end subroutine l_vect_free

  subroutine l_vect_ins(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    integer(psb_ipk_), intent(in)               :: irl(:)
    integer(psb_lpk_), intent(in)        :: val(:,:)
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i

    info = 0
    if (.not.allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if

    call  x%v%ins(n,irl,val,dupl,info)

  end subroutine l_vect_ins


  subroutine l_vect_cnv(x,mold)
    class(psb_l_multivect_type), intent(inout) :: x
    class(psb_l_base_multivect_type), intent(in), optional :: mold
    class(psb_l_base_multivect_type), allocatable :: tmp
    integer(psb_ipk_) :: info

    if (present(mold)) then
      allocate(tmp,stat=info,mold=mold)
    else
      allocate(tmp,stat=info, mold=psb_l_get_base_multivect_default())
    endif
    if (allocated(x%v)) then
      call x%v%sync()
      if (info == psb_success_) call tmp%bld(x%v%v)
      call x%v%free(info)
    end if
    call move_alloc(tmp,x%v)
  end subroutine l_vect_cnv


end module psb_l_multivect_mod
