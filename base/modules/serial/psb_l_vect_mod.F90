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

  use psb_realloc_mod
  use psb_l_base_vect_mod
  use psb_i_vect_mod

  type psb_l_vect_type
    class(psb_l_base_vect_type), allocatable :: v
    integer(psb_ipk_) :: nrmv = 0
    integer(psb_ipk_) :: remote_build=psb_matbld_noremote_
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
       & l_vect_all, l_vect_reall, l_vect_zero,  l_vect_asb, &
       & l_vect_gthab, l_vect_gthzv, l_vect_sctb, &
       & l_vect_free, l_vect_ins_a, l_vect_ins_v, l_vect_bld_x, &
       & l_vect_bld_mn, l_vect_bld_en, l_vect_get_vect, &
       & l_vect_cnv, l_vect_set_scal, &
       & l_vect_set_vect, l_vect_clone, l_vect_sync, l_vect_is_host, &
       & l_vect_is_dev, l_vect_is_sync, l_vect_set_host, &
       & l_vect_set_dev, l_vect_set_sync, &
       & l_vect_set_remote_build, l_is_remote_build, &
       & l_vect_set_dupl, l_get_dupl, l_vect_set_nrmv, l_get_nrmv


  class(psb_l_base_vect_type), allocatable, target,&
       & save, private :: psb_l_base_vect_default

  
  interface 
    module function l_vect_get_dupl(x) result(res)
      implicit none
      class(psb_l_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function l_vect_get_dupl
  end interface

  interface 
    module subroutine l_vect_set_dupl(x,val)
      class(psb_l_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine l_vect_set_dupl
  end interface

  interface 
    module function l_vect_get_ncfs(x) result(res)
      class(psb_l_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function l_vect_get_ncfs
  end interface

  interface 
    module subroutine l_vect_set_ncfs(x,val)
      class(psb_l_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine l_vect_set_ncfs
  end interface

  interface 
    module function l_vect_get_state(x) result(res)
      class(psb_l_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function l_vect_get_state
  end interface

  interface 
    module function l_vect_is_null(x) result(res)
      class(psb_l_vect_type), intent(in) :: x
      logical :: res
    end function l_vect_is_null
  end interface

  interface 
    module function l_vect_is_bld(x) result(res)
      class(psb_l_vect_type), intent(in) :: x
      logical :: res
    end function l_vect_is_bld
  end interface

  interface 
    module function l_vect_is_upd(x) result(res)
      class(psb_l_vect_type), intent(in) :: x
      logical :: res
    end function l_vect_is_upd
  end interface

  interface 
    module function l_vect_is_asb(x) result(res)
      class(psb_l_vect_type), intent(in) :: x
      logical :: res
    end function l_vect_is_asb
  end interface

  interface 
    module subroutine  l_vect_set_state(n,x)
      class(psb_l_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine l_vect_set_state
  end interface

  interface 
    module subroutine  l_vect_set_null(x)
      class(psb_l_vect_type), intent(inout) :: x
    end subroutine l_vect_set_null
  end interface

  interface 
    module subroutine  l_vect_set_bld(x)
      class(psb_l_vect_type), intent(inout) :: x
    end subroutine l_vect_set_bld
  end interface

  interface 
    module subroutine  l_vect_set_upd(x)
      class(psb_l_vect_type), intent(inout) :: x
    end subroutine l_vect_set_upd
  end interface

  interface 
    module subroutine  l_vect_set_asb(x)
      class(psb_l_vect_type), intent(inout) :: x
    end subroutine l_vect_set_asb
  end interface

  interface 
    module function l_vect_get_nrmv(x) result(res)
      class(psb_l_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function l_vect_get_nrmv
  end interface

  interface 
    module subroutine l_vect_set_nrmv(x,val)
      class(psb_l_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: val
    end subroutine l_vect_set_nrmv
  end interface

  interface 
    module function l_vect_is_remote_build(x) result(res)
      class(psb_l_vect_type), intent(in) :: x
      logical :: res
    end function l_vect_is_remote_build
  end interface
  
  interface 
    module subroutine l_vect_set_remote_build(x,val)
      class(psb_l_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine l_vect_set_remote_build
  end interface
        
  interface psb_set_vect_default
    module subroutine  psb_l_set_vect_default(v)
      class(psb_l_base_vect_type), intent(in) :: v
    end subroutine psb_l_set_vect_default
  end interface

  interface psb_get_vect_default
    module function psb_l_get_vect_default(v) result(res)
      class(psb_l_vect_type), intent(in) :: v
      class(psb_l_base_vect_type), pointer :: res
    end function psb_l_get_vect_default
  end interface

  interface 
    module subroutine  psb_l_clear_vect_default()
    end subroutine psb_l_clear_vect_default
  end interface

  interface 
    module function psb_l_get_base_vect_default() result(res)
      class(psb_l_base_vect_type), pointer :: res
    end function psb_l_get_base_vect_default
  end interface

  interface 
    module subroutine l_vect_clone(x,y,info)
      class(psb_l_vect_type), intent(inout) :: x
      class(psb_l_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)        :: info
    end subroutine l_vect_clone
  end interface

  interface 
    module subroutine l_vect_bld_x(x,invect,mold,scratch)
      integer(psb_lpk_), intent(in)          :: invect(:)
      class(psb_l_vect_type), intent(inout) :: x
      class(psb_l_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine l_vect_bld_x
  end interface

  interface 
    module subroutine l_vect_bld_mn(x,n,mold,scratch)
      integer(psb_mpk_), intent(in) :: n
      class(psb_l_vect_type), intent(inout) :: x
      class(psb_l_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine l_vect_bld_mn
  end interface

  interface 
    module subroutine l_vect_bld_en(x,n,mold,scratch)
      integer(psb_epk_), intent(in) :: n
      class(psb_l_vect_type), intent(inout) :: x
      class(psb_l_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine l_vect_bld_en
  end interface

  interface 
    module function  l_vect_get_vect(x,n) result(res)
      class(psb_l_vect_type), intent(inout)  :: x
      integer(psb_lpk_), allocatable                 :: res(:)
      integer(psb_ipk_), optional :: n
    end function l_vect_get_vect
  end interface

  interface 
    module subroutine l_vect_set_scal(x,val,first,last)
      class(psb_l_vect_type), intent(inout)  :: x
      integer(psb_lpk_), intent(in) :: val
      integer(psb_ipk_), optional :: first, last
    end subroutine l_vect_set_scal
  end interface

  interface 
    module subroutine l_vect_set_vect(x,val,first,last)
      class(psb_l_vect_type), intent(inout) :: x
      integer(psb_lpk_), intent(in)         :: val(:)
      integer(psb_ipk_), optional :: first, last
    end subroutine l_vect_set_vect
  end interface

  interface 
    module subroutine l_vect_check_addr(x)
      class(psb_l_vect_type), intent(inout) :: x
    end subroutine l_vect_check_addr
  end interface

  interface 
    module function l_vect_get_nrows(x) result(res)
      class(psb_l_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function l_vect_get_nrows
  end interface

  interface 
    module function l_vect_sizeof(x) result(res)
      class(psb_l_vect_type), intent(in) :: x
      integer(psb_epk_) :: res
    end function l_vect_sizeof
  end interface

  interface 
    module function l_vect_get_fmt(x) result(res)
      class(psb_l_vect_type), intent(in) :: x
      character(len=5) :: res
    end function l_vect_get_fmt
  end interface

  interface 
    module subroutine l_vect_all(n, x, info, mold)
      integer(psb_ipk_), intent(in)           :: n
      class(psb_l_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)      :: info
      class(psb_l_base_vect_type), intent(in), optional :: mold
    end subroutine l_vect_all
  end interface

  interface 
    module subroutine l_vect_reinit(x, info, clear)
      class(psb_l_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: clear
    end subroutine l_vect_reinit
  end interface
 
  interface 
    module subroutine l_vect_reall(n, x, info)
      integer(psb_ipk_), intent(in)         :: n
      class(psb_l_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)        :: info
    end subroutine l_vect_reall
  end interface

  interface 
    module subroutine l_vect_zero(x)
      class(psb_l_vect_type), intent(inout)    :: x
    end subroutine l_vect_zero
  end interface

  interface 
    module subroutine l_vect_asb(n, x, info, scratch)
      integer(psb_ipk_), intent(in)         :: n
      class(psb_l_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional        :: scratch
    end subroutine l_vect_asb
  end interface

  interface 
    module subroutine l_vect_gthab(n,idx,alpha,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_lpk_) :: alpha, beta, y(:)
      class(psb_l_vect_type) :: x
    end subroutine l_vect_gthab
  end interface

  interface 
    module subroutine l_vect_gthzv(n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_lpk_) ::  y(:)
      class(psb_l_vect_type) :: x
    end subroutine l_vect_gthzv
  end interface

  interface 
    module subroutine l_vect_sctb(n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_lpk_) :: beta, x(:)
      class(psb_l_vect_type) :: y
    end subroutine l_vect_sctb
  end interface

  interface 
    module subroutine l_vect_free(x, info)
      class(psb_l_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine l_vect_free
  end interface
  
  interface 
    module subroutine l_vect_ins_a(n,irl,val,x,maxr,info)
      class(psb_l_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n, maxr
      integer(psb_ipk_), intent(in)               :: irl(:)
      integer(psb_lpk_), intent(in)        :: val(:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine l_vect_ins_a
  end interface

  interface 
    module subroutine l_vect_ins_v(n,irl,val,x,maxr,info)
      class(psb_l_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n, maxr
      class(psb_i_vect_type), intent(inout)       :: irl
      class(psb_l_vect_type), intent(inout)       :: val
      integer(psb_ipk_), intent(out)              :: info
    end subroutine l_vect_ins_v
  end interface

  interface 
    module subroutine l_vect_cnv(x,mold)
      class(psb_l_vect_type), intent(inout) :: x
      class(psb_l_base_vect_type), intent(in), optional :: mold
      class(psb_l_base_vect_type), allocatable :: tmp
    end subroutine l_vect_cnv
  end interface
  
  interface 
    module subroutine l_vect_sync(x)
      class(psb_l_vect_type), intent(inout) :: x
    end subroutine l_vect_sync
  end interface

  interface 
    module subroutine l_vect_set_sync(x)
      class(psb_l_vect_type), intent(inout) :: x
    end subroutine l_vect_set_sync
  end interface

  interface 
    module subroutine l_vect_set_host(x)
      class(psb_l_vect_type), intent(inout) :: x
    end subroutine l_vect_set_host
  end interface

  interface 
    module subroutine l_vect_set_dev(x)
      class(psb_l_vect_type), intent(inout) :: x
    end subroutine l_vect_set_dev
  end interface

  interface 
    module function l_vect_is_sync(x) result(res)
      logical :: res
      class(psb_l_vect_type), intent(inout) :: x
    end function l_vect_is_sync
  end interface

  interface 
    module function l_vect_is_host(x) result(res)
      logical :: res
      class(psb_l_vect_type), intent(inout) :: x
    end function l_vect_is_host
  end interface

  interface 
    module function l_vect_is_dev(x) result(res)
      logical :: res
      class(psb_l_vect_type), intent(inout) :: x
    end function l_vect_is_dev
  end interface




contains

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

!!$  function l_vect_get_nrows(x) result(res)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(in) :: x
!!$    integer(psb_ipk_) :: res
!!$    res = 0
!!$    if (allocated(x%v)) res = x%v%get_nrows()
!!$  end function l_vect_get_nrows
!!$
!!$  function l_vect_sizeof(x) result(res)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(in) :: x
!!$    integer(psb_epk_) :: res
!!$    res = 0
!!$    if (allocated(x%v)) res = x%v%sizeof()
!!$  end function l_vect_sizeof
!!$
!!$  function l_vect_get_fmt(x) result(res)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(in) :: x
!!$    character(len=5) :: res
!!$    res = 'NULL'
!!$    if (allocated(x%v)) res = x%v%get_fmt()
!!$  end function l_vect_get_fmt
!!$
!!$  subroutine l_vect_all(n, x, info, mold)
!!$
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(out)      :: info
!!$    class(psb_l_base_vect_type), intent(in), optional :: mold
!!$
!!$    if (allocated(x%v)) &
!!$         & call x%free(info)
!!$
!!$    if (present(mold)) then
!!$      allocate(x%v,stat=info,mold=mold)
!!$    else
!!$      allocate(psb_l_base_vect_type :: x%v,stat=info)
!!$    endif
!!$    if (info == 0) then
!!$      call x%v%all(n,info)
!!$    else
!!$      info = psb_err_alloc_dealloc_
!!$    end if
!!$    call x%set_bld()
!!$  end subroutine l_vect_all
!!$
!!$  subroutine l_vect_reinit(x, info, clear)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(out)      :: info
!!$    logical, intent(in), optional       :: clear
!!$
!!$    if (allocated(x%v)) call x%v%reinit(info,clear)
!!$    call x%set_upd()
!!$
!!$  end subroutine l_vect_reinit
!!$ 
!!$  subroutine l_vect_reall(n, x, info)
!!$
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)         :: n
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(out)        :: info
!!$
!!$    info = 0
!!$    if (.not.allocated(x%v)) &
!!$         & call x%all(n,info)
!!$    if (info == 0) &
!!$         & call x%asb(n,info)
!!$
!!$  end subroutine l_vect_reall
!!$
!!$  subroutine l_vect_zero(x)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout)    :: x
!!$
!!$    if (allocated(x%v)) call x%v%zero()
!!$
!!$  end subroutine l_vect_zero
!!$
!!$  subroutine l_vect_asb(n, x, info, scratch)
!!$    use psi_serial_mod
!!$    use psb_realloc_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)         :: n
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(out)        :: info
!!$    logical, intent(in), optional        :: scratch
!!$
!!$    if (allocated(x%v)) then
!!$      call x%v%asb(n,info,scratch=scratch)
!!$      call x%set_asb()
!!$    end if
!!$  end subroutine l_vect_asb
!!$
!!$  subroutine l_vect_gthab(n,idx,alpha,x,beta,y)
!!$    use psi_serial_mod
!!$    integer(psb_mpk_) :: n
!!$    integer(psb_ipk_) :: idx(:)
!!$    integer(psb_lpk_) :: alpha, beta, y(:)
!!$    class(psb_l_vect_type) :: x
!!$
!!$    if (allocated(x%v)) &
!!$         &  call x%v%gth(n,idx,alpha,beta,y)
!!$
!!$  end subroutine l_vect_gthab
!!$
!!$  subroutine l_vect_gthzv(n,idx,x,y)
!!$    use psi_serial_mod
!!$    integer(psb_mpk_) :: n
!!$    integer(psb_ipk_) :: idx(:)
!!$    integer(psb_lpk_) ::  y(:)
!!$    class(psb_l_vect_type) :: x
!!$
!!$    if (allocated(x%v)) &
!!$         &  call x%v%gth(n,idx,y)
!!$
!!$  end subroutine l_vect_gthzv
!!$
!!$  subroutine l_vect_sctb(n,idx,x,beta,y)
!!$    use psi_serial_mod
!!$    integer(psb_mpk_) :: n
!!$    integer(psb_ipk_) :: idx(:)
!!$    integer(psb_lpk_) :: beta, x(:)
!!$    class(psb_l_vect_type) :: y
!!$
!!$    if (allocated(y%v)) &
!!$         &  call y%v%sct(n,idx,x,beta)
!!$
!!$  end subroutine l_vect_sctb
!!$
!!$  subroutine l_vect_free(x, info)
!!$    use psi_serial_mod
!!$    use psb_realloc_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    info = 0
!!$    if (allocated(x%v)) then
!!$      call x%v%free(info)
!!$      if (info == 0) deallocate(x%v,stat=info)
!!$    end if
!!$
!!$  end subroutine l_vect_free
!!$
!!$  subroutine l_vect_ins_a(n,irl,val,x,maxr,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    integer(psb_ipk_), intent(in)               :: n, maxr
!!$    integer(psb_ipk_), intent(in)               :: irl(:)
!!$    integer(psb_lpk_), intent(in)        :: val(:)
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    integer(psb_ipk_) :: i, dupl
!!$
!!$    info = 0
!!$    if (.not.allocated(x%v)) then
!!$      info = psb_err_invalid_vect_state_
!!$      return
!!$    end if
!!$    dupl = x%get_dupl()
!!$    call  x%v%ins(n,irl,val,dupl,maxr,info)
!!$
!!$  end subroutine l_vect_ins_a
!!$
!!$  subroutine l_vect_ins_v(n,irl,val,x,maxr,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    integer(psb_ipk_), intent(in)               :: n, maxr
!!$    class(psb_i_vect_type), intent(inout)       :: irl
!!$    class(psb_l_vect_type), intent(inout)       :: val
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    integer(psb_ipk_) :: i, dupl
!!$
!!$    info = 0
!!$    if (.not.(allocated(x%v).and.allocated(irl%v).and.allocated(val%v))) then
!!$      info = psb_err_invalid_vect_state_
!!$      return
!!$    end if
!!$    dupl = x%get_dupl()
!!$    call  x%v%ins(n,irl%v,val%v,dupl,maxr,info)
!!$
!!$  end subroutine l_vect_ins_v
!!$
!!$
!!$  subroutine l_vect_cnv(x,mold)
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$    class(psb_l_base_vect_type), intent(in), optional :: mold
!!$    class(psb_l_base_vect_type), allocatable :: tmp
!!$
!!$    integer(psb_ipk_) :: info
!!$
!!$    info = psb_success_
!!$    if (present(mold)) then
!!$      allocate(tmp,stat=info,mold=mold)
!!$    else
!!$      allocate(tmp,stat=info,mold=psb_l_get_base_vect_default())
!!$    end if
!!$    if (allocated(x%v)) then
!!$      if (allocated(x%v%v)) then 
!!$        call x%v%sync()
!!$        if (info == psb_success_) call tmp%bld(x%v%v)
!!$        call x%v%base_cpy(tmp)
!!$        call x%v%free(info)
!!$      endif
!!$    end if
!!$    call move_alloc(tmp,x%v)
!!$
!!$  end subroutine l_vect_cnv
!!$
!!$
!!$  subroutine l_vect_sync(x)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$
!!$    if (allocated(x%v)) &
!!$         & call x%v%sync()
!!$
!!$  end subroutine l_vect_sync
!!$
!!$  subroutine l_vect_set_sync(x)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$
!!$    if (allocated(x%v)) &
!!$         & call x%v%set_sync()
!!$
!!$  end subroutine l_vect_set_sync
!!$
!!$  subroutine l_vect_set_host(x)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$
!!$    if (allocated(x%v)) &
!!$         & call x%v%set_host()
!!$
!!$  end subroutine l_vect_set_host
!!$
!!$  subroutine l_vect_set_dev(x)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$
!!$    if (allocated(x%v)) &
!!$         & call x%v%set_dev()
!!$
!!$  end subroutine l_vect_set_dev
!!$
!!$  function l_vect_is_sync(x) result(res)
!!$    implicit none
!!$    logical :: res
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$
!!$    res = .true.
!!$    if (allocated(x%v)) &
!!$         & res = x%v%is_sync()
!!$
!!$  end function l_vect_is_sync
!!$
!!$  function l_vect_is_host(x) result(res)
!!$    implicit none
!!$    logical :: res
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$
!!$    res = .true.
!!$    if (allocated(x%v)) &
!!$         & res = x%v%is_host()
!!$
!!$  end function l_vect_is_host
!!$
!!$  function l_vect_is_dev(x) result(res)
!!$    implicit none
!!$    logical :: res
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$
!!$    res = .false.
!!$    if (allocated(x%v)) &
!!$         & res =  x%v%is_dev()
!!$
!!$  end function l_vect_is_dev
!!$
!!$  @NOTINTS@
!!$
!!$  function l_vect_get_entry(x,index) result(res)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)        :: index
!!$    integer(psb_lpk_) :: res
!!$    res = lzero
!!$    if (allocated(x%v)) res = x%v%get_entry(index)
!!$  end function l_vect_get_entry
!!$
!!$  subroutine l_vect_set_entry(x,index,val)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)        :: index
!!$    integer(psb_lpk_) :: val
!!$    if (allocated(x%v)) call x%v%set_entry(index,val)
!!$  end subroutine l_vect_set_entry
!!$
!!$  function l_vect_dot_v(n,x,y) result(res)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout) :: x, y
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_lpk_)                :: res
!!$
!!$    res = lzero
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & res = x%v%dot(n,y%v)
!!$
!!$  end function l_vect_dot_v
!!$
!!$  function l_vect_dot_a(n,x,y) result(res)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$    integer(psb_lpk_), intent(in)    :: y(:)
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_lpk_)                :: res
!!$
!!$    res = lzero
!!$    if (allocated(x%v)) &
!!$         & res = x%v%dot_a(n,y)
!!$
!!$  end function l_vect_dot_a
!!$
!!$  subroutine l_vect_axpby_v(m,alpha, x, beta, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    integer(psb_lpk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    if (allocated(x%v).and.allocated(y%v)) then
!!$      call y%v%axpby(m,alpha,x%v,beta,info)
!!$    else
!!$      info = psb_err_invalid_vect_state_
!!$    end if
!!$
!!$  end subroutine l_vect_axpby_v
!!$
!!$  subroutine l_vect_axpby_v2(m,alpha, x, beta, y, z, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)            :: m
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    class(psb_l_vect_type), intent(inout)  :: z
!!$    integer(psb_lpk_), intent (in)             :: alpha, beta
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$    if (allocated(x%v).and.allocated(y%v)) then
!!$      call z%v%axpby(m,alpha,x%v,beta,y%v,info)
!!$    else
!!$      info = psb_err_invalid_vect_state_
!!$    end if
!!$
!!$  end subroutine l_vect_axpby_v2
!!$
!!$  subroutine l_vect_axpby_a(m,alpha, x, beta, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    integer(psb_lpk_), intent(in)        :: x(:)
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    integer(psb_lpk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    if (allocated(y%v)) &
!!$         & call y%v%axpby(m,alpha,x,beta,info)
!!$
!!$  end subroutine l_vect_axpby_a
!!$
!!$  subroutine l_vect_axpby_a2(m,alpha, x, beta, y, z, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    integer(psb_lpk_), intent(in)                 :: x(:)
!!$    integer(psb_lpk_), intent(in)                 :: y(:)
!!$    class(psb_l_vect_type), intent(inout)     :: z
!!$    integer(psb_lpk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    if (allocated(z%v)) &
!!$         & call z%v%axpby(m,alpha,x,beta,y,info)
!!$
!!$  end subroutine l_vect_axpby_a2
!!$
!!$  subroutine l_vect_upd_xyz(m,alpha,beta,gamma,delta,x, y, z, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)            :: m
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    class(psb_l_vect_type), intent(inout)  :: z
!!$    integer(psb_lpk_), intent (in)     :: alpha, beta, gamma, delta
!!$    integer(psb_ipk_), intent(out)   :: info
!!$
!!$    if (allocated(z%v)) &
!!$         call z%v%upd_xyz(m,alpha,beta,gamma,delta,x%v,y%v,info)
!!$    
!!$  end subroutine l_vect_upd_xyz
!!$
!!$  subroutine l_vect_xyzw(m,a,b,c,d,e,f,x, y, z, w, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)            :: m
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    class(psb_l_vect_type), intent(inout)  :: z
!!$    class(psb_l_vect_type), intent(inout)  :: w
!!$    integer(psb_lpk_), intent (in)     :: a, b, c, d, e, f
!!$    integer(psb_ipk_), intent(out)   :: info
!!$
!!$    if (allocated(w%v)) &
!!$         call w%v%xyzw(m,a,b,c,d,e,f,x%v,y%v,z%v,info)
!!$    
!!$  end subroutine l_vect_xyzw
!!$
!!$
!!$  subroutine l_vect_mlt_v(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & call y%v%mlt(x%v,info)
!!$
!!$  end subroutine l_vect_mlt_v
!!$
!!$  subroutine l_vect_mlt_a(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_lpk_), intent(in)        :: x(:)
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$
!!$    info = 0
!!$    if (allocated(y%v)) &
!!$         & call y%v%mlt(x,info)
!!$
!!$  end subroutine l_vect_mlt_a
!!$
!!$
!!$  subroutine l_vect_mlt_a_2(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_lpk_), intent(in)         :: alpha,beta
!!$    integer(psb_lpk_), intent(in)         :: y(:)
!!$    integer(psb_lpk_), intent(in)         :: x(:)
!!$    class(psb_l_vect_type), intent(inout) :: z
!!$    integer(psb_ipk_), intent(out)                  :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(z%v)) &
!!$         & call z%v%mlt(alpha,x,y,beta,info)
!!$
!!$  end subroutine l_vect_mlt_a_2
!!$
!!$  subroutine l_vect_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_lpk_), intent(in)          :: alpha,beta
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    class(psb_l_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)                   :: info
!!$    character(len=1), intent(in), optional :: conjgx, conjgy
!!$
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v).and.&
!!$         & allocated(z%v)) &
!!$         & call z%v%mlt(alpha,x%v,y%v,beta,info,conjgx,conjgy)
!!$
!!$  end subroutine l_vect_mlt_v_2
!!$
!!$  subroutine l_vect_mlt_av(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_lpk_), intent(in)        :: alpha,beta
!!$    integer(psb_lpk_), intent(in)        :: x(:)
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    class(psb_l_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(z%v).and.allocated(y%v)) &
!!$         & call z%v%mlt(alpha,x,y%v,beta,info)
!!$
!!$  end subroutine l_vect_mlt_av
!!$
!!$  subroutine l_vect_mlt_va(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_lpk_), intent(in)        :: alpha,beta
!!$    integer(psb_lpk_), intent(in)        :: y(:)
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$
!!$    if (allocated(z%v).and.allocated(x%v)) &
!!$         & call z%v%mlt(alpha,x%v,y,beta,info)
!!$
!!$  end subroutine l_vect_mlt_va
!!$
!!$  subroutine l_vect_div_v(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & call x%v%div(y%v,info)
!!$
!!$  end subroutine l_vect_div_v
!!$
!!$  subroutine l_vect_div_v2( x, y, z, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    class(psb_l_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v).and.allocated(z%v)) &
!!$         & call z%v%div(x%v,y%v,info)
!!$
!!$  end subroutine l_vect_div_v2
!!$
!!$  subroutine l_vect_div_v_check(x, y, info, flag)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$    logical, intent(in) :: flag
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & call x%v%div(y%v,info,flag)
!!$
!!$  end subroutine l_vect_div_v_check
!!$
!!$  subroutine l_vect_div_v2_check(x, y, z, info, flag)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    class(psb_l_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$    logical, intent(in) :: flag
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v).and.allocated(z%v)) &
!!$         & call z%v%div(x%v,y%v,info,flag)
!!$
!!$  end subroutine l_vect_div_v2_check
!!$
!!$  subroutine l_vect_div_a2(x, y, z, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_lpk_), intent(in) :: x(:)
!!$    integer(psb_lpk_), intent(in)    :: y(:)
!!$    class(psb_l_vect_type), intent(inout) :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(z%v)) &
!!$         & call z%v%div(x,y,info)
!!$
!!$  end subroutine l_vect_div_a2
!!$
!!$  subroutine l_vect_div_a2_check(x, y, z, info,flag)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_lpk_), intent(in) :: x(:)
!!$    integer(psb_lpk_), intent(in)    :: y(:)
!!$    class(psb_l_vect_type), intent(inout) :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$    logical, intent(in) :: flag
!!$
!!$    info = 0
!!$    if (allocated(z%v)) &
!!$         & call z%v%div(x,y,info,flag)
!!$
!!$  end subroutine l_vect_div_a2_check
!!$
!!$  subroutine l_vect_inv_v(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & call y%v%inv(x%v,info)
!!$
!!$  end subroutine l_vect_inv_v
!!$
!!$  subroutine l_vect_inv_v_check(x, y, info, flag)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$    logical, intent(in) :: flag
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & call y%v%inv(x%v,info,flag)
!!$
!!$  end subroutine l_vect_inv_v_check
!!$
!!$  subroutine l_vect_inv_a2(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_lpk_), intent(inout) :: x(:)
!!$    class(psb_l_vect_type), intent(inout) :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(y%v)) &
!!$         & call y%v%inv(x,info)
!!$
!!$  end subroutine l_vect_inv_a2
!!$
!!$  subroutine l_vect_inv_a2_check(x, y, info,flag)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_lpk_), intent(inout) :: x(:)
!!$    class(psb_l_vect_type), intent(inout) :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$    logical, intent(in) :: flag
!!$
!!$    info = 0
!!$    if (allocated(y%v)) &
!!$         & call y%v%inv(x,info,flag)
!!$
!!$  end subroutine l_vect_inv_a2_check
!!$
!!$  subroutine l_vect_acmp_a2(x,c,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_lpk_), intent(in)              :: c
!!$    integer(psb_lpk_), intent(inout)           :: x(:)
!!$    class(psb_l_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$    info = 0
!!$    if (allocated(z%v)) &
!!$         & call z%acmp(x,c,info)
!!$
!!$  end subroutine l_vect_acmp_a2
!!$
!!$  subroutine l_vect_acmp_v2(x,c,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_lpk_), intent(in)              :: c
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(z%v)) &
!!$         & call z%v%acmp(x%v,c,info)
!!$
!!$  end subroutine l_vect_acmp_v2
!!$
!!$  subroutine l_vect_scal(alpha, x)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    integer(psb_lpk_), intent (in)       :: alpha
!!$
!!$    if (allocated(x%v)) call x%v%scal(alpha)
!!$
!!$  end subroutine l_vect_scal
!!$
!!$  subroutine l_vect_absval1(x)
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$
!!$    if (allocated(x%v)) &
!!$         &  call x%v%absval()
!!$
!!$  end subroutine l_vect_absval1
!!$
!!$  subroutine l_vect_absval2(x,y)
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$
!!$    if (allocated(x%v)) then
!!$      if (.not.allocated(y%v))  call y%bld(psb_size(x%v%v))
!!$      call x%v%absval(y%v)
!!$    end if
!!$  end subroutine l_vect_absval2
!!$
!!$  function l_vect_nrm2(n,x) result(res)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_lpk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%nrm2(n)
!!$    else
!!$      res = lzero
!!$    end if
!!$
!!$  end function l_vect_nrm2
!!$
!!$  function l_vect_nrm2_weight(n,x,w,aux) result(res)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$    class(psb_l_vect_type), intent(inout) :: w
!!$    class(psb_l_vect_type), intent(inout), optional :: aux
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_lpk_)                        :: res
!!$    integer(psb_ipk_)                       :: info
!!$
!!$    ! Temp vectors
!!$    type(psb_l_vect_type) :: wtemp
!!$
!!$    info = 0
!!$    if( allocated(w%v) ) then
!!$      if (.not.present(aux)) then
!!$        allocate(wtemp%v, mold=w%v)
!!$        call wtemp%v%bld(w%get_vect())
!!$      else
!!$        call psb_geaxpby(n,lone,w%v%v,lzero,aux%v%v,info)
!!$      end if
!!$    else
!!$      info = -1
!!$    end if
!!$    if (info /= 0 ) then
!!$      res = -lone
!!$      return
!!$    end if
!!$
!!$    if (allocated(x%v)) then
!!$      if (.not.present(aux)) then
!!$        call wtemp%v%mlt(x%v,info)
!!$        res = wtemp%v%nrm2(n)
!!$      else
!!$        call aux%v%mlt(x%v,info)
!!$        res = aux%v%nrm2(n)
!!$      end if
!!$    else
!!$      res = lzero
!!$    end if
!!$
!!$    if (.not.present(aux)) then
!!$      call wtemp%free(info)
!!$    end if
!!$
!!$  end function l_vect_nrm2_weight
!!$ 
!!$  function l_vect_nrm2_weight_mask(n,x,w,id,info,aux) result(res)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$    class(psb_l_vect_type), intent(inout) :: w
!!$    class(psb_l_vect_type), intent(inout) :: id
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_lpk_)                        :: res
!!$    integer(psb_ipk_), intent(out)          :: info
!!$    class(psb_l_vect_type), intent(inout), optional :: aux
!!$
!!$    ! Temp vectors
!!$    type(psb_l_vect_type) :: wtemp
!!$
!!$    info = 0 
!!$    if( allocated(w%v) ) then
!!$      if (.not.present(aux)) then
!!$        allocate(wtemp%v, mold=w%v)
!!$        call wtemp%v%bld(w%get_vect())
!!$      else
!!$        call psb_geaxpby(n,lone,w%v%v,lzero,aux%v%v,info)
!!$      end if
!!$    else
!!$      info = -1
!!$    end if
!!$    if (info /= 0 ) then
!!$      res = -lone
!!$      return
!!$    end if
!!$
!!$
!!$    if (allocated(x%v).and.allocated(id%v)) then
!!$      if (.not.present(aux)) then
!!$        where( abs(id%v%v) <= lzero) wtemp%v%v = lzero
!!$        call wtemp%set_host()
!!$        call wtemp%v%mlt(x%v,info)
!!$        res = wtemp%v%nrm2(n)
!!$      else
!!$        where( abs(id%v%v) <= lzero) aux%v%v = lzero
!!$        call aux%set_host()
!!$        call aux%v%mlt(x%v,info)
!!$        res = aux%v%nrm2(n)
!!$      end if
!!$    else
!!$      res = lzero
!!$    end if
!!$
!!$    if (.not.present(aux)) then
!!$      call wtemp%free(info)
!!$    end if
!!$
!!$  end function l_vect_nrm2_weight_mask
!!$
!!$  function l_vect_amax(n,x) result(res)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_lpk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%amax(n)
!!$    else
!!$      res = lzero
!!$    end if
!!$
!!$  end function l_vect_amax
!!$
!!$  @NOTCPLXS@
!!$  function l_vect_min(n,x) result(res)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_lpk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%minreal(n)
!!$    else
!!$      res = HUGE(lone)
!!$    end if
!!$
!!$  end function l_vect_min
!!$  @NOTCPLXE@
!!$
!!$  function l_vect_asum(n,x) result(res)
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_lpk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%asum(n)
!!$    else
!!$      res = lzero
!!$    end if
!!$
!!$  end function l_vect_asum
!!$
!!$  @NOTINTE@
!!$
!!$  @NOTINTS@
!!$  @NOTCPLXS@
!!$  subroutine l_vect_mask_a(c,x,m,t,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_lpk_), intent(inout)          :: c(:)
!!$    integer(psb_lpk_), intent(inout)          :: x(:)
!!$    logical, intent(out)                     :: t;
!!$    class(psb_l_vect_type), intent(inout)  :: m
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$    info = 0
!!$    if (allocated(m%v)) &
!!$         & call m%mask(c,x,t,info)
!!$
!!$  end subroutine l_vect_mask_a
!!$
!!$  subroutine l_vect_mask_v(c,x,m,t,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout)  :: c
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: m
!!$    logical, intent(out)                     :: t;
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(c%v)) &
!!$         & call m%v%mask(x%v,c%v,t,info)
!!$
!!$  end subroutine l_vect_mask_v
!!$
!!$  function l_vect_minquotient_v(x, y, info) result(z)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: y
!!$    integer(psb_lpk_)                         :: z
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & z = x%v%minquotient(y%v,info)
!!$
!!$  end function l_vect_minquotient_v
!!$
!!$  function l_vect_minquotient_a2(x, y, info) result(z)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    integer(psb_lpk_), intent(inout)           :: y(:)
!!$    integer(psb_ipk_), intent(out)           :: info
!!$    integer(psb_lpk_)                         :: z
!!$
!!$    info = 0
!!$    z =  x%v%minquotient(y,info)
!!$
!!$  end function l_vect_minquotient_a2
!!$
!!$
!!$  @NOTINTE@
!!$  @NOTCPLXE@
!!$
!!$  @NOTINTS@
!!$  subroutine l_vect_addconst_a2(x,b,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_lpk_), intent(in)             :: b
!!$    integer(psb_lpk_), intent(inout)           :: x(:)
!!$    class(psb_l_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$    info = 0
!!$    if (allocated(z%v)) &
!!$         & call z%addconst(x,b,info)
!!$
!!$  end subroutine l_vect_addconst_a2
!!$
!!$  subroutine l_vect_addconst_v2(x,b,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_lpk_), intent(in)             :: b
!!$    class(psb_l_vect_type), intent(inout)  :: x
!!$    class(psb_l_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(z%v)) &
!!$         & call z%v%addconst(x%v,b,info)
!!$
!!$  end subroutine l_vect_addconst_v2
!!$  @NOTINTE@

end module psb_l_vect_mod

module psb_l_multivect_mod

  use psb_l_base_multivect_mod
  use psb_const_mod
  use psb_i_vect_mod


  !private

  type psb_l_multivect_type
    class(psb_l_base_multivect_type), allocatable :: v
    integer(psb_ipk_) :: nrmv = 0
    integer(psb_ipk_) :: remote_build=psb_matbld_noremote_
    integer(psb_ipk_) :: dupl = psb_dupl_add_
    integer(psb_lpk_), allocatable :: rmtv(:,:)
  contains
    procedure, pass(x) :: get_nrows => l_mvect_get_nrows
    procedure, pass(x) :: get_ncols => l_mvect_get_ncols
    procedure, pass(x) :: sizeof   => l_mvect_sizeof
    procedure, pass(x) :: get_fmt  => l_mvect_get_fmt
    procedure, pass(x) :: is_remote_build => l_mvect_is_remote_build
    procedure, pass(x) :: set_remote_build => l_mvect_set_remote_build
    procedure, pass(x) :: get_dupl => l_mvect_get_dupl
    procedure, pass(x) :: set_dupl => l_mvect_set_dupl

    procedure, pass(x) :: all      => l_mvect_all
    procedure, pass(x) :: reall    => l_mvect_reall
    procedure, pass(x) :: zero     => l_mvect_zero
    procedure, pass(x) :: asb      => l_mvect_asb
    procedure, pass(x) :: sync     => l_mvect_sync
    procedure, pass(x) :: free     => l_mvect_free
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

  public  :: psb_l_multivect, psb_l_multivect_type,&
       & psb_set_multivect_default, psb_get_multivect_default, &
       & psb_l_base_multivect_type

  private
  interface psb_l_multivect
    module procedure constructor, size_const
  end interface psb_l_multivect

  class(psb_l_base_multivect_type), allocatable, target,&
       & save, private :: psb_l_base_multivect_default

  
  interface 
    module function l_mvect_get_dupl(x) result(res)
      class(psb_l_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function l_mvect_get_dupl
  end interface

  interface 
    module subroutine l_mvect_set_dupl(x,val)
      class(psb_l_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine l_mvect_set_dupl
  end interface

  interface 
    module function l_mvect_is_remote_build(x) result(res)
      class(psb_l_multivect_type), intent(in) :: x
      logical :: res
    end function l_mvect_is_remote_build
  end interface

  interface 
    module subroutine l_mvect_set_remote_build(x,val)
      class(psb_l_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine l_mvect_set_remote_build
  end interface

  interface 
    module subroutine  psb_l_set_multivect_default(v)
      class(psb_l_base_multivect_type), intent(in) :: v
    end subroutine psb_l_set_multivect_default
  end interface

  interface 
    module function psb_l_get_multivect_default(v) result(res)
      class(psb_l_multivect_type), intent(in) :: v
      class(psb_l_base_multivect_type), pointer :: res
    end function psb_l_get_multivect_default
  end interface

  interface 
    module function psb_l_get_base_multivect_default() result(res)
      class(psb_l_base_multivect_type), pointer :: res
    end function psb_l_get_base_multivect_default
  end interface

  interface 
    module subroutine l_mvect_clone(x,y,info)
      class(psb_l_multivect_type), intent(inout) :: x
      class(psb_l_multivect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)        :: info
    end subroutine l_mvect_clone
  end interface

  interface 
    module subroutine l_mvect_bld_x(x,invect,mold)
      integer(psb_lpk_), intent(in)          :: invect(:,:)
      class(psb_l_multivect_type), intent(out) :: x
      class(psb_l_base_multivect_type), intent(in), optional :: mold
    end subroutine l_mvect_bld_x
  end interface


  interface 
    module subroutine l_mvect_bld_n(x,m,n,mold,scratch)
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_l_multivect_type), intent(out) :: x
      class(psb_l_base_multivect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine l_mvect_bld_n
  end interface

  interface 
    module function  l_mvect_get_vect(x) result(res)
      class(psb_l_multivect_type), intent(inout)  :: x
      integer(psb_lpk_), allocatable                 :: res(:,:)
    end function l_mvect_get_vect
  end interface

  interface 
    module subroutine l_mvect_set_scal(x,val)
      class(psb_l_multivect_type), intent(inout)  :: x
      integer(psb_lpk_), intent(in) :: val
    end subroutine l_mvect_set_scal
  end interface

  interface 
    module subroutine l_mvect_set_vect(x,val)
      class(psb_l_multivect_type), intent(inout) :: x
      integer(psb_lpk_), intent(in)         :: val(:,:)
    end subroutine l_mvect_set_vect
  end interface

  interface 
    module function l_mvect_get_nrows(x) result(res)
      class(psb_l_multivect_type), intent(in) :: x
      integer(psb_ipk_)  :: res
    end function l_mvect_get_nrows
  end interface

  interface 
    module function l_mvect_get_ncols(x) result(res)
      class(psb_l_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function l_mvect_get_ncols
  end interface

  interface 
    module function l_mvect_sizeof(x) result(res)
      class(psb_l_multivect_type), intent(in) :: x
      integer(psb_epk_) :: res
    end function l_mvect_sizeof
  end interface

  interface 
    module function l_mvect_get_fmt(x) result(res)
      class(psb_l_multivect_type), intent(in) :: x
      character(len=5) :: res
    end function l_mvect_get_fmt
  end interface

  interface 
    module subroutine l_mvect_all(m,n, x, info, mold)
      integer(psb_ipk_), intent(in)       :: m,n
      class(psb_l_multivect_type), intent(out) :: x
      class(psb_l_base_multivect_type), intent(in), optional :: mold
      integer(psb_ipk_), intent(out)      :: info
    end subroutine l_mvect_all
  end interface

  interface 
    module subroutine l_mvect_reall(m,n, x, info)
      integer(psb_ipk_), intent(in)         :: m,n
      class(psb_l_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)        :: info
    end subroutine l_mvect_reall
  end interface

  interface 
    module subroutine l_mvect_zero(x)
      class(psb_l_multivect_type), intent(inout)    :: x
    end subroutine l_mvect_zero
  end interface

  interface 
    module subroutine l_mvect_asb(m,n, x, info)
      integer(psb_ipk_), intent(in)              :: m,n
      class(psb_l_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine l_mvect_asb
  end interface

  interface 
    module subroutine l_mvect_sync(x)
      class(psb_l_multivect_type), intent(inout) :: x
    end subroutine l_mvect_sync
  end interface

  interface 
    module subroutine l_mvect_gthab(n,idx,alpha,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_lpk_) :: alpha, beta, y(:)
      class(psb_l_multivect_type) :: x
    end subroutine l_mvect_gthab
  end interface

  interface 
    module subroutine l_mvect_gthzv(n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_lpk_) ::  y(:)
      class(psb_l_multivect_type) :: x
    end subroutine l_mvect_gthzv
  end interface

  interface 
    module subroutine l_mvect_gthzv_x(i,n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      integer(psb_lpk_) ::  y(:)
      class(psb_l_multivect_type) :: x
    end subroutine l_mvect_gthzv_x
  end interface

  interface 
    module subroutine l_mvect_sctb(n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_lpk_) :: beta, x(:)
      class(psb_l_multivect_type) :: y
    end subroutine l_mvect_sctb
  end interface

  interface 
    module subroutine l_mvect_sctb_x(i,n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      integer(psb_lpk_) :: beta, x(:)
      class(psb_l_multivect_type) :: y
    end subroutine l_mvect_sctb_x
  end interface

  interface 
    module subroutine l_mvect_free(x, info)
      class(psb_l_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine l_mvect_free
  end interface

  interface 
    module subroutine l_mvect_ins(n,irl,val,x,maxr,info)
      class(psb_l_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n,maxr
      integer(psb_ipk_), intent(in)               :: irl(:)
      integer(psb_lpk_), intent(in)        :: val(:,:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine l_mvect_ins
  end interface

  interface 
    module subroutine l_mvect_cnv(x,mold)
      class(psb_l_multivect_type), intent(inout) :: x
      class(psb_l_base_multivect_type), intent(in), optional :: mold
    end subroutine l_mvect_cnv
  end interface


contains
  
 function constructor(x) result(this)
    integer(psb_lpk_)   :: x(:,:)
    type(psb_l_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld_x(x)
    call this%asb(size(x,dim=1,kind=psb_ipk_),size(x,dim=2,kind=psb_ipk_),info)

  end function constructor

  function size_const(m,n) result(this)
    integer(psb_ipk_), intent(in) :: m,n
    type(psb_l_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld_n(m,n)
    call this%asb(m,n,info)

  end function size_const

end module psb_l_multivect_mod
