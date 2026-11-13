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
! package: psb_i_vect_mod
!
! This module contains the definition of the psb_i_vect type which
! is the outer container for dense vectors.
! Therefore all methods simply invoke the corresponding methods of the
! inner component.
!
module psb_i_vect_mod

  use psb_realloc_mod
  use psb_i_base_vect_mod

  type psb_i_vect_type
    class(psb_i_base_vect_type), allocatable :: v
    integer(psb_ipk_) :: nrmv = 0
    integer(psb_ipk_) :: remote_build=psb_matbld_noremote_
    integer(psb_ipk_) :: dupl = psb_dupl_add_
    integer(psb_ipk_), allocatable :: rmtv(:)
    integer(psb_lpk_), allocatable :: rmidx(:)
  contains
    procedure, pass(x) :: get_nrows => i_vect_get_nrows
    procedure, pass(x) :: sizeof   => i_vect_sizeof
    procedure, pass(x) :: get_fmt  => i_vect_get_fmt
    procedure, pass(x) :: is_remote_build => i_vect_is_remote_build
    procedure, pass(x) :: set_remote_build => i_vect_set_remote_build
    procedure, pass(x) :: get_nrmv => i_vect_get_nrmv
    procedure, pass(x) :: set_nrmv => i_vect_set_nrmv
    procedure, pass(x) :: all      => i_vect_all
    procedure, pass(x) :: reall    => i_vect_reall
    procedure, pass(x) :: zero     => i_vect_zero
    procedure, pass(x) :: asb      => i_vect_asb
    procedure, pass(x) :: set_dupl => i_vect_set_dupl 
    procedure, pass(x) :: get_dupl => i_vect_get_dupl
    procedure, pass(x) :: set_ncfs => i_vect_set_ncfs 
    procedure, pass(x) :: get_ncfs => i_vect_get_ncfs
    procedure, pass(x) :: set_state => i_vect_set_state
    procedure, pass(x) :: set_null  => i_vect_set_null
    procedure, pass(x) :: set_bld   => i_vect_set_bld
    procedure, pass(x) :: set_upd   => i_vect_set_upd
    procedure, pass(x) :: set_asb   => i_vect_set_asb
    procedure, pass(x) :: get_state => i_vect_get_state
    procedure, pass(x) :: is_null   => i_vect_is_null
    procedure, pass(x) :: is_bld    => i_vect_is_bld
    procedure, pass(x) :: is_upd    => i_vect_is_upd
    procedure, pass(x) :: is_asb    => i_vect_is_asb
    procedure, pass(x) :: reinit    => i_vect_reinit

    procedure, pass(x) :: gthab    => i_vect_gthab
    procedure, pass(x) :: gthzv    => i_vect_gthzv
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => i_vect_sctb
    generic, public    :: sct      => sctb
    procedure, pass(x) :: free     => i_vect_free
    procedure, pass(x) :: ins_a    => i_vect_ins_a
    procedure, pass(x) :: ins_v    => i_vect_ins_v
    generic, public    :: ins      => ins_v, ins_a
    procedure, pass(x) :: bld_x    => i_vect_bld_x
    procedure, pass(x) :: bld_mn   => i_vect_bld_mn
    procedure, pass(x) :: bld_en   => i_vect_bld_en
    generic, public    :: bld      => bld_x, bld_mn, bld_en
    procedure, pass(x) :: get_vect => i_vect_get_vect
    procedure, pass(x) :: cnv      => i_vect_cnv
    procedure, pass(x) :: set_scal => i_vect_set_scal
    procedure, pass(x) :: set_vect => i_vect_set_vect
    generic, public    :: set      => set_vect, set_scal
    procedure, pass(x) :: clone    => i_vect_clone

    procedure, pass(x) :: sync     => i_vect_sync
    procedure, pass(x) :: is_host  => i_vect_is_host
    procedure, pass(x) :: is_dev   => i_vect_is_dev
    procedure, pass(x) :: is_sync  => i_vect_is_sync
    procedure, pass(x) :: set_host => i_vect_set_host
    procedure, pass(x) :: set_dev  => i_vect_set_dev
    procedure, pass(x) :: set_sync => i_vect_set_sync
    procedure, pass(x) :: check_addr => i_vect_check_addr



  end type psb_i_vect_type

  public  :: psb_i_vect
  private :: constructor, size_const
  interface psb_i_vect
    module procedure constructor, size_const
  end interface psb_i_vect

  private :: i_vect_get_nrows, i_vect_sizeof, i_vect_get_fmt, &
       & i_vect_all, i_vect_reall, i_vect_zero,  i_vect_asb, &
       & i_vect_gthab, i_vect_gthzv, i_vect_sctb, &
       & i_vect_free, i_vect_ins_a, i_vect_ins_v, i_vect_bld_x, &
       & i_vect_bld_mn, i_vect_bld_en, i_vect_get_vect, &
       & i_vect_cnv, i_vect_set_scal, &
       & i_vect_set_vect, i_vect_clone, i_vect_sync, i_vect_is_host, &
       & i_vect_is_dev, i_vect_is_sync, i_vect_set_host, &
       & i_vect_set_dev, i_vect_set_sync, &
       & i_vect_set_remote_build, i_is_remote_build, &
       & i_vect_set_dupl, i_get_dupl, i_vect_set_nrmv, i_get_nrmv


  class(psb_i_base_vect_type), allocatable, target,&
       & save, private :: psb_i_base_vect_default

  
  interface 
    module function i_vect_get_dupl(x) result(res)
      implicit none
      class(psb_i_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_vect_get_dupl
  end interface

  interface 
    module subroutine i_vect_set_dupl(x,val)
      class(psb_i_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine i_vect_set_dupl
  end interface

  interface 
    module function i_vect_get_ncfs(x) result(res)
      class(psb_i_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_vect_get_ncfs
  end interface

  interface 
    module subroutine i_vect_set_ncfs(x,val)
      class(psb_i_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine i_vect_set_ncfs
  end interface

  interface 
    module function i_vect_get_state(x) result(res)
      class(psb_i_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_vect_get_state
  end interface

  interface 
    module function i_vect_is_null(x) result(res)
      class(psb_i_vect_type), intent(in) :: x
      logical :: res
    end function i_vect_is_null
  end interface

  interface 
    module function i_vect_is_bld(x) result(res)
      class(psb_i_vect_type), intent(in) :: x
      logical :: res
    end function i_vect_is_bld
  end interface

  interface 
    module function i_vect_is_upd(x) result(res)
      class(psb_i_vect_type), intent(in) :: x
      logical :: res
    end function i_vect_is_upd
  end interface

  interface 
    module function i_vect_is_asb(x) result(res)
      class(psb_i_vect_type), intent(in) :: x
      logical :: res
    end function i_vect_is_asb
  end interface

  interface 
    module subroutine  i_vect_set_state(n,x)
      class(psb_i_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine i_vect_set_state
  end interface

  interface 
    module subroutine  i_vect_set_null(x)
      class(psb_i_vect_type), intent(inout) :: x
    end subroutine i_vect_set_null
  end interface

  interface 
    module subroutine  i_vect_set_bld(x)
      class(psb_i_vect_type), intent(inout) :: x
    end subroutine i_vect_set_bld
  end interface

  interface 
    module subroutine  i_vect_set_upd(x)
      class(psb_i_vect_type), intent(inout) :: x
    end subroutine i_vect_set_upd
  end interface

  interface 
    module subroutine  i_vect_set_asb(x)
      class(psb_i_vect_type), intent(inout) :: x
    end subroutine i_vect_set_asb
  end interface

  interface 
    module function i_vect_get_nrmv(x) result(res)
      class(psb_i_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_vect_get_nrmv
  end interface

  interface 
    module subroutine i_vect_set_nrmv(x,val)
      class(psb_i_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: val
    end subroutine i_vect_set_nrmv
  end interface

  interface 
    module function i_vect_is_remote_build(x) result(res)
      class(psb_i_vect_type), intent(in) :: x
      logical :: res
    end function i_vect_is_remote_build
  end interface
  
  interface 
    module subroutine i_vect_set_remote_build(x,val)
      class(psb_i_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine i_vect_set_remote_build
  end interface
        
  interface psb_set_vect_default
    module subroutine  psb_i_set_vect_default(v)
      class(psb_i_base_vect_type), intent(in) :: v
    end subroutine psb_i_set_vect_default
  end interface

  interface psb_get_vect_default
    module function psb_i_get_vect_default(v) result(res)
      class(psb_i_vect_type), intent(in) :: v
      class(psb_i_base_vect_type), pointer :: res
    end function psb_i_get_vect_default
  end interface

  interface 
    module subroutine  psb_i_clear_vect_default()
    end subroutine psb_i_clear_vect_default
  end interface

  interface 
    module function psb_i_get_base_vect_default() result(res)
      class(psb_i_base_vect_type), pointer :: res
    end function psb_i_get_base_vect_default
  end interface

  interface 
    module subroutine i_vect_clone(x,y,info)
      class(psb_i_vect_type), intent(inout) :: x
      class(psb_i_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)        :: info
    end subroutine i_vect_clone
  end interface

  interface 
    module subroutine i_vect_bld_x(x,invect,mold,scratch)
      integer(psb_ipk_), intent(in)          :: invect(:)
      class(psb_i_vect_type), intent(inout) :: x
      class(psb_i_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine i_vect_bld_x
  end interface

  interface 
    module subroutine i_vect_bld_mn(x,n,mold,scratch)
      integer(psb_mpk_), intent(in) :: n
      class(psb_i_vect_type), intent(inout) :: x
      class(psb_i_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine i_vect_bld_mn
  end interface

  interface 
    module subroutine i_vect_bld_en(x,n,mold,scratch)
      integer(psb_epk_), intent(in) :: n
      class(psb_i_vect_type), intent(inout) :: x
      class(psb_i_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine i_vect_bld_en
  end interface

  interface 
    module function  i_vect_get_vect(x,n) result(res)
      class(psb_i_vect_type), intent(inout)  :: x
      integer(psb_ipk_), allocatable                 :: res(:)
      integer(psb_ipk_), optional :: n
    end function i_vect_get_vect
  end interface

  interface 
    module subroutine i_vect_set_scal(x,val,first,last)
      class(psb_i_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in) :: val
      integer(psb_ipk_), optional :: first, last
    end subroutine i_vect_set_scal
  end interface

  interface 
    module subroutine i_vect_set_vect(x,val,first,last)
      class(psb_i_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)         :: val(:)
      integer(psb_ipk_), optional :: first, last
    end subroutine i_vect_set_vect
  end interface

  interface 
    module subroutine i_vect_check_addr(x)
      class(psb_i_vect_type), intent(inout) :: x
    end subroutine i_vect_check_addr
  end interface

  interface 
    module function i_vect_get_nrows(x) result(res)
      class(psb_i_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_vect_get_nrows
  end interface

  interface 
    module function i_vect_sizeof(x) result(res)
      class(psb_i_vect_type), intent(in) :: x
      integer(psb_epk_) :: res
    end function i_vect_sizeof
  end interface

  interface 
    module function i_vect_get_fmt(x) result(res)
      class(psb_i_vect_type), intent(in) :: x
      character(len=5) :: res
    end function i_vect_get_fmt
  end interface

  interface 
    module subroutine i_vect_all(n, x, info, mold)
      integer(psb_ipk_), intent(in)           :: n
      class(psb_i_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)      :: info
      class(psb_i_base_vect_type), intent(in), optional :: mold
    end subroutine i_vect_all
  end interface

  interface 
    module subroutine i_vect_reinit(x, info, clear)
      class(psb_i_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: clear
    end subroutine i_vect_reinit
  end interface
 
  interface 
    module subroutine i_vect_reall(n, x, info)
      integer(psb_ipk_), intent(in)         :: n
      class(psb_i_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)        :: info
    end subroutine i_vect_reall
  end interface

  interface 
    module subroutine i_vect_zero(x)
      class(psb_i_vect_type), intent(inout)    :: x
    end subroutine i_vect_zero
  end interface

  interface 
    module subroutine i_vect_asb(n, x, info, scratch)
      integer(psb_ipk_), intent(in)         :: n
      class(psb_i_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional        :: scratch
    end subroutine i_vect_asb
  end interface

  interface 
    module subroutine i_vect_gthab(n,idx,alpha,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_ipk_) :: alpha, beta, y(:)
      class(psb_i_vect_type) :: x
    end subroutine i_vect_gthab
  end interface

  interface 
    module subroutine i_vect_gthzv(n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_ipk_) ::  y(:)
      class(psb_i_vect_type) :: x
    end subroutine i_vect_gthzv
  end interface

  interface 
    module subroutine i_vect_sctb(n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_ipk_) :: beta, x(:)
      class(psb_i_vect_type) :: y
    end subroutine i_vect_sctb
  end interface

  interface 
    module subroutine i_vect_free(x, info)
      class(psb_i_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine i_vect_free
  end interface
  
  interface 
    module subroutine i_vect_ins_a(n,irl,val,x,maxr,info)
      class(psb_i_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n, maxr
      integer(psb_ipk_), intent(in)               :: irl(:)
      integer(psb_ipk_), intent(in)        :: val(:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine i_vect_ins_a
  end interface

  interface 
    module subroutine i_vect_ins_v(n,irl,val,x,maxr,info)
      class(psb_i_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n, maxr
      class(psb_i_vect_type), intent(inout)       :: irl
      class(psb_i_vect_type), intent(inout)       :: val
      integer(psb_ipk_), intent(out)              :: info
    end subroutine i_vect_ins_v
  end interface

  interface 
    module subroutine i_vect_cnv(x,mold)
      class(psb_i_vect_type), intent(inout) :: x
      class(psb_i_base_vect_type), intent(in), optional :: mold
      class(psb_i_base_vect_type), allocatable :: tmp
    end subroutine i_vect_cnv
  end interface
  
  interface 
    module subroutine i_vect_sync(x)
      class(psb_i_vect_type), intent(inout) :: x
    end subroutine i_vect_sync
  end interface

  interface 
    module subroutine i_vect_set_sync(x)
      class(psb_i_vect_type), intent(inout) :: x
    end subroutine i_vect_set_sync
  end interface

  interface 
    module subroutine i_vect_set_host(x)
      class(psb_i_vect_type), intent(inout) :: x
    end subroutine i_vect_set_host
  end interface

  interface 
    module subroutine i_vect_set_dev(x)
      class(psb_i_vect_type), intent(inout) :: x
    end subroutine i_vect_set_dev
  end interface

  interface 
    module function i_vect_is_sync(x) result(res)
      logical :: res
      class(psb_i_vect_type), intent(inout) :: x
    end function i_vect_is_sync
  end interface

  interface 
    module function i_vect_is_host(x) result(res)
      logical :: res
      class(psb_i_vect_type), intent(inout) :: x
    end function i_vect_is_host
  end interface

  interface 
    module function i_vect_is_dev(x) result(res)
      logical :: res
      class(psb_i_vect_type), intent(inout) :: x
    end function i_vect_is_dev
  end interface




contains

  function constructor(x) result(this)
    integer(psb_ipk_)   :: x(:)
    type(psb_i_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(x)
    call this%asb(size(x,kind=psb_ipk_),info)

  end function constructor


  function size_const(n) result(this)
    integer(psb_ipk_), intent(in) :: n
    type(psb_i_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(n)
    call this%asb(n,info)

  end function size_const

!!$  function i_vect_get_nrows(x) result(res)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(in) :: x
!!$    integer(psb_ipk_) :: res
!!$    res = 0
!!$    if (allocated(x%v)) res = x%v%get_nrows()
!!$  end function i_vect_get_nrows
!!$
!!$  function i_vect_sizeof(x) result(res)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(in) :: x
!!$    integer(psb_epk_) :: res
!!$    res = 0
!!$    if (allocated(x%v)) res = x%v%sizeof()
!!$  end function i_vect_sizeof
!!$
!!$  function i_vect_get_fmt(x) result(res)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(in) :: x
!!$    character(len=5) :: res
!!$    res = 'NULL'
!!$    if (allocated(x%v)) res = x%v%get_fmt()
!!$  end function i_vect_get_fmt
!!$
!!$  subroutine i_vect_all(n, x, info, mold)
!!$
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(out)      :: info
!!$    class(psb_i_base_vect_type), intent(in), optional :: mold
!!$
!!$    if (allocated(x%v)) &
!!$         & call x%free(info)
!!$
!!$    if (present(mold)) then
!!$      allocate(x%v,stat=info,mold=mold)
!!$    else
!!$      allocate(psb_i_base_vect_type :: x%v,stat=info)
!!$    endif
!!$    if (info == 0) then
!!$      call x%v%all(n,info)
!!$    else
!!$      info = psb_err_alloc_dealloc_
!!$    end if
!!$    call x%set_bld()
!!$  end subroutine i_vect_all
!!$
!!$  subroutine i_vect_reinit(x, info, clear)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(out)      :: info
!!$    logical, intent(in), optional       :: clear
!!$
!!$    if (allocated(x%v)) call x%v%reinit(info,clear)
!!$    call x%set_upd()
!!$
!!$  end subroutine i_vect_reinit
!!$ 
!!$  subroutine i_vect_reall(n, x, info)
!!$
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)         :: n
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(out)        :: info
!!$
!!$    info = 0
!!$    if (.not.allocated(x%v)) &
!!$         & call x%all(n,info)
!!$    if (info == 0) &
!!$         & call x%asb(n,info)
!!$
!!$  end subroutine i_vect_reall
!!$
!!$  subroutine i_vect_zero(x)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout)    :: x
!!$
!!$    if (allocated(x%v)) call x%v%zero()
!!$
!!$  end subroutine i_vect_zero
!!$
!!$  subroutine i_vect_asb(n, x, info, scratch)
!!$    use psi_serial_mod
!!$    use psb_realloc_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)         :: n
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(out)        :: info
!!$    logical, intent(in), optional        :: scratch
!!$
!!$    if (allocated(x%v)) then
!!$      call x%v%asb(n,info,scratch=scratch)
!!$      call x%set_asb()
!!$    end if
!!$  end subroutine i_vect_asb
!!$
!!$  subroutine i_vect_gthab(n,idx,alpha,x,beta,y)
!!$    use psi_serial_mod
!!$    integer(psb_mpk_) :: n
!!$    integer(psb_ipk_) :: idx(:)
!!$    integer(psb_ipk_) :: alpha, beta, y(:)
!!$    class(psb_i_vect_type) :: x
!!$
!!$    if (allocated(x%v)) &
!!$         &  call x%v%gth(n,idx,alpha,beta,y)
!!$
!!$  end subroutine i_vect_gthab
!!$
!!$  subroutine i_vect_gthzv(n,idx,x,y)
!!$    use psi_serial_mod
!!$    integer(psb_mpk_) :: n
!!$    integer(psb_ipk_) :: idx(:)
!!$    integer(psb_ipk_) ::  y(:)
!!$    class(psb_i_vect_type) :: x
!!$
!!$    if (allocated(x%v)) &
!!$         &  call x%v%gth(n,idx,y)
!!$
!!$  end subroutine i_vect_gthzv
!!$
!!$  subroutine i_vect_sctb(n,idx,x,beta,y)
!!$    use psi_serial_mod
!!$    integer(psb_mpk_) :: n
!!$    integer(psb_ipk_) :: idx(:)
!!$    integer(psb_ipk_) :: beta, x(:)
!!$    class(psb_i_vect_type) :: y
!!$
!!$    if (allocated(y%v)) &
!!$         &  call y%v%sct(n,idx,x,beta)
!!$
!!$  end subroutine i_vect_sctb
!!$
!!$  subroutine i_vect_free(x, info)
!!$    use psi_serial_mod
!!$    use psb_realloc_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    info = 0
!!$    if (allocated(x%v)) then
!!$      call x%v%free(info)
!!$      if (info == 0) deallocate(x%v,stat=info)
!!$    end if
!!$
!!$  end subroutine i_vect_free
!!$
!!$  subroutine i_vect_ins_a(n,irl,val,x,maxr,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    integer(psb_ipk_), intent(in)               :: n, maxr
!!$    integer(psb_ipk_), intent(in)               :: irl(:)
!!$    integer(psb_ipk_), intent(in)        :: val(:)
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
!!$  end subroutine i_vect_ins_a
!!$
!!$  subroutine i_vect_ins_v(n,irl,val,x,maxr,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    integer(psb_ipk_), intent(in)               :: n, maxr
!!$    class(psb_i_vect_type), intent(inout)       :: irl
!!$    class(psb_i_vect_type), intent(inout)       :: val
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
!!$  end subroutine i_vect_ins_v
!!$
!!$
!!$  subroutine i_vect_cnv(x,mold)
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$    class(psb_i_base_vect_type), intent(in), optional :: mold
!!$    class(psb_i_base_vect_type), allocatable :: tmp
!!$
!!$    integer(psb_ipk_) :: info
!!$
!!$    info = psb_success_
!!$    if (present(mold)) then
!!$      allocate(tmp,stat=info,mold=mold)
!!$    else
!!$      allocate(tmp,stat=info,mold=psb_i_get_base_vect_default())
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
!!$  end subroutine i_vect_cnv
!!$
!!$
!!$  subroutine i_vect_sync(x)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$
!!$    if (allocated(x%v)) &
!!$         & call x%v%sync()
!!$
!!$  end subroutine i_vect_sync
!!$
!!$  subroutine i_vect_set_sync(x)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$
!!$    if (allocated(x%v)) &
!!$         & call x%v%set_sync()
!!$
!!$  end subroutine i_vect_set_sync
!!$
!!$  subroutine i_vect_set_host(x)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$
!!$    if (allocated(x%v)) &
!!$         & call x%v%set_host()
!!$
!!$  end subroutine i_vect_set_host
!!$
!!$  subroutine i_vect_set_dev(x)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$
!!$    if (allocated(x%v)) &
!!$         & call x%v%set_dev()
!!$
!!$  end subroutine i_vect_set_dev
!!$
!!$  function i_vect_is_sync(x) result(res)
!!$    implicit none
!!$    logical :: res
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$
!!$    res = .true.
!!$    if (allocated(x%v)) &
!!$         & res = x%v%is_sync()
!!$
!!$  end function i_vect_is_sync
!!$
!!$  function i_vect_is_host(x) result(res)
!!$    implicit none
!!$    logical :: res
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$
!!$    res = .true.
!!$    if (allocated(x%v)) &
!!$         & res = x%v%is_host()
!!$
!!$  end function i_vect_is_host
!!$
!!$  function i_vect_is_dev(x) result(res)
!!$    implicit none
!!$    logical :: res
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$
!!$    res = .false.
!!$    if (allocated(x%v)) &
!!$         & res =  x%v%is_dev()
!!$
!!$  end function i_vect_is_dev
!!$
!!$  @NOTINTS@
!!$
!!$  function i_vect_get_entry(x,index) result(res)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)        :: index
!!$    integer(psb_ipk_) :: res
!!$    res = izero
!!$    if (allocated(x%v)) res = x%v%get_entry(index)
!!$  end function i_vect_get_entry
!!$
!!$  subroutine i_vect_set_entry(x,index,val)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)        :: index
!!$    integer(psb_ipk_) :: val
!!$    if (allocated(x%v)) call x%v%set_entry(index,val)
!!$  end subroutine i_vect_set_entry
!!$
!!$  function i_vect_dot_v(n,x,y) result(res)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout) :: x, y
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_ipk_)                :: res
!!$
!!$    res = izero
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & res = x%v%dot(n,y%v)
!!$
!!$  end function i_vect_dot_v
!!$
!!$  function i_vect_dot_a(n,x,y) result(res)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)    :: y(:)
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_ipk_)                :: res
!!$
!!$    res = izero
!!$    if (allocated(x%v)) &
!!$         & res = x%v%dot_a(n,y)
!!$
!!$  end function i_vect_dot_a
!!$
!!$  subroutine i_vect_axpby_v(m,alpha, x, beta, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    if (allocated(x%v).and.allocated(y%v)) then
!!$      call y%v%axpby(m,alpha,x%v,beta,info)
!!$    else
!!$      info = psb_err_invalid_vect_state_
!!$    end if
!!$
!!$  end subroutine i_vect_axpby_v
!!$
!!$  subroutine i_vect_axpby_v2(m,alpha, x, beta, y, z, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)            :: m
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    class(psb_i_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent (in)             :: alpha, beta
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$    if (allocated(x%v).and.allocated(y%v)) then
!!$      call z%v%axpby(m,alpha,x%v,beta,y%v,info)
!!$    else
!!$      info = psb_err_invalid_vect_state_
!!$    end if
!!$
!!$  end subroutine i_vect_axpby_v2
!!$
!!$  subroutine i_vect_axpby_a(m,alpha, x, beta, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    integer(psb_ipk_), intent(in)        :: x(:)
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    if (allocated(y%v)) &
!!$         & call y%v%axpby(m,alpha,x,beta,info)
!!$
!!$  end subroutine i_vect_axpby_a
!!$
!!$  subroutine i_vect_axpby_a2(m,alpha, x, beta, y, z, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    integer(psb_ipk_), intent(in)                 :: x(:)
!!$    integer(psb_ipk_), intent(in)                 :: y(:)
!!$    class(psb_i_vect_type), intent(inout)     :: z
!!$    integer(psb_ipk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    if (allocated(z%v)) &
!!$         & call z%v%axpby(m,alpha,x,beta,y,info)
!!$
!!$  end subroutine i_vect_axpby_a2
!!$
!!$  subroutine i_vect_upd_xyz(m,alpha,beta,gamma,delta,x, y, z, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)            :: m
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    class(psb_i_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent (in)     :: alpha, beta, gamma, delta
!!$    integer(psb_ipk_), intent(out)   :: info
!!$
!!$    if (allocated(z%v)) &
!!$         call z%v%upd_xyz(m,alpha,beta,gamma,delta,x%v,y%v,info)
!!$    
!!$  end subroutine i_vect_upd_xyz
!!$
!!$  subroutine i_vect_xyzw(m,a,b,c,d,e,f,x, y, z, w, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)            :: m
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    class(psb_i_vect_type), intent(inout)  :: z
!!$    class(psb_i_vect_type), intent(inout)  :: w
!!$    integer(psb_ipk_), intent (in)     :: a, b, c, d, e, f
!!$    integer(psb_ipk_), intent(out)   :: info
!!$
!!$    if (allocated(w%v)) &
!!$         call w%v%xyzw(m,a,b,c,d,e,f,x%v,y%v,z%v,info)
!!$    
!!$  end subroutine i_vect_xyzw
!!$
!!$
!!$  subroutine i_vect_mlt_v(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & call y%v%mlt(x%v,info)
!!$
!!$  end subroutine i_vect_mlt_v
!!$
!!$  subroutine i_vect_mlt_a(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)        :: x(:)
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$
!!$    info = 0
!!$    if (allocated(y%v)) &
!!$         & call y%v%mlt(x,info)
!!$
!!$  end subroutine i_vect_mlt_a
!!$
!!$
!!$  subroutine i_vect_mlt_a_2(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)         :: alpha,beta
!!$    integer(psb_ipk_), intent(in)         :: y(:)
!!$    integer(psb_ipk_), intent(in)         :: x(:)
!!$    class(psb_i_vect_type), intent(inout) :: z
!!$    integer(psb_ipk_), intent(out)                  :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(z%v)) &
!!$         & call z%v%mlt(alpha,x,y,beta,info)
!!$
!!$  end subroutine i_vect_mlt_a_2
!!$
!!$  subroutine i_vect_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)          :: alpha,beta
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    class(psb_i_vect_type), intent(inout)  :: z
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
!!$  end subroutine i_vect_mlt_v_2
!!$
!!$  subroutine i_vect_mlt_av(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)        :: alpha,beta
!!$    integer(psb_ipk_), intent(in)        :: x(:)
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    class(psb_i_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(z%v).and.allocated(y%v)) &
!!$         & call z%v%mlt(alpha,x,y%v,beta,info)
!!$
!!$  end subroutine i_vect_mlt_av
!!$
!!$  subroutine i_vect_mlt_va(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)        :: alpha,beta
!!$    integer(psb_ipk_), intent(in)        :: y(:)
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$
!!$    if (allocated(z%v).and.allocated(x%v)) &
!!$         & call z%v%mlt(alpha,x%v,y,beta,info)
!!$
!!$  end subroutine i_vect_mlt_va
!!$
!!$  subroutine i_vect_div_v(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & call x%v%div(y%v,info)
!!$
!!$  end subroutine i_vect_div_v
!!$
!!$  subroutine i_vect_div_v2( x, y, z, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    class(psb_i_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v).and.allocated(z%v)) &
!!$         & call z%v%div(x%v,y%v,info)
!!$
!!$  end subroutine i_vect_div_v2
!!$
!!$  subroutine i_vect_div_v_check(x, y, info, flag)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$    logical, intent(in) :: flag
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & call x%v%div(y%v,info,flag)
!!$
!!$  end subroutine i_vect_div_v_check
!!$
!!$  subroutine i_vect_div_v2_check(x, y, z, info, flag)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    class(psb_i_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$    logical, intent(in) :: flag
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v).and.allocated(z%v)) &
!!$         & call z%v%div(x%v,y%v,info,flag)
!!$
!!$  end subroutine i_vect_div_v2_check
!!$
!!$  subroutine i_vect_div_a2(x, y, z, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in) :: x(:)
!!$    integer(psb_ipk_), intent(in)    :: y(:)
!!$    class(psb_i_vect_type), intent(inout) :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(z%v)) &
!!$         & call z%v%div(x,y,info)
!!$
!!$  end subroutine i_vect_div_a2
!!$
!!$  subroutine i_vect_div_a2_check(x, y, z, info,flag)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in) :: x(:)
!!$    integer(psb_ipk_), intent(in)    :: y(:)
!!$    class(psb_i_vect_type), intent(inout) :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$    logical, intent(in) :: flag
!!$
!!$    info = 0
!!$    if (allocated(z%v)) &
!!$         & call z%v%div(x,y,info,flag)
!!$
!!$  end subroutine i_vect_div_a2_check
!!$
!!$  subroutine i_vect_inv_v(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & call y%v%inv(x%v,info)
!!$
!!$  end subroutine i_vect_inv_v
!!$
!!$  subroutine i_vect_inv_v_check(x, y, info, flag)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$    logical, intent(in) :: flag
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & call y%v%inv(x%v,info,flag)
!!$
!!$  end subroutine i_vect_inv_v_check
!!$
!!$  subroutine i_vect_inv_a2(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(inout) :: x(:)
!!$    class(psb_i_vect_type), intent(inout) :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(y%v)) &
!!$         & call y%v%inv(x,info)
!!$
!!$  end subroutine i_vect_inv_a2
!!$
!!$  subroutine i_vect_inv_a2_check(x, y, info,flag)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(inout) :: x(:)
!!$    class(psb_i_vect_type), intent(inout) :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$    logical, intent(in) :: flag
!!$
!!$    info = 0
!!$    if (allocated(y%v)) &
!!$         & call y%v%inv(x,info,flag)
!!$
!!$  end subroutine i_vect_inv_a2_check
!!$
!!$  subroutine i_vect_acmp_a2(x,c,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)              :: c
!!$    integer(psb_ipk_), intent(inout)           :: x(:)
!!$    class(psb_i_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$    info = 0
!!$    if (allocated(z%v)) &
!!$         & call z%acmp(x,c,info)
!!$
!!$  end subroutine i_vect_acmp_a2
!!$
!!$  subroutine i_vect_acmp_v2(x,c,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)              :: c
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(z%v)) &
!!$         & call z%v%acmp(x%v,c,info)
!!$
!!$  end subroutine i_vect_acmp_v2
!!$
!!$  subroutine i_vect_scal(alpha, x)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    integer(psb_ipk_), intent (in)       :: alpha
!!$
!!$    if (allocated(x%v)) call x%v%scal(alpha)
!!$
!!$  end subroutine i_vect_scal
!!$
!!$  subroutine i_vect_absval1(x)
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$
!!$    if (allocated(x%v)) &
!!$         &  call x%v%absval()
!!$
!!$  end subroutine i_vect_absval1
!!$
!!$  subroutine i_vect_absval2(x,y)
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$
!!$    if (allocated(x%v)) then
!!$      if (.not.allocated(y%v))  call y%bld(psb_size(x%v%v))
!!$      call x%v%absval(y%v)
!!$    end if
!!$  end subroutine i_vect_absval2
!!$
!!$  function i_vect_nrm2(n,x) result(res)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_ipk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%nrm2(n)
!!$    else
!!$      res = izero
!!$    end if
!!$
!!$  end function i_vect_nrm2
!!$
!!$  function i_vect_nrm2_weight(n,x,w,aux) result(res)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$    class(psb_i_vect_type), intent(inout) :: w
!!$    class(psb_i_vect_type), intent(inout), optional :: aux
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_ipk_)                        :: res
!!$    integer(psb_ipk_)                       :: info
!!$
!!$    ! Temp vectors
!!$    type(psb_i_vect_type) :: wtemp
!!$
!!$    info = 0
!!$    if( allocated(w%v) ) then
!!$      if (.not.present(aux)) then
!!$        allocate(wtemp%v, mold=w%v)
!!$        call wtemp%v%bld(w%get_vect())
!!$      else
!!$        call psb_geaxpby(n,ione,w%v%v,izero,aux%v%v,info)
!!$      end if
!!$    else
!!$      info = -1
!!$    end if
!!$    if (info /= 0 ) then
!!$      res = -ione
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
!!$      res = izero
!!$    end if
!!$
!!$    if (.not.present(aux)) then
!!$      call wtemp%free(info)
!!$    end if
!!$
!!$  end function i_vect_nrm2_weight
!!$ 
!!$  function i_vect_nrm2_weight_mask(n,x,w,id,info,aux) result(res)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$    class(psb_i_vect_type), intent(inout) :: w
!!$    class(psb_i_vect_type), intent(inout) :: id
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_ipk_)                        :: res
!!$    integer(psb_ipk_), intent(out)          :: info
!!$    class(psb_i_vect_type), intent(inout), optional :: aux
!!$
!!$    ! Temp vectors
!!$    type(psb_i_vect_type) :: wtemp
!!$
!!$    info = 0 
!!$    if( allocated(w%v) ) then
!!$      if (.not.present(aux)) then
!!$        allocate(wtemp%v, mold=w%v)
!!$        call wtemp%v%bld(w%get_vect())
!!$      else
!!$        call psb_geaxpby(n,ione,w%v%v,izero,aux%v%v,info)
!!$      end if
!!$    else
!!$      info = -1
!!$    end if
!!$    if (info /= 0 ) then
!!$      res = -ione
!!$      return
!!$    end if
!!$
!!$
!!$    if (allocated(x%v).and.allocated(id%v)) then
!!$      if (.not.present(aux)) then
!!$        where( abs(id%v%v) <= izero) wtemp%v%v = izero
!!$        call wtemp%set_host()
!!$        call wtemp%v%mlt(x%v,info)
!!$        res = wtemp%v%nrm2(n)
!!$      else
!!$        where( abs(id%v%v) <= izero) aux%v%v = izero
!!$        call aux%set_host()
!!$        call aux%v%mlt(x%v,info)
!!$        res = aux%v%nrm2(n)
!!$      end if
!!$    else
!!$      res = izero
!!$    end if
!!$
!!$    if (.not.present(aux)) then
!!$      call wtemp%free(info)
!!$    end if
!!$
!!$  end function i_vect_nrm2_weight_mask
!!$
!!$  function i_vect_amax(n,x) result(res)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_ipk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%amax(n)
!!$    else
!!$      res = izero
!!$    end if
!!$
!!$  end function i_vect_amax
!!$
!!$  @NOTCPLXS@
!!$  function i_vect_min(n,x) result(res)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_ipk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%minreal(n)
!!$    else
!!$      res = HUGE(ione)
!!$    end if
!!$
!!$  end function i_vect_min
!!$  @NOTCPLXE@
!!$
!!$  function i_vect_asum(n,x) result(res)
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    integer(psb_ipk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%asum(n)
!!$    else
!!$      res = izero
!!$    end if
!!$
!!$  end function i_vect_asum
!!$
!!$  @NOTINTE@
!!$
!!$  @NOTINTS@
!!$  @NOTCPLXS@
!!$  subroutine i_vect_mask_a(c,x,m,t,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(inout)          :: c(:)
!!$    integer(psb_ipk_), intent(inout)          :: x(:)
!!$    logical, intent(out)                     :: t;
!!$    class(psb_i_vect_type), intent(inout)  :: m
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$    info = 0
!!$    if (allocated(m%v)) &
!!$         & call m%mask(c,x,t,info)
!!$
!!$  end subroutine i_vect_mask_a
!!$
!!$  subroutine i_vect_mask_v(c,x,m,t,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout)  :: c
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: m
!!$    logical, intent(out)                     :: t;
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(c%v)) &
!!$         & call m%v%mask(x%v,c%v,t,info)
!!$
!!$  end subroutine i_vect_mask_v
!!$
!!$  function i_vect_minquotient_v(x, y, info) result(z)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: y
!!$    integer(psb_ipk_)                         :: z
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & z = x%v%minquotient(y%v,info)
!!$
!!$  end function i_vect_minquotient_v
!!$
!!$  function i_vect_minquotient_a2(x, y, info) result(z)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    integer(psb_ipk_), intent(inout)           :: y(:)
!!$    integer(psb_ipk_), intent(out)           :: info
!!$    integer(psb_ipk_)                         :: z
!!$
!!$    info = 0
!!$    z =  x%v%minquotient(y,info)
!!$
!!$  end function i_vect_minquotient_a2
!!$
!!$
!!$  @NOTINTE@
!!$  @NOTCPLXE@
!!$
!!$  @NOTINTS@
!!$  subroutine i_vect_addconst_a2(x,b,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)             :: b
!!$    integer(psb_ipk_), intent(inout)           :: x(:)
!!$    class(psb_i_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$    info = 0
!!$    if (allocated(z%v)) &
!!$         & call z%addconst(x,b,info)
!!$
!!$  end subroutine i_vect_addconst_a2
!!$
!!$  subroutine i_vect_addconst_v2(x,b,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)             :: b
!!$    class(psb_i_vect_type), intent(inout)  :: x
!!$    class(psb_i_vect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)           :: info
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(z%v)) &
!!$         & call z%v%addconst(x%v,b,info)
!!$
!!$  end subroutine i_vect_addconst_v2
!!$  @NOTINTE@

end module psb_i_vect_mod

module psb_i_multivect_mod

  use psb_i_base_multivect_mod
  use psb_const_mod
  use psb_i_vect_mod


  !private

  type psb_i_multivect_type
    class(psb_i_base_multivect_type), allocatable :: v
    integer(psb_ipk_) :: nrmv = 0
    integer(psb_ipk_) :: remote_build=psb_matbld_noremote_
    integer(psb_ipk_) :: dupl = psb_dupl_add_
    integer(psb_ipk_), allocatable :: rmtv(:,:)
  contains
    procedure, pass(x) :: get_nrows => i_mvect_get_nrows
    procedure, pass(x) :: get_ncols => i_mvect_get_ncols
    procedure, pass(x) :: sizeof   => i_mvect_sizeof
    procedure, pass(x) :: get_fmt  => i_mvect_get_fmt
    procedure, pass(x) :: is_remote_build => i_mvect_is_remote_build
    procedure, pass(x) :: set_remote_build => i_mvect_set_remote_build
    procedure, pass(x) :: get_dupl => i_mvect_get_dupl
    procedure, pass(x) :: set_dupl => i_mvect_set_dupl

    procedure, pass(x) :: all      => i_mvect_all
    procedure, pass(x) :: reall    => i_mvect_reall
    procedure, pass(x) :: zero     => i_mvect_zero
    procedure, pass(x) :: asb      => i_mvect_asb
    procedure, pass(x) :: sync     => i_mvect_sync
    procedure, pass(x) :: free     => i_mvect_free
    procedure, pass(x) :: ins      => i_mvect_ins
    procedure, pass(x) :: bld_x    => i_mvect_bld_x
    procedure, pass(x) :: bld_n    => i_mvect_bld_n
    generic, public    :: bld      => bld_x, bld_n
    procedure, pass(x) :: get_vect => i_mvect_get_vect
    procedure, pass(x) :: cnv      => i_mvect_cnv
    procedure, pass(x) :: set_scal => i_mvect_set_scal
    procedure, pass(x) :: set_vect => i_mvect_set_vect
    generic, public    :: set      => set_vect, set_scal
    procedure, pass(x) :: clone    => i_mvect_clone
    procedure, pass(x) :: gthab    => i_mvect_gthab
    procedure, pass(x) :: gthzv    => i_mvect_gthzv
    procedure, pass(x) :: gthzv_x  => i_mvect_gthzv_x
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => i_mvect_sctb
    procedure, pass(y) :: sctb_x   => i_mvect_sctb_x
    generic, public    :: sct      => sctb, sctb_x
  end type psb_i_multivect_type

  public  :: psb_i_multivect, psb_i_multivect_type,&
       & psb_set_multivect_default, psb_get_multivect_default, &
       & psb_i_base_multivect_type

  private
  interface psb_i_multivect
    module procedure constructor, size_const
  end interface psb_i_multivect

  class(psb_i_base_multivect_type), allocatable, target,&
       & save, private :: psb_i_base_multivect_default

  
  interface 
    module function i_mvect_get_dupl(x) result(res)
      class(psb_i_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_mvect_get_dupl
  end interface

  interface 
    module subroutine i_mvect_set_dupl(x,val)
      class(psb_i_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine i_mvect_set_dupl
  end interface

  interface 
    module function i_mvect_is_remote_build(x) result(res)
      class(psb_i_multivect_type), intent(in) :: x
      logical :: res
    end function i_mvect_is_remote_build
  end interface

  interface 
    module subroutine i_mvect_set_remote_build(x,val)
      class(psb_i_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine i_mvect_set_remote_build
  end interface

  interface 
    module subroutine  psb_i_set_multivect_default(v)
      class(psb_i_base_multivect_type), intent(in) :: v
    end subroutine psb_i_set_multivect_default
  end interface

  interface 
    module function psb_i_get_multivect_default(v) result(res)
      class(psb_i_multivect_type), intent(in) :: v
      class(psb_i_base_multivect_type), pointer :: res
    end function psb_i_get_multivect_default
  end interface

  interface 
    module function psb_i_get_base_multivect_default() result(res)
      class(psb_i_base_multivect_type), pointer :: res
    end function psb_i_get_base_multivect_default
  end interface

  interface 
    module subroutine i_mvect_clone(x,y,info)
      class(psb_i_multivect_type), intent(inout) :: x
      class(psb_i_multivect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)        :: info
    end subroutine i_mvect_clone
  end interface

  interface 
    module subroutine i_mvect_bld_x(x,invect,mold)
      integer(psb_ipk_), intent(in)          :: invect(:,:)
      class(psb_i_multivect_type), intent(out) :: x
      class(psb_i_base_multivect_type), intent(in), optional :: mold
    end subroutine i_mvect_bld_x
  end interface


  interface 
    module subroutine i_mvect_bld_n(x,m,n,mold,scratch)
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_i_multivect_type), intent(out) :: x
      class(psb_i_base_multivect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine i_mvect_bld_n
  end interface

  interface 
    module function  i_mvect_get_vect(x) result(res)
      class(psb_i_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), allocatable                 :: res(:,:)
    end function i_mvect_get_vect
  end interface

  interface 
    module subroutine i_mvect_set_scal(x,val)
      class(psb_i_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in) :: val
    end subroutine i_mvect_set_scal
  end interface

  interface 
    module subroutine i_mvect_set_vect(x,val)
      class(psb_i_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)         :: val(:,:)
    end subroutine i_mvect_set_vect
  end interface

  interface 
    module function i_mvect_get_nrows(x) result(res)
      class(psb_i_multivect_type), intent(in) :: x
      integer(psb_ipk_)  :: res
    end function i_mvect_get_nrows
  end interface

  interface 
    module function i_mvect_get_ncols(x) result(res)
      class(psb_i_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_mvect_get_ncols
  end interface

  interface 
    module function i_mvect_sizeof(x) result(res)
      class(psb_i_multivect_type), intent(in) :: x
      integer(psb_epk_) :: res
    end function i_mvect_sizeof
  end interface

  interface 
    module function i_mvect_get_fmt(x) result(res)
      class(psb_i_multivect_type), intent(in) :: x
      character(len=5) :: res
    end function i_mvect_get_fmt
  end interface

  interface 
    module subroutine i_mvect_all(m,n, x, info, mold)
      integer(psb_ipk_), intent(in)       :: m,n
      class(psb_i_multivect_type), intent(out) :: x
      class(psb_i_base_multivect_type), intent(in), optional :: mold
      integer(psb_ipk_), intent(out)      :: info
    end subroutine i_mvect_all
  end interface

  interface 
    module subroutine i_mvect_reall(m,n, x, info)
      integer(psb_ipk_), intent(in)         :: m,n
      class(psb_i_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)        :: info
    end subroutine i_mvect_reall
  end interface

  interface 
    module subroutine i_mvect_zero(x)
      class(psb_i_multivect_type), intent(inout)    :: x
    end subroutine i_mvect_zero
  end interface

  interface 
    module subroutine i_mvect_asb(m,n, x, info)
      integer(psb_ipk_), intent(in)              :: m,n
      class(psb_i_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine i_mvect_asb
  end interface

  interface 
    module subroutine i_mvect_sync(x)
      class(psb_i_multivect_type), intent(inout) :: x
    end subroutine i_mvect_sync
  end interface

  interface 
    module subroutine i_mvect_gthab(n,idx,alpha,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_ipk_) :: alpha, beta, y(:)
      class(psb_i_multivect_type) :: x
    end subroutine i_mvect_gthab
  end interface

  interface 
    module subroutine i_mvect_gthzv(n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_ipk_) ::  y(:)
      class(psb_i_multivect_type) :: x
    end subroutine i_mvect_gthzv
  end interface

  interface 
    module subroutine i_mvect_gthzv_x(i,n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      integer(psb_ipk_) ::  y(:)
      class(psb_i_multivect_type) :: x
    end subroutine i_mvect_gthzv_x
  end interface

  interface 
    module subroutine i_mvect_sctb(n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_ipk_) :: beta, x(:)
      class(psb_i_multivect_type) :: y
    end subroutine i_mvect_sctb
  end interface

  interface 
    module subroutine i_mvect_sctb_x(i,n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      integer(psb_ipk_) :: beta, x(:)
      class(psb_i_multivect_type) :: y
    end subroutine i_mvect_sctb_x
  end interface

  interface 
    module subroutine i_mvect_free(x, info)
      class(psb_i_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine i_mvect_free
  end interface

  interface 
    module subroutine i_mvect_ins(n,irl,val,x,maxr,info)
      class(psb_i_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n,maxr
      integer(psb_ipk_), intent(in)               :: irl(:)
      integer(psb_ipk_), intent(in)        :: val(:,:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine i_mvect_ins
  end interface

  interface 
    module subroutine i_mvect_cnv(x,mold)
      class(psb_i_multivect_type), intent(inout) :: x
      class(psb_i_base_multivect_type), intent(in), optional :: mold
    end subroutine i_mvect_cnv
  end interface


contains
  
 function constructor(x) result(this)
    integer(psb_ipk_)   :: x(:,:)
    type(psb_i_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld_x(x)
    call this%asb(size(x,dim=1,kind=psb_ipk_),size(x,dim=2,kind=psb_ipk_),info)

  end function constructor

  function size_const(m,n) result(this)
    integer(psb_ipk_), intent(in) :: m,n
    type(psb_i_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld_n(m,n)
    call this%asb(m,n,info)

  end function size_const

end module psb_i_multivect_mod
