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

  public  :: psb_l_vect, psb_l_vect_type,&
       & psb_l_set_vect_default, psb_l_get_vect_default, &
       & psb_l_clear_vect_default, psb_l_base_vect_type
  
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
  
  subroutine  psb_l_set_vect_default(v)
    class(psb_l_base_vect_type), intent(in) :: v

    if (allocated(psb_l_base_vect_default)) then
      deallocate(psb_l_base_vect_default)
    end if
    allocate(psb_l_base_vect_default, mold=v)

  end subroutine psb_l_set_vect_default

  function psb_l_get_vect_default(v) result(res)
    class(psb_l_vect_type), intent(in) :: v
    class(psb_l_base_vect_type), pointer :: res

    res => psb_l_get_base_vect_default()

  end function psb_l_get_vect_default

  subroutine  psb_l_clear_vect_default()

    if (allocated(psb_l_base_vect_default)) then
      deallocate(psb_l_base_vect_default)
    end if

  end subroutine psb_l_clear_vect_default

  function psb_l_get_base_vect_default() result(res)
    class(psb_l_base_vect_type), pointer :: res

    if (.not.allocated(psb_l_base_vect_default)) then
      allocate(psb_l_base_vect_type :: psb_l_base_vect_default)
    end if

    res => psb_l_base_vect_default

  end function psb_l_get_base_vect_default

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
       & psb_l_set_multivect_default, psb_l_get_base_multivect_default, &
       & psb_l_clear_multivect_default, psb_l_base_multivect_type

  interface psb_l_multivect
    module procedure constructor, size_const
  end interface psb_l_multivect

  private

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
      integer(psb_ipk_) :: info
    end subroutine l_mvect_cnv
  end interface


contains

  subroutine  psb_l_set_multivect_default(v)
    class(psb_l_base_multivect_type), intent(in) :: v

    if (allocated(psb_l_base_multivect_default)) then
      deallocate(psb_l_base_multivect_default)
    end if
    allocate(psb_l_base_multivect_default, mold=v)

  end subroutine psb_l_set_multivect_default

!!$  function psb_l_get_multivect_default(v) result(res)
!!$    class(psb_l_multivect_type), intent(in) :: v
!!$    class(psb_l_base_multivect_type), pointer :: res
!!$
!!$    res => psb_l_get_base_multivect_default()
!!$
!!$  end function psb_l_get_multivect_default
!!$

  function psb_l_get_base_multivect_default() result(res)
    class(psb_l_base_multivect_type), pointer :: res

    if (.not.allocated(psb_l_base_multivect_default)) then
      allocate(psb_l_base_multivect_type :: psb_l_base_multivect_default)
    end if

    res => psb_l_base_multivect_default

  end function psb_l_get_base_multivect_default

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

  
  subroutine  psb_l_clear_multivect_default()

    if (allocated(psb_l_base_multivect_default)) then
      deallocate(psb_l_base_multivect_default)
    end if

  end subroutine psb_l_clear_multivect_default

end module psb_l_multivect_mod
