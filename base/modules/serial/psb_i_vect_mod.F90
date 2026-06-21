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

  public  :: psb_i_vect, psb_i_vect_type,&
       & psb_i_set_vect_default, psb_i_get_vect_default, &
       & psb_i_clear_vect_default, psb_i_base_vect_type
  
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
  
  subroutine  psb_i_set_vect_default(v)
    class(psb_i_base_vect_type), intent(in) :: v

    if (allocated(psb_i_base_vect_default)) then
      deallocate(psb_i_base_vect_default)
    end if
    allocate(psb_i_base_vect_default, mold=v)

  end subroutine psb_i_set_vect_default

  function psb_i_get_vect_default(v) result(res)
    class(psb_i_vect_type), intent(in) :: v
    class(psb_i_base_vect_type), pointer :: res

    res => psb_i_get_base_vect_default()

  end function psb_i_get_vect_default

  subroutine  psb_i_clear_vect_default()

    if (allocated(psb_i_base_vect_default)) then
      deallocate(psb_i_base_vect_default)
    end if

  end subroutine psb_i_clear_vect_default

  function psb_i_get_base_vect_default() result(res)
    class(psb_i_base_vect_type), pointer :: res

    if (.not.allocated(psb_i_base_vect_default)) then
      allocate(psb_i_base_vect_type :: psb_i_base_vect_default)
    end if

    res => psb_i_base_vect_default

  end function psb_i_get_base_vect_default

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
       & psb_i_set_multivect_default, psb_i_get_base_multivect_default, &
       & psb_i_clear_multivect_default, psb_i_base_multivect_type

  interface psb_i_multivect
    module procedure constructor, size_const
  end interface psb_i_multivect

  private

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
      integer(psb_ipk_) :: info
    end subroutine i_mvect_cnv
  end interface


contains

  subroutine  psb_i_set_multivect_default(v)
    class(psb_i_base_multivect_type), intent(in) :: v

    if (allocated(psb_i_base_multivect_default)) then
      deallocate(psb_i_base_multivect_default)
    end if
    allocate(psb_i_base_multivect_default, mold=v)

  end subroutine psb_i_set_multivect_default

!!$  function psb_i_get_multivect_default(v) result(res)
!!$    class(psb_i_multivect_type), intent(in) :: v
!!$    class(psb_i_base_multivect_type), pointer :: res
!!$
!!$    res => psb_i_get_base_multivect_default()
!!$
!!$  end function psb_i_get_multivect_default
!!$

  function psb_i_get_base_multivect_default() result(res)
    class(psb_i_base_multivect_type), pointer :: res

    if (.not.allocated(psb_i_base_multivect_default)) then
      allocate(psb_i_base_multivect_type :: psb_i_base_multivect_default)
    end if

    res => psb_i_base_multivect_default

  end function psb_i_get_base_multivect_default

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

  
  subroutine  psb_i_clear_multivect_default()

    if (allocated(psb_i_base_multivect_default)) then
      deallocate(psb_i_base_multivect_default)
    end if

  end subroutine psb_i_clear_multivect_default

end module psb_i_multivect_mod
