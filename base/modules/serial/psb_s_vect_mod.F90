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
! package: psb_s_vect_mod
!
! This module contains the definition of the psb_s_vect type which
! is the outer container for dense vectors.
! Therefore all methods simply invoke the corresponding methods of the
! inner component.
!
module psb_s_vect_mod

  use psb_s_base_vect_mod
  use psb_i_vect_mod

  type psb_s_vect_type
    class(psb_s_base_vect_type), allocatable :: v
  contains
    procedure, pass(x) :: get_nrows => s_vect_get_nrows
    procedure, pass(x) :: sizeof   => s_vect_sizeof
    procedure, pass(x) :: get_fmt  => s_vect_get_fmt
    procedure, pass(x) :: all      => s_vect_all
    procedure, pass(x) :: reall    => s_vect_reall
    procedure, pass(x) :: zero     => s_vect_zero
    procedure, pass(x) :: asb      => s_vect_asb
    procedure, pass(x) :: gthab    => s_vect_gthab
    procedure, pass(x) :: gthzv    => s_vect_gthzv
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => s_vect_sctb
    generic, public    :: sct      => sctb
    procedure, pass(x) :: free     => s_vect_free
    procedure, pass(x) :: ins_a    => s_vect_ins_a
    procedure, pass(x) :: ins_v    => s_vect_ins_v
    generic, public    :: ins      => ins_v, ins_a
    procedure, pass(x) :: bld_x    => s_vect_bld_x
    procedure, pass(x) :: bld_mn   => s_vect_bld_mn
    procedure, pass(x) :: bld_en   => s_vect_bld_en
    generic, public    :: bld      => bld_x, bld_mn, bld_en
    procedure, pass(x) :: get_vect => s_vect_get_vect
    procedure, pass(x) :: cnv      => s_vect_cnv
    procedure, pass(x) :: set_scal => s_vect_set_scal
    procedure, pass(x) :: set_vect => s_vect_set_vect
    generic, public    :: set      => set_vect, set_scal
    procedure, pass(x) :: clone    => s_vect_clone

    procedure, pass(x) :: sync     => s_vect_sync
    procedure, pass(x) :: is_host  => s_vect_is_host
    procedure, pass(x) :: is_dev   => s_vect_is_dev
    procedure, pass(x) :: is_sync  => s_vect_is_sync
    procedure, pass(x) :: set_host => s_vect_set_host
    procedure, pass(x) :: set_dev  => s_vect_set_dev
    procedure, pass(x) :: set_sync => s_vect_set_sync

    procedure, pass(x) :: dot_v    => s_vect_dot_v
    procedure, pass(x) :: dot_a    => s_vect_dot_a
    generic, public    :: dot      => dot_v, dot_a
    procedure, pass(y) :: axpby_v  => s_vect_axpby_v
    procedure, pass(y) :: axpby_a  => s_vect_axpby_a
    generic, public    :: axpby    => axpby_v, axpby_a
    procedure, pass(y) :: mlt_v    => s_vect_mlt_v
    procedure, pass(y) :: mlt_a    => s_vect_mlt_a
    procedure, pass(z) :: mlt_a_2  => s_vect_mlt_a_2
    procedure, pass(z) :: mlt_v_2  => s_vect_mlt_v_2
    procedure, pass(z) :: mlt_va   => s_vect_mlt_va
    procedure, pass(z) :: mlt_av   => s_vect_mlt_av
    generic, public    :: mlt      => mlt_v, mlt_a, mlt_a_2,&
         & mlt_v_2, mlt_av, mlt_va
    procedure, pass(x) :: div_v    => s_vect_div_v
    procedure, pass(x) :: div_v_check => s_vect_div_v_check
    procedure, pass(z) :: div_a2   => s_vect_div_a2
    procedure, pass(z) :: div_a2_check => s_vect_div_a2_check
    generic, public    :: div      => div_v, div_v_check, div_a2, div_a2_check
    procedure, pass(y) :: inv_v    => s_vect_inv_v
    procedure, pass(y) :: inv_v_check => s_vect_inv_v_check
    procedure, pass(y) :: inv_a2   => s_vect_inv_a2
    procedure, pass(y) :: inv_a2_check => s_vect_inv_a2_check
    generic, public    :: inv      => inv_v, inv_v_check, inv_a2, inv_a2_check
    procedure, pass(x) :: scal     => s_vect_scal
    procedure, pass(x) :: absval1  => s_vect_absval1
    procedure, pass(x) :: absval2  => s_vect_absval2
    generic, public    :: absval   => absval1, absval2
    procedure, pass(x) :: nrm2std  => s_vect_nrm2
    procedure, pass(x) :: nrm2weight => s_vect_nrm2_weight
    procedure, pass(x) :: nrm2weightmask => s_vect_nrm2_weight_mask
    generic, public    :: nrm2     => nrm2std, nrm2weight, nrm2weightmask
    procedure, pass(x) :: amax     => s_vect_amax
    procedure, pass(x) :: asum     => s_vect_asum
    procedure, pass(z) :: acmp_a2   => s_vect_acmp_a2
    procedure, pass(z) :: acmp_v2   => s_vect_acmp_v2
    generic, public    :: acmp      => acmp_a2, acmp_v2

    procedure, pass(m) :: mask_v => s_vect_mask_v
    procedure, pass(m) :: mask_a => s_vect_mask_a
    generic, public    :: mask => mask_a, mask_v

  end type psb_s_vect_type

  public  :: psb_s_vect
  private :: constructor, size_const
  interface psb_s_vect
    module procedure constructor, size_const
  end interface psb_s_vect

  private :: s_vect_get_nrows, s_vect_sizeof, s_vect_get_fmt, &
       & s_vect_all, s_vect_reall, s_vect_zero,  s_vect_asb, &
       & s_vect_gthab, s_vect_gthzv, s_vect_sctb, &
       & s_vect_free, s_vect_ins_a, s_vect_ins_v, s_vect_bld_x, &
       & s_vect_bld_mn, s_vect_bld_en, s_vect_get_vect, &
       & s_vect_cnv, s_vect_set_scal, &
       & s_vect_set_vect, s_vect_clone, s_vect_sync, s_vect_is_host, &
       & s_vect_is_dev, s_vect_is_sync, s_vect_set_host, &
       & s_vect_set_dev, s_vect_set_sync

  private ::  s_vect_dot_v, s_vect_dot_a, s_vect_axpby_v, s_vect_axpby_a, &
       & s_vect_mlt_v, s_vect_mlt_a, s_vect_mlt_a_2, s_vect_mlt_v_2, &
       & s_vect_mlt_va, s_vect_mlt_av, s_vect_scal, s_vect_absval1, &
       & s_vect_absval2, s_vect_nrm2, s_vect_amax, s_vect_asum


!  @NOTCPLXS@
!  @NOTINTS@
!  private :: s_vect_acmp_a2, s_vect_acmp_v2
!  @NOTINTE@
!  @NOTCPLXE@

  class(psb_s_base_vect_type), allocatable, target,&
       & save, private :: psb_s_base_vect_default

  interface psb_set_vect_default
    module procedure psb_s_set_vect_default
  end interface psb_set_vect_default

  interface psb_get_vect_default
    module procedure psb_s_get_vect_default
  end interface psb_get_vect_default


contains


  subroutine  psb_s_set_vect_default(v)
    implicit none
    class(psb_s_base_vect_type), intent(in) :: v

    if (allocated(psb_s_base_vect_default)) then
      deallocate(psb_s_base_vect_default)
    end if
    allocate(psb_s_base_vect_default, mold=v)

  end subroutine psb_s_set_vect_default

  function psb_s_get_vect_default(v) result(res)
    implicit none
    class(psb_s_vect_type), intent(in) :: v
    class(psb_s_base_vect_type), pointer :: res

    res => psb_s_get_base_vect_default()

  end function psb_s_get_vect_default


  function psb_s_get_base_vect_default() result(res)
    implicit none
    class(psb_s_base_vect_type), pointer :: res

    if (.not.allocated(psb_s_base_vect_default)) then
      allocate(psb_s_base_vect_type :: psb_s_base_vect_default)
    end if

    res => psb_s_base_vect_default

  end function psb_s_get_base_vect_default


  subroutine s_vect_clone(x,y,info)
    implicit none
    class(psb_s_vect_type), intent(inout) :: x
    class(psb_s_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    call y%free(info)
    if ((info==0).and.allocated(x%v)) then
      call y%bld(x%get_vect(),mold=x%v)
    end if
  end subroutine s_vect_clone

  subroutine s_vect_bld_x(x,invect,mold)
    real(psb_spk_), intent(in)          :: invect(:)
    class(psb_s_vect_type), intent(inout) :: x
    class(psb_s_base_vect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info

    info = psb_success_
    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_s_get_base_vect_default())
    endif

    if (info == psb_success_) call x%v%bld(invect)

  end subroutine s_vect_bld_x


  subroutine s_vect_bld_mn(x,n,mold)
    integer(psb_mpk_), intent(in) :: n
    class(psb_s_vect_type), intent(inout) :: x
    class(psb_s_base_vect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info
    class(psb_s_base_vect_type), pointer :: mld

    info = psb_success_
    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_s_get_base_vect_default())
    endif
    if (info == psb_success_) call x%v%bld(n)

  end subroutine s_vect_bld_mn


  subroutine s_vect_bld_en(x,n,mold)
    integer(psb_epk_), intent(in) :: n
    class(psb_s_vect_type), intent(inout) :: x
    class(psb_s_base_vect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info

    info = psb_success_

    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_s_get_base_vect_default())
    endif
    if (info == psb_success_) call x%v%bld(n)

  end subroutine s_vect_bld_en

  function  s_vect_get_vect(x,n) result(res)
    class(psb_s_vect_type), intent(inout)  :: x
    real(psb_spk_), allocatable                 :: res(:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional :: n

    if (allocated(x%v)) then
      res = x%v%get_vect(n)
    end if
  end function s_vect_get_vect

  subroutine s_vect_set_scal(x,val,first,last)
    class(psb_s_vect_type), intent(inout)  :: x
    real(psb_spk_), intent(in) :: val
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val,first,last)

  end subroutine s_vect_set_scal

  subroutine s_vect_set_vect(x,val,first,last)
    class(psb_s_vect_type), intent(inout) :: x
    real(psb_spk_), intent(in)         :: val(:)
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val,first,last)

  end subroutine s_vect_set_vect


  function constructor(x) result(this)
    real(psb_spk_)   :: x(:)
    type(psb_s_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(x)
    call this%asb(size(x,kind=psb_ipk_),info)

  end function constructor


  function size_const(n) result(this)
    integer(psb_ipk_), intent(in) :: n
    type(psb_s_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(n)
    call this%asb(n,info)

  end function size_const

  function s_vect_get_nrows(x) result(res)
    implicit none
    class(psb_s_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_nrows()
  end function s_vect_get_nrows

  function s_vect_sizeof(x) result(res)
    implicit none
    class(psb_s_vect_type), intent(in) :: x
    integer(psb_epk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%sizeof()
  end function s_vect_sizeof

  function s_vect_get_fmt(x) result(res)
    implicit none
    class(psb_s_vect_type), intent(in) :: x
    character(len=5) :: res
    res = 'NULL'
    if (allocated(x%v)) res = x%v%get_fmt()
  end function s_vect_get_fmt

  subroutine s_vect_all(n, x, info, mold)

    implicit none
    integer(psb_ipk_), intent(in)           :: n
    class(psb_s_vect_type), intent(inout) :: x
    class(psb_s_base_vect_type), intent(in), optional :: mold
    integer(psb_ipk_), intent(out)      :: info

    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(psb_s_base_vect_type :: x%v,stat=info)
    endif
    if (info == 0) then
      call x%v%all(n,info)
    else
      info = psb_err_alloc_dealloc_
    end if

  end subroutine s_vect_all

  subroutine s_vect_reall(n, x, info)

    implicit none
    integer(psb_ipk_), intent(in)         :: n
    class(psb_s_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    info = 0
    if (.not.allocated(x%v)) &
         & call x%all(n,info)
    if (info == 0) &
         & call x%asb(n,info)

  end subroutine s_vect_reall

  subroutine s_vect_zero(x)
    use psi_serial_mod
    implicit none
    class(psb_s_vect_type), intent(inout)    :: x

    if (allocated(x%v)) call x%v%zero()

  end subroutine s_vect_zero

  subroutine s_vect_asb(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    integer(psb_ipk_), intent(in)              :: n
    class(psb_s_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%v)) &
         & call x%v%asb(n,info)

  end subroutine s_vect_asb

  subroutine s_vect_gthab(n,idx,alpha,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    real(psb_spk_) :: alpha, beta, y(:)
    class(psb_s_vect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,alpha,beta,y)

  end subroutine s_vect_gthab

  subroutine s_vect_gthzv(n,idx,x,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    real(psb_spk_) ::  y(:)
    class(psb_s_vect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,y)

  end subroutine s_vect_gthzv

  subroutine s_vect_sctb(n,idx,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    real(psb_spk_) :: beta, x(:)
    class(psb_s_vect_type) :: y

    if (allocated(y%v)) &
         &  call y%v%sct(n,idx,x,beta)

  end subroutine s_vect_sctb

  subroutine s_vect_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    class(psb_s_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (allocated(x%v)) then
      call x%v%free(info)
      if (info == 0) deallocate(x%v,stat=info)
    end if

  end subroutine s_vect_free

  subroutine s_vect_ins_a(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none
    class(psb_s_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    integer(psb_ipk_), intent(in)               :: irl(:)
    real(psb_spk_), intent(in)        :: val(:)
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i

    info = 0
    if (.not.allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if

    call  x%v%ins(n,irl,val,dupl,info)

  end subroutine s_vect_ins_a

  subroutine s_vect_ins_v(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none
    class(psb_s_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    class(psb_i_vect_type), intent(inout)       :: irl
    class(psb_s_vect_type), intent(inout)       :: val
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i

    info = 0
    if (.not.(allocated(x%v).and.allocated(irl%v).and.allocated(val%v))) then
      info = psb_err_invalid_vect_state_
      return
    end if

    call  x%v%ins(n,irl%v,val%v,dupl,info)

  end subroutine s_vect_ins_v


  subroutine s_vect_cnv(x,mold)
    class(psb_s_vect_type), intent(inout) :: x
    class(psb_s_base_vect_type), intent(in), optional :: mold
    class(psb_s_base_vect_type), allocatable :: tmp

    integer(psb_ipk_) :: info

    info = psb_success_
    if (present(mold)) then
      allocate(tmp,stat=info,mold=mold)
    else
      allocate(tmp,stat=info,mold=psb_s_get_base_vect_default())
    end if
    if (allocated(x%v)) then
      call x%v%sync()
      if (info == psb_success_) call tmp%bld(x%v%v)
      call x%v%free(info)
    end if
    call move_alloc(tmp,x%v)

  end subroutine s_vect_cnv


  subroutine s_vect_sync(x)
    implicit none
    class(psb_s_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%sync()

  end subroutine s_vect_sync

  subroutine s_vect_set_sync(x)
    implicit none
    class(psb_s_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%set_sync()

  end subroutine s_vect_set_sync

  subroutine s_vect_set_host(x)
    implicit none
    class(psb_s_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%set_host()

  end subroutine s_vect_set_host

  subroutine s_vect_set_dev(x)
    implicit none
    class(psb_s_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%set_dev()

  end subroutine s_vect_set_dev

  function s_vect_is_sync(x) result(res)
    implicit none
    logical :: res
    class(psb_s_vect_type), intent(inout) :: x

    res = .true.
    if (allocated(x%v)) &
         & res = x%v%is_sync()

  end function s_vect_is_sync

  function s_vect_is_host(x) result(res)
    implicit none
    logical :: res
    class(psb_s_vect_type), intent(inout) :: x

    res = .true.
    if (allocated(x%v)) &
         & res = x%v%is_host()

  end function s_vect_is_host

  function s_vect_is_dev(x) result(res)
    implicit none
    logical :: res
    class(psb_s_vect_type), intent(inout) :: x

    res = .false.
    if (allocated(x%v)) &
         & res =  x%v%is_dev()

  end function s_vect_is_dev


  function s_vect_dot_v(n,x,y) result(res)
    implicit none
    class(psb_s_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)                :: res

    res = szero
    if (allocated(x%v).and.allocated(y%v)) &
         & res = x%v%dot(n,y%v)

  end function s_vect_dot_v

  function s_vect_dot_a(n,x,y) result(res)
    implicit none
    class(psb_s_vect_type), intent(inout) :: x
    real(psb_spk_), intent(in)    :: y(:)
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)                :: res

    res = szero
    if (allocated(x%v)) &
         & res = x%v%dot(n,y)

  end function s_vect_dot_a

  subroutine s_vect_axpby_v(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    class(psb_s_vect_type), intent(inout)  :: x
    class(psb_s_vect_type), intent(inout)  :: y
    real(psb_spk_), intent (in)       :: alpha, beta
    integer(psb_ipk_), intent(out)              :: info

    if (allocated(x%v).and.allocated(y%v)) then
      call y%v%axpby(m,alpha,x%v,beta,info)
    else
      info = psb_err_invalid_vect_state_
    end if

  end subroutine s_vect_axpby_v

  subroutine s_vect_axpby_a(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    real(psb_spk_), intent(in)        :: x(:)
    class(psb_s_vect_type), intent(inout)  :: y
    real(psb_spk_), intent (in)       :: alpha, beta
    integer(psb_ipk_), intent(out)              :: info

    if (allocated(y%v)) &
         & call y%v%axpby(m,alpha,x,beta,info)

  end subroutine s_vect_axpby_a


  subroutine s_vect_mlt_v(x, y, info)
    use psi_serial_mod
    implicit none
    class(psb_s_vect_type), intent(inout)  :: x
    class(psb_s_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(x%v).and.allocated(y%v)) &
         & call y%v%mlt(x%v,info)

  end subroutine s_vect_mlt_v

  subroutine s_vect_mlt_a(x, y, info)
    use psi_serial_mod
    implicit none
    real(psb_spk_), intent(in)        :: x(:)
    class(psb_s_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n


    info = 0
    if (allocated(y%v)) &
         & call y%v%mlt(x,info)

  end subroutine s_vect_mlt_a


  subroutine s_vect_mlt_a_2(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none
    real(psb_spk_), intent(in)         :: alpha,beta
    real(psb_spk_), intent(in)         :: y(:)
    real(psb_spk_), intent(in)         :: x(:)
    class(psb_s_vect_type), intent(inout) :: z
    integer(psb_ipk_), intent(out)                  :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(z%v)) &
         & call z%v%mlt(alpha,x,y,beta,info)

  end subroutine s_vect_mlt_a_2

  subroutine s_vect_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
    use psi_serial_mod
    implicit none
    real(psb_spk_), intent(in)          :: alpha,beta
    class(psb_s_vect_type), intent(inout)  :: x
    class(psb_s_vect_type), intent(inout)  :: y
    class(psb_s_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)                   :: info
    character(len=1), intent(in), optional :: conjgx, conjgy

    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(x%v).and.allocated(y%v).and.&
         & allocated(z%v)) &
         & call z%v%mlt(alpha,x%v,y%v,beta,info,conjgx,conjgy)

  end subroutine s_vect_mlt_v_2

  subroutine s_vect_mlt_av(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none
    real(psb_spk_), intent(in)        :: alpha,beta
    real(psb_spk_), intent(in)        :: x(:)
    class(psb_s_vect_type), intent(inout)  :: y
    class(psb_s_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(z%v).and.allocated(y%v)) &
         & call z%v%mlt(alpha,x,y%v,beta,info)

  end subroutine s_vect_mlt_av

  subroutine s_vect_mlt_va(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none
    real(psb_spk_), intent(in)        :: alpha,beta
    real(psb_spk_), intent(in)        :: y(:)
    class(psb_s_vect_type), intent(inout)  :: x
    class(psb_s_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0

    if (allocated(z%v).and.allocated(x%v)) &
         & call z%v%mlt(alpha,x%v,y,beta,info)

  end subroutine s_vect_mlt_va

  subroutine s_vect_div_v(x, y, info)
    use psi_serial_mod
    implicit none
    class(psb_s_vect_type), intent(inout)  :: x
    class(psb_s_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(x%v).and.allocated(y%v)) &
         & call x%v%div(y%v,info)

  end subroutine s_vect_div_v

  subroutine s_vect_div_v_check(x, y, info, flag)
    use psi_serial_mod
    implicit none
    class(psb_s_vect_type), intent(inout)  :: x
    class(psb_s_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n
    logical, intent(in) :: flag

    info = 0
    if (allocated(x%v).and.allocated(y%v)) &
         & call x%v%div(y%v,info,flag)

  end subroutine s_vect_div_v_check

  subroutine s_vect_div_a2(x, y, z, info)
    use psi_serial_mod
    implicit none
    real(psb_spk_), intent(in) :: x(:)
    real(psb_spk_), intent(in)    :: y(:)
    class(psb_s_vect_type), intent(inout) :: z
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(z%v)) &
         & call z%v%div(x,y,info)

  end subroutine s_vect_div_a2

  subroutine s_vect_div_a2_check(x, y, z, info,flag)
    use psi_serial_mod
    implicit none
    real(psb_spk_), intent(in) :: x(:)
    real(psb_spk_), intent(in)    :: y(:)
    class(psb_s_vect_type), intent(inout) :: z
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n
    logical, intent(in) :: flag

    info = 0
    if (allocated(z%v)) &
         & call z%v%div(x,y,info,flag)

  end subroutine s_vect_div_a2_check

  subroutine s_vect_inv_v(x, y, info)
    use psi_serial_mod
    implicit none
    class(psb_s_vect_type), intent(inout)  :: x
    class(psb_s_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(x%v).and.allocated(y%v)) &
         & call y%v%inv(x%v,info)

  end subroutine s_vect_inv_v

  subroutine s_vect_inv_v_check(x, y, info, flag)
    use psi_serial_mod
    implicit none
    class(psb_s_vect_type), intent(inout)  :: x
    class(psb_s_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n
    logical, intent(in) :: flag

    info = 0
    if (allocated(x%v).and.allocated(y%v)) &
         & call y%v%inv(x%v,info,flag)

  end subroutine s_vect_inv_v_check

  subroutine s_vect_inv_a2(x, y, info)
    use psi_serial_mod
    implicit none
    real(psb_spk_), intent(inout) :: x(:)
    class(psb_s_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(y%v)) &
         & call y%v%inv(x,info)

  end subroutine s_vect_inv_a2

  subroutine s_vect_inv_a2_check(x, y, info,flag)
    use psi_serial_mod
    implicit none
    real(psb_spk_), intent(inout) :: x(:)
    class(psb_s_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n
    logical, intent(in) :: flag

    info = 0
    if (allocated(y%v)) &
         & call y%v%inv(x,info,flag)

  end subroutine s_vect_inv_a2_check

  subroutine s_vect_acmp_a2(x,c,z,info)
    use psi_serial_mod
    implicit none
    real(psb_spk_), intent(in)              :: c
    real(psb_spk_), intent(inout)           :: x(:)
    class(psb_s_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)           :: info

    info = 0
    if (allocated(z%v)) &
         & call z%acmp(x,c,info)

  end subroutine s_vect_acmp_a2

  subroutine s_vect_acmp_v2(x,c,z,info)
    use psi_serial_mod
    implicit none
    real(psb_spk_), intent(in)              :: c
    class(psb_s_vect_type), intent(inout)  :: x
    class(psb_s_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)           :: info

    info = 0
    if (allocated(x%v).and.allocated(z%v)) &
         & call z%v%acmp(x%v,c,info)

  end subroutine s_vect_acmp_v2

  subroutine s_vect_scal(alpha, x)
    use psi_serial_mod
    implicit none
    class(psb_s_vect_type), intent(inout)  :: x
    real(psb_spk_), intent (in)       :: alpha

    if (allocated(x%v)) call x%v%scal(alpha)

  end subroutine s_vect_scal

  subroutine s_vect_absval1(x)
    class(psb_s_vect_type), intent(inout)  :: x

    if (allocated(x%v)) &
         &  call x%v%absval()

  end subroutine s_vect_absval1

  subroutine s_vect_absval2(x,y)
    class(psb_s_vect_type), intent(inout)  :: x
    class(psb_s_vect_type), intent(inout)  :: y

    if (allocated(x%v)) then
      if (.not.allocated(y%v))  call y%bld(psb_size(x%v%v))
      call x%v%absval(y%v)
    end if
  end subroutine s_vect_absval2

  function s_vect_nrm2(n,x) result(res)
    implicit none
    class(psb_s_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)                :: res

    if (allocated(x%v)) then
      res = x%v%nrm2(n)
    else
      res = szero
    end if

  end function s_vect_nrm2

  function s_vect_nrm2_weight(n,x,w) result(res)
    implicit none
    class(psb_s_vect_type), intent(inout) :: x
    class(psb_s_vect_type), intent(inout) :: w
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)                        :: res
    integer(psb_ipk_)                       :: info

    if (allocated(x%v).and.allocated(w%v)) then
      call w%v%mlt(x%v,info)
      res = w%v%nrm2(n)
    else
      res = szero
    end if

  end function s_vect_nrm2_weight

  function s_vect_nrm2_weight_mask(n,x,w,id) result(res)
    implicit none
    class(psb_s_vect_type), intent(inout) :: x
    class(psb_s_vect_type), intent(inout) :: w
    class(psb_s_vect_type), intent(inout) :: id
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)                        :: res
    integer(psb_ipk_)                       :: info

    if (allocated(x%v).and.allocated(w%v).and.allocated(id%v)) then
      call w%v%acmp(id%v,szero,info)
      call w%v%mlt(x%v,info)
      res = w%v%nrm2(n)
    else
      res = szero
    end if

  end function s_vect_nrm2_weight_mask

  function s_vect_amax(n,x) result(res)
    implicit none
    class(psb_s_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)                :: res

    if (allocated(x%v)) then
      res = x%v%amax(n)
    else
      res = szero
    end if

  end function s_vect_amax

  function s_vect_asum(n,x) result(res)
    implicit none
    class(psb_s_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)                :: res

    if (allocated(x%v)) then
      res = x%v%asum(n)
    else
      res = szero
    end if

  end function s_vect_asum


  subroutine s_vect_mask_a(c,x,m,t,info)
    use psi_serial_mod
    implicit none
    real(psb_spk_), intent(inout)          :: c(:)
    real(psb_spk_), intent(inout)          :: x(:)
    logical, intent(out)                     :: t;
    class(psb_s_vect_type), intent(inout)  :: m
    integer(psb_ipk_), intent(out)           :: info

    info = 0
    if (allocated(m%v)) &
         & call m%mask(c,x,t,info)

  end subroutine s_vect_mask_a

  subroutine s_vect_mask_v(c,x,m,t,info)
    use psi_serial_mod
    implicit none
    class(psb_s_vect_type), intent(inout)  :: c
    class(psb_s_vect_type), intent(inout)  :: x
    class(psb_s_vect_type), intent(inout)  :: m
    logical, intent(out)                     :: t;
    integer(psb_ipk_), intent(out)           :: info

    info = 0
    if (allocated(x%v).and.allocated(c%v)) &
         & call m%v%mask(x%v,c%v,t,info)

  end subroutine s_vect_mask_v

end module psb_s_vect_mod



module psb_s_multivect_mod

  use psb_s_base_multivect_mod
  use psb_const_mod
  use psb_i_vect_mod


  !private

  type psb_s_multivect_type
    class(psb_s_base_multivect_type), allocatable :: v
  contains
    procedure, pass(x) :: get_nrows => s_vect_get_nrows
    procedure, pass(x) :: get_ncols => s_vect_get_ncols
    procedure, pass(x) :: sizeof   => s_vect_sizeof
    procedure, pass(x) :: get_fmt  => s_vect_get_fmt

    procedure, pass(x) :: all      => s_vect_all
    procedure, pass(x) :: reall    => s_vect_reall
    procedure, pass(x) :: zero     => s_vect_zero
    procedure, pass(x) :: asb      => s_vect_asb
    procedure, pass(x) :: sync     => s_vect_sync
    procedure, pass(x) :: free     => s_vect_free
    procedure, pass(x) :: ins      => s_vect_ins
    procedure, pass(x) :: bld_x    => s_vect_bld_x
    procedure, pass(x) :: bld_n    => s_vect_bld_n
    generic, public    :: bld      => bld_x, bld_n
    procedure, pass(x) :: get_vect => s_vect_get_vect
    procedure, pass(x) :: cnv      => s_vect_cnv
    procedure, pass(x) :: set_scal => s_vect_set_scal
    procedure, pass(x) :: set_vect => s_vect_set_vect
    generic, public    :: set      => set_vect, set_scal
    procedure, pass(x) :: clone    => s_vect_clone
    procedure, pass(x) :: gthab    => s_vect_gthab
    procedure, pass(x) :: gthzv    => s_vect_gthzv
    procedure, pass(x) :: gthzv_x  => s_vect_gthzv_x
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => s_vect_sctb
    procedure, pass(y) :: sctb_x   => s_vect_sctb_x
    generic, public    :: sct      => sctb, sctb_x
!!$    procedure, pass(x) :: dot_v    => s_vect_dot_v
!!$    procedure, pass(x) :: dot_a    => s_vect_dot_a
!!$    generic, public    :: dot      => dot_v, dot_a
!!$    procedure, pass(y) :: axpby_v  => s_vect_axpby_v
!!$    procedure, pass(y) :: axpby_a  => s_vect_axpby_a
!!$    generic, public    :: axpby    => axpby_v, axpby_a
!!$    procedure, pass(y) :: mlt_v    => s_vect_mlt_v
!!$    procedure, pass(y) :: mlt_a    => s_vect_mlt_a
!!$    procedure, pass(z) :: mlt_a_2  => s_vect_mlt_a_2
!!$    procedure, pass(z) :: mlt_v_2  => s_vect_mlt_v_2
!!$    procedure, pass(z) :: mlt_va   => s_vect_mlt_va
!!$    procedure, pass(z) :: mlt_av   => s_vect_mlt_av
!!$    generic, public    :: mlt      => mlt_v, mlt_a, mlt_a_2,&
!!$         & mlt_v_2, mlt_av, mlt_va
!!$    procedure, pass(x) :: scal     => s_vect_scal
!!$    procedure, pass(x) :: nrm2     => s_vect_nrm2
!!$    procedure, pass(x) :: amax     => s_vect_amax
!!$    procedure, pass(x) :: asum     => s_vect_asum
  end type psb_s_multivect_type

  public  :: psb_s_multivect, psb_s_multivect_type,&
       & psb_set_multivect_default, psb_get_multivect_default, &
       & psb_s_base_multivect_type

  private
  interface psb_s_multivect
    module procedure constructor, size_const
  end interface psb_s_multivect

  class(psb_s_base_multivect_type), allocatable, target,&
       & save, private :: psb_s_base_multivect_default

  interface psb_set_multivect_default
    module procedure psb_s_set_multivect_default
  end interface psb_set_multivect_default

  interface psb_get_multivect_default
    module procedure psb_s_get_multivect_default
  end interface psb_get_multivect_default


contains


  subroutine  psb_s_set_multivect_default(v)
    implicit none
    class(psb_s_base_multivect_type), intent(in) :: v

    if (allocated(psb_s_base_multivect_default)) then
      deallocate(psb_s_base_multivect_default)
    end if
    allocate(psb_s_base_multivect_default, mold=v)

  end subroutine psb_s_set_multivect_default

  function psb_s_get_multivect_default(v) result(res)
    implicit none
    class(psb_s_multivect_type), intent(in) :: v
    class(psb_s_base_multivect_type), pointer :: res

    res => psb_s_get_base_multivect_default()

  end function psb_s_get_multivect_default


  function psb_s_get_base_multivect_default() result(res)
    implicit none
    class(psb_s_base_multivect_type), pointer :: res

    if (.not.allocated(psb_s_base_multivect_default)) then
      allocate(psb_s_base_multivect_type :: psb_s_base_multivect_default)
    end if

    res => psb_s_base_multivect_default

  end function psb_s_get_base_multivect_default


  subroutine s_vect_clone(x,y,info)
    implicit none
    class(psb_s_multivect_type), intent(inout) :: x
    class(psb_s_multivect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    call y%free(info)
    if ((info==0).and.allocated(x%v)) then
      call y%bld(x%get_vect(),mold=x%v)
    end if
  end subroutine s_vect_clone

  subroutine s_vect_bld_x(x,invect,mold)
    real(psb_spk_), intent(in)          :: invect(:,:)
    class(psb_s_multivect_type), intent(out) :: x
    class(psb_s_base_multivect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info
    class(psb_s_base_multivect_type), pointer :: mld

    info = psb_success_
    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_s_get_base_multivect_default())
    endif

    if (info == psb_success_) call x%v%bld(invect)

  end subroutine s_vect_bld_x


  subroutine s_vect_bld_n(x,m,n,mold)
    integer(psb_ipk_), intent(in) :: m,n
    class(psb_s_multivect_type), intent(out) :: x
    class(psb_s_base_multivect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info

    info = psb_success_
    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_s_get_base_multivect_default())
    endif
    if (info == psb_success_) call x%v%bld(m,n)

  end subroutine s_vect_bld_n

  function  s_vect_get_vect(x) result(res)
    class(psb_s_multivect_type), intent(inout)  :: x
    real(psb_spk_), allocatable                 :: res(:,:)
    integer(psb_ipk_) :: info

    if (allocated(x%v)) then
      res = x%v%get_vect()
    end if
  end function s_vect_get_vect

  subroutine s_vect_set_scal(x,val)
    class(psb_s_multivect_type), intent(inout)  :: x
    real(psb_spk_), intent(in) :: val

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val)

  end subroutine s_vect_set_scal

  subroutine s_vect_set_vect(x,val)
    class(psb_s_multivect_type), intent(inout) :: x
    real(psb_spk_), intent(in)         :: val(:,:)

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val)

  end subroutine s_vect_set_vect


  function constructor(x) result(this)
    real(psb_spk_)   :: x(:,:)
    type(psb_s_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(x)
    call this%asb(size(x,dim=1,kind=psb_ipk_),size(x,dim=2,kind=psb_ipk_),info)

  end function constructor


  function size_const(m,n) result(this)
    integer(psb_ipk_), intent(in) :: m,n
    type(psb_s_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(m,n)
    call this%asb(m,n,info)

  end function size_const

  function s_vect_get_nrows(x) result(res)
    implicit none
    class(psb_s_multivect_type), intent(in) :: x
    integer(psb_ipk_)  :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_nrows()
  end function s_vect_get_nrows

  function s_vect_get_ncols(x) result(res)
    implicit none
    class(psb_s_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_ncols()
  end function s_vect_get_ncols

  function s_vect_sizeof(x) result(res)
    implicit none
    class(psb_s_multivect_type), intent(in) :: x
    integer(psb_epk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%sizeof()
  end function s_vect_sizeof

  function s_vect_get_fmt(x) result(res)
    implicit none
    class(psb_s_multivect_type), intent(in) :: x
    character(len=5) :: res
    res = 'NULL'
    if (allocated(x%v)) res = x%v%get_fmt()
  end function s_vect_get_fmt

  subroutine s_vect_all(m,n, x, info, mold)

    implicit none
    integer(psb_ipk_), intent(in)       :: m,n
    class(psb_s_multivect_type), intent(out) :: x
    class(psb_s_base_multivect_type), intent(in), optional :: mold
    integer(psb_ipk_), intent(out)      :: info

    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(psb_s_base_multivect_type :: x%v,stat=info)
    endif
    if (info == 0) then
      call x%v%all(m,n,info)
    else
      info = psb_err_alloc_dealloc_
    end if

  end subroutine s_vect_all

  subroutine s_vect_reall(m,n, x, info)

    implicit none
    integer(psb_ipk_), intent(in)         :: m,n
    class(psb_s_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    info = 0
    if (.not.allocated(x%v)) &
         & call x%all(m,n,info)
    if (info == 0) &
         & call x%asb(m,n,info)

  end subroutine s_vect_reall

  subroutine s_vect_zero(x)
    use psi_serial_mod
    implicit none
    class(psb_s_multivect_type), intent(inout)    :: x

    if (allocated(x%v)) call x%v%zero()

  end subroutine s_vect_zero

  subroutine s_vect_asb(m,n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    integer(psb_ipk_), intent(in)              :: m,n
    class(psb_s_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%v)) &
         & call x%v%asb(m,n,info)

  end subroutine s_vect_asb

  subroutine s_vect_sync(x)
    implicit none
    class(psb_s_multivect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%sync()

  end subroutine s_vect_sync

  subroutine s_vect_gthab(n,idx,alpha,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    real(psb_spk_) :: alpha, beta, y(:)
    class(psb_s_multivect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,alpha,beta,y)

  end subroutine s_vect_gthab

  subroutine s_vect_gthzv(n,idx,x,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    real(psb_spk_) ::  y(:)
    class(psb_s_multivect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,y)

  end subroutine s_vect_gthzv

  subroutine s_vect_gthzv_x(i,n,idx,x,y)
    use psi_serial_mod
    integer(psb_ipk_) :: i,n
    class(psb_i_base_vect_type) :: idx
    real(psb_spk_) ::  y(:)
    class(psb_s_multivect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(i,n,idx,y)

  end subroutine s_vect_gthzv_x

  subroutine s_vect_sctb(n,idx,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    real(psb_spk_) :: beta, x(:)
    class(psb_s_multivect_type) :: y

    if (allocated(y%v)) &
         &  call y%v%sct(n,idx,x,beta)

  end subroutine s_vect_sctb

  subroutine s_vect_sctb_x(i,n,idx,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    real(psb_spk_) :: beta, x(:)
    class(psb_s_multivect_type) :: y

    if (allocated(y%v)) &
         &  call y%v%sct(i,n,idx,x,beta)

  end subroutine s_vect_sctb_x

  subroutine s_vect_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    class(psb_s_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (allocated(x%v)) then
      call x%v%free(info)
      if (info == 0) deallocate(x%v,stat=info)
    end if

  end subroutine s_vect_free

  subroutine s_vect_ins(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none
    class(psb_s_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    integer(psb_ipk_), intent(in)               :: irl(:)
    real(psb_spk_), intent(in)        :: val(:,:)
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i

    info = 0
    if (.not.allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if

    call  x%v%ins(n,irl,val,dupl,info)

  end subroutine s_vect_ins


  subroutine s_vect_cnv(x,mold)
    class(psb_s_multivect_type), intent(inout) :: x
    class(psb_s_base_multivect_type), intent(in), optional :: mold
    class(psb_s_base_multivect_type), allocatable :: tmp
    integer(psb_ipk_) :: info

    if (present(mold)) then
      allocate(tmp,stat=info,mold=mold)
    else
      allocate(tmp,stat=info, mold=psb_s_get_base_multivect_default())
    endif
    if (allocated(x%v)) then
      call x%v%sync()
      if (info == psb_success_) call tmp%bld(x%v%v)
      call x%v%free(info)
    end if
    call move_alloc(tmp,x%v)
  end subroutine s_vect_cnv


!!$  function s_vect_dot_v(n,x,y) result(res)
!!$    implicit none
!!$    class(psb_s_multivect_type), intent(inout) :: x, y
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_spk_)                :: res
!!$
!!$    res = szero
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & res = x%v%dot(n,y%v)
!!$
!!$  end function s_vect_dot_v
!!$
!!$  function s_vect_dot_a(n,x,y) result(res)
!!$    implicit none
!!$    class(psb_s_multivect_type), intent(inout) :: x
!!$    real(psb_spk_), intent(in)    :: y(:)
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_spk_)                :: res
!!$
!!$    res = szero
!!$    if (allocated(x%v)) &
!!$         & res = x%v%dot(n,y)
!!$
!!$  end function s_vect_dot_a
!!$
!!$  subroutine s_vect_axpby_v(m,alpha, x, beta, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    class(psb_s_multivect_type), intent(inout)  :: x
!!$    class(psb_s_multivect_type), intent(inout)  :: y
!!$    real(psb_spk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    if (allocated(x%v).and.allocated(y%v)) then
!!$      call y%v%axpby(m,alpha,x%v,beta,info)
!!$    else
!!$      info = psb_err_invalid_vect_state_
!!$    end if
!!$
!!$  end subroutine s_vect_axpby_v
!!$
!!$  subroutine s_vect_axpby_a(m,alpha, x, beta, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    real(psb_spk_), intent(in)        :: x(:)
!!$    class(psb_s_multivect_type), intent(inout)  :: y
!!$    real(psb_spk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    if (allocated(y%v)) &
!!$         & call y%v%axpby(m,alpha,x,beta,info)
!!$
!!$  end subroutine s_vect_axpby_a
!!$
!!$
!!$  subroutine s_vect_mlt_v(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_s_multivect_type), intent(inout)  :: x
!!$    class(psb_s_multivect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & call y%v%mlt(x%v,info)
!!$
!!$  end subroutine s_vect_mlt_v
!!$
!!$  subroutine s_vect_mlt_a(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    real(psb_spk_), intent(in)        :: x(:)
!!$    class(psb_s_multivect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$
!!$    info = 0
!!$    if (allocated(y%v)) &
!!$         & call y%v%mlt(x,info)
!!$
!!$  end subroutine s_vect_mlt_a
!!$
!!$
!!$  subroutine s_vect_mlt_a_2(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    real(psb_spk_), intent(in)         :: alpha,beta
!!$    real(psb_spk_), intent(in)         :: y(:)
!!$    real(psb_spk_), intent(in)         :: x(:)
!!$    class(psb_s_multivect_type), intent(inout) :: z
!!$    integer(psb_ipk_), intent(out)                  :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(z%v)) &
!!$         & call z%v%mlt(alpha,x,y,beta,info)
!!$
!!$  end subroutine s_vect_mlt_a_2
!!$
!!$  subroutine s_vect_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
!!$    use psi_serial_mod
!!$    implicit none
!!$    real(psb_spk_), intent(in)          :: alpha,beta
!!$    class(psb_s_multivect_type), intent(inout)  :: x
!!$    class(psb_s_multivect_type), intent(inout)  :: y
!!$    class(psb_s_multivect_type), intent(inout)  :: z
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
!!$  end subroutine s_vect_mlt_v_2
!!$
!!$  subroutine s_vect_mlt_av(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    real(psb_spk_), intent(in)        :: alpha,beta
!!$    real(psb_spk_), intent(in)        :: x(:)
!!$    class(psb_s_multivect_type), intent(inout)  :: y
!!$    class(psb_s_multivect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(z%v).and.allocated(y%v)) &
!!$         & call z%v%mlt(alpha,x,y%v,beta,info)
!!$
!!$  end subroutine s_vect_mlt_av
!!$
!!$  subroutine s_vect_mlt_va(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    real(psb_spk_), intent(in)        :: alpha,beta
!!$    real(psb_spk_), intent(in)        :: y(:)
!!$    class(psb_s_multivect_type), intent(inout)  :: x
!!$    class(psb_s_multivect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$
!!$    if (allocated(z%v).and.allocated(x%v)) &
!!$         & call z%v%mlt(alpha,x%v,y,beta,info)
!!$
!!$  end subroutine s_vect_mlt_va
!!$
!!$  subroutine s_vect_scal(alpha, x)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_s_multivect_type), intent(inout)  :: x
!!$    real(psb_spk_), intent (in)       :: alpha
!!$
!!$    if (allocated(x%v)) call x%v%scal(alpha)
!!$
!!$  end subroutine s_vect_scal
!!$
!!$
!!$  function s_vect_nrm2(n,x) result(res)
!!$    implicit none
!!$    class(psb_s_multivect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_spk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%nrm2(n)
!!$    else
!!$      res = szero
!!$    end if
!!$
!!$  end function s_vect_nrm2
!!$
!!$  function s_vect_amax(n,x) result(res)
!!$    implicit none
!!$    class(psb_s_multivect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_spk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%amax(n)
!!$    else
!!$      res = szero
!!$    end if
!!$
!!$  end function s_vect_amax
!!$
!!$  function s_vect_asum(n,x) result(res)
!!$    implicit none
!!$    class(psb_s_multivect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_spk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%asum(n)
!!$    else
!!$      res = szero
!!$    end if
!!$
!!$  end function s_vect_asum

end module psb_s_multivect_mod
