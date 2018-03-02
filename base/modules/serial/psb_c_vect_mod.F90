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
! package: psb_c_vect_mod
!
! This module contains the definition of the psb_c_vect type which
! is the outer container for dense vectors.
! Therefore all methods simply invoke the corresponding methods of the
! inner component. 
!
module psb_c_vect_mod

  use psb_c_base_vect_mod
  use psb_i_vect_mod

  type psb_c_vect_type
    class(psb_c_base_vect_type), allocatable :: v 
  contains
    procedure, pass(x) :: get_nrows => c_vect_get_nrows
    procedure, pass(x) :: sizeof   => c_vect_sizeof
    procedure, pass(x) :: get_fmt  => c_vect_get_fmt
    procedure, pass(x) :: all      => c_vect_all
    procedure, pass(x) :: reall    => c_vect_reall
    procedure, pass(x) :: zero     => c_vect_zero
    procedure, pass(x) :: asb      => c_vect_asb
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

    procedure, pass(x) :: dot_v    => c_vect_dot_v
    procedure, pass(x) :: dot_a    => c_vect_dot_a
    generic, public    :: dot      => dot_v, dot_a
    procedure, pass(y) :: axpby_v  => c_vect_axpby_v
    procedure, pass(y) :: axpby_a  => c_vect_axpby_a
    generic, public    :: axpby    => axpby_v, axpby_a
    procedure, pass(y) :: mlt_v    => c_vect_mlt_v
    procedure, pass(y) :: mlt_a    => c_vect_mlt_a
    procedure, pass(z) :: mlt_a_2  => c_vect_mlt_a_2
    procedure, pass(z) :: mlt_v_2  => c_vect_mlt_v_2
    procedure, pass(z) :: mlt_va   => c_vect_mlt_va
    procedure, pass(z) :: mlt_av   => c_vect_mlt_av
    generic, public    :: mlt      => mlt_v, mlt_a, mlt_a_2,&
         & mlt_v_2, mlt_av, mlt_va
    procedure, pass(x) :: scal     => c_vect_scal
    procedure, pass(x) :: absval1  => c_vect_absval1
    procedure, pass(x) :: absval2  => c_vect_absval2
    generic, public    :: absval   => absval1, absval2
    procedure, pass(x) :: nrm2     => c_vect_nrm2
    procedure, pass(x) :: amax     => c_vect_amax
    procedure, pass(x) :: asum     => c_vect_asum                  
  end type psb_c_vect_type

  public  :: psb_c_vect
  private :: constructor, size_const
  interface psb_c_vect
    module procedure constructor, size_const
  end interface psb_c_vect

  private :: c_vect_get_nrows, c_vect_sizeof, c_vect_get_fmt, &
       & c_vect_all, c_vect_reall, c_vect_zero,  c_vect_asb, &
       & c_vect_gthab, c_vect_gthzv, c_vect_sctb, &
       & c_vect_free, c_vect_ins_a, c_vect_ins_v, c_vect_bld_x, &
       & c_vect_bld_mn, c_vect_bld_en, c_vect_get_vect, &
       & c_vect_cnv, c_vect_set_scal, &
       & c_vect_set_vect, c_vect_clone, c_vect_sync, c_vect_is_host, &
       & c_vect_is_dev, c_vect_is_sync, c_vect_set_host, &
       & c_vect_set_dev, c_vect_set_sync

  private ::  c_vect_dot_v, c_vect_dot_a, c_vect_axpby_v, c_vect_axpby_a, &
       & c_vect_mlt_v, c_vect_mlt_a, c_vect_mlt_a_2, c_vect_mlt_v_2, &
       & c_vect_mlt_va, c_vect_mlt_av, c_vect_scal, c_vect_absval1, &
       & c_vect_absval2, c_vect_nrm2, c_vect_amax, c_vect_asum                  



  class(psb_c_base_vect_type), allocatable, target,&
       & save, private :: psb_c_base_vect_default

  interface psb_set_vect_default
    module procedure psb_c_set_vect_default
  end interface psb_set_vect_default

  interface psb_get_vect_default
    module procedure psb_c_get_vect_default
  end interface psb_get_vect_default


contains


  subroutine  psb_c_set_vect_default(v) 
    implicit none 
    class(psb_c_base_vect_type), intent(in) :: v

    if (allocated(psb_c_base_vect_default)) then 
      deallocate(psb_c_base_vect_default)
    end if
    allocate(psb_c_base_vect_default, mold=v)

  end subroutine psb_c_set_vect_default

  function psb_c_get_vect_default(v) result(res)
    implicit none 
    class(psb_c_vect_type), intent(in) :: v
    class(psb_c_base_vect_type), pointer :: res

    res => psb_c_get_base_vect_default()

  end function psb_c_get_vect_default


  function psb_c_get_base_vect_default() result(res)
    implicit none 
    class(psb_c_base_vect_type), pointer :: res

    if (.not.allocated(psb_c_base_vect_default)) then 
      allocate(psb_c_base_vect_type :: psb_c_base_vect_default)
    end if

    res => psb_c_base_vect_default

  end function psb_c_get_base_vect_default


  subroutine c_vect_clone(x,y,info)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    call y%free(info)
    if ((info==0).and.allocated(x%v)) then 
      call y%bld(x%get_vect(),mold=x%v)
    end if
  end subroutine c_vect_clone

  subroutine c_vect_bld_x(x,invect,mold)
    complex(psb_spk_), intent(in)          :: invect(:)
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_base_vect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info
    class(psb_c_base_vect_type), pointer :: mld

    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then 
#ifdef HAVE_MOLD
      allocate(x%v,stat=info,mold=mold)
#else
      call mold%mold(x%v,info)
#endif
    else
#ifdef HAVE_MOLD
      allocate(x%v,stat=info, mold=psb_c_get_base_vect_default())
#else 
      mld = psb_c_get_base_vect_default()
      call mld%mold(x%v,info)
#endif
    endif

    if (info == psb_success_) call x%v%bld(invect)

  end subroutine c_vect_bld_x


  subroutine c_vect_bld_mn(x,n,mold)
    integer(psb_mpk_), intent(in) :: n
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_base_vect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info
    class(psb_c_base_vect_type), pointer :: mld


    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then 
#ifdef HAVE_MOLD
      allocate(x%v,stat=info,mold=mold)
#else
      call mold%mold(x%v,info)
#endif
    else
#ifdef HAVE_MOLD
      allocate(x%v,stat=info, mold=psb_c_get_base_vect_default())
#else 
      mld = psb_c_get_base_vect_default()
      call mld%mold(x%v,info)
#endif
    endif
    if (info == psb_success_) call x%v%bld(n)

  end subroutine c_vect_bld_mn


  subroutine c_vect_bld_en(x,n,mold)
    integer(psb_epk_), intent(in) :: n
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_base_vect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info
    class(psb_c_base_vect_type), pointer :: mld


    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then 
#ifdef HAVE_MOLD
      allocate(x%v,stat=info,mold=mold)
#else
      call mold%mold(x%v,info)
#endif
    else
#ifdef HAVE_MOLD
      allocate(x%v,stat=info, mold=psb_c_get_base_vect_default())
#else 
      mld = psb_c_get_base_vect_default()
      call mld%mold(x%v,info)
#endif
    endif
    if (info == psb_success_) call x%v%bld(n)

  end subroutine c_vect_bld_en

  function  c_vect_get_vect(x) result(res)
    class(psb_c_vect_type), intent(inout)  :: x
    complex(psb_spk_), allocatable                 :: res(:)
    integer(psb_ipk_) :: info

    if (allocated(x%v)) then
      res = x%v%get_vect()
    end if
  end function c_vect_get_vect

  subroutine c_vect_set_scal(x,val,first,last)
    class(psb_c_vect_type), intent(inout)  :: x
    complex(psb_spk_), intent(in) :: val
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val,first,last)

  end subroutine c_vect_set_scal

  subroutine c_vect_set_vect(x,val,first,last)
    class(psb_c_vect_type), intent(inout) :: x
    complex(psb_spk_), intent(in)         :: val(:)
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val,first,last)

  end subroutine c_vect_set_vect


  function constructor(x) result(this)
    complex(psb_spk_)   :: x(:)
    type(psb_c_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(x)
    call this%asb(size(x,kind=psb_ipk_),info)

  end function constructor


  function size_const(n) result(this)
    integer(psb_ipk_), intent(in) :: n
    type(psb_c_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(n)
    call this%asb(n,info)

  end function size_const

  function c_vect_get_nrows(x) result(res)
    implicit none 
    class(psb_c_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_nrows()
  end function c_vect_get_nrows

  function c_vect_sizeof(x) result(res)
    implicit none 
    class(psb_c_vect_type), intent(in) :: x
    integer(psb_epk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%sizeof()
  end function c_vect_sizeof

  function c_vect_get_fmt(x) result(res)
    implicit none 
    class(psb_c_vect_type), intent(in) :: x
    character(len=5) :: res
    res = 'NULL'
    if (allocated(x%v)) res = x%v%get_fmt()
  end function c_vect_get_fmt
  subroutine c_vect_all(n, x, info, mold)

    implicit none 
    integer(psb_ipk_), intent(in)           :: n
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_base_vect_type), intent(in), optional :: mold
    integer(psb_ipk_), intent(out)      :: info

    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then 
#ifdef HAVE_MOLD
      allocate(x%v,stat=info,mold=mold)
#else
      call mold%mold(x%v,info)
#endif
    else
      allocate(psb_c_base_vect_type :: x%v,stat=info)
    endif
    if (info == 0) then 
      call x%v%all(n,info)
    else
      info = psb_err_alloc_dealloc_
    end if

  end subroutine c_vect_all

  subroutine c_vect_reall(n, x, info)

    implicit none 
    integer(psb_ipk_), intent(in)         :: n
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    info = 0 
    if (.not.allocated(x%v)) &
         & call x%all(n,info)
    if (info == 0) &
         & call x%asb(n,info)

  end subroutine c_vect_reall

  subroutine c_vect_zero(x)
    use psi_serial_mod
    implicit none 
    class(psb_c_vect_type), intent(inout)    :: x

    if (allocated(x%v)) call x%v%zero()

  end subroutine c_vect_zero

  subroutine c_vect_asb(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer(psb_ipk_), intent(in)              :: n
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%v)) &
         & call x%v%asb(n,info)

  end subroutine c_vect_asb

  subroutine c_vect_gthab(n,idx,alpha,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    complex(psb_spk_) :: alpha, beta, y(:)
    class(psb_c_vect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,alpha,beta,y)

  end subroutine c_vect_gthab

  subroutine c_vect_gthzv(n,idx,x,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    complex(psb_spk_) ::  y(:)
    class(psb_c_vect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,y)

  end subroutine c_vect_gthzv

  subroutine c_vect_sctb(n,idx,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    complex(psb_spk_) :: beta, x(:)
    class(psb_c_vect_type) :: y

    if (allocated(y%v)) &
         &  call y%v%sct(n,idx,x,beta)

  end subroutine c_vect_sctb

  subroutine c_vect_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    class(psb_c_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (allocated(x%v)) then 
      call x%v%free(info)
      if (info == 0) deallocate(x%v,stat=info)
    end if

  end subroutine c_vect_free

  subroutine c_vect_ins_a(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none 
    class(psb_c_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    integer(psb_ipk_), intent(in)               :: irl(:)
    complex(psb_spk_), intent(in)        :: val(:)
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i

    info = 0
    if (.not.allocated(x%v)) then 
      info = psb_err_invalid_vect_state_
      return
    end if

    call  x%v%ins(n,irl,val,dupl,info)

  end subroutine c_vect_ins_a

  subroutine c_vect_ins_v(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none 
    class(psb_c_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    class(psb_i_vect_type), intent(inout)       :: irl
    class(psb_c_vect_type), intent(inout)       :: val
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i

    info = 0
    if (.not.(allocated(x%v).and.allocated(irl%v).and.allocated(val%v))) then 
      info = psb_err_invalid_vect_state_
      return
    end if

    call  x%v%ins(n,irl%v,val%v,dupl,info)

  end subroutine c_vect_ins_v


  subroutine c_vect_cnv(x,mold)
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_base_vect_type), intent(in), optional :: mold
    class(psb_c_base_vect_type), allocatable :: tmp
    integer(psb_ipk_) :: info

    if (present(mold)) then 
#ifdef HAVE_MOLD
      allocate(tmp,stat=info,mold=mold)
#else
      call mold%mold(tmp,info)
#endif
      if (allocated(x%v)) then 
        call x%v%sync()
        if (info == psb_success_) call tmp%bld(x%v%v)
        call x%v%free(info)
      end if
      call move_alloc(tmp,x%v)
    end if
  end subroutine c_vect_cnv


  subroutine c_vect_sync(x)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%sync()

  end subroutine c_vect_sync

  subroutine c_vect_set_sync(x)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%set_sync()

  end subroutine c_vect_set_sync

  subroutine c_vect_set_host(x)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%set_host()

  end subroutine c_vect_set_host

  subroutine c_vect_set_dev(x)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%set_dev()

  end subroutine c_vect_set_dev

  function c_vect_is_sync(x) result(res)
    implicit none 
    logical :: res
    class(psb_c_vect_type), intent(inout) :: x

    res = .true.
    if (allocated(x%v)) &
         & res = x%v%is_sync()

  end function c_vect_is_sync

  function c_vect_is_host(x) result(res)
    implicit none 
    logical :: res
    class(psb_c_vect_type), intent(inout) :: x

    res = .true.
    if (allocated(x%v)) &
         & res = x%v%is_host()

  end function c_vect_is_host

  function c_vect_is_dev(x) result(res)
    implicit none 
    logical :: res
    class(psb_c_vect_type), intent(inout) :: x

    res = .false. 
    if (allocated(x%v)) &
         & res =  x%v%is_dev()

  end function c_vect_is_dev


  function c_vect_dot_v(n,x,y) result(res)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(in)           :: n
    complex(psb_spk_)                :: res

    res = czero
    if (allocated(x%v).and.allocated(y%v)) &
         & res = x%v%dot(n,y%v)

  end function c_vect_dot_v

  function c_vect_dot_a(n,x,y) result(res)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x
    complex(psb_spk_), intent(in)    :: y(:)
    integer(psb_ipk_), intent(in)           :: n
    complex(psb_spk_)                :: res

    res = czero
    if (allocated(x%v)) &
         & res = x%v%dot(n,y)

  end function c_vect_dot_a

  subroutine c_vect_axpby_v(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none 
    integer(psb_ipk_), intent(in)               :: m
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    complex(psb_spk_), intent (in)       :: alpha, beta
    integer(psb_ipk_), intent(out)              :: info

    if (allocated(x%v).and.allocated(y%v)) then 
      call y%v%axpby(m,alpha,x%v,beta,info)
    else
      info = psb_err_invalid_vect_state_
    end if

  end subroutine c_vect_axpby_v

  subroutine c_vect_axpby_a(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none 
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)        :: x(:)
    class(psb_c_vect_type), intent(inout)  :: y
    complex(psb_spk_), intent (in)       :: alpha, beta
    integer(psb_ipk_), intent(out)              :: info

    if (allocated(y%v)) &
         & call y%v%axpby(m,alpha,x,beta,info)

  end subroutine c_vect_axpby_a


  subroutine c_vect_mlt_v(x, y, info)
    use psi_serial_mod
    implicit none 
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info    
    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(x%v).and.allocated(y%v)) &
         & call y%v%mlt(x%v,info)

  end subroutine c_vect_mlt_v

  subroutine c_vect_mlt_a(x, y, info)
    use psi_serial_mod
    implicit none 
    complex(psb_spk_), intent(in)        :: x(:)
    class(psb_c_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n


    info = 0
    if (allocated(y%v)) &
         & call y%v%mlt(x,info)

  end subroutine c_vect_mlt_a


  subroutine c_vect_mlt_a_2(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    complex(psb_spk_), intent(in)         :: alpha,beta
    complex(psb_spk_), intent(in)         :: y(:)
    complex(psb_spk_), intent(in)         :: x(:)
    class(psb_c_vect_type), intent(inout) :: z
    integer(psb_ipk_), intent(out)                  :: info
    integer(psb_ipk_) :: i, n

    info = 0    
    if (allocated(z%v)) &
         & call z%v%mlt(alpha,x,y,beta,info)

  end subroutine c_vect_mlt_a_2

  subroutine c_vect_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
    use psi_serial_mod
    implicit none 
    complex(psb_spk_), intent(in)          :: alpha,beta
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    class(psb_c_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)                   :: info    
    character(len=1), intent(in), optional :: conjgx, conjgy

    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(x%v).and.allocated(y%v).and.&
         & allocated(z%v)) &
         & call z%v%mlt(alpha,x%v,y%v,beta,info,conjgx,conjgy)

  end subroutine c_vect_mlt_v_2

  subroutine c_vect_mlt_av(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    complex(psb_spk_), intent(in)        :: alpha,beta
    complex(psb_spk_), intent(in)        :: x(:)
    class(psb_c_vect_type), intent(inout)  :: y
    class(psb_c_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info    
    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(z%v).and.allocated(y%v)) &
         & call z%v%mlt(alpha,x,y%v,beta,info)

  end subroutine c_vect_mlt_av

  subroutine c_vect_mlt_va(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    complex(psb_spk_), intent(in)        :: alpha,beta
    complex(psb_spk_), intent(in)        :: y(:)
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info    
    integer(psb_ipk_) :: i, n

    info = 0

    if (allocated(z%v).and.allocated(x%v)) &
         & call z%v%mlt(alpha,x%v,y,beta,info)

  end subroutine c_vect_mlt_va

  subroutine c_vect_scal(alpha, x)
    use psi_serial_mod
    implicit none 
    class(psb_c_vect_type), intent(inout)  :: x
    complex(psb_spk_), intent (in)       :: alpha

    if (allocated(x%v)) call x%v%scal(alpha)

  end subroutine c_vect_scal

  subroutine c_vect_absval1(x)
    class(psb_c_vect_type), intent(inout)  :: x

    if (allocated(x%v)) &
         &  call x%v%absval()

  end subroutine c_vect_absval1

  subroutine c_vect_absval2(x,y)
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y

    if (allocated(x%v)) then 
      if (.not.allocated(y%v))  call y%bld(psb_size(x%v%v))
      call x%v%absval(y%v)
    end if
  end subroutine c_vect_absval2

  function c_vect_nrm2(n,x) result(res)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)                :: res

    if (allocated(x%v)) then 
      res = x%v%nrm2(n)
    else
      res = szero
    end if

  end function c_vect_nrm2

  function c_vect_amax(n,x) result(res)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)                :: res

    if (allocated(x%v)) then 
      res = x%v%amax(n)
    else
      res = szero
    end if

  end function c_vect_amax

  function c_vect_asum(n,x) result(res)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)                :: res

    if (allocated(x%v)) then 
      res = x%v%asum(n)
    else
      res = szero
    end if

  end function c_vect_asum


end module psb_c_vect_mod



module psb_c_multivect_mod

  use psb_c_base_multivect_mod
  use psb_const_mod
  use psb_i_vect_mod


  !private

  type psb_c_multivect_type
    class(psb_c_base_multivect_type), allocatable :: v 
  contains
    procedure, pass(x) :: get_nrows => c_vect_get_nrows
    procedure, pass(x) :: get_ncols => c_vect_get_ncols
    procedure, pass(x) :: sizeof   => c_vect_sizeof
    procedure, pass(x) :: get_fmt  => c_vect_get_fmt

    procedure, pass(x) :: all      => c_vect_all
    procedure, pass(x) :: reall    => c_vect_reall
    procedure, pass(x) :: zero     => c_vect_zero
    procedure, pass(x) :: asb      => c_vect_asb
    procedure, pass(x) :: sync     => c_vect_sync
    procedure, pass(x) :: free     => c_vect_free
    procedure, pass(x) :: ins      => c_vect_ins
    procedure, pass(x) :: bld_x    => c_vect_bld_x
    procedure, pass(x) :: bld_n    => c_vect_bld_n
    generic, public    :: bld      => bld_x, bld_n
    procedure, pass(x) :: get_vect => c_vect_get_vect
    procedure, pass(x) :: cnv      => c_vect_cnv
    procedure, pass(x) :: set_scal => c_vect_set_scal
    procedure, pass(x) :: set_vect => c_vect_set_vect
    generic, public    :: set      => set_vect, set_scal
    procedure, pass(x) :: clone    => c_vect_clone
    procedure, pass(x) :: gthab    => c_vect_gthab
    procedure, pass(x) :: gthzv    => c_vect_gthzv
    procedure, pass(x) :: gthzv_x  => c_vect_gthzv_x
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => c_vect_sctb
    procedure, pass(y) :: sctb_x   => c_vect_sctb_x
    generic, public    :: sct      => sctb, sctb_x
!!$    procedure, pass(x) :: dot_v    => c_vect_dot_v
!!$    procedure, pass(x) :: dot_a    => c_vect_dot_a
!!$    generic, public    :: dot      => dot_v, dot_a
!!$    procedure, pass(y) :: axpby_v  => c_vect_axpby_v
!!$    procedure, pass(y) :: axpby_a  => c_vect_axpby_a
!!$    generic, public    :: axpby    => axpby_v, axpby_a
!!$    procedure, pass(y) :: mlt_v    => c_vect_mlt_v
!!$    procedure, pass(y) :: mlt_a    => c_vect_mlt_a
!!$    procedure, pass(z) :: mlt_a_2  => c_vect_mlt_a_2
!!$    procedure, pass(z) :: mlt_v_2  => c_vect_mlt_v_2
!!$    procedure, pass(z) :: mlt_va   => c_vect_mlt_va
!!$    procedure, pass(z) :: mlt_av   => c_vect_mlt_av
!!$    generic, public    :: mlt      => mlt_v, mlt_a, mlt_a_2,&
!!$         & mlt_v_2, mlt_av, mlt_va
!!$    procedure, pass(x) :: scal     => c_vect_scal
!!$    procedure, pass(x) :: nrm2     => c_vect_nrm2
!!$    procedure, pass(x) :: amax     => c_vect_amax
!!$    procedure, pass(x) :: asum     => c_vect_asum
  end type psb_c_multivect_type

  public  :: psb_c_multivect, psb_c_multivect_type,&
       & psb_set_multivect_default, psb_get_multivect_default, &
       & psb_c_base_multivect_type

  private
  interface psb_c_multivect
    module procedure constructor, size_const
  end interface psb_c_multivect

  class(psb_c_base_multivect_type), allocatable, target,&
       & save, private :: psb_c_base_multivect_default

  interface psb_set_multivect_default
    module procedure psb_c_set_multivect_default
  end interface psb_set_multivect_default

  interface psb_get_vect_default
    module procedure psb_c_get_multivect_default
  end interface psb_get_vect_default


contains


  subroutine  psb_c_set_multivect_default(v) 
    implicit none 
    class(psb_c_base_multivect_type), intent(in) :: v

    if (allocated(psb_c_base_multivect_default)) then 
      deallocate(psb_c_base_multivect_default)
    end if
    allocate(psb_c_base_multivect_default, mold=v)

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

    if (.not.allocated(psb_c_base_multivect_default)) then 
      allocate(psb_c_base_multivect_type :: psb_c_base_multivect_default)
    end if

    res => psb_c_base_multivect_default

  end function psb_c_get_base_multivect_default


  subroutine c_vect_clone(x,y,info)
    implicit none 
    class(psb_c_multivect_type), intent(inout) :: x
    class(psb_c_multivect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    call y%free(info)
    if ((info==0).and.allocated(x%v)) then 
      call y%bld(x%get_vect(),mold=x%v)
    end if
  end subroutine c_vect_clone

  subroutine c_vect_bld_x(x,invect,mold)
    complex(psb_spk_), intent(in)          :: invect(:,:)
    class(psb_c_multivect_type), intent(out) :: x
    class(psb_c_base_multivect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info
    class(psb_c_base_multivect_type), pointer :: mld

    if (present(mold)) then 
#ifdef HAVE_MOLD
      allocate(x%v,stat=info,mold=mold)
#else
      call mold%mold(x%v,info)
#endif
    else
#ifdef HAVE_MOLD
      allocate(x%v,stat=info, mold=psb_c_get_base_multivect_default())
#else 
      mld = psb_c_get_base_multivect_default()
      call mld%mold(x%v,info)
#endif
    endif

    if (info == psb_success_) call x%v%bld(invect)

  end subroutine c_vect_bld_x


  subroutine c_vect_bld_n(x,m,n,mold)
    integer(psb_ipk_), intent(in) :: m,n
    class(psb_c_multivect_type), intent(out) :: x
    class(psb_c_base_multivect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info
    class(psb_c_base_multivect_type), pointer :: mld

    if (present(mold)) then 
#ifdef HAVE_MOLD
      allocate(x%v,stat=info,mold=mold)
#else
      call mold%mold(x%v,info)
#endif
    else
#ifdef HAVE_MOLD
      allocate(x%v,stat=info, mold=psb_c_get_base_multivect_default())
#else 
      mld = psb_c_get_base_multivect_default()
      call mld%mold(x%v,info)
#endif
    endif
    if (info == psb_success_) call x%v%bld(m,n)

  end subroutine c_vect_bld_n

  function  c_vect_get_vect(x) result(res)
    class(psb_c_multivect_type), intent(inout)  :: x
    complex(psb_spk_), allocatable                 :: res(:,:)
    integer(psb_ipk_) :: info

    if (allocated(x%v)) then
      res = x%v%get_vect()
    end if
  end function c_vect_get_vect

  subroutine c_vect_set_scal(x,val)
    class(psb_c_multivect_type), intent(inout)  :: x
    complex(psb_spk_), intent(in) :: val

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val)

  end subroutine c_vect_set_scal

  subroutine c_vect_set_vect(x,val)
    class(psb_c_multivect_type), intent(inout) :: x
    complex(psb_spk_), intent(in)         :: val(:,:)

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val)

  end subroutine c_vect_set_vect


  function constructor(x) result(this)
    complex(psb_spk_)   :: x(:,:)
    type(psb_c_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(x)
    call this%asb(size(x,dim=1,kind=psb_ipk_),size(x,dim=2,kind=psb_ipk_),info)

  end function constructor


  function size_const(m,n) result(this)
    integer(psb_ipk_), intent(in) :: m,n
    type(psb_c_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(m,n)
    call this%asb(m,n,info)

  end function size_const

  function c_vect_get_nrows(x) result(res)
    implicit none 
    class(psb_c_multivect_type), intent(in) :: x
    integer(psb_ipk_)  :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_nrows()
  end function c_vect_get_nrows

  function c_vect_get_ncols(x) result(res)
    implicit none 
    class(psb_c_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_ncols()
  end function c_vect_get_ncols

  function c_vect_sizeof(x) result(res)
    implicit none 
    class(psb_c_multivect_type), intent(in) :: x
    integer(psb_epk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%sizeof()
  end function c_vect_sizeof

  function c_vect_get_fmt(x) result(res)
    implicit none 
    class(psb_c_multivect_type), intent(in) :: x
    character(len=5) :: res
    res = 'NULL'
    if (allocated(x%v)) res = x%v%get_fmt()
  end function c_vect_get_fmt

  subroutine c_vect_all(m,n, x, info, mold)

    implicit none 
    integer(psb_ipk_), intent(in)       :: m,n
    class(psb_c_multivect_type), intent(out) :: x
    class(psb_c_base_multivect_type), intent(in), optional :: mold
    integer(psb_ipk_), intent(out)      :: info

    if (present(mold)) then 
#ifdef HAVE_MOLD
      allocate(x%v,stat=info,mold=mold)
#else
      call mold%mold(x%v,info)
#endif
    else
      allocate(psb_c_base_multivect_type :: x%v,stat=info)
    endif
    if (info == 0) then 
      call x%v%all(m,n,info)
    else
      info = psb_err_alloc_dealloc_
    end if

  end subroutine c_vect_all

  subroutine c_vect_reall(m,n, x, info)

    implicit none 
    integer(psb_ipk_), intent(in)         :: m,n
    class(psb_c_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    info = 0 
    if (.not.allocated(x%v)) &
         & call x%all(m,n,info)
    if (info == 0) &
         & call x%asb(m,n,info)

  end subroutine c_vect_reall

  subroutine c_vect_zero(x)
    use psi_serial_mod
    implicit none 
    class(psb_c_multivect_type), intent(inout)    :: x

    if (allocated(x%v)) call x%v%zero()

  end subroutine c_vect_zero

  subroutine c_vect_asb(m,n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer(psb_ipk_), intent(in)              :: m,n
    class(psb_c_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%v)) &
         & call x%v%asb(m,n,info)

  end subroutine c_vect_asb

  subroutine c_vect_sync(x)
    implicit none 
    class(psb_c_multivect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%sync()

  end subroutine c_vect_sync

  subroutine c_vect_gthab(n,idx,alpha,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    complex(psb_spk_) :: alpha, beta, y(:)
    class(psb_c_multivect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,alpha,beta,y)

  end subroutine c_vect_gthab

  subroutine c_vect_gthzv(n,idx,x,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    complex(psb_spk_) ::  y(:)
    class(psb_c_multivect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,y)

  end subroutine c_vect_gthzv

  subroutine c_vect_gthzv_x(i,n,idx,x,y)
    use psi_serial_mod
    integer(psb_ipk_) :: i,n
    class(psb_i_base_vect_type) :: idx
    complex(psb_spk_) ::  y(:)
    class(psb_c_multivect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(i,n,idx,y)

  end subroutine c_vect_gthzv_x

  subroutine c_vect_sctb(n,idx,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    complex(psb_spk_) :: beta, x(:)
    class(psb_c_multivect_type) :: y

    if (allocated(y%v)) &
         &  call y%v%sct(n,idx,x,beta)

  end subroutine c_vect_sctb

  subroutine c_vect_sctb_x(i,n,idx,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    complex(psb_spk_) :: beta, x(:)
    class(psb_c_multivect_type) :: y

    if (allocated(y%v)) &
         &  call y%v%sct(i,n,idx,x,beta)

  end subroutine c_vect_sctb_x

  subroutine c_vect_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    class(psb_c_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (allocated(x%v)) then 
      call x%v%free(info)
      if (info == 0) deallocate(x%v,stat=info)
    end if

  end subroutine c_vect_free

  subroutine c_vect_ins(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none 
    class(psb_c_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    integer(psb_ipk_), intent(in)               :: irl(:)
    complex(psb_spk_), intent(in)        :: val(:,:)
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i

    info = 0
    if (.not.allocated(x%v)) then 
      info = psb_err_invalid_vect_state_
      return
    end if

    call  x%v%ins(n,irl,val,dupl,info)

  end subroutine c_vect_ins


  subroutine c_vect_cnv(x,mold)
    class(psb_c_multivect_type), intent(inout) :: x
    class(psb_c_base_multivect_type), intent(in), optional :: mold
    class(psb_c_base_multivect_type), allocatable :: tmp
    integer(psb_ipk_) :: info

    if (present(mold)) then 
#ifdef HAVE_MOLD
      allocate(tmp,stat=info,mold=mold)
#else
      call mold%mold(tmp,info)
#endif
      if (allocated(x%v)) then 
        call x%v%sync()
        if (info == psb_success_) call tmp%bld(x%v%v)
        call x%v%free(info)
      end if
      call move_alloc(tmp,x%v)
    end if
  end subroutine c_vect_cnv


!!$  function c_vect_dot_v(n,x,y) result(res)
!!$    implicit none 
!!$    class(psb_c_multivect_type), intent(inout) :: x, y
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    complex(psb_spk_)                :: res
!!$
!!$    res = czero
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & res = x%v%dot(n,y%v)
!!$
!!$  end function c_vect_dot_v
!!$
!!$  function c_vect_dot_a(n,x,y) result(res)
!!$    implicit none 
!!$    class(psb_c_multivect_type), intent(inout) :: x
!!$    complex(psb_spk_), intent(in)    :: y(:)
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    complex(psb_spk_)                :: res
!!$    
!!$    res = czero
!!$    if (allocated(x%v)) &
!!$         & res = x%v%dot(n,y)
!!$    
!!$  end function c_vect_dot_a
!!$    
!!$  subroutine c_vect_axpby_v(m,alpha, x, beta, y, info)
!!$    use psi_serial_mod
!!$    implicit none 
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    class(psb_c_multivect_type), intent(inout)  :: x
!!$    class(psb_c_multivect_type), intent(inout)  :: y
!!$    complex(psb_spk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    
!!$    if (allocated(x%v).and.allocated(y%v)) then 
!!$      call y%v%axpby(m,alpha,x%v,beta,info)
!!$    else
!!$      info = psb_err_invalid_vect_state_
!!$    end if
!!$
!!$  end subroutine c_vect_axpby_v
!!$
!!$  subroutine c_vect_axpby_a(m,alpha, x, beta, y, info)
!!$    use psi_serial_mod
!!$    implicit none 
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    complex(psb_spk_), intent(in)        :: x(:)
!!$    class(psb_c_multivect_type), intent(inout)  :: y
!!$    complex(psb_spk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    
!!$    if (allocated(y%v)) &
!!$         & call y%v%axpby(m,alpha,x,beta,info)
!!$    
!!$  end subroutine c_vect_axpby_a
!!$
!!$    
!!$  subroutine c_vect_mlt_v(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none 
!!$    class(psb_c_multivect_type), intent(inout)  :: x
!!$    class(psb_c_multivect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info    
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & call y%v%mlt(x%v,info)
!!$
!!$  end subroutine c_vect_mlt_v
!!$
!!$  subroutine c_vect_mlt_a(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none 
!!$    complex(psb_spk_), intent(in)        :: x(:)
!!$    class(psb_c_multivect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$
!!$    info = 0
!!$    if (allocated(y%v)) &
!!$         & call y%v%mlt(x,info)
!!$    
!!$  end subroutine c_vect_mlt_a
!!$
!!$
!!$  subroutine c_vect_mlt_a_2(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none 
!!$    complex(psb_spk_), intent(in)         :: alpha,beta
!!$    complex(psb_spk_), intent(in)         :: y(:)
!!$    complex(psb_spk_), intent(in)         :: x(:)
!!$    class(psb_c_multivect_type), intent(inout) :: z
!!$    integer(psb_ipk_), intent(out)                  :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0    
!!$    if (allocated(z%v)) &
!!$         & call z%v%mlt(alpha,x,y,beta,info)
!!$    
!!$  end subroutine c_vect_mlt_a_2
!!$
!!$  subroutine c_vect_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
!!$    use psi_serial_mod
!!$    implicit none 
!!$    complex(psb_spk_), intent(in)          :: alpha,beta
!!$    class(psb_c_multivect_type), intent(inout)  :: x
!!$    class(psb_c_multivect_type), intent(inout)  :: y
!!$    class(psb_c_multivect_type), intent(inout)  :: z
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
!!$  end subroutine c_vect_mlt_v_2
!!$
!!$  subroutine c_vect_mlt_av(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none 
!!$    complex(psb_spk_), intent(in)        :: alpha,beta
!!$    complex(psb_spk_), intent(in)        :: x(:)
!!$    class(psb_c_multivect_type), intent(inout)  :: y
!!$    class(psb_c_multivect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info    
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(z%v).and.allocated(y%v)) &
!!$         & call z%v%mlt(alpha,x,y%v,beta,info)
!!$
!!$  end subroutine c_vect_mlt_av
!!$
!!$  subroutine c_vect_mlt_va(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none 
!!$    complex(psb_spk_), intent(in)        :: alpha,beta
!!$    complex(psb_spk_), intent(in)        :: y(:)
!!$    class(psb_c_multivect_type), intent(inout)  :: x
!!$    class(psb_c_multivect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info    
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    
!!$    if (allocated(z%v).and.allocated(x%v)) &
!!$         & call z%v%mlt(alpha,x%v,y,beta,info)
!!$
!!$  end subroutine c_vect_mlt_va
!!$
!!$  subroutine c_vect_scal(alpha, x)
!!$    use psi_serial_mod
!!$    implicit none 
!!$    class(psb_c_multivect_type), intent(inout)  :: x
!!$    complex(psb_spk_), intent (in)       :: alpha
!!$    
!!$    if (allocated(x%v)) call x%v%scal(alpha)
!!$
!!$  end subroutine c_vect_scal
!!$
!!$
!!$  function c_vect_nrm2(n,x) result(res)
!!$    implicit none 
!!$    class(psb_c_multivect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_spk_)                :: res
!!$    
!!$    if (allocated(x%v)) then 
!!$      res = x%v%nrm2(n)
!!$    else
!!$      res = szero
!!$    end if
!!$
!!$  end function c_vect_nrm2
!!$  
!!$  function c_vect_amax(n,x) result(res)
!!$    implicit none 
!!$    class(psb_c_multivect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_spk_)                :: res
!!$
!!$    if (allocated(x%v)) then 
!!$      res = x%v%amax(n)
!!$    else
!!$      res = szero
!!$    end if
!!$
!!$  end function c_vect_amax
!!$
!!$  function c_vect_asum(n,x) result(res)
!!$    implicit none 
!!$    class(psb_c_multivect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_spk_)                :: res
!!$
!!$    if (allocated(x%v)) then 
!!$      res = x%v%asum(n)
!!$    else
!!$      res = szero
!!$    end if
!!$
!!$  end function c_vect_asum

end module psb_c_multivect_mod
