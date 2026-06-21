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
! package: psb_c_vect_mod
!
! This module contains the definition of the psb_c_vect type which
! is the outer container for dense vectors.
! Therefore all methods simply invoke the corresponding methods of the
! inner component.
!
submodule (psb_c_vect_mod) psb_c_vect_impl
  use psb_base_mod
  use psi_serial_mod
  implicit none 

contains

  module function c_vect_get_dupl(x) result(res)
    class(psb_c_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    if (allocated(x%v)) then 
      res = x%v%get_dupl()
    else
      res = psb_dupl_null_
    end if
  end function c_vect_get_dupl

  module subroutine c_vect_set_dupl(x,val)
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if (allocated(x%v)) then 
      if (present(val)) then
        call x%v%set_dupl(val)
      else
        call x%v%set_dupl(psb_dupl_def_)
      end if
    end if
  end subroutine c_vect_set_dupl

  module function c_vect_get_ncfs(x) result(res)
    class(psb_c_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    if (allocated(x%v)) then 
      res = x%v%get_ncfs()
    else
      res = 0
    end if
  end function c_vect_get_ncfs

  module subroutine c_vect_set_ncfs(x,val)
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if (allocated(x%v)) then 
      if (present(val)) then
        call x%v%set_ncfs(val)
      else
        call x%v%set_ncfs(0)
      end if
    end if
  end subroutine c_vect_set_ncfs

  module function c_vect_get_state(x) result(res)
    class(psb_c_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    if (allocated(x%v)) then 
      res = x%v%get_state()
    else
      res = psb_vect_null_
    end if
  end function c_vect_get_state

  module function c_vect_is_null(x) result(res)
    class(psb_c_vect_type), intent(in) :: x
    logical :: res
    res = (x%get_state() == psb_vect_null_)
  end function c_vect_is_null

  module function c_vect_is_bld(x) result(res)
    class(psb_c_vect_type), intent(in) :: x
    logical :: res
    res = (x%get_state() == psb_vect_bld_)
  end function c_vect_is_bld

  module function c_vect_is_upd(x) result(res)
    class(psb_c_vect_type), intent(in) :: x
    logical :: res
    res = (x%get_state() == psb_vect_upd_)
  end function c_vect_is_upd

  module function c_vect_is_asb(x) result(res)
    class(psb_c_vect_type), intent(in) :: x
    logical :: res
    res = (x%get_state() == psb_vect_asb_)
  end function c_vect_is_asb

  module subroutine  c_vect_set_state(n,x)
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n
    if (allocated(x%v)) then 
      call x%v%set_state(n)
    end if
  end subroutine c_vect_set_state

  module subroutine  c_vect_set_null(x)
    class(psb_c_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_null_)
  end subroutine c_vect_set_null

  module subroutine  c_vect_set_bld(x)
    class(psb_c_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_bld_)
  end subroutine c_vect_set_bld

  module subroutine  c_vect_set_upd(x)
    class(psb_c_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_upd_)
  end subroutine c_vect_set_upd

  module subroutine  c_vect_set_asb(x)
    class(psb_c_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_asb_)
  end subroutine c_vect_set_asb

  module function c_vect_get_nrmv(x) result(res)
    class(psb_c_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = x%nrmv
  end function c_vect_get_nrmv

  module subroutine c_vect_set_nrmv(x,val)
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: val

    x%nrmv = val
  end subroutine c_vect_set_nrmv

  module function c_vect_is_remote_build(x) result(res)
    class(psb_c_vect_type), intent(in) :: x
    logical :: res
    res = (x%remote_build == psb_matbld_remote_)
  end function c_vect_is_remote_build

  module subroutine c_vect_set_remote_build(x,val)
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if (present(val)) then
      x%remote_build = val
    else
      x%remote_build = psb_matbld_remote_
    end if
  end subroutine c_vect_set_remote_build

  module subroutine c_vect_clone(x,y,info)
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    call y%free(info)
    if ((info==0).and.allocated(x%v)) then
      !
      ! Using sourced allocation here creates
      ! problems with handling of memory allocated
      ! elsewhere (e.g. accelerators), hence delegation
      ! to %bld method
      ! 
      call y%bld(x%get_vect(),mold=x%v)
    end if
  end subroutine c_vect_clone

  module subroutine c_vect_bld_x(x,invect,mold,scratch)
    complex(psb_spk_), intent(in)          :: invect(:)
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_base_vect_type), intent(in), optional :: mold
    logical, intent(in), optional        :: scratch

    logical :: scratch_
    integer(psb_ipk_) :: info

    if (present(scratch)) then
      scratch_ = scratch
    else
      scratch_ = .false.
    end if

    info = psb_success_
    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_c_get_base_vect_default())
    endif

    if (info == psb_success_) call x%v%bld(invect,scratch=scratch_)

  end subroutine c_vect_bld_x

  module subroutine c_vect_bld_mn(x,n,mold,scratch)
    integer(psb_mpk_), intent(in) :: n
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_base_vect_type), intent(in), optional :: mold
    logical, intent(in), optional        :: scratch

    logical :: scratch_
    integer(psb_ipk_) :: info
    class(psb_c_base_vect_type), pointer :: mld
    if (present(scratch)) then
      scratch_ = scratch
    else
      scratch_ = .false.
    end if

    info = psb_success_
    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_c_get_base_vect_default())
    endif
    if (info == psb_success_) call x%v%bld(n,scratch=scratch_)

  end subroutine c_vect_bld_mn

  module subroutine c_vect_bld_en(x,n,mold,scratch)
    integer(psb_epk_), intent(in) :: n
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_base_vect_type), intent(in), optional :: mold
    logical, intent(in), optional        :: scratch

    logical :: scratch_
    integer(psb_ipk_) :: info

    info = psb_success_
    if (present(scratch)) then
      scratch_ = scratch
    else
      scratch_ = .false.
    end if

    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_c_get_base_vect_default())
    endif
    if (info == psb_success_) call x%v%bld(n,scratch=scratch_)

  end subroutine c_vect_bld_en

  module function  c_vect_get_vect(x,n) result(res)
    class(psb_c_vect_type), intent(inout)  :: x
    complex(psb_spk_), allocatable                 :: res(:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional :: n

    if (allocated(x%v)) then
      res = x%v%get_vect(n)
    end if
  end function c_vect_get_vect

  module subroutine c_vect_set_scal(x,val,first,last)
    class(psb_c_vect_type), intent(inout)  :: x
    complex(psb_spk_), intent(in) :: val
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val,first,last)

  end subroutine c_vect_set_scal

  module subroutine c_vect_set_vect(x,val,first,last)
    class(psb_c_vect_type), intent(inout) :: x
    complex(psb_spk_), intent(in)         :: val(:)
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val,first,last)

  end subroutine c_vect_set_vect

  module subroutine c_vect_check_addr(x)
    class(psb_c_vect_type), intent(inout) :: x

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%check_addr()

  end subroutine c_vect_check_addr

  module function c_vect_get_nrows(x) result(res)
    class(psb_c_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_nrows()
  end function c_vect_get_nrows

  module function c_vect_sizeof(x) result(res)
    class(psb_c_vect_type), intent(in) :: x
    integer(psb_epk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%sizeof()
  end function c_vect_sizeof

  module function c_vect_get_fmt(x) result(res)
    class(psb_c_vect_type), intent(in) :: x
    character(len=5) :: res
    res = 'NULL'
    if (allocated(x%v)) res = x%v%get_fmt()
  end function c_vect_get_fmt

  module subroutine c_vect_all(n, x, info, mold)

    integer(psb_ipk_), intent(in)           :: n
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)      :: info
    class(psb_c_base_vect_type), intent(in), optional :: mold

    if (allocated(x%v)) &
         & call x%free(info)

    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(psb_c_base_vect_type :: x%v,stat=info)
    endif
    if (info == 0) then
      call x%v%all(n,info)
    else
      info = psb_err_alloc_dealloc_
    end if
    call x%set_bld()
  end subroutine c_vect_all

  module subroutine c_vect_reinit(x, info, clear)
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)      :: info
    logical, intent(in), optional       :: clear

    if (allocated(x%v)) call x%v%reinit(info,clear)
    call x%set_upd()

  end subroutine c_vect_reinit
 
  module subroutine c_vect_reall(n, x, info)

    integer(psb_ipk_), intent(in)         :: n
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    info = 0
    if (.not.allocated(x%v)) &
         & call x%all(n,info)
    if (info == 0) &
         & call x%asb(n,info)

  end subroutine c_vect_reall

  module subroutine c_vect_zero(x)
    class(psb_c_vect_type), intent(inout)    :: x

    if (allocated(x%v)) call x%v%zero()

  end subroutine c_vect_zero

  module subroutine c_vect_asb(n, x, info, scratch)
    integer(psb_ipk_), intent(in)         :: n
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info
    logical, intent(in), optional        :: scratch

    if (allocated(x%v)) then
      call x%v%asb(n,info,scratch=scratch)
      call x%set_asb()
    end if
  end subroutine c_vect_asb

  module subroutine c_vect_gthab(n,idx,alpha,x,beta,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    complex(psb_spk_) :: alpha, beta, y(:)
    class(psb_c_vect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,alpha,beta,y)

  end subroutine c_vect_gthab

  module subroutine c_vect_gthzv(n,idx,x,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    complex(psb_spk_) ::  y(:)
    class(psb_c_vect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,y)

  end subroutine c_vect_gthzv

  module subroutine c_vect_sctb(n,idx,x,beta,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    complex(psb_spk_) :: beta, x(:)
    class(psb_c_vect_type) :: y

    if (allocated(y%v)) &
         &  call y%v%sct(n,idx,x,beta)

  end subroutine c_vect_sctb

  module subroutine c_vect_free(x, info)
    class(psb_c_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (allocated(x%v)) then
      call x%v%free(info)
      if (info == 0) deallocate(x%v,stat=info)
    end if

  end subroutine c_vect_free

  module subroutine c_vect_ins_a(n,irl,val,x,maxr,info)
    class(psb_c_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, maxr
    integer(psb_ipk_), intent(in)               :: irl(:)
    complex(psb_spk_), intent(in)        :: val(:)
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i, dupl

    info = 0
    if (.not.allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if
    dupl = x%get_dupl()
    call  x%v%ins(n,irl,val,dupl,maxr,info)

  end subroutine c_vect_ins_a

  module subroutine c_vect_ins_v(n,irl,val,x,maxr,info)
    class(psb_c_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, maxr
    class(psb_i_vect_type), intent(inout)       :: irl
    class(psb_c_vect_type), intent(inout)       :: val
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i, dupl

    info = 0
    if (.not.(allocated(x%v).and.allocated(irl%v).and.allocated(val%v))) then
      info = psb_err_invalid_vect_state_
      return
    end if
    dupl = x%get_dupl()
    call  x%v%ins(n,irl%v,val%v,dupl,maxr,info)

  end subroutine c_vect_ins_v

  module subroutine c_vect_cnv(x,mold)
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_base_vect_type), intent(in), optional :: mold
    class(psb_c_base_vect_type), allocatable :: tmp

    integer(psb_ipk_) :: info

    info = psb_success_
    if (present(mold)) then
      allocate(tmp,stat=info,mold=mold)
    else
      allocate(tmp,stat=info,mold=psb_c_get_base_vect_default())
    end if
    if (allocated(x%v)) then
      if (allocated(x%v%v)) then 
        call x%v%sync()
        if (info == psb_success_) call tmp%bld(x%v%v)
        call x%v%base_cpy(tmp)
        call x%v%free(info)
      endif
    end if
    call move_alloc(tmp,x%v)

  end subroutine c_vect_cnv

  module subroutine c_vect_sync(x)
    class(psb_c_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%sync()

  end subroutine c_vect_sync

  module subroutine c_vect_set_sync(x)
    class(psb_c_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%set_sync()

  end subroutine c_vect_set_sync

  module subroutine c_vect_set_host(x)
    class(psb_c_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%set_host()

  end subroutine c_vect_set_host

  module subroutine c_vect_set_dev(x)
    class(psb_c_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%set_dev()

  end subroutine c_vect_set_dev

  module function c_vect_is_sync(x) result(res)
    logical :: res
    class(psb_c_vect_type), intent(inout) :: x

    res = .true.
    if (allocated(x%v)) &
         & res = x%v%is_sync()

  end function c_vect_is_sync

  module function c_vect_is_host(x) result(res)
    logical :: res
    class(psb_c_vect_type), intent(inout) :: x

    res = .true.
    if (allocated(x%v)) &
         & res = x%v%is_host()

  end function c_vect_is_host

  module function c_vect_is_dev(x) result(res)
    logical :: res
    class(psb_c_vect_type), intent(inout) :: x

    res = .false.
    if (allocated(x%v)) &
         & res =  x%v%is_dev()

  end function c_vect_is_dev

  module function c_vect_get_entry(x,index) result(res)
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: index
    complex(psb_spk_) :: res
    res = czero
    if (allocated(x%v)) res = x%v%get_entry(index)
  end function c_vect_get_entry

  module subroutine c_vect_set_entry(x,index,val)
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: index
    complex(psb_spk_) :: val
    if (allocated(x%v)) call x%v%set_entry(index,val)
  end subroutine c_vect_set_entry

  module function c_vect_dot_v(n,x,y) result(res)
    class(psb_c_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(in)           :: n
    complex(psb_spk_)                :: res

    res = czero
    if (allocated(x%v).and.allocated(y%v)) &
         & res = x%v%dot(n,y%v)

  end function c_vect_dot_v

  module function c_vect_dot_a(n,x,y) result(res)
    class(psb_c_vect_type), intent(inout) :: x
    complex(psb_spk_), intent(in)    :: y(:)
    integer(psb_ipk_), intent(in)           :: n
    complex(psb_spk_)                :: res

    res = czero
    if (allocated(x%v)) &
         & res = x%v%dot_a(n,y)

  end function c_vect_dot_a

  module subroutine c_vect_axpby_v(m,alpha, x, beta, y, info)
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

  module subroutine c_vect_axpby_v2(m,alpha, x, beta, y, z, info)
    integer(psb_ipk_), intent(in)            :: m
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    class(psb_c_vect_type), intent(inout)  :: z
    complex(psb_spk_), intent (in)             :: alpha, beta
    integer(psb_ipk_), intent(out)           :: info

    if (allocated(x%v).and.allocated(y%v)) then
      call z%v%axpby(m,alpha,x%v,beta,y%v,info)
    else
      info = psb_err_invalid_vect_state_
    end if

  end subroutine c_vect_axpby_v2

  module subroutine c_vect_axpby_a(m,alpha, x, beta, y, info)
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)        :: x(:)
    class(psb_c_vect_type), intent(inout)  :: y
    complex(psb_spk_), intent (in)       :: alpha, beta
    integer(psb_ipk_), intent(out)              :: info

    if (allocated(y%v)) &
         & call y%v%axpby(m,alpha,x,beta,info)

  end subroutine c_vect_axpby_a

  module subroutine c_vect_axpby_a2(m,alpha, x, beta, y, z, info)
    integer(psb_ipk_), intent(in)               :: m
    complex(psb_spk_), intent(in)                 :: x(:)
    complex(psb_spk_), intent(in)                 :: y(:)
    class(psb_c_vect_type), intent(inout)     :: z
    complex(psb_spk_), intent (in)       :: alpha, beta
    integer(psb_ipk_), intent(out)              :: info

    if (allocated(z%v)) &
         & call z%v%axpby(m,alpha,x,beta,y,info)

  end subroutine c_vect_axpby_a2

  module subroutine c_vect_upd_xyz(m,alpha,beta,gamma,delta,x, y, z, info)
    integer(psb_ipk_), intent(in)            :: m
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    class(psb_c_vect_type), intent(inout)  :: z
    complex(psb_spk_), intent (in)     :: alpha, beta, gamma, delta
    integer(psb_ipk_), intent(out)   :: info

    if (allocated(z%v)) &
         call z%v%upd_xyz(m,alpha,beta,gamma,delta,x%v,y%v,info)
    
  end subroutine c_vect_upd_xyz

  module subroutine c_vect_xyzw(m,a,b,c,d,e,f,x, y, z, w, info)
    integer(psb_ipk_), intent(in)            :: m
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    class(psb_c_vect_type), intent(inout)  :: z
    class(psb_c_vect_type), intent(inout)  :: w
    complex(psb_spk_), intent (in)     :: a, b, c, d, e, f
    integer(psb_ipk_), intent(out)   :: info

    if (allocated(w%v)) &
         call w%v%xyzw(m,a,b,c,d,e,f,x%v,y%v,z%v,info)
    
  end subroutine c_vect_xyzw

  module subroutine c_vect_mlt_v(x, y, info)
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(x%v).and.allocated(y%v)) &
         & call y%v%mlt(x%v,info)

  end subroutine c_vect_mlt_v

  module subroutine c_vect_mlt_a(x, y, info)
    complex(psb_spk_), intent(in)        :: x(:)
    class(psb_c_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n


    info = 0
    if (allocated(y%v)) &
         & call y%v%mlt(x,info)

  end subroutine c_vect_mlt_a

  module subroutine c_vect_mlt_a_2(alpha,x,y,beta,z,info)
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

  module subroutine c_vect_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
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

  module subroutine c_vect_mlt_av(alpha,x,y,beta,z,info)
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

  module subroutine c_vect_mlt_va(alpha,x,y,beta,z,info)
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

  module subroutine c_vect_div_v(x, y, info)
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(x%v).and.allocated(y%v)) &
         & call y%v%div(x%v,info)

  end subroutine c_vect_div_v

  module subroutine c_vect_div_v2( x, y, z, info)
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    class(psb_c_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(x%v).and.allocated(y%v).and.allocated(z%v)) &
         & call z%v%div(x%v,y%v,info)

  end subroutine c_vect_div_v2

  module subroutine c_vect_div_v_check(x, y, info, flag)
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n
    logical, intent(in) :: flag

    info = 0
    if (allocated(x%v).and.allocated(y%v)) &
         & call y%v%div(x%v,info,flag)

  end subroutine c_vect_div_v_check

  module subroutine c_vect_div_v2_check(x, y, z, info, flag)
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    class(psb_c_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n
    logical, intent(in) :: flag

    info = 0
    if (allocated(x%v).and.allocated(y%v).and.allocated(z%v)) &
         & call z%v%div(x%v,y%v,info,flag)

  end subroutine c_vect_div_v2_check

  module subroutine c_vect_div_a2(x, y, z, info)
    complex(psb_spk_), intent(in) :: x(:)
    complex(psb_spk_), intent(in)    :: y(:)
    class(psb_c_vect_type), intent(inout) :: z
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(z%v)) &
         & call z%v%div(x,y,info)

  end subroutine c_vect_div_a2

  module subroutine c_vect_div_a2_check(x, y, z, info,flag)
    complex(psb_spk_), intent(in) :: x(:)
    complex(psb_spk_), intent(in)    :: y(:)
    class(psb_c_vect_type), intent(inout) :: z
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n
    logical, intent(in) :: flag

    info = 0
    if (allocated(z%v)) &
         & call z%v%div(x,y,info,flag)

  end subroutine c_vect_div_a2_check

  module subroutine c_vect_inv_v(x, y, info)
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(x%v).and.allocated(y%v)) &
         & call y%v%inv(x%v,info)

  end subroutine c_vect_inv_v

  module subroutine c_vect_inv_v_check(x, y, info, flag)
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n
    logical, intent(in) :: flag

    info = 0
    if (allocated(x%v).and.allocated(y%v)) &
         & call y%v%inv(x%v,info,flag)

  end subroutine c_vect_inv_v_check

  module subroutine c_vect_inv_a2(x, y, info)
    complex(psb_spk_), intent(inout) :: x(:)
    class(psb_c_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (allocated(y%v)) &
         & call y%v%inv(x,info)

  end subroutine c_vect_inv_a2

  module subroutine c_vect_inv_a2_check(x, y, info,flag)
    complex(psb_spk_), intent(inout) :: x(:)
    class(psb_c_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n
    logical, intent(in) :: flag

    info = 0
    if (allocated(y%v)) &
         & call y%v%inv(x,info,flag)

  end subroutine c_vect_inv_a2_check

  module subroutine c_vect_acmp_a2(x,c,z,info)
    real(psb_spk_), intent(in)              :: c
    complex(psb_spk_), intent(inout)           :: x(:)
    class(psb_c_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)           :: info

    info = 0
    if (allocated(z%v)) &
         & call z%acmp(x,c,info)

  end subroutine c_vect_acmp_a2

  module subroutine c_vect_acmp_v2(x,c,z,info)
    real(psb_spk_), intent(in)              :: c
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)           :: info

    info = 0
    if (allocated(x%v).and.allocated(z%v)) &
         & call z%v%acmp(x%v,c,info)

  end subroutine c_vect_acmp_v2

  module subroutine c_vect_scal(alpha, x)
    class(psb_c_vect_type), intent(inout)  :: x
    complex(psb_spk_), intent (in)       :: alpha

    if (allocated(x%v)) call x%v%scal(alpha)

  end subroutine c_vect_scal

  module subroutine c_vect_absval1(x)
    class(psb_c_vect_type), intent(inout)  :: x

    if (allocated(x%v)) &
         &  call x%v%absval()

  end subroutine c_vect_absval1

  module subroutine c_vect_absval2(x,y)
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y

    if (allocated(x%v)) then
      if (.not.allocated(y%v))  call y%bld(psb_size(x%v%v))
      call x%v%absval(y%v)
    end if
  end subroutine c_vect_absval2

  module function c_vect_nrm2(n,x) result(res)
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)                :: res

    if (allocated(x%v)) then
      res = x%v%nrm2(n)
    else
      res = szero
    end if

  end function c_vect_nrm2

  module function c_vect_nrm2_weight(n,x,w,aux) result(res)
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_vect_type), intent(inout) :: w
    class(psb_c_vect_type), intent(inout), optional :: aux
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)                        :: res
    integer(psb_ipk_)                       :: info

    ! Temp vectors
    type(psb_c_vect_type) :: wtemp

    info = 0
    if( allocated(w%v) ) then
      if (.not.present(aux)) then
        allocate(wtemp%v, mold=w%v)
        call wtemp%v%bld(w%get_vect())
      else
        call psb_geaxpby(n,cone,w%v%v,czero,aux%v%v,info)
      end if
    else
      info = -1
    end if
    if (info /= 0 ) then
      res = -sone
      return
    end if

    if (allocated(x%v)) then
      if (.not.present(aux)) then
        call wtemp%v%mlt(x%v,info)
        res = wtemp%v%nrm2(n)
      else
        call aux%v%mlt(x%v,info)
        res = aux%v%nrm2(n)
      end if
    else
      res = szero
    end if

    if (.not.present(aux)) then
      call wtemp%free(info)
    end if

  end function c_vect_nrm2_weight
 
  module function c_vect_nrm2_weight_mask(n,x,w,id,info,aux) result(res)
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_vect_type), intent(inout) :: w
    class(psb_c_vect_type), intent(inout) :: id
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)                        :: res
    integer(psb_ipk_), intent(out)          :: info
    class(psb_c_vect_type), intent(inout), optional :: aux

    ! Temp vectors
    type(psb_c_vect_type) :: wtemp

    info = 0 
    if( allocated(w%v) ) then
      if (.not.present(aux)) then
        allocate(wtemp%v, mold=w%v)
        call wtemp%v%bld(w%get_vect())
      else
        call psb_geaxpby(n,cone,w%v%v,czero,aux%v%v,info)
      end if
    else
      info = -1
    end if
    if (info /= 0 ) then
      res = -sone
      return
    end if

    if (allocated(x%v).and.allocated(id%v)) then
      if (.not.present(aux)) then
        where( abs(id%v%v) <= szero) wtemp%v%v = szero
        call wtemp%set_host()
        call wtemp%v%mlt(x%v,info)
        res = wtemp%v%nrm2(n)
      else
        where( abs(id%v%v) <= szero) aux%v%v = szero
        call aux%set_host()
        call aux%v%mlt(x%v,info)
        res = aux%v%nrm2(n)
      end if
    else
      res = szero
    end if

    if (.not.present(aux)) then
      call wtemp%free(info)
    end if

  end function c_vect_nrm2_weight_mask

  module function c_vect_amax(n,x) result(res)
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)                :: res

    if (allocated(x%v)) then
      res = x%v%amax(n)
    else
      res = szero
    end if

  end function c_vect_amax


  module function c_vect_asum(n,x) result(res)
    class(psb_c_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_spk_)                :: res

    if (allocated(x%v)) then
      res = x%v%asum(n)
    else
      res = szero
    end if

  end function c_vect_asum

  module subroutine c_vect_addconst_a2(x,b,z,info)
    real(psb_spk_), intent(in)             :: b
    complex(psb_spk_), intent(inout)           :: x(:)
    class(psb_c_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)           :: info

    info = 0
    if (allocated(z%v)) &
         & call z%addconst(x,b,info)

  end subroutine c_vect_addconst_a2

  module subroutine c_vect_addconst_v2(x,b,z,info)
    real(psb_spk_), intent(in)             :: b
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)           :: info

    info = 0
    if (allocated(x%v).and.allocated(z%v)) &
         & call z%v%addconst(x%v,b,info)

  end subroutine c_vect_addconst_v2

end submodule psb_c_vect_impl


submodule (psb_c_multivect_mod) psb_c_multivect_impl
  use psb_base_mod
  use psi_serial_mod

contains
  
  module function c_mvect_get_dupl(x) result(res)
    class(psb_c_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = x%dupl
  end function c_mvect_get_dupl

  module subroutine c_mvect_set_dupl(x,val)
    class(psb_c_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if (present(val)) then
      x%dupl = val
    else
      x%dupl = psb_dupl_def_
    end if
  end subroutine c_mvect_set_dupl
        
  module function c_mvect_is_remote_build(x) result(res)
    class(psb_c_multivect_type), intent(in) :: x
    logical :: res
    res = (x%remote_build == psb_matbld_remote_)
  end function c_mvect_is_remote_build

  module subroutine c_mvect_set_remote_build(x,val)
    class(psb_c_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if (present(val)) then
      x%remote_build = val
    else
      x%remote_build = psb_matbld_remote_
    end if
  end subroutine c_mvect_set_remote_build
       
  module subroutine c_mvect_clone(x,y,info)
    class(psb_c_multivect_type), intent(inout) :: x
    class(psb_c_multivect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    call y%free(info)
    if ((info==0).and.allocated(x%v)) then
      call y%bld_x(x%get_vect(),mold=x%v)
    end if
  end subroutine c_mvect_clone

  module subroutine c_mvect_bld_x(x,invect,mold)
    complex(psb_spk_), intent(in)          :: invect(:,:)
    class(psb_c_multivect_type), intent(out) :: x
    class(psb_c_base_multivect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info
    class(psb_c_base_multivect_type), pointer :: mld

    info = psb_success_
    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_c_get_base_multivect_default())
    endif

    if (info == psb_success_) call x%v%bld(invect)

  end subroutine c_mvect_bld_x

  module subroutine c_mvect_bld_n(x,m,n,mold,scratch)
    integer(psb_ipk_), intent(in) :: m,n
    class(psb_c_multivect_type), intent(out) :: x
    class(psb_c_base_multivect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info
    logical, intent(in), optional        :: scratch
    
    info = psb_success_
    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_c_get_base_multivect_default())
    endif
    if (info == psb_success_) call x%v%bld(m,n,scratch=scratch)

  end subroutine c_mvect_bld_n

  module function  c_mvect_get_vect(x) result(res)
    class(psb_c_multivect_type), intent(inout)  :: x
    complex(psb_spk_), allocatable                 :: res(:,:)
    integer(psb_ipk_) :: info

    if (allocated(x%v)) then
      res = x%v%get_vect()
    end if
  end function c_mvect_get_vect

  module subroutine c_mvect_set_scal(x,val)
    class(psb_c_multivect_type), intent(inout)  :: x
    complex(psb_spk_), intent(in) :: val

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val)

  end subroutine c_mvect_set_scal

  module subroutine c_mvect_set_vect(x,val)
    class(psb_c_multivect_type), intent(inout) :: x
    complex(psb_spk_), intent(in)         :: val(:,:)

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val)

  end subroutine c_mvect_set_vect

  module function c_mvect_get_nrows(x) result(res)
    class(psb_c_multivect_type), intent(in) :: x
    integer(psb_ipk_)  :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_nrows()
  end function c_mvect_get_nrows

  module function c_mvect_get_ncols(x) result(res)
    class(psb_c_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_ncols()
  end function c_mvect_get_ncols

  module function c_mvect_sizeof(x) result(res)
    class(psb_c_multivect_type), intent(in) :: x
    integer(psb_epk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%sizeof()
  end function c_mvect_sizeof

  module function c_mvect_get_fmt(x) result(res)
    class(psb_c_multivect_type), intent(in) :: x
    character(len=5) :: res
    res = 'NULL'
    if (allocated(x%v)) res = x%v%get_fmt()
  end function c_mvect_get_fmt

  module subroutine c_mvect_all(m,n, x, info, mold)

    integer(psb_ipk_), intent(in)       :: m,n
    class(psb_c_multivect_type), intent(out) :: x
    class(psb_c_base_multivect_type), intent(in), optional :: mold
    integer(psb_ipk_), intent(out)      :: info

    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(psb_c_base_multivect_type :: x%v,stat=info)
    endif
    if (info == 0) then
      call x%v%all(m,n,info)
    else
      info = psb_err_alloc_dealloc_
    end if

  end subroutine c_mvect_all

  module subroutine c_mvect_reall(m,n, x, info)

    integer(psb_ipk_), intent(in)         :: m,n
    class(psb_c_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    info = 0
    if (.not.allocated(x%v)) &
         & call x%all(m,n,info)
    if (info == 0) &
         & call x%asb(m,n,info)

  end subroutine c_mvect_reall

  module subroutine c_mvect_zero(x)
    class(psb_c_multivect_type), intent(inout)    :: x

    if (allocated(x%v)) call x%v%zero()

  end subroutine c_mvect_zero

  module subroutine c_mvect_asb(m,n, x, info)
    integer(psb_ipk_), intent(in)              :: m,n
    class(psb_c_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%v)) &
         & call x%v%asb(m,n,info)

  end subroutine c_mvect_asb

  module subroutine c_mvect_sync(x)
    class(psb_c_multivect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%sync()

  end subroutine c_mvect_sync

  module subroutine c_mvect_gthab(n,idx,alpha,x,beta,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    complex(psb_spk_) :: alpha, beta, y(:)
    class(psb_c_multivect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,alpha,beta,y)

  end subroutine c_mvect_gthab

  module subroutine c_mvect_gthzv(n,idx,x,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    complex(psb_spk_) ::  y(:)
    class(psb_c_multivect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,y)

  end subroutine c_mvect_gthzv

  module subroutine c_mvect_gthzv_x(i,n,idx,x,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: i
    class(psb_i_base_vect_type) :: idx
    complex(psb_spk_) ::  y(:)
    class(psb_c_multivect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(i,n,idx,y)

  end subroutine c_mvect_gthzv_x

  module subroutine c_mvect_sctb(n,idx,x,beta,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    complex(psb_spk_) :: beta, x(:)
    class(psb_c_multivect_type) :: y

    if (allocated(y%v)) &
         &  call y%v%sct(n,idx,x,beta)

  end subroutine c_mvect_sctb

  module subroutine c_mvect_sctb_x(i,n,idx,x,beta,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: i
    class(psb_i_base_vect_type) :: idx
    complex(psb_spk_) :: beta, x(:)
    class(psb_c_multivect_type) :: y

    if (allocated(y%v)) &
         &  call y%v%sct(i,n,idx,x,beta)

  end subroutine c_mvect_sctb_x

  module subroutine c_mvect_free(x, info)
    class(psb_c_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (allocated(x%v)) then
      call x%v%free(info)
      if (info == 0) deallocate(x%v,stat=info)
    end if

  end subroutine c_mvect_free

  module subroutine c_mvect_ins(n,irl,val,x,maxr,info)
    class(psb_c_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n,maxr
    integer(psb_ipk_), intent(in)               :: irl(:)
    complex(psb_spk_), intent(in)        :: val(:,:)
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i, dupl

    info = 0
    if (.not.allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if
    dupl = x%get_dupl()
    call  x%v%ins(n,irl,val,dupl,maxr,info)

  end subroutine c_mvect_ins

  module subroutine c_mvect_cnv(x,mold)
    class(psb_c_multivect_type), intent(inout) :: x
    class(psb_c_base_multivect_type), intent(in), optional :: mold
    class(psb_c_base_multivect_type), allocatable :: tmp
    integer(psb_ipk_) :: info

    if (present(mold)) then
      allocate(tmp,stat=info,mold=mold)
    else
      allocate(tmp,stat=info, mold=psb_c_get_base_multivect_default())
    endif
    if (allocated(x%v)) then
      call x%v%sync()
      if (info == psb_success_) call tmp%bld(x%v%v)
      call x%v%free(info)
    end if
    call move_alloc(tmp,x%v)
  end subroutine c_mvect_cnv

!!$  module function c_mvect_dot_v(n,x,y) result(res)
!!$    class(psb_c_multivect_type), intent(inout) :: x, y
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    complex(psb_spk_)                :: res
!!$
!!$    res = czero
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & res = x%v%dot(n,y%v)
!!$
!!$  end function c_mvect_dot_v
!!$
!!$  function c_mvect_dot_a(n,x,y) result(res)
!!$    class(psb_c_multivect_type), intent(inout) :: x
!!$    complex(psb_spk_), intent(in)    :: y(:)
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    complex(psb_spk_)                :: res
!!$
!!$    res = czero
!!$    if (allocated(x%v)) &
!!$         & res = x%v%dot(n,y)
!!$
!!$  end function c_mvect_dot_a
!!$
!!$  module subroutine c_mvect_axpby_v(m,alpha, x, beta, y, info)
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    class(psb_c_multivect_type), intent(inout)  :: x
!!$    class(psb_c_multivect_type), intent(inout)  :: y
!!$    complex(psb_spk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    if (allocated(x%v).and.allocated(y%v)) then
!!$      call y%v%axpby(m,alpha,x%v,beta,info)
!!$    else
!!$      info = psb_err_invalid_mvect_state_
!!$    end if
!!$
!!$  end subroutine c_mvect_axpby_v
!!$
!!$  subroutine c_mvect_axpby_a(m,alpha, x, beta, y, info)
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    complex(psb_spk_), intent(in)        :: x(:)
!!$    class(psb_c_multivect_type), intent(inout)  :: y
!!$    complex(psb_spk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    if (allocated(y%v)) &
!!$         & call y%v%axpby(m,alpha,x,beta,info)
!!$
!!$  end subroutine c_mvect_axpby_a
!!$
!!$
!!$  subroutine c_mvect_mlt_v(x, y, info)
!!$    class(psb_c_multivect_type), intent(inout)  :: x
!!$    class(psb_c_multivect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & call y%v%mlt(x%v,info)
!!$
!!$  end subroutine c_mvect_mlt_v
!!$
!!$  subroutine c_mvect_mlt_a(x, y, info)
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
!!$  end subroutine c_mvect_mlt_a
!!$
!!$
!!$  subroutine c_mvect_mlt_a_2(alpha,x,y,beta,z,info)
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
!!$  end subroutine c_mvect_mlt_a_2
!!$
!!$  subroutine c_mvect_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
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
!!$  end subroutine c_mvect_mlt_v_2
!!$
!!$  subroutine c_mvect_mlt_av(alpha,x,y,beta,z,info)
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
!!$  end subroutine c_mvect_mlt_av
!!$
!!$  subroutine c_mvect_mlt_va(alpha,x,y,beta,z,info)
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
!!$  end subroutine c_mvect_mlt_va
!!$
!!$  subroutine c_mvect_scal(alpha, x)
!!$    class(psb_c_multivect_type), intent(inout)  :: x
!!$    complex(psb_spk_), intent (in)       :: alpha
!!$
!!$    if (allocated(x%v)) call x%v%scal(alpha)
!!$
!!$  end subroutine c_mvect_scal
!!$
!!$
!!$  function c_mvect_nrm2(n,x) result(res)
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
!!$  end function c_mvect_nrm2
!!$
!!$  function c_mvect_amax(n,x) result(res)
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
!!$  end function c_mvect_amax
!!$
!!$  function c_mvect_asum(n,x) result(res)
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
!!$  end function c_mvect_asum

end submodule psb_c_multivect_impl
