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
submodule (psb_l_vect_mod) psb_l_vect_impl
  use psb_base_mod
  use psi_serial_mod
  implicit none 

contains

  module function l_vect_get_dupl(x) result(res)
    class(psb_l_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    if (allocated(x%v)) then 
      res = x%v%get_dupl()
    else
      res = psb_dupl_null_
    end if
  end function l_vect_get_dupl

  module subroutine l_vect_set_dupl(x,val)
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if (allocated(x%v)) then 
      if (present(val)) then
        call x%v%set_dupl(val)
      else
        call x%v%set_dupl(psb_dupl_def_)
      end if
    end if
  end subroutine l_vect_set_dupl

  module function l_vect_get_ncfs(x) result(res)
    class(psb_l_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    if (allocated(x%v)) then 
      res = x%v%get_ncfs()
    else
      res = 0
    end if
  end function l_vect_get_ncfs

  module subroutine l_vect_set_ncfs(x,val)
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if (allocated(x%v)) then 
      if (present(val)) then
        call x%v%set_ncfs(val)
      else
        call x%v%set_ncfs(0)
      end if
    end if
  end subroutine l_vect_set_ncfs

  module function l_vect_get_state(x) result(res)
    class(psb_l_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    if (allocated(x%v)) then 
      res = x%v%get_state()
    else
      res = psb_vect_null_
    end if
  end function l_vect_get_state

  module function l_vect_is_null(x) result(res)
    class(psb_l_vect_type), intent(in) :: x
    logical :: res
    res = (x%get_state() == psb_vect_null_)
  end function l_vect_is_null

  module function l_vect_is_bld(x) result(res)
    class(psb_l_vect_type), intent(in) :: x
    logical :: res
    res = (x%get_state() == psb_vect_bld_)
  end function l_vect_is_bld

  module function l_vect_is_upd(x) result(res)
    class(psb_l_vect_type), intent(in) :: x
    logical :: res
    res = (x%get_state() == psb_vect_upd_)
  end function l_vect_is_upd

  module function l_vect_is_asb(x) result(res)
    class(psb_l_vect_type), intent(in) :: x
    logical :: res
    res = (x%get_state() == psb_vect_asb_)
  end function l_vect_is_asb

  module subroutine  l_vect_set_state(n,x)
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n
    if (allocated(x%v)) then 
      call x%v%set_state(n)
    end if
  end subroutine l_vect_set_state

  module subroutine  l_vect_set_null(x)
    class(psb_l_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_null_)
  end subroutine l_vect_set_null

  module subroutine  l_vect_set_bld(x)
    class(psb_l_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_bld_)
  end subroutine l_vect_set_bld

  module subroutine  l_vect_set_upd(x)
    class(psb_l_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_upd_)
  end subroutine l_vect_set_upd

  module subroutine  l_vect_set_asb(x)
    class(psb_l_vect_type), intent(inout) :: x

    call x%set_state(psb_vect_asb_)
  end subroutine l_vect_set_asb

  module function l_vect_get_nrmv(x) result(res)
    class(psb_l_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = x%nrmv
  end function l_vect_get_nrmv

  module subroutine l_vect_set_nrmv(x,val)
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: val

    x%nrmv = val
  end subroutine l_vect_set_nrmv

  module function l_vect_is_remote_build(x) result(res)
    class(psb_l_vect_type), intent(in) :: x
    logical :: res
    res = (x%remote_build == psb_matbld_remote_)
  end function l_vect_is_remote_build

  module subroutine l_vect_set_remote_build(x,val)
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if (present(val)) then
      x%remote_build = val
    else
      x%remote_build = psb_matbld_remote_
    end if
  end subroutine l_vect_set_remote_build

  module subroutine l_vect_clone(x,y,info)
    class(psb_l_vect_type), intent(inout) :: x
    class(psb_l_vect_type), intent(inout) :: y
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
  end subroutine l_vect_clone

  module subroutine l_vect_bld_x(x,invect,mold,scratch)
    integer(psb_lpk_), intent(in)          :: invect(:)
    class(psb_l_vect_type), intent(inout) :: x
    class(psb_l_base_vect_type), intent(in), optional :: mold
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
      allocate(x%v,stat=info, mold=psb_l_get_base_vect_default())
    endif

    if (info == psb_success_) call x%v%bld(invect,scratch=scratch_)

  end subroutine l_vect_bld_x

  module subroutine l_vect_bld_mn(x,n,mold,scratch)
    integer(psb_mpk_), intent(in) :: n
    class(psb_l_vect_type), intent(inout) :: x
    class(psb_l_base_vect_type), intent(in), optional :: mold
    logical, intent(in), optional        :: scratch

    logical :: scratch_
    integer(psb_ipk_) :: info
    class(psb_l_base_vect_type), pointer :: mld
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
      allocate(x%v,stat=info, mold=psb_l_get_base_vect_default())
    endif
    if (info == psb_success_) call x%v%bld(n,scratch=scratch_)

  end subroutine l_vect_bld_mn

  module subroutine l_vect_bld_en(x,n,mold,scratch)
    integer(psb_epk_), intent(in) :: n
    class(psb_l_vect_type), intent(inout) :: x
    class(psb_l_base_vect_type), intent(in), optional :: mold
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
      allocate(x%v,stat=info, mold=psb_l_get_base_vect_default())
    endif
    if (info == psb_success_) call x%v%bld(n,scratch=scratch_)

  end subroutine l_vect_bld_en

  module function  l_vect_get_vect(x,n) result(res)
    class(psb_l_vect_type), intent(inout)  :: x
    integer(psb_lpk_), allocatable                 :: res(:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional :: n

    if (allocated(x%v)) then
      res = x%v%get_vect(n)
    end if
  end function l_vect_get_vect

  module subroutine l_vect_set_scal(x,val,first,last)
    class(psb_l_vect_type), intent(inout)  :: x
    integer(psb_lpk_), intent(in) :: val
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val,first,last)

  end subroutine l_vect_set_scal

  module subroutine l_vect_set_vect(x,val,first,last)
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_lpk_), intent(in)         :: val(:)
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val,first,last)

  end subroutine l_vect_set_vect

  module subroutine l_vect_check_addr(x)
    class(psb_l_vect_type), intent(inout) :: x

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%check_addr()

  end subroutine l_vect_check_addr

  module function l_vect_get_nrows(x) result(res)
    class(psb_l_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_nrows()
  end function l_vect_get_nrows

  module function l_vect_sizeof(x) result(res)
    class(psb_l_vect_type), intent(in) :: x
    integer(psb_epk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%sizeof()
  end function l_vect_sizeof

  module function l_vect_get_fmt(x) result(res)
    class(psb_l_vect_type), intent(in) :: x
    character(len=5) :: res
    res = 'NULL'
    if (allocated(x%v)) res = x%v%get_fmt()
  end function l_vect_get_fmt

  module subroutine l_vect_all(n, x, info, mold)

    integer(psb_ipk_), intent(in)           :: n
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)      :: info
    class(psb_l_base_vect_type), intent(in), optional :: mold

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
    call x%set_bld()
  end subroutine l_vect_all

  module subroutine l_vect_reinit(x, info, clear)
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)      :: info
    logical, intent(in), optional       :: clear

    if (allocated(x%v)) call x%v%reinit(info,clear)
    call x%set_upd()

  end subroutine l_vect_reinit
 
  module subroutine l_vect_reall(n, x, info)

    integer(psb_ipk_), intent(in)         :: n
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    info = 0
    if (.not.allocated(x%v)) &
         & call x%all(n,info)
    if (info == 0) &
         & call x%asb(n,info)

  end subroutine l_vect_reall

  module subroutine l_vect_zero(x)
    class(psb_l_vect_type), intent(inout)    :: x

    if (allocated(x%v)) call x%v%zero()

  end subroutine l_vect_zero

  module subroutine l_vect_asb(n, x, info, scratch)
    integer(psb_ipk_), intent(in)         :: n
    class(psb_l_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info
    logical, intent(in), optional        :: scratch

    if (allocated(x%v)) then
      call x%v%asb(n,info,scratch=scratch)
      call x%set_asb()
    end if
  end subroutine l_vect_asb

  module subroutine l_vect_gthab(n,idx,alpha,x,beta,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_lpk_) :: alpha, beta, y(:)
    class(psb_l_vect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,alpha,beta,y)

  end subroutine l_vect_gthab

  module subroutine l_vect_gthzv(n,idx,x,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_lpk_) ::  y(:)
    class(psb_l_vect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,y)

  end subroutine l_vect_gthzv

  module subroutine l_vect_sctb(n,idx,x,beta,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_lpk_) :: beta, x(:)
    class(psb_l_vect_type) :: y

    if (allocated(y%v)) &
         &  call y%v%sct(n,idx,x,beta)

  end subroutine l_vect_sctb

  module subroutine l_vect_free(x, info)
    class(psb_l_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (allocated(x%v)) then
      call x%v%free(info)
      if (info == 0) deallocate(x%v,stat=info)
    end if

  end subroutine l_vect_free

  module subroutine l_vect_ins_a(n,irl,val,x,maxr,info)
    class(psb_l_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, maxr
    integer(psb_ipk_), intent(in)               :: irl(:)
    integer(psb_lpk_), intent(in)        :: val(:)
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i, dupl

    info = 0
    if (.not.allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if
    dupl = x%get_dupl()
    call  x%v%ins(n,irl,val,dupl,maxr,info)

  end subroutine l_vect_ins_a

  module subroutine l_vect_ins_v(n,irl,val,x,maxr,info)
    class(psb_l_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, maxr
    class(psb_i_vect_type), intent(inout)       :: irl
    class(psb_l_vect_type), intent(inout)       :: val
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i, dupl

    info = 0
    if (.not.(allocated(x%v).and.allocated(irl%v).and.allocated(val%v))) then
      info = psb_err_invalid_vect_state_
      return
    end if
    dupl = x%get_dupl()
    call  x%v%ins(n,irl%v,val%v,dupl,maxr,info)

  end subroutine l_vect_ins_v

  module subroutine l_vect_cnv(x,mold)
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
      if (allocated(x%v%v)) then 
        call x%v%sync()
        if (info == psb_success_) call tmp%bld(x%v%v)
        call x%v%base_cpy(tmp)
        call x%v%free(info)
      endif
    end if
    call move_alloc(tmp,x%v)

  end subroutine l_vect_cnv

  module subroutine l_vect_sync(x)
    class(psb_l_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%sync()

  end subroutine l_vect_sync

  module subroutine l_vect_set_sync(x)
    class(psb_l_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%set_sync()

  end subroutine l_vect_set_sync

  module subroutine l_vect_set_host(x)
    class(psb_l_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%set_host()

  end subroutine l_vect_set_host

  module subroutine l_vect_set_dev(x)
    class(psb_l_vect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%set_dev()

  end subroutine l_vect_set_dev

  module function l_vect_is_sync(x) result(res)
    logical :: res
    class(psb_l_vect_type), intent(inout) :: x

    res = .true.
    if (allocated(x%v)) &
         & res = x%v%is_sync()

  end function l_vect_is_sync

  module function l_vect_is_host(x) result(res)
    logical :: res
    class(psb_l_vect_type), intent(inout) :: x

    res = .true.
    if (allocated(x%v)) &
         & res = x%v%is_host()

  end function l_vect_is_host

  module function l_vect_is_dev(x) result(res)
    logical :: res
    class(psb_l_vect_type), intent(inout) :: x

    res = .false.
    if (allocated(x%v)) &
         & res =  x%v%is_dev()

  end function l_vect_is_dev


end submodule psb_l_vect_impl


submodule (psb_l_multivect_mod) psb_l_multivect_impl
  use psb_base_mod
  use psi_serial_mod

contains
  
  module function l_mvect_get_dupl(x) result(res)
    class(psb_l_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = x%dupl
  end function l_mvect_get_dupl

  module subroutine l_mvect_set_dupl(x,val)
    class(psb_l_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if (present(val)) then
      x%dupl = val
    else
      x%dupl = psb_dupl_def_
    end if
  end subroutine l_mvect_set_dupl
        
  module function l_mvect_is_remote_build(x) result(res)
    class(psb_l_multivect_type), intent(in) :: x
    logical :: res
    res = (x%remote_build == psb_matbld_remote_)
  end function l_mvect_is_remote_build

  module subroutine l_mvect_set_remote_build(x,val)
    class(psb_l_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in), optional :: val

    if (present(val)) then
      x%remote_build = val
    else
      x%remote_build = psb_matbld_remote_
    end if
  end subroutine l_mvect_set_remote_build
       
  module subroutine l_mvect_clone(x,y,info)
    class(psb_l_multivect_type), intent(inout) :: x
    class(psb_l_multivect_type), intent(inout) :: y
    integer(psb_ipk_), intent(out)        :: info

    info = psb_success_
    call y%free(info)
    if ((info==0).and.allocated(x%v)) then
      call y%bld_x(x%get_vect(),mold=x%v)
    end if
  end subroutine l_mvect_clone

  module subroutine l_mvect_bld_x(x,invect,mold)
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

  end subroutine l_mvect_bld_x

  module subroutine l_mvect_bld_n(x,m,n,mold,scratch)
    integer(psb_ipk_), intent(in) :: m,n
    class(psb_l_multivect_type), intent(out) :: x
    class(psb_l_base_multivect_type), intent(in), optional :: mold
    integer(psb_ipk_) :: info
    logical, intent(in), optional        :: scratch
    
    info = psb_success_
    if (present(mold)) then
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(x%v,stat=info, mold=psb_l_get_base_multivect_default())
    endif
    if (info == psb_success_) call x%v%bld(m,n,scratch=scratch)

  end subroutine l_mvect_bld_n

  module function  l_mvect_get_vect(x) result(res)
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_lpk_), allocatable                 :: res(:,:)
    integer(psb_ipk_) :: info

    if (allocated(x%v)) then
      res = x%v%get_vect()
    end if
  end function l_mvect_get_vect

  module subroutine l_mvect_set_scal(x,val)
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_lpk_), intent(in) :: val

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val)

  end subroutine l_mvect_set_scal

  module subroutine l_mvect_set_vect(x,val)
    class(psb_l_multivect_type), intent(inout) :: x
    integer(psb_lpk_), intent(in)         :: val(:,:)

    integer(psb_ipk_) :: info
    if (allocated(x%v)) call x%v%set(val)

  end subroutine l_mvect_set_vect

  module function l_mvect_get_nrows(x) result(res)
    class(psb_l_multivect_type), intent(in) :: x
    integer(psb_ipk_)  :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_nrows()
  end function l_mvect_get_nrows

  module function l_mvect_get_ncols(x) result(res)
    class(psb_l_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_ncols()
  end function l_mvect_get_ncols

  module function l_mvect_sizeof(x) result(res)
    class(psb_l_multivect_type), intent(in) :: x
    integer(psb_epk_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%sizeof()
  end function l_mvect_sizeof

  module function l_mvect_get_fmt(x) result(res)
    class(psb_l_multivect_type), intent(in) :: x
    character(len=5) :: res
    res = 'NULL'
    if (allocated(x%v)) res = x%v%get_fmt()
  end function l_mvect_get_fmt

  module subroutine l_mvect_all(m,n, x, info, mold)

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

  end subroutine l_mvect_all

  module subroutine l_mvect_reall(m,n, x, info)

    integer(psb_ipk_), intent(in)         :: m,n
    class(psb_l_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)        :: info

    info = 0
    if (.not.allocated(x%v)) &
         & call x%all(m,n,info)
    if (info == 0) &
         & call x%asb(m,n,info)

  end subroutine l_mvect_reall

  module subroutine l_mvect_zero(x)
    class(psb_l_multivect_type), intent(inout)    :: x

    if (allocated(x%v)) call x%v%zero()

  end subroutine l_mvect_zero

  module subroutine l_mvect_asb(m,n, x, info)
    integer(psb_ipk_), intent(in)              :: m,n
    class(psb_l_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%v)) &
         & call x%v%asb(m,n,info)

  end subroutine l_mvect_asb

  module subroutine l_mvect_sync(x)
    class(psb_l_multivect_type), intent(inout) :: x

    if (allocated(x%v)) &
         & call x%v%sync()

  end subroutine l_mvect_sync

  module subroutine l_mvect_gthab(n,idx,alpha,x,beta,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_lpk_) :: alpha, beta, y(:)
    class(psb_l_multivect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,alpha,beta,y)

  end subroutine l_mvect_gthab

  module subroutine l_mvect_gthzv(n,idx,x,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_lpk_) ::  y(:)
    class(psb_l_multivect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,y)

  end subroutine l_mvect_gthzv

  module subroutine l_mvect_gthzv_x(i,n,idx,x,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: i
    class(psb_i_base_vect_type) :: idx
    integer(psb_lpk_) ::  y(:)
    class(psb_l_multivect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(i,n,idx,y)

  end subroutine l_mvect_gthzv_x

  module subroutine l_mvect_sctb(n,idx,x,beta,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_lpk_) :: beta, x(:)
    class(psb_l_multivect_type) :: y

    if (allocated(y%v)) &
         &  call y%v%sct(n,idx,x,beta)

  end subroutine l_mvect_sctb

  module subroutine l_mvect_sctb_x(i,n,idx,x,beta,y)
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: i
    class(psb_i_base_vect_type) :: idx
    integer(psb_lpk_) :: beta, x(:)
    class(psb_l_multivect_type) :: y

    if (allocated(y%v)) &
         &  call y%v%sct(i,n,idx,x,beta)

  end subroutine l_mvect_sctb_x

  module subroutine l_mvect_free(x, info)
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (allocated(x%v)) then
      call x%v%free(info)
      if (info == 0) deallocate(x%v,stat=info)
    end if

  end subroutine l_mvect_free

  module subroutine l_mvect_ins(n,irl,val,x,maxr,info)
    class(psb_l_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n,maxr
    integer(psb_ipk_), intent(in)               :: irl(:)
    integer(psb_lpk_), intent(in)        :: val(:,:)
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i, dupl

    info = 0
    if (.not.allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      return
    end if
    dupl = x%get_dupl()
    call  x%v%ins(n,irl,val,dupl,maxr,info)

  end subroutine l_mvect_ins

  module subroutine l_mvect_cnv(x,mold)
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
  end subroutine l_mvect_cnv


end submodule psb_l_multivect_impl
