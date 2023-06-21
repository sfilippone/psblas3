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
module psb_m_realloc_mod
  use psb_const_mod

  implicit none

  !
  ! psb_realloc will reallocate the input array to have exactly 
  ! the size specified, possibly shortening it. 
  !
  Interface psb_realloc
    module procedure psb_r_m_s
    module procedure psb_r_m_m_rk1
    module procedure psb_r_m_m_rk2
    module procedure psb_r_e_m_rk1
    module procedure psb_r_e_m_rk2
    module procedure psb_r_me_m_rk2
    module procedure psb_r_em_m_rk2

    module procedure psb_r_m_2_m_rk1
    module procedure psb_r_e_2_m_rk1
    
  end Interface psb_realloc

  interface psb_move_alloc
    module procedure psb_move_alloc_m_rk1, psb_move_alloc_m_rk2
  end interface psb_move_alloc

  Interface psb_safe_ab_cpy
    module procedure psb_ab_cpy_m_s, psb_ab_cpy_m_rk1, psb_ab_cpy_m_rk2
  end Interface psb_safe_ab_cpy

  Interface psb_safe_cpy
    module procedure psb_cpy_m_rk1, psb_cpy_m_rk2
  end Interface psb_safe_cpy

  !
  ! psb_ensure_size will reallocate the input array if necessary
  ! to guarantee that its size is at least as large as the 
  ! value required, usually with some room to spare.
  !
  interface psb_ensure_size
    module procedure psb_ensure_m_sz_m_rk1, psb_ensure_e_sz_m_rk1
  end Interface psb_ensure_size

  !
  ! psb_size returns 0 if argument is not allocated. 
  !
  interface psb_size
    module procedure psb_size_m_rk1, psb_size_m_rk2
  end interface psb_size


Contains

  Subroutine psb_r_m_s(rrax,info)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_mpk_), allocatable, intent(inout) :: rrax
    integer(psb_ipk_) :: info

    ! ...Local Variables
    integer(psb_ipk_) :: err_act,err
    character(len=30)  :: name
    logical, parameter :: debug=.false.

    name='psb_r_m_s'
    call psb_erractionsave(err_act)
    info=psb_success_ 

    if (.not.allocated(rrax)) then 
      Allocate(rrax,stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, l_err=(/1_psb_lpk_/), &
             & a_err='integer(psb_mpk_)')
        goto 9999
      end if
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_r_m_s

  Subroutine psb_r_m_m_rk1(len,rrax,info,pad,lb)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_mpk_),Intent(in) :: len
    integer(psb_mpk_), allocatable, intent(inout) :: rrax(:)
    integer(psb_ipk_) :: info
    integer(psb_mpk_), optional, intent(in) :: pad
    integer(psb_mpk_), optional, intent(in) :: lb

    ! ...Local Variables
    integer(psb_mpk_),allocatable  :: tmp(:)
    integer(psb_mpk_) :: dim, lb_, lbi,ub_, i
    integer(psb_ipk_) :: err_act,err
    character(len=30)  :: name
    logical, parameter :: debug=.false.

    name='psb_r_m_m_rk1'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (debug) write(psb_err_unit,*) 'reallocate D',len

    if (present(lb)) then
      lb_ = lb
    else
      lb_ = 1
    endif
    if ((len<0)) then 
      err=4025
      call psb_errpush(err,name, l_err=(/len*1_psb_lpk_/), &
           & a_err='integer(psb_mpk_)')
      goto 9999
    end if
    ub_ = lb_ + len-1

    if (allocated(rrax)) then 
      dim = size(rrax)
      lbi = lbound(rrax,1)
      If ((dim /= len).or.(lbi /= lb_))  Then
        Allocate(tmp(lb_:ub_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, l_err=(/len*1_psb_lpk_/), &
               & a_err='integer(psb_mpk_)')
          goto 9999
        end if
        tmp(lb_:lb_-1+min(len,dim))=rrax(lbi:lbi-1+min(len,dim))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim = 0
      Allocate(rrax(lb_:ub_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, l_err=(/len*1_psb_lpk_/), &
             & a_err='integer(psb_mpk_)')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      !$omp parallel do private(i) shared(dim,len)
      do i=lb_-1+dim+1,lb_-1+len
        rrax(i) = pad
      end do
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_r_m_m_rk1

  Subroutine psb_r_m_m_rk2(len1,len2,rrax,info,pad,lb1,lb2)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_mpk_),Intent(in) :: len1,len2
    integer(psb_mpk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    integer(psb_mpk_), optional, intent(in) :: pad
    integer(psb_mpk_),Intent(in), optional  :: lb1,lb2

    ! ...Local Variables

    integer(psb_mpk_),allocatable  :: tmp(:,:)
    integer(psb_ipk_) :: err_act,err
    integer(psb_mpk_) :: dim,dim2,lb1_, lb2_, ub1_, ub2_,lbi1, lbi2, i
    character(len=30)  :: name

    name='psb_r_m_m_rk2'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (present(lb1)) then 
      lb1_ = lb1
    else
      lb1_ = 1
    endif
    if (present(lb2)) then 
      lb2_ = lb2
    else
      lb2_ = 1
    endif
    ub1_ = lb1_ + len1 -1
    ub2_ = lb2_ + len2 -1

    if (len1 < 0) then
      err=4025 
      call psb_errpush(err,name, l_err=(/len1*1_psb_lpk_/), &
           & a_err='integer(psb_mpk_)')
      goto 9999
    end if
    if (len2 < 0) then
      err=4025
      call psb_errpush(err,name, l_err=(/len2*1_psb_lpk_/), &
           & a_err='integer(psb_mpk_)')
      goto 9999
    end if


    if (allocated(rrax)) then 
      dim  = size(rrax,1)
      lbi1 = lbound(rrax,1)
      dim2 = size(rrax,2)
      lbi2 = lbound(rrax,2)
      If ((dim /= len1).or.(dim2 /= len2).or.(lbi1 /= lb1_)&
           &  .or.(lbi2 /= lb2_)) Then
        Allocate(tmp(lb1_:ub1_,lb2_:ub2_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, l_err=(/len1*1_psb_lpk_*len2/), &
               & a_err='integer(psb_mpk_)')
          goto 9999
        end if
        tmp(lb1_:lb1_-1+min(len1,dim),lb2_:lb2_-1+min(len2,dim2)) = &
             & rrax(lbi1:lbi1-1+min(len1,dim),lbi2:lbi2-1+min(len2,dim2))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(lb1_:ub1_,lb2_:ub2_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, l_err=(/len1*1_psb_lpk_*len2/), &
             & a_err='integer(psb_mpk_)')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      !$omp parallel do private(i) shared(lb1_,dim,len1)
      do i=lb1_-1+dim+1,lb1_-1+len1
        rrax(i,:) = pad
      end do
      !$omp parallel do private(i) shared(lb1_,dim,len1,lb2_,dim2,len2)
      do i=lb1_,lb1_-1+len1
        rrax(i,lb2_-1+dim2+1:lb2_-1+len2) = pad
      end do
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_r_m_m_rk2


  Subroutine psb_r_e_m_rk1(len,rrax,info,pad,lb)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_epk_),Intent(in) :: len
    integer(psb_mpk_), allocatable, intent(inout) :: rrax(:)
    integer(psb_ipk_) :: info
    integer(psb_mpk_), optional, intent(in) :: pad
    integer(psb_epk_), optional, intent(in) :: lb

    ! ...Local Variables
    integer(psb_mpk_),allocatable  :: tmp(:)
    integer(psb_epk_) :: dim, lb_, lbi,ub_
    integer(psb_ipk_) :: err_act,err
    character(len=30)  :: name
    logical, parameter :: debug=.false.

    name='psb_r_m_m_rk1'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (debug) write(psb_err_unit,*) 'reallocate D',len

    if (present(lb)) then
      lb_ = lb
    else
      lb_ = 1
    endif
    if ((len<0)) then 
      err=4025
      call psb_errpush(err,name, e_err=(/len/), &
           & a_err='integer(psb_mpk_)')
      goto 9999
    end if
    ub_ = lb_ + len-1

    if (allocated(rrax)) then 
      dim = size(rrax)
      lbi = lbound(rrax,1)
      If ((dim /= len).or.(lbi /= lb_))  Then
        Allocate(tmp(lb_:ub_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, e_err=(/len/), &
               & a_err='integer(psb_mpk_)')
          goto 9999
        end if
        tmp(lb_:lb_-1+min(len,dim))=rrax(lbi:lbi-1+min(len,dim))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim = 0
      Allocate(rrax(lb_:ub_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, e_err=(/len/), &
             & a_err='integer(psb_mpk_)')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb_-1+dim+1:lb_-1+len) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_r_e_m_rk1

  Subroutine psb_r_e_m_rk2(len1,len2,rrax,info,pad,lb1,lb2)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_epk_),Intent(in) :: len1,len2
    integer(psb_mpk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    integer(psb_mpk_), optional, intent(in) :: pad
    integer(psb_epk_),Intent(in), optional  :: lb1,lb2

    ! ...Local Variables

    integer(psb_mpk_),allocatable  :: tmp(:,:)
    integer(psb_ipk_) :: err_act,err
    integer(psb_epk_) :: dim,dim2,lb1_, lb2_, ub1_, ub2_,lbi1, lbi2
    character(len=30)  :: name

    name='psb_r_e_m_rk2'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (present(lb1)) then 
      lb1_ = lb1
    else
      lb1_ = 1
    endif
    if (present(lb2)) then 
      lb2_ = lb2
    else
      lb2_ = 1
    endif
    ub1_ = lb1_ + len1 -1
    ub2_ = lb2_ + len2 -1

    if (len1 < 0) then
      err=4025
      call psb_errpush(err,name, e_err=(/len1/), &
           & a_err='integer(psb_mpk_)')
      goto 9999
    end if
    if (len2 < 0) then
      err=4025
      call psb_errpush(err,name, e_err=(/len2/), &
           & a_err='integer(psb_mpk_)')
      goto 9999
    end if


    if (allocated(rrax)) then 
      dim  = size(rrax,1)
      lbi1 = lbound(rrax,1)
      dim2 = size(rrax,2)
      lbi2 = lbound(rrax,2)
      If ((dim /= len1).or.(dim2 /= len2).or.(lbi1 /= lb1_)&
           &  .or.(lbi2 /= lb2_)) Then
        Allocate(tmp(lb1_:ub1_,lb2_:ub2_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, e_err=(/(len1*len2)/), &
               & a_err='integer(psb_mpk_)')
          goto 9999
        end if
        tmp(lb1_:lb1_-1+min(len1,dim),lb2_:lb2_-1+min(len2,dim2)) = &
             & rrax(lbi1:lbi1-1+min(len1,dim),lbi2:lbi2-1+min(len2,dim2))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(lb1_:ub1_,lb2_:ub2_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, e_err=(/(len1*len2)/), &
             & a_err='integer(psb_mpk_)')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb1_-1+dim+1:lb1_-1+len1,:) = pad
      rrax(lb1_:lb1_-1+dim,lb2_-1+dim2+1:lb2_-1+len2) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_r_e_m_rk2

  Subroutine psb_r_me_m_rk2(len1,len2,rrax,info,pad,lb1,lb2)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_mpk_),Intent(in) :: len1
    integer(psb_epk_),Intent(in) :: len2
    integer(psb_mpk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    integer(psb_mpk_), optional, intent(in) :: pad
    integer(psb_mpk_),Intent(in), optional  :: lb1,lb2

    ! ...Local Variables

    integer(psb_mpk_),allocatable  :: tmp(:,:)
    integer(psb_ipk_) :: err_act,err
    integer(psb_epk_) :: dim,lb1_, lb2_, ub1_, ub2_,lbi1, lbi2
    integer(psb_epk_) :: dim2
    character(len=30)  :: name

    name='psb_r_me_m_rk2'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (present(lb1)) then 
      lb1_ = lb1
    else
      lb1_ = 1
    endif
    if (present(lb2)) then 
      lb2_ = lb2
    else
      lb2_ = 1
    endif
    ub1_ = lb1_ + len1 -1
    ub2_ = lb2_ + len2 -1

    if (len1 < 0) then
      err=4025 
      call psb_errpush(err,name, m_err=(/len1/), &
           & a_err='integer(psb_mpk_)')
      goto 9999
    end if
    if (len2 < 0) then
      err=4025
      call psb_errpush(err,name, e_err=(/len2/), &
           & a_err='integer(psb_mpk_)')
      goto 9999
    end if


    if (allocated(rrax)) then 
      dim  = size(rrax,1)
      lbi1 = lbound(rrax,1)
      dim2 = size(rrax,2)
      lbi2 = lbound(rrax,2)
      If ((dim /= len1).or.(dim2 /= len2).or.(lbi1 /= lb1_)&
           &  .or.(lbi2 /= lb2_)) Then
        Allocate(tmp(lb1_:ub1_,lb2_:ub2_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, e_err=(/len1*len2/), &
               & a_err='integer(psb_mpk_)')
          goto 9999
        end if
        tmp(lb1_:lb1_-1+min(len1,dim),lb2_:lb2_-1+min(len2,dim2)) = &
             & rrax(lbi1:lbi1-1+min(len1,dim),lbi2:lbi2-1+min(len2,dim2))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(lb1_:ub1_,lb2_:ub2_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name,e_err=(/len1*len2/),&
             &  a_err='integer(psb_mpk_)')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb1_-1+dim+1:lb1_-1+len1,:) = pad
      rrax(lb1_:lb1_-1+dim,lb2_-1+dim2+1:lb2_-1+len2) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_r_me_m_rk2

  Subroutine psb_r_em_m_rk2(len1,len2,rrax,info,pad,lb1,lb2)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_epk_),Intent(in) :: len1
    integer(psb_mpk_),Intent(in) :: len2
    integer(psb_mpk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    integer(psb_mpk_), optional, intent(in) :: pad
    integer(psb_mpk_),Intent(in), optional  :: lb1,lb2

    ! ...Local Variables

    integer(psb_mpk_),allocatable  :: tmp(:,:)
    integer(psb_ipk_) :: err_act,err
    integer(psb_epk_) :: dim2,lb1_, lb2_, ub1_, ub2_,lbi1, lbi2
    integer(psb_epk_) :: dim
    character(len=30)  :: name

    name='psb_r_me_m_rk2'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (present(lb1)) then 
      lb1_ = lb1
    else
      lb1_ = 1
    endif
    if (present(lb2)) then 
      lb2_ = lb2
    else
      lb2_ = 1
    endif
    ub1_ = lb1_ + len1 -1
    ub2_ = lb2_ + len2 -1

    if (len1 < 0) then
      err=4025 
      call psb_errpush(err,name, e_err=(/len1/), &
           & a_err='integer(psb_mpk_)')
      goto 9999
    end if
    if (len2 < 0) then
      err=4025
      call psb_errpush(err,name, m_err=(/len2/), &
           & a_err='integer(psb_mpk_)')
      goto 9999
    end if


    if (allocated(rrax)) then 
      dim  = size(rrax,1)
      lbi1 = lbound(rrax,1)
      dim2 = size(rrax,2)
      lbi2 = lbound(rrax,2)
      If ((dim /= len1).or.(dim2 /= len2).or.(lbi1 /= lb1_)&
           &  .or.(lbi2 /= lb2_)) Then
        Allocate(tmp(lb1_:ub1_,lb2_:ub2_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, e_err=(/len1*len2/), &
               & a_err='integer(psb_mpk_)')
          goto 9999
        end if
        tmp(lb1_:lb1_-1+min(len1,dim),lb2_:lb2_-1+min(len2,dim2)) = &
             & rrax(lbi1:lbi1-1+min(len1,dim),lbi2:lbi2-1+min(len2,dim2))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(lb1_:ub1_,lb2_:ub2_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, e_err=(/len1*len2/), &
             & a_err='integer(psb_mpk_)')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb1_-1+dim+1:lb1_-1+len1,:) = pad
      rrax(lb1_:lb1_-1+dim,lb2_-1+dim2+1:lb2_-1+len2) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_r_em_m_rk2

  Subroutine psb_r_m_2_m_rk1(len,rrax,y,info,pad)
    use psb_error_mod
    ! ...Subroutine Arguments  

    integer(psb_mpk_),Intent(in) :: len
    integer(psb_mpk_),allocatable, intent(inout) :: rrax(:),y(:)
    integer(psb_ipk_) :: info
    integer(psb_mpk_), optional, intent(in) :: pad
    character(len=30)  :: name
    integer(psb_ipk_) :: err_act, err

    name='psb_r_m_2_m_rk1'
    call psb_erractionsave(err_act)
    info=psb_success_

    if(psb_errstatus_fatal()) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    call psb_realloc(len,rrax,info,pad=pad)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,y,info,pad=pad)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  End Subroutine psb_r_m_2_m_rk1

  Subroutine psb_r_e_2_m_rk1(len,rrax,y,info,pad)
    use psb_error_mod
    ! ...Subroutine Arguments  

    integer(psb_epk_),Intent(in) :: len
    integer(psb_mpk_),allocatable, intent(inout) :: rrax(:),y(:)
    integer(psb_ipk_) :: info
    integer(psb_mpk_), optional, intent(in) :: pad
    character(len=30)  :: name
    integer(psb_ipk_) :: err_act, err

    name='psb_r_m_2_m_rk1'
    call psb_erractionsave(err_act)
    info=psb_success_

    if(psb_errstatus_fatal()) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    call psb_realloc(len,rrax,info,pad=pad)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,y,info,pad=pad)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  End Subroutine psb_r_e_2_m_rk1



  subroutine psb_ab_cpy_m_s(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_mpk_), allocatable, intent(in)  :: vin
    integer(psb_mpk_), allocatable, intent(out) :: vout
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: err_act
    character(len=30)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_ab_cpy_m_s'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_errstatus_fatal()) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    if (allocated(vin)) then 
      call psb_realloc(vout,info)
      if (info /= psb_success_) then     
        info=psb_err_from_subroutine_
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout = vin
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_ab_cpy_m_s

  subroutine psb_ab_cpy_m_rk1(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_mpk_), allocatable, intent(in)  :: vin(:)
    integer(psb_mpk_), allocatable, intent(out) :: vout(:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz,err_act,lb
    character(len=30)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_ab_cpy_m_rk1'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_errstatus_fatal()) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    if (allocated(vin)) then 
      isz = size(vin)
      lb  = lbound(vin,1)
      call psb_realloc(isz,vout,info,lb=lb)
      if (info /= psb_success_) then     
        info=psb_err_from_subroutine_
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        !$omp workshare
        vout(:) = vin(:)
        !$omp end workshare
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_ab_cpy_m_rk1

  subroutine psb_ab_cpy_m_rk2(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_mpk_), allocatable, intent(in)  :: vin(:,:)
    integer(psb_mpk_), allocatable, intent(out) :: vout(:,:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz1, isz2,err_act, lb1, lb2 
    character(len=30)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_ab_cpy_m_rk2'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_errstatus_fatal()) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    if (allocated(vin)) then 
      isz1 = size(vin,1)
      isz2 = size(vin,2)
      lb1  = lbound(vin,1)
      lb2  = lbound(vin,2)
      call psb_realloc(isz1,isz2,vout,info,lb1=lb1,lb2=lb2)
      if (info /= psb_success_) then     
        info=psb_err_from_subroutine_
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        !$omp workshare
        vout(:,:) = vin(:,:)
        !$omp end workshare
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_ab_cpy_m_rk2


  subroutine psb_cpy_m_rk1(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_mpk_), intent(in)               :: vin(:)
    integer(psb_mpk_), allocatable, intent(out) :: vout(:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz,err_act,lb
    character(len=30)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_cpy_m_rk1'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_errstatus_fatal()) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if
    isz = size(vin)
    lb  = lbound(vin,1)
    call psb_realloc(isz,vout,info,lb=lb)
    if (info /= psb_success_) then     
      info=psb_err_from_subroutine_
      char_err='psb_realloc'
      call psb_errpush(info,name,a_err=char_err)
      goto 9999
    else
      vout(:) = vin(:)
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_cpy_m_rk1

  subroutine psb_cpy_m_rk2(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_mpk_), intent(in)               :: vin(:,:)
    integer(psb_mpk_), allocatable, intent(out) :: vout(:,:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz1, isz2,err_act, lb1, lb2 
    character(len=30)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_errstatus_fatal()) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    isz1 = size(vin,1)
    isz2 = size(vin,2)
    lb1  = lbound(vin,1)
    lb2  = lbound(vin,2)
    call psb_realloc(isz1,isz2,vout,info,lb1=lb1,lb2=lb2)
    if (info /= psb_success_) then     
      info=psb_err_from_subroutine_
      char_err='psb_realloc'
      call psb_errpush(info,name,a_err=char_err)
      goto 9999
    else
      vout(:,:) = vin(:,:)
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_cpy_m_rk2


  function psb_size_m_rk1(vin) result(val)
    integer(psb_epk_) :: val
    integer(psb_mpk_), allocatable, intent(in) :: vin(:)

    if (.not.allocated(vin)) then 
      val = 0
    else
      val = size(vin)
    end if
  end function psb_size_m_rk1


  function psb_size_m_rk2(vin,dim) result(val)
    integer(psb_epk_) :: val
    integer(psb_mpk_), allocatable, intent(in) :: vin(:,:)
    integer(psb_ipk_), optional :: dim
    integer(psb_ipk_) :: dim_


    if (.not.allocated(vin)) then 
      val = 0
    else
      if (present(dim)) then 
        dim_= dim
        val = size(vin,dim=dim_)
      else
        val = size(vin)
      end if
    end if
  end function psb_size_m_rk2

  Subroutine psb_ensure_m_sz_m_rk1(len,v,info,pad,addsz,newsz)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_mpk_),Intent(in)                 :: len
    integer(psb_mpk_),allocatable, intent(inout) :: v(:)
    integer(psb_ipk_) :: info
    integer(psb_mpk_), optional, intent(in)          :: addsz,newsz
    integer(psb_mpk_), optional, intent(in) :: pad
    ! ...Local Variables
    character(len=30)  :: name
    logical, parameter :: debug=.false.
    integer(psb_ipk_) :: err_act
    integer(psb_mpk_) :: isz

    name='psb_ensure_m_sz_m_rk1'
    call psb_erractionsave(err_act)
    info = psb_success_

    if (psb_errstatus_fatal()) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

!!$    If (len > psb_size(v)) Then
!!$      if (present(newsz)) then 
!!$        isz = (max(len+1,newsz))
!!$      else
!!$        if (present(addsz)) then 
!!$          isz = len+max(1,addsz)
!!$        else
!!$          isz = max(len+10, int(1.25*len))
!!$        endif
!!$      endif
!!$
!!$      call psb_realloc(isz,v,info,pad=pad)
!!$      if (info /= psb_success_) then
!!$        info=psb_err_from_subroutine_
!!$        call psb_errpush(info,name,a_err='psb_realloc')
!!$        goto 9999
!!$      End If
!!$    end If
    isz = psb_size(v)
    If (len > isz) Then
#if defined(OPENMP)
      !$OMP CRITICAL
      if (len > isz) then
        if (present(newsz)) then
          isz = max(len+1,1,newsz)
        else if (present(addsz)) then
          isz = max(len,1,isz+addsz)
        else
          isz = max(len,1,int(1.25*isz))
        endif

        call psb_realloc(isz,v,info,pad=pad)
      end if
      !$OMP END CRITICAL

      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_realloc')
	goto 9999
      end if
#else
      if (len > isz) then
        if (present(newsz)) then
          isz = max(len+1,1,newsz)
        else if (present(addsz)) then
          isz = max(len,1,isz+addsz)
        else
          isz = max(len,1,int(1.25*isz))
        endif

        call psb_realloc(isz,v,info,pad=pad)
      end if

      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_realloc')
        goto 9999
      End If
#endif
    end If

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return


  End Subroutine psb_ensure_m_sz_m_rk1

  Subroutine psb_ensure_e_sz_m_rk1(len,v,info,pad,addsz,newsz)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_epk_),Intent(in)                 :: len
    integer(psb_mpk_),allocatable, intent(inout) :: v(:)
    integer(psb_ipk_) :: info
    integer(psb_epk_), optional, intent(in)          :: addsz,newsz
    integer(psb_mpk_), optional, intent(in) :: pad
    ! ...Local Variables
    character(len=30)  :: name
    logical, parameter :: debug=.false.
    integer(psb_ipk_) :: err_act
    integer(psb_epk_) :: isz

    name='psb_ensure_m_sz_m_rk1'
    call psb_erractionsave(err_act)
    info = psb_success_

    if (psb_errstatus_fatal()) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if
    isz = psb_size(v)
    If (len > isz) Then
      if (present(newsz)) then
        isz = max(len+1,1,newsz)
      else if (present(addsz)) then
        isz = max(len,1,isz+addsz)
      else
        isz = max(len,1,int(1.25*isz))
      endif

      call psb_realloc(isz,v,info,pad=pad)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_realloc')
        goto 9999
      End If
    end If

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return


  End Subroutine psb_ensure_e_sz_m_rk1

  Subroutine psb_move_alloc_m_rk1(vin,vout,info)
    use psb_error_mod
    integer(psb_mpk_), allocatable, intent(inout) :: vin(:),vout(:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_
    call move_alloc(vin,vout)

  end Subroutine psb_move_alloc_m_rk1

  Subroutine psb_move_alloc_m_rk2(vin,vout,info)
    use psb_error_mod
    integer(psb_mpk_), allocatable, intent(inout) :: vin(:,:),vout(:,:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_

    call move_alloc(vin,vout)

  end Subroutine psb_move_alloc_m_rk2

end module psb_m_realloc_mod
