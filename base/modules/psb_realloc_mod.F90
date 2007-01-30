!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
module psb_realloc_mod
  use psb_const_mod
  implicit none
  
  Interface psb_realloc
    module procedure psb_dreallocate1i
    module procedure psb_dreallocate2i
    module procedure psb_dreallocate2i1d
    module procedure psb_dreallocate1d
    module procedure psb_dreallocated2
    module procedure psb_dreallocatei2
    module procedure psb_dreallocate2i1z
    module procedure psb_dreallocate1z
    module procedure psb_dreallocatez2
  end Interface

  interface psb_transfer
    module procedure psb_dtransfer1d
    module procedure psb_dtransfer2d
    module procedure psb_itransfer1d
    module procedure psb_itransfer2d
    module procedure psb_ztransfer1d
    module procedure psb_ztransfer2d
  end interface

  Interface psb_safe_cpy
    module procedure psb_icpy1d,psb_icpy2d, &
         & psb_dcpy1d, psb_dcpy2d, psb_zcpy1d, psb_zcpy2d
  end Interface

  Interface psb_check_size
    module procedure psb_icksz1d, psb_dcksz1d, psb_zcksz1d
  end Interface

  interface psb_size
    module procedure psb_isize1d, psb_isize2d,&
         & psb_dsize1d, psb_dsize2d,&
         & psb_zsize1d, psb_zsize2d
  end interface
  
  
Contains

  subroutine psb_icpy1d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    Integer,allocatable, intent(in)  :: vin(:)
    Integer,allocatable, intent(out) :: vout(:)
    integer         :: info
    ! ...Local Variables

    Integer :: isz,err_act,lb
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_cpy1d'
    call psb_erractionsave(err_act)

    if(psb_get_errstatus().ne.0) return 
    info = 0
    if (allocated(vin)) then 
      isz = size(vin)
      lb  = lbound(vin,1)
      call psb_realloc(isz,vout,info,lb=lb)
      if (info /= 0) then     
        info=4010
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:) = vin(:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  end subroutine psb_icpy1d

  subroutine psb_icpy2d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    Integer, allocatable, intent(in)  :: vin(:,:)
    Integer, allocatable, intent(out) :: vout(:,:)
    integer         :: info
    ! ...Local Variables

    Integer :: isz1, isz2,err_act, lb1, lb2
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_cpy1d'
    call psb_erractionsave(err_act)

    if(psb_get_errstatus().ne.0) return 
    info = 0
    if (allocated(vin)) then 
      isz1 = size(vin,1)
      isz2 = size(vin,2)
      lb1  = lbound(vin,1)
      lb2  = lbound(vin,2)
      call psb_realloc(isz1,isz2,vout,info,lb1=lb1,lb2=lb2)
      if (info /= 0) then     
        info=4010
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:,:) = vin(:,:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  end subroutine psb_icpy2d
  
  subroutine psb_dcpy1d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    real(kind(1.d0)), allocatable, intent(in)  :: vin(:)
    real(kind(1.d0)), allocatable, intent(out) :: vout(:)
    integer         :: info
    ! ...Local Variables

    Integer :: isz,err_act,lb
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_cpy1d'
    call psb_erractionsave(err_act)

    if(psb_get_errstatus().ne.0) return 
    info = 0
    if (allocated(vin)) then 
      isz = size(vin)
      lb  = lbound(vin,1)
      call psb_realloc(isz,vout,info,lb=lb)
      if (info /= 0) then     
        info=4010
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:) = vin(:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  end subroutine psb_dcpy1d
  
  subroutine psb_dcpy2d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    real(kind(1.d0)), allocatable, intent(in)  :: vin(:,:)
    real(kind(1.d0)), allocatable, intent(out) :: vout(:,:)
    integer         :: info
    ! ...Local Variables

    Integer :: isz1, isz2,err_act, lb1, lb2 
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_cpy1d'
    call psb_erractionsave(err_act)

    if(psb_get_errstatus().ne.0) return 
    info = 0
    if (allocated(vin)) then 
      isz1 = size(vin,1)
      isz2 = size(vin,2)
      lb1  = lbound(vin,1)
      lb2  = lbound(vin,2)
      call psb_realloc(isz1,isz2,vout,info,lb1=lb1,lb2=lb2)
      if (info /= 0) then     
        info=4010
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:,:) = vin(:,:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  end subroutine psb_dcpy2d
  
  subroutine psb_zcpy1d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    complex(kind(1.d0)), allocatable, intent(in)  :: vin(:)
    complex(kind(1.d0)), allocatable, intent(out) :: vout(:)
    integer         :: info
    ! ...Local Variables

    Integer :: isz,err_act,lb
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_cpy1d'
    call psb_erractionsave(err_act)

    if(psb_get_errstatus().ne.0) return 
    info = 0
    if (allocated(vin)) then 
      isz = size(vin)
      lb  = lbound(vin,1)
      call psb_realloc(isz,vout,info,lb=lb)
      if (info /= 0) then     
        info=4010
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:) = vin(:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  end subroutine psb_zcpy1d
  
  subroutine psb_zcpy2d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    complex(kind(1.d0)), allocatable, intent(in)  :: vin(:,:)
    complex(kind(1.d0)), allocatable, intent(out) :: vout(:,:)
    integer         :: info
    ! ...Local Variables

    Integer :: isz1, isz2,err_act, lb1, lb2 
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_cpy1d'
    call psb_erractionsave(err_act)

    if(psb_get_errstatus().ne.0) return 
    info = 0
    if (allocated(vin)) then 
      isz1 = size(vin,1)
      isz2 = size(vin,2)
      lb1  = lbound(vin,1)
      lb2  = lbound(vin,2)
      call psb_realloc(isz1,isz2,vout,info,lb1=lb1,lb2=lb2)
      if (info /= 0) then     
        info=4010
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:,:) = vin(:,:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  end subroutine psb_zcpy2d
  
  function psb_isize1d(vin)
    integer :: psb_isize1d
    integer, allocatable, intent(in) :: vin(:)
    
    if (.not.allocated(vin)) then 
      psb_isize1d = 0
    else
      psb_isize1d = size(vin)
    end if
  end function psb_isize1d
  function psb_isize2d(vin,dim)
    integer :: psb_isize2d
    integer, allocatable, intent(in) :: vin(:,:)
    integer, optional :: dim

    if (.not.allocated(vin)) then 
      psb_isize2d = 0
    else
      if (present(dim)) then 
        psb_isize2d = size(vin,dim=dim)
      else
        psb_isize2d = size(vin)
      end if
    end if
  end function psb_isize2d
  
  function psb_dsize1d(vin)
    integer :: psb_dsize1d
    real(kind(1.d0)), allocatable, intent(in) :: vin(:)
    
    if (.not.allocated(vin)) then 
      psb_dsize1d = 0
    else
      psb_dsize1d = size(vin)
    end if
  end function psb_dsize1d
  function psb_dsize2d(vin,dim)
    integer :: psb_dsize2d
    real(kind(1.d0)), allocatable, intent(in) :: vin(:,:)
    integer, optional :: dim

    if (.not.allocated(vin)) then 
      psb_dsize2d = 0
    else
      if (present(dim)) then 
        psb_dsize2d = size(vin,dim=dim)
      else
        psb_dsize2d = size(vin)
      end if
    end if
  end function psb_dsize2d

  
  function psb_zsize1d(vin)
    integer :: psb_zsize1d
    complex(kind(1.d0)), allocatable, intent(in) :: vin(:)
    
    if (.not.allocated(vin)) then 
      psb_zsize1d = 0
    else
      psb_zsize1d = size(vin)
    end if
  end function psb_zsize1d

  function psb_zsize2d(vin,dim)
    integer :: psb_zsize2d
    complex(kind(1.d0)), allocatable, intent(in) :: vin(:,:)
    integer, optional :: dim

    if (.not.allocated(vin)) then 
      psb_zsize2d = 0
    else
      if (present(dim)) then 
        psb_zsize2d = size(vin,dim=dim)
      else
        psb_zsize2d = size(vin)
      end if
    end if
  end function psb_zsize2d


  Subroutine psb_icksz1d(len,v,info,pad)
    use psb_error_mod

    ! ...Subroutine Arguments  
    Integer,Intent(in)                 :: len
    Integer,allocatable, intent(inout) :: v(:)
    integer         :: info
    integer, optional, intent(in) :: pad
    ! ...Local Variables
    character(len=20)  :: name
    logical, parameter :: debug=.false.
    integer :: isz, err_act

    name='psb_check_size'
    call psb_erractionsave(err_act)

    if(psb_get_errstatus().ne.0) return 
    info=0
    
    If (len > psb_size(v)) Then
      isz = max((3*psb_size(v))/2,(len+1))
      if (present(pad)) then
        call psb_realloc(isz,v,info,pad=pad)
      else
        call psb_realloc(isz,v,info)
        if (info /= 0) then
          info=4010
          call psb_errpush(info,name,a_err='psb_realloc')
          goto 9999
        end if
      End If
    end If

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return


  End Subroutine psb_icksz1d


  Subroutine psb_dcksz1d(len,v,info,pad)
    use psb_error_mod

    ! ...Subroutine Arguments  
    Integer,Intent(in)                 :: len
    real(kind(1.d0)),allocatable, intent(inout) :: v(:)
    integer         :: info
    real(kind(1.d0)), optional, intent(in) :: pad
    ! ...Local Variables
    character(len=20)  :: name
    logical, parameter :: debug=.false.
    integer :: isz, err_act

    name='psb_check_size'
    call psb_erractionsave(err_act)

    if(psb_get_errstatus().ne.0) return 
    info=0
    
    If (len > psb_size(v)) Then
      isz = max((3*psb_size(v))/2,(len+1))
      if (present(pad)) then
        call psb_realloc(isz,v,info,pad=pad)
      else
        call psb_realloc(isz,v,info)
        if (info /= 0) then
          info=4010
          call psb_errpush(info,name,a_err='psb_realloc')
          goto 9999
        end if
      End If
    end If

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return


  End Subroutine psb_dcksz1d


  Subroutine psb_zcksz1d(len,v,info,pad)
    use psb_error_mod

    ! ...Subroutine Arguments  
    Integer,Intent(in)                 :: len
    complex(kind(1.d0)),allocatable, intent(inout) :: v(:)
    integer         :: info
    complex(kind(1.d0)), optional, intent(in) :: pad
    ! ...Local Variables
    character(len=20)  :: name
    logical, parameter :: debug=.false.
    integer :: isz, err_act

    name='psb_check_size'
    call psb_erractionsave(err_act)

    if(psb_get_errstatus().ne.0) return 
    info=0
    
    If (len > psb_size(v)) Then
      isz = max((3*psb_size(v))/2,(len+1))
      if (present(pad)) then
        call psb_realloc(isz,v,info,pad=pad)
      else
        call psb_realloc(isz,v,info)
        if (info /= 0) then
          info=4010
          call psb_errpush(info,name,a_err='psb_realloc')
          goto 9999
        end if
      End If
    end If

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return


  End Subroutine psb_zcksz1d


  Subroutine psb_dreallocate1i(len,rrax,info,pad,lb)
    use psb_error_mod

    ! ...Subroutine Arguments  
    Integer,Intent(in)                 :: len
    Integer,allocatable, intent(inout) :: rrax(:)
    integer         :: info
    integer, optional, intent(in) :: pad
    integer, optional, intent(in) :: lb
    ! ...Local Variables
    Integer,allocatable  :: tmp(:)
    Integer :: dim, err_act, err,i,lb_, lbi, ub_
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_dreallocate1i' 
    call psb_erractionsave(err_act)

    if (debug) write(0,*) 'reallocate I',len
    if (psb_get_errstatus().ne.0) return 
    info=0
    if (present(lb)) then
      lb_ = lb
    else
      lb_ = 1
    endif
    if ((len<0).or.(len>25*1024*1024)) then 
      err=2025
      call psb_errpush(err,name,i_err=(/len,0,0,0,0/))
      goto 9999
    end if
    ub_ = lb_+len-1

    if (allocated(rrax)) then 
      dim = size(rrax)
      lbi = lbound(rrax,1) 
      If ((dim /= len).or.(lbi /= lb_))  Then
        Allocate(tmp(lb_:ub_),stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        tmp(lb_:lb_-1+min(len,dim))=rrax(lbi:lbi-1+min(len,dim))
        call psb_transfer(tmp,rrax,info)
      end if
    else
      dim = 0
      allocate(rrax(lb_:ub_),stat=info)
      if (info /= 0) then
        err=4000
        call psb_errpush(err,name)
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb_-1+dim+1:lb_-1+len) = pad
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return


  End Subroutine psb_dreallocate1i


  Subroutine psb_dreallocate1d(len,rrax,info,pad,lb)
    use psb_error_mod

    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len
    Real(kind(1.d0)),allocatable, intent(inout) :: rrax(:)
    integer :: info
    real(kind(1.d0)), optional, intent(in) :: pad
    integer, optional, intent(in) :: lb

    ! ...Local Variables
    Real(kind(1.d0)),allocatable  :: tmp(:)
    Integer :: dim,err_act,err,m, lb_, lbi,ub_
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_dreallocate1d'
    call psb_erractionsave(err_act)
    info = 0 
    if (debug) write(0,*) 'reallocate D',len

    if (present(lb)) then
      lb_ = lb
    else
      lb_ = 1
    endif
    if ((len<0).or.(len>25*1024*1024)) then 
      err=2025
      call psb_errpush(err,name,i_err=(/len,0,0,0,0/))
      goto 9999
    end if
    ub_ = lb_ + len-1

    if (allocated(rrax)) then 
      dim = size(rrax)
      lbi = lbound(rrax,1)
      If ((dim /= len).or.(lbi /= lb_))  Then
        Allocate(tmp(lb_:ub_),stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        tmp(lb_:lb_-1+min(len,dim))=rrax(lbi:lbi-1+min(len,dim))
        call psb_transfer(tmp,rrax,info)
      End If
    else
      dim = 0
      Allocate(rrax(lb_:ub_),stat=info)
      if (info /= 0) then
        err=4000
        call psb_errpush(err,name)
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb_-1+dim+1:lb_-1+len) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  End Subroutine psb_dreallocate1d


  Subroutine psb_dreallocate1z(len,rrax,info,pad,lb)
    use psb_error_mod

    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len
    complex(kind(1.d0)),allocatable, intent(inout):: rrax(:)
    integer :: info
    complex(kind(1.d0)), optional, intent(in) :: pad
    integer, optional, intent(in) :: lb

    ! ...Local Variables
    complex(kind(1.d0)),allocatable  :: tmp(:)
    Integer :: dim,err_act,err,i,lb_,ub_,lbi
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_dreallocate1z'
    call psb_erractionsave(err_act)
    info = 0
    if (debug) write(0,*) 'reallocate Z',len    
    if (present(lb)) then
      lb_ = lb
    else
      lb_ = 1
    endif
    if ((len<0).or.(len>25*1024*1024)) then 
      err=2025
      call psb_errpush(err,name,i_err=(/len,0,0,0,0/))
      goto 9999
    end if
    ub_ = lb_+len-1

    if (allocated(rrax)) then 
      dim = size(rrax)
      lbi = lbound(rrax,1) 
      If ((dim /= len).or.(lbi /= lb_))  Then
        Allocate(tmp(lb_:ub_),stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        tmp(lb_:lb_-1+min(len,dim))=rrax(lbi:lbi-1+min(len,dim))
        call psb_transfer(tmp,rrax,info)
      end if
    else
      dim = 0
      Allocate(rrax(lb_:ub_),stat=info)
      if (info /= 0) then
        err=4000
        call psb_errpush(err,name)
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb_-1+dim+1:lb_-1+len) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  End Subroutine psb_dreallocate1z



  Subroutine psb_dreallocated2(len1,len2,rrax,info,pad,lb1,lb2)
    use psb_error_mod
    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len1,len2
    Real(kind(1.d0)),allocatable :: rrax(:,:)
    integer :: info
    real(kind(1.d0)), optional, intent(in) :: pad
    Integer,Intent(in), optional  :: lb1,lb2

    ! ...Local Variables

    Real(kind(1.d0)),allocatable  :: tmp(:,:)
    Integer :: dim,err_act,err,i, m, dim2,lb1_, lb2_, ub1_, ub2_,&
         & lbi1, lbi2
    character(len=20)  :: name

    name='psb_dreallocated2'
    call psb_erractionsave(err_act)
    info = 0 
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

    if (allocated(rrax)) then 
      dim  = size(rrax,1)
      lbi1 = lbound(rrax,1)
      dim2 = size(rrax,2)
      lbi2 = lbound(rrax,2)
      If ((dim /= len1).or.(dim2 /= len2).or.(lbi1 /= lb1_)&
           &  .or.(lbi2 /= lb2_)) Then
        Allocate(tmp(lb1_:ub1_,lb2_:ub2_),stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        tmp(lb1_:lb1_-1+min(len1,dim),lb2_:lb2_-1+min(len2,dim2)) = &
             & rrax(lbi1:lbi1-1+min(len1,dim),lbi2:lbi2-1+min(len2,dim2))
        call psb_transfer(tmp,rrax,info)
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(lb1_:ub1_,lb2_:ub2_),stat=info)
      if (info /= 0) then
        err=4000
        call psb_errpush(err,name)
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
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  End Subroutine psb_dreallocated2



  Subroutine psb_dreallocatez2(len1,len2,rrax,info,pad,lb1,lb2)
    use psb_error_mod
    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len1,len2
    complex(kind(1.d0)),allocatable :: rrax(:,:)
    integer :: info
    complex(kind(1.d0)), optional, intent(in) :: pad
    Integer,Intent(in), optional  :: lb1,lb2

    ! ...Local Variables

    complex(kind(1.d0)),allocatable  :: tmp(:,:)
    Integer :: dim,err_act,err,i, m, dim2,lb1_, lb2_, ub1_, ub2_,&
         & lbi1, lbi2
    character(len=20)  :: name

    name='psb_dreallocatez2'
    call psb_erractionsave(err_act)
    info = 0 
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

    if (allocated(rrax)) then 
      dim  = size(rrax,1)
      lbi1 = lbound(rrax,1)
      dim2 = size(rrax,2)
      lbi2 = lbound(rrax,2)
      If ((dim /= len1).or.(dim2 /= len2).or.(lbi1 /= lb1_)&
           &  .or.(lbi2 /= lb2_)) Then
        Allocate(tmp(lb1_:ub1_,lb2_:ub2_),stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        tmp(lb1_:lb1_-1+min(len1,dim),lb2_:lb2_-1+min(len2,dim2)) = &
             & rrax(lbi1:lbi1-1+min(len1,dim),lbi2:lbi2-1+min(len2,dim2))
        call psb_transfer(tmp,rrax,info)
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(lb1_:ub1_,lb2_:ub2_),stat=info)
      if (info /= 0) then
        err=4000
        call psb_errpush(err,name)
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
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  End Subroutine psb_dreallocatez2


  Subroutine psb_dreallocatei2(len1,len2,rrax,info,pad,lb1,lb2)
    use psb_error_mod
    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len1,len2
    integer,allocatable :: rrax(:,:)
    integer :: info
    integer, optional, intent(in) :: pad
    Integer,Intent(in), optional  :: lb1,lb2

    ! ...Local Variables
    integer,allocatable  :: tmp(:,:)
    Integer :: dim,err_act,err,i, m, dim2,lb1_, lb2_, ub1_, ub2_,&
         & lbi1, lbi2
    character(len=20)  :: name

    name='psb_dreallocatei2'
    call psb_erractionsave(err_act)
    info = 0 
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

    if (allocated(rrax)) then 
      dim  = size(rrax,1)
      lbi1 = lbound(rrax,1)
      dim2 = size(rrax,2)
      lbi2 = lbound(rrax,2)
      If ((dim /= len1).or.(dim2 /= len2).or.(lbi1 /= lb1_)&
           &  .or.(lbi2 /= lb2_)) Then
        Allocate(tmp(lb1_:ub1_,lb2_:ub2_),stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        tmp(lb1_:lb1_-1+min(len1,dim),lb2_:lb2_-1+min(len2,dim2)) = &
             & rrax(lbi1:lbi1-1+min(len1,dim),lbi2:lbi2-1+min(len2,dim2))
        call psb_transfer(tmp,rrax,info)
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(lb1_:ub1_,lb2_:ub2_),stat=info)
      if (info /= 0) then
        err=4000
        call psb_errpush(err,name)
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
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  End Subroutine psb_dreallocatei2

  Subroutine psb_dreallocate2i(len,rrax,y,info,pad)
    use psb_error_mod
    ! ...Subroutine Arguments  

    Integer,Intent(in) :: len
    Integer,allocatable, intent(inout) :: rrax(:),y(:)
    integer :: info
    integer, optional, intent(in) :: pad
    character(len=20)  :: name
    integer :: err_act, err

    name='psb_dreallocate2i'
    call psb_erractionsave(err_act)

    if(psb_get_errstatus().ne.0) return 
    info=0
    call psb_dreallocate1i(len,rrax,info,pad=pad)
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_dreallocate1i(len,y,info,pad=pad)
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  End Subroutine psb_dreallocate2i




  Subroutine psb_dreallocate2i1d(len,rrax,y,z,info)
    use psb_error_mod
    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len
    Integer,allocatable, intent(inout)  :: rrax(:),y(:)
    Real(Kind(1.d0)),allocatable, intent(inout) :: z(:)
    integer :: info
    character(len=20)  :: name
    integer :: err_act, err

    name='psb_dreallocate2i1d'
    call psb_erractionsave(err_act)


    info = 0
    call psb_dreallocate1i(len,rrax,info)
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_dreallocate1i(len,y,info)    
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_dreallocate1d(len,z,info)
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return
  End Subroutine psb_dreallocate2i1d



  Subroutine psb_dreallocate2i1z(len,rrax,y,z,info)
    use psb_error_mod
    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len
    Integer,allocatable, intent(inout) :: rrax(:),y(:)
    complex(Kind(1.d0)),allocatable, intent(inout) :: z(:)
    integer :: info
    character(len=20)  :: name
    integer :: err_act, err

    name='psb_dreallocate2i1d'
    call psb_erractionsave(err_act)


    info = 0
    call psb_dreallocate1i(len,rrax,info)
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_dreallocate1i(len,y,info)    
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_dreallocate1z(len,z,info)
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return
  End Subroutine psb_dreallocate2i1z

  Subroutine psb_dtransfer1d(vin,vout,info)
    use psb_error_mod
    real(kind(1.d0)), allocatable, intent(inout) :: vin(:),vout(:)
    integer, intent(out) :: info 
    !
    ! 
    info = 0
#ifdef HAVE_MOVE_ALLOC    
    
    if (allocated(vin)) then 
      call move_alloc(vin,vout)
    else if (allocated(vout)) then 
      write(0,*) 'transfer: Clearing output'
      deallocate(vout)
    end if

#else      
    
    if (.not.allocated(vin) ) then 
      if (allocated(vout)) then 
        deallocate(vout,stat=info)
      end if
    else if (allocated(vin)) then 
      if (.not.allocated(vout)) then 
        allocate(vout(size(vin)),stat=info)
        if (info /= 0) return
      else
        if (size(vout) /= size(vin)) then 
          deallocate(vout,stat=info)
          if (info /= 0) return
          allocate(vout(size(vin)),stat=info)
          if (info /= 0) return
        end if
      end if
      vout = vin
      deallocate(vin,stat=info)
    end if
#endif
  end Subroutine psb_dtransfer1d

  Subroutine psb_dtransfer2d(vin,vout,info)
    use psb_error_mod
    real(kind(1.d0)), allocatable, intent(inout) :: vin(:,:),vout(:,:)
    integer, intent(out) :: info 
    !
    ! 
    info = 0
#ifdef HAVE_MOVE_ALLOC
    if (allocated(vin)) then 
      call move_alloc(vin,vout)
    else if (allocated(vout)) then 
      deallocate(vout)
    end if
#else
    
    if (.not.allocated(vin) ) then 
      if (allocated(vout)) then 
        deallocate(vout,stat=info)
      end if
    else if (allocated(vin)) then 
      if (.not.allocated(vout)) then 
        allocate(vout(size(vin,1),size(vin,2)),stat=info)
        if (info /= 0) return
      else
        if (size(vout) /= size(vin)) then 
          deallocate(vout,stat=info)
          if (info /= 0) return
          allocate(vout(size(vin,1),size(vin,2)),stat=info)
          if (info /= 0) return
        end if
      end if
      vout = vin
      deallocate(vin,stat=info)
    end if
#endif
  end Subroutine psb_dtransfer2d

  Subroutine psb_ztransfer1d(vin,vout,info)
    use psb_error_mod
    complex(kind(1.d0)), allocatable, intent(inout) :: vin(:),vout(:)
    integer, intent(out) :: info 
    !
    ! 
    info = 0
#ifdef HAVE_MOVE_ALLOC
    if (allocated(vin)) then 
      call move_alloc(vin,vout)
    else if (allocated(vout)) then 
      deallocate(vout)
    end if
#else
    if (.not.allocated(vin) ) then 
      if (allocated(vout)) then 
        deallocate(vout,stat=info)
      end if
    else if (allocated(vin)) then 
      if (.not.allocated(vout)) then 
        allocate(vout(size(vin)),stat=info)
        if (info /= 0) return
      else
        if (size(vout) /= size(vin)) then 
          deallocate(vout,stat=info)
          if (info /= 0) return
          allocate(vout(size(vin)),stat=info)
          if (info /= 0) return
        end if
      end if
      vout = vin
      deallocate(vin,stat=info)
    end if
#endif
  end Subroutine psb_ztransfer1d

  Subroutine psb_ztransfer2d(vin,vout,info)
    use psb_error_mod
    complex(kind(1.d0)), allocatable, intent(inout) :: vin(:,:),vout(:,:)
    integer, intent(out) :: info 
    !
    ! 
    info = 0
#ifdef HAVE_MOVE_ALLOC
    if (allocated(vin)) then 
      call move_alloc(vin,vout)
    else if (allocated(vout)) then 
      deallocate(vout)
    end if
#else
    if (.not.allocated(vin) ) then 
      if (allocated(vout)) then 
        deallocate(vout,stat=info)
      end if
    else if (allocated(vin)) then 
      if (.not.allocated(vout)) then 
        allocate(vout(size(vin,1),size(vin,2)),stat=info)
        if (info /= 0) return
      else
        if (size(vout) /= size(vin)) then 
          deallocate(vout,stat=info)
          if (info /= 0) return
          allocate(vout(size(vin,1),size(vin,2)),stat=info)
          if (info /= 0) return
        end if
      end if
      vout = vin
      deallocate(vin,stat=info)
    end if
#endif
  end Subroutine psb_ztransfer2d

  Subroutine psb_itransfer1d(vin,vout,info)
    use psb_error_mod
    integer, allocatable, intent(inout) :: vin(:),vout(:)
    integer, intent(out) :: info 
    !
    ! 
    info = 0
#ifdef HAVE_MOVE_ALLOC
    if (allocated(vin)) then 
      call move_alloc(vin,vout)
    else if (allocated(vout)) then 
      write(0,*) 'transfer: Clearing output'
      deallocate(vout)
    end if
#else
    if (.not.allocated(vin) ) then 
      if (allocated(vout)) then 
        deallocate(vout,stat=info)
      end if
    else if (allocated(vin)) then 
      if (.not.allocated(vout)) then 
        allocate(vout(size(vin)),stat=info)
        if (info /= 0) return
      else
        if (size(vout) /= size(vin)) then 
          deallocate(vout,stat=info)
          if (info /= 0) return
          allocate(vout(size(vin)),stat=info)
          if (info /= 0) return
        end if
      end if
      vout = vin
      deallocate(vin,stat=info)
    end if
#endif
  end Subroutine psb_itransfer1d

  Subroutine psb_itransfer2d(vin,vout,info)
    use psb_error_mod
    integer, allocatable, intent(inout) :: vin(:,:),vout(:,:)
    integer, intent(out) :: info 
    !
    ! 
    info = 0
#ifdef HAVE_MOVE_ALLOC
    if (allocated(vin)) then 
      call move_alloc(vin,vout)
    else if (allocated(vout)) then 
      deallocate(vout)
    end if
#else
    if (.not.allocated(vin) ) then 
      if (allocated(vout)) then 
        deallocate(vout,stat=info)
      end if
    else if (allocated(vin)) then 
      if (.not.allocated(vout)) then 
        allocate(vout(size(vin,1),size(vin,2)),stat=info)
        if (info /= 0) return
      else
        if (size(vout) /= size(vin)) then 
          deallocate(vout,stat=info)
          if (info /= 0) return
          allocate(vout(size(vin,1),size(vin,2)),stat=info)
          if (info /= 0) return
        end if
      end if
      vout = vin
      deallocate(vin,stat=info)
    end if
#endif
  end Subroutine psb_itransfer2d

end module psb_realloc_mod
