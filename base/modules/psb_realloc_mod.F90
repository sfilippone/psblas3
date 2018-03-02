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
module psb_realloc_mod
  use psb_const_mod
  use psb_m_realloc_mod
  use psb_e_realloc_mod
  use psb_s_realloc_mod
  use psb_d_realloc_mod
  use psb_c_realloc_mod
  use psb_z_realloc_mod
  
  implicit none
 
  !
  ! psb_realloc will reallocate the input array to have exactly 
  ! the size specified, possibly shortening it. 
  !
  Interface psb_realloc
    module procedure psb_reallocate2i1d
    module procedure psb_reallocate2i1s
    module procedure psb_reallocate2i1z
    module procedure psb_reallocate2i1c
  end Interface psb_realloc


  logical, private :: do_maybe_free_buffer = .true.

Contains
  
  function psb_get_maybe_free_buffer() result(res)
    logical :: res
    res = do_maybe_free_buffer
  end function psb_get_maybe_free_buffer

  subroutine psb_set_maybe_free_buffer(val)
    logical, intent(in) :: val
    do_maybe_free_buffer = val
  end subroutine psb_set_maybe_free_buffer
  


  Subroutine psb_reallocate2i1s(len,rrax,y,z,info)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len
    integer(psb_ipk_),allocatable, intent(inout)  :: rrax(:),y(:)
    Real(psb_spk_),allocatable, intent(inout) :: z(:)
    integer(psb_ipk_) :: info
    character(len=20)  :: name
    integer(psb_ipk_) :: err_act, err
    logical, parameter :: debug=.false.

    name='psb_reallocate2i1s'
    call psb_erractionsave(err_act)


    info=psb_success_
    call psb_realloc(len,rrax,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,y,info)    
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,z,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  End Subroutine psb_reallocate2i1s


  Subroutine psb_reallocate2i1d(len,rrax,y,z,info)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len
    integer(psb_ipk_),allocatable, intent(inout)  :: rrax(:),y(:)
    Real(psb_dpk_),allocatable, intent(inout) :: z(:)
    integer(psb_ipk_) :: info
    character(len=20)  :: name
    integer(psb_ipk_) :: err_act, err

    name='psb_reallocate2i1d'
    call psb_erractionsave(err_act)

    info=psb_success_

    call psb_realloc(len,rrax,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,y,info)    
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,z,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  End Subroutine psb_reallocate2i1d



  Subroutine psb_reallocate2i1c(len,rrax,y,z,info)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len
    integer(psb_ipk_),allocatable, intent(inout) :: rrax(:),y(:)
    complex(psb_spk_),allocatable, intent(inout) :: z(:)
    integer(psb_ipk_) :: info
    character(len=20)  :: name
    integer(psb_ipk_) :: err_act, err

    name='psb_reallocate2i1c'
    call psb_erractionsave(err_act)


    info=psb_success_
    call psb_realloc(len,rrax,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,y,info)    
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,z,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  End Subroutine psb_reallocate2i1c

  Subroutine psb_reallocate2i1z(len,rrax,y,z,info)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len
    integer(psb_ipk_),allocatable, intent(inout) :: rrax(:),y(:)
    complex(psb_dpk_),allocatable, intent(inout) :: z(:)
    integer(psb_ipk_) :: info
    character(len=20)  :: name
    integer(psb_ipk_) :: err_act, err

    name='psb_reallocate2i1z'
    call psb_erractionsave(err_act)

    info=psb_success_
    call psb_realloc(len,rrax,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,y,info)    
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,z,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  End Subroutine psb_reallocate2i1z

end module psb_realloc_mod
