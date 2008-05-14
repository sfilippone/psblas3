!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
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
! File:  psb_dspgtblk.f90 
! Subroutine: psb_dspgtblk
!    Gets one or more rows from a sparse matrix. 
! Arguments:
!*****************************************************************************
!*                                                                           *
!* Takes a specified row from matrix A and copies into matrix B (possibly    *
!*  appending to B). Output is always COO. Input might be anything,          *
!*                                                                           *
!*****************************************************************************
subroutine psb_dspsetbld1(a,info)
  ! Output is always in  COO format  into B, irrespective of 
  ! the input format 
  use psb_spmat_type
  use psb_const_mod
  use psb_string_mod
  use psb_serial_mod, psb_protect_name => psb_dspsetbld1

  implicit none

  type(psb_dspmat_type), intent(inout) :: a
  integer, intent(out)                 :: info      

  logical            :: append_,srt_
  integer            :: i,j, err_act
  character(len=20)  :: name

  name='psb_sp_setbld'
  info  = 0
  call psb_erractionsave(err_act)
  
  select case(psb_sp_getifld(psb_state_,a,info))
  case (psb_spmat_asb_,psb_spmat_upd_) 
    call psb_spcnv(a,info,afmt='COO')
    if (info == 0) call psb_sp_setifld(psb_spmat_bld_,psb_state_,a,info)
  case (psb_spmat_bld_) 
    ! do nothing
    
  case default
    info=591     
    call psb_errpush(info,name)
    goto 9999
  end select


  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Internal call')
    goto 9999
  end if

  call psb_erractionrestore(err_act)

  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return


end subroutine psb_dspsetbld1
subroutine psb_dspsetbld2(a,b,info)
  ! Output is always in  COO format  into B, irrespective of 
  ! the input format 
  use psb_spmat_type
  use psb_const_mod
  use psb_string_mod
  use psb_serial_mod, psb_protect_name => psb_dspsetbld2

  implicit none

  type(psb_dspmat_type), intent(in)  :: a
  type(psb_dspmat_type), intent(out) :: b
  integer, intent(out)               :: info      

  logical            :: append_,srt_
  integer            :: i,j, err_act
  character(len=20)  :: name

  name='psb_sp_setbld'
  info  = 0
  call psb_erractionsave(err_act)


  select case(psb_sp_getifld(psb_state_,a,info))
  case (psb_spmat_asb_,psb_spmat_upd_) 
    call psb_spcnv(a,b,info,afmt='COO')
    if (info == 0) call psb_sp_setifld(psb_spmat_bld_,psb_state_,b,info)

  case (psb_spmat_bld_) 
    call psb_sp_clone(a,b,info) 

  case default
    info=591     
    call psb_errpush(info,name)
    goto 9999
  end select

  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Internal call')
    goto 9999
  end if

  call psb_erractionrestore(err_act)

  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return


end subroutine psb_dspsetbld2

