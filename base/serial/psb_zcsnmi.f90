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
! File:  psb_zcsnmi.f90 
! Subroutine: 
! Parameters:

real(kind(1.d0)) function psb_zcsnmi(a,info,trans)

  use psb_spmat_type
  use psb_error_mod
  implicit none

  type(psb_zspmat_type), intent(in)  :: a
  integer, intent(out)               :: info
  character, optional                :: trans

  interface
     real(kind(1.d0)) function zcsnmi(trans,m,n,fida,descra,a,ia1,ia2,&
          &                 infoa,ierror)
       integer          ::  m,n, ierror
       character        ::  trans
       integer          ::  ia1(*),ia2(*),infoa(*)
       character        ::  descra*11, fida*5
       complex(kind(1.d0)) ::  a(*)
     end function zcsnmi
  end interface

  integer             :: err_act
  character           :: itrans
  character(len=20)   :: name, ch_err

  name='psb_zcsnmi'
  call psb_erractionsave(err_act)

  if(present(trans)) then
     itrans=trans
  else
     itrans='N'
  end if

  psb_zcsnmi = zcsnmi(itrans,a%m,a%k,a%fida,a%descra,a%aspk,a%ia1,a%ia2,a%infoa,info)
  if(info/=0) then
     psb_zcsnmi = -1
     info=4010
     ch_err='psb_zcsnmi'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if
  

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
     call psb_error()
     return
  end if
  return

end function psb_zcsnmi
