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
! File:  psb_zspscal.f90 
! Subroutine: 
! Parameters:

!*****************************************************************************
!*                                                                           *
!*                                                                           *
!*****************************************************************************
subroutine psb_zspscal(a,d,info)
  ! the input format 
  use psb_spmat_type
  use psb_error_mod
  use psb_const_mod
  use psb_string_mod
  implicit none

  type(psb_zspmat_type), intent(inout) :: a
  integer, intent(out)                 :: info
  complex(kind(1.d0)), intent(in)         :: d(:) 

  integer :: i,j,k,nr, nz,err_act
  character(len=20)                 :: name, ch_err

  name='psb_zspscal'
  info  = 0
  call psb_erractionsave(err_act)

  select case(toupper(a%fida(1:3)))
  case  ('CSR')

    do i=1, a%m
      do j=a%ia2(i),a%ia2(i+1)-1
        a%aspk(j) = a%aspk(j) * d(i)
      end do
    end do

  case ('COO') 

    do i=1,a%infoa(psb_nnz_)
      j=a%ia1(i)
      a%aspk(i) = a%aspk(i) * d(j)
    enddo

  case ('JAD') 
    info=135
    ch_err=a%fida(1:3)
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  case default
    info=136
    ch_err=a%fida(1:3)
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
     call psb_error()
     return
  end if
  return

end subroutine psb_zspscal

