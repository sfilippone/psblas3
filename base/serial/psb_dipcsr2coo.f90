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
! File:  psb_dipcsr2coo.f90 
! Subroutine: 
! Parameters:

Subroutine psb_dipcsr2coo(a,info)
  use psb_spmat_type
  use psb_const_mod
  use psb_error_mod
  use psb_string_mod
  use psb_realloc_mod
  implicit none

  !....Parameters...
  Type(psb_dspmat_type), intent(inout) :: A
  Integer, intent(out)                 :: info

  !locals
  Integer              :: nza, nr
  integer              :: i,j,err_act
  logical, parameter   :: debug=.false.
  integer, allocatable :: iaux(:), itemp(:)
  character(len=20)    :: name, ch_err

  name='psb_dipcsr2coo'
  info  = 0
  call psb_erractionsave(err_act)

  if (toupper(a%fida) /= 'CSR') then 
    info = 5
    call psb_errpush(info,name)
    goto 9999
  end if

  nr  = a%m 
  nza = a%ia2(nr+1) - 1
  allocate(iaux(nza),stat=info)
  if (info /=0) then 
    write(0,*) 'Failed allocation ',info, nza
    return
  end if
!!$  write(0,*) 'ipcsr2coo ',a%m      
  call psb_transfer(a%ia2,itemp,info)
  call psb_transfer(a%ia1,a%ia2,info)
  call psb_transfer(iaux,a%ia1,info)
  
  do i=1, nr
    do j=itemp(i),itemp(i+1)-1
      a%ia1(j) = i
    end do
  end do
  
  a%fida='COO'
  a%infoa(psb_nnz_) = nza
  a%infoa(psb_srtd_) = psb_isrtdcoo_
  a%infoa(psb_upd_) = psb_upd_srch_

  deallocate(itemp)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
     call psb_error()
     return
  end if
  return

end Subroutine psb_dipcsr2coo
