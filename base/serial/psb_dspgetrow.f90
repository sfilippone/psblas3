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
! File:  psb_dspgetrow.f90 
! Subroutine: psb_dspgetrow
!    Gets one or more rows from a sparse matrix. 
! Parameters:
!*****************************************************************************
!*                                                                           *
!*                                                                           *
!*****************************************************************************
subroutine psb_dspgetrow(irw,a,nz,ia,ja,val,info,iren,lrw,append,nzin)
  ! Output is always in  COO format 
  use psb_spmat_type
  use psb_const_mod
  use psb_getrow_mod
  use psb_string_mod
  use psb_serial_mod, psb_protect_name => psb_dspgetrow
  implicit none

  type(psb_dspmat_type), intent(in)    :: a
  integer, intent(in)                  :: irw
  integer, intent(out)                 :: nz
  integer, allocatable, intent(inout)  :: ia(:), ja(:)
  real(kind(1.d0)), allocatable,  intent(inout)    :: val(:)
  integer,intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer, intent(in), optional        :: iren(:)
  integer, intent(in), optional        :: lrw, nzin

  logical :: append_ 
  integer :: i,j,k,ip,jp,nr,idx,iret,nzin_, nza, lrw_, irw_, err_act
  character(len=20)                 :: name, ch_err

  name='psb_spgetrow'
  info  = 0

  call psb_erractionsave(err_act)

  irw_ = irw 
  if (present(lrw)) then
    lrw_ = lrw
  else
    lrw_ = irw
  endif
  if (lrw_ < irw) then
    write(0,*) 'SPGTBLK input error: fixing lrw',irw,lrw_
    lrw_ = irw
  end if
  if (present(append)) then
    append_=append
  else
    append_=.false.
  endif


  if ((append_).and.(present(nzin))) then 
    nzin_ = nzin
  else
    nzin_ = 0
  endif

  if (toupper(a%fida) == 'CSR') then 
     call csr_getrow(irw_,a,nz,ia,ja,val,nzin_,append_,lrw_,info,iren)
     
  else if (toupper(a%fida) == 'COO') then 
     call coo_getrow(irw_,a,nz,ia,ja,val,nzin_,append_,lrw_,info,iren)
     
  else if (toupper(a%fida) == 'JAD') then 
     call jad_getrow(irw_,a,nz,ia,ja,val,nzin_,append_,lrw_,info,iren)

  else
     info=136
     ch_err=a%fida(1:3)
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if
  
  if (info /= 0) goto 9999
!!$  call psb_erractionrestore(err_act)
  return
  
9999 continue
!!$  call psb_erractionrestore(err_act)
  call psb_erractionsave(err_act)
  if (err_act.eq.psb_act_abort_) then
     call psb_error()
     return
  end if
  return
  

end subroutine psb_dspgetrow

