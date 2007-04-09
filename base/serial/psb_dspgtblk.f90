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
! File:  psb_dspgtblk.f90 
! Subroutine: psb_dspgtblk
!    Gets one or more rows from a sparse matrix. 
! Parameters:
!*****************************************************************************
!*                                                                           *
!* Takes a specified row from matrix A and copies into matrix B (possibly    *
!*  appending to B). Output is always COO. Input might be anything,          *
!*                                                                           *
!*****************************************************************************
subroutine psb_dspgtblk(irw,a,b,info,append,iren,lrw,srt)
  ! Output is always in  COO format  into B, irrespective of 
  ! the input format 
  use psb_spmat_type
  use psb_const_mod
  use psb_string_mod
  use psb_serial_mod, psb_protect_name => psb_dspgtblk
  implicit none

  type(psb_dspmat_type), intent(in)     :: a
  integer, intent(in)                   :: irw
  type(psb_dspmat_type), intent(inout)  :: b
  integer,intent(out)                   :: info
  logical, intent(in), optional         :: append
  integer, intent(in), target, optional :: iren(:)
  integer, intent(in), optional         :: lrw
  logical, intent(in), optional         :: srt

  logical            :: append_,srt_
  integer            :: i,j,k,ip,jp,nr,idx, nz,iret,nzb, nza, lrw_, irw_, err_act
  character(len=20)  :: name, ch_err

  name='psb_spgtblk'
  info  = 0
!!$  call psb_erractionsave(err_act)

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
    append_ = append
  else
    append_ = .false.
  endif

  if (present(srt)) then
    srt_ = srt
  else
    srt_ = .true.
  endif


  if (append_) then 
    nzb = b%infoa(psb_nnz_)
  else
    nzb = 0
    b%m = 0 
    b%k = 0
    b%descra = a%descra
  endif
  b%fida = 'COO'
 
  call psb_sp_getrow(irw,a,nz,b%ia1,b%ia2,b%aspk,info,iren=iren,&
       & lrw=lrw_,append=append_,nzin=nzb)
  if (.not.allocated(b%pl)) then 
    allocate(b%pl(1),stat=info)
    b%pl = 0
  endif
  if (.not.allocated(b%pr)) then 
    allocate(b%pr(1),stat=info)
    b%pr = 0
  endif
  b%infoa(psb_nnz_) = nzb+nz
  b%m = b%m+lrw_-irw+1
  b%k = max(b%k,a%k)
  if (srt_) call psb_fixcoo(b,info)
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


end subroutine psb_dspgtblk

