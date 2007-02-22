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
! File:  psb_zspgetrow.f90 
! Subroutine: psb_zspgetrow
!    Gets one or more rows from a sparse matrix. 
! Parameters:

!*****************************************************************************
!*                                                                           *
!* Takes a specified row from matrix A and copies into NZ,IA,JA,VAL  in COO  *
!* format.                                                                   *
!*                                                                           *
!*****************************************************************************
subroutine psb_zspgetrow(irw,a,nz,ia,ja,val,info,iren,lrw)
  use psb_spmat_type
  use psb_string_mod
  use psb_serial_mod, psb_protect_name => psb_zspgetrow
  type(psb_zspmat_type), intent(in) :: a
  integer, intent(in)       :: irw
  integer, intent(out)      :: nz
  integer, intent(inout)    :: ia(:), ja(:)
  complex(kind(1.d0)),  intent(inout)    :: val(:)
  integer, intent(in), target, optional :: iren(:)
  integer, intent(in), optional :: lrw
  integer, intent(out)  :: info

  integer               :: lrw_, ierr(5), err_act
  type(psb_zspmat_type) :: b
  integer, pointer      :: iren_(:)
  character(len=20)     :: name, ch_err


  name='psb_sp_getrow'
  info  = 0
  call psb_erractionsave(err_act)
  call psb_set_erraction(0)  

  call psb_nullify_sp(b)

  if (present(lrw)) then
    lrw_ = lrw
  else
    lrw_ = irw
  endif
  if (lrw_ < irw) then
    write(0,*) 'SPGETROW input error: fixing lrw',irw,lrw_
    lrw_ = irw
  end if
  call psb_sp_all(lrw_-irw+1,lrw_-irw+1,b,info)

  if (present(iren)) then
    call psb_sp_getblk(irw,a,b,info,iren=iren,lrw=lrw_)
  else 
    call psb_sp_getblk(irw,a,b,info,lrw=lrw_)
  end if
  if (info /= 0) then     
    info=136
    ch_err=a%fida(1:3)
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (toupper(b%fida) /= 'COO') then 
    info=4010
    ch_err=a%fida(1:3)
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  endif

  nz = b%infoa(psb_nnz_)

  if (size(ia)>= nz) then 
    ia(1:nz) = b%ia1(1:nz)
  else
    info    = 135
    ierr(1) = 4
    ierr(2) = size(ia)
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  endif

  if (size(ja)>= nz) then 
    ja(1:nz) = b%ia2(1:nz)
  else
    info    = 135
    ierr(1) = 5
    ierr(2) = size(ja)
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  endif

  if (size(val)>= nz) then 
    val(1:nz) = b%aspk(1:nz)
  else
    info    = 135
    ierr(1) = 6
    ierr(2) = size(val)
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  endif


  call psb_sp_free(b,info)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return


end subroutine psb_zspgetrow

