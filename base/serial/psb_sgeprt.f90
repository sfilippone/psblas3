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
! File:  psb_scsprt.f90 
! Subroutine: 
! Arguments:

!*****************************************************************************
!*                                                                           *
!* Print out a matrix.                                                       *
!*  Should really align with the F77 version under the SERIAL dir, which     *
!*  does a nice printout in MatrixMarket format; this would be a quick job.  *
!*                                                                           *
!*  Handles both a shift in the row/col indices and a fuctional transform    *
!*  on the indices.                                                          *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*****************************************************************************
subroutine psb_sgeprtn2(fname,a,head)
  use psb_serial_mod, psb_protect_name => psb_sgeprtn2
  implicit none 
  
  character(len=*), intent(in)  :: fname   
  real(psb_spk_), intent(in)    :: a(:,:)
  character(len=*), optional    :: head

  !
  integer(psb_ipk_) :: iout, info
  logical        :: isopen
  
  ! Search for an unused unit to write
  iout = 7
  do 
    inquire(unit=iout, opened=isopen)
    if (.not.isopen) exit
    iout = iout + 1
    if (iout > 99) exit
  end do
  if (iout > 99) then 
    write(psb_err_unit,*) 'Error: could not find a free unit for I/O'
    return
  end if
  open(iout,file=fname,iostat=info)
  if (info == psb_success_) then 
    call psb_geprt(iout,a,head)
    close(iout)
  else
    write(psb_err_unit,*) 'Error: could not open ',fname,' for output'
  end if

end subroutine psb_sgeprtn2

subroutine psb_sgeprtn1(fname,a,head)
  use psb_serial_mod, psb_protect_name => psb_sgeprtn1
  implicit none 
  
  character(len=*), intent(in)  :: fname   
  real(psb_spk_), intent(in)    :: a(:)
  character(len=*), optional    :: head

  !
  integer(psb_ipk_) :: iout, info
  logical        :: isopen
  
  ! Search for an unused unit to write
  iout = 7
  do 
    inquire(unit=iout, opened=isopen)
    if (.not.isopen) exit
    iout = iout + 1
    if (iout > 99) exit
  end do
  if (iout > 99) then 
    write(psb_err_unit,*) 'Error: could not find a free unit for I/O'
    return
  end if
  open(iout,file=fname,iostat=info)
  if (info == psb_success_) then 
    call psb_geprt(iout,a,head)
    close(iout)
  else
    write(psb_err_unit,*) 'Error: could not open ',fname,' for output'
  end if

end subroutine psb_sgeprtn1

subroutine psb_sgeprt2(iout,a,head)
  use psb_serial_mod, psb_protect_name => psb_sgeprt2
  implicit none 

  integer(psb_ipk_), intent(in)            :: iout
  real(psb_spk_), intent(in)     :: a(:,:)
  character(len=*), optional     :: head
  character(len=80)              :: frmtv 
  integer(psb_ipk_) :: irs,ics,i,j, nmx, ni, nrow, ncol

  write(iout,'(a)') '%%MatrixMarket matrix array real general'
  write(iout,'(a)') '% '//trim(head)
  write(iout,'(a)') '% '
  nrow = size(a,1) 
  ncol = size(a,2) 
  write(iout,*) nrow,ncol

  write(frmtv,'(a,i3.3,a)') '(',ncol,'(es26.18,1x))'

  do i=1,nrow
    write(iout,frmtv) a(i,1:ncol)
  end do

  if (iout /= 6) close(iout)

  return 
  ! open failed
901 write(psb_err_unit,*) 'geprt: could not open file ',&
       & iout,' for output'
  return
end subroutine psb_sgeprt2

subroutine psb_sgeprt1(iout,a,head)
  use psb_serial_mod, psb_protect_name => psb_sgeprt1
  implicit none 

  integer(psb_ipk_), intent(in)            :: iout
  real(psb_spk_), intent(in)     :: a(:)
  character(len=*), optional     :: head
  character(len=80)              :: frmtv 
  integer(psb_ipk_) :: irs,ics,i,j, nmx, ni, nrow, ncol

  write(iout,'(a)') '%%MatrixMarket matrix array real general'
  write(iout,'(a)') '% '//trim(head)
  write(iout,'(a)') '% '
  nrow = size(a,1) 
  ncol = 1
  write(iout,*) nrow

  write(frmtv,'(a,i3.3,a)') '(',ncol,'(es26.18,1x))'

  do i=1,nrow
    write(iout,frmtv) a(i)
  end do

  if (iout /= 6) close(iout)

  return 
  ! open failed
901 write(psb_err_unit,*) 'geprt: could not open file ',&
       & iout,' for output'
  return
end subroutine psb_sgeprt1
