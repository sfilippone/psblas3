!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
!
!  Warning: MM does not define a format for an array with integer entries.
!  Hence we hijack the REAL format, but this could lead to errors when
!  used with non-integer files. 
!
subroutine mm_ivet_read(b, info, iunit, filename)   
  use psb_base_mod
  implicit none
  integer(psb_ipk_), allocatable, intent(out)  :: b(:)
  integer(psb_ipk_), intent(out)        :: info
  integer(psb_ipk_), optional, intent(in)          :: iunit
  character(len=*), optional, intent(in) :: filename
  integer(psb_ipk_) :: nrow, ncol, i,root, np,  me,  ircode, j, infile
  character            :: mmheader*15, fmt*15, object*10, type*10, sym*15,&
       & line*1024

  info = psb_success_
  if (present(filename)) then
    if (filename == '-') then 
      infile=5
    else
      if (present(iunit)) then 
        infile=iunit
      else
        infile=99
      endif
      open(infile,file=filename, status='OLD', err=901, action='READ')
    endif
  else 
    if (present(iunit)) then 
      infile=iunit
    else
      infile=5
    endif
  endif

  read(infile,fmt=*, end=902) mmheader, object, fmt, type, sym

  if ( (object /= 'matrix').or.(fmt /= 'array')) then
    write(psb_err_unit,*) 'read_rhs: input file type not yet supported'
    info = -3
    return
  end if

  do 
    read(infile,fmt='(a)') line
    if (line(1:1) /= '%')  exit
  end do

  read(line,fmt=*)nrow,ncol

  if ((psb_tolower(type) == 'real').and.(psb_tolower(sym) == 'general')) then
    allocate(b(nrow),stat = ircode)
    if (ircode /= 0)   goto 993
    do i=1,nrow
      read(infile,fmt=*,end=902) b(i)
    end do

  end if      ! read right hand sides
  if (infile /= 5) close(infile)

  return 
  ! open failed
901 write(psb_err_unit,*) 'mm_vet_read: could not open file ',&
       & infile,' for input'
  info = -1
  return

902 write(psb_err_unit,*) 'mmv_vet_read: unexpected end of file ',infile,&
       & ' during input'
  info = -2
  return
993 write(psb_err_unit,*) 'mm_vet_read: memory allocation failure'
  info = -3
  return
end subroutine mm_ivet_read

subroutine mm_ivet2_read(b, info, iunit, filename)   
  use psb_base_mod
  implicit none
  integer(psb_ipk_), allocatable, intent(out)  :: b(:,:)
  integer(psb_ipk_), intent(out)        :: info
  integer(psb_ipk_), optional, intent(in)          :: iunit
  character(len=*), optional, intent(in) :: filename
  integer(psb_ipk_) :: nrow, ncol, i,root, np,  me,  ircode, j, infile
  character            :: mmheader*15, fmt*15, object*10, type*10, sym*15,&
       & line*1024

  info = psb_success_
  if (present(filename)) then
    if (filename == '-') then 
      infile=5
    else
      if (present(iunit)) then 
        infile=iunit
      else
        infile=99
      endif
      open(infile,file=filename, status='OLD', err=901, action='READ')
    endif
  else 
    if (present(iunit)) then 
      infile=iunit
    else
      infile=5
    endif
  endif

  read(infile,fmt=*, end=902) mmheader, object, fmt, type, sym

  if ( (object /= 'matrix').or.(fmt /= 'array')) then
    write(psb_err_unit,*) 'read_rhs: input file type not yet supported'
    info = -3
    return
  end if

  do 
    read(infile,fmt='(a)') line
    if (line(1:1) /= '%')  exit
  end do

  read(line,fmt=*)nrow,ncol
  
  if ((psb_tolower(type) == 'real').and.(psb_tolower(sym) == 'general')) then
    allocate(b(nrow,ncol),stat = ircode)
    if (ircode /= 0)   goto 993
    read(infile,fmt=*,end=902) ((b(i,j), i=1,nrow),j=1,ncol)

  end if      ! read right hand sides
  if (infile /= 5) close(infile)

  return 
  ! open failed
901 write(psb_err_unit,*) 'mm_vet_read: could not open file ',&
       & infile,' for input'
  info = -1
  return

902 write(psb_err_unit,*) 'mmv_vet_read: unexpected end of file ',infile,&
       & ' during input'
  info = -2
  return
993 write(psb_err_unit,*) 'mm_vet_read: memory allocation failure'
  info = -3
  return
end subroutine mm_ivet2_read

subroutine mm_ivet2_write(b, header, info, iunit, filename)   
  use psb_base_mod
  implicit none
  integer(psb_ipk_), intent(in)  :: b(:,:)
  character(len=*), intent(in) :: header
  integer(psb_ipk_), intent(out)        :: info
  integer(psb_ipk_), optional, intent(in)          :: iunit
  character(len=*), optional, intent(in) :: filename
  integer(psb_ipk_) :: nrow, ncol, i,root, np,  me,  ircode, j, outfile

  character(len=80)                 :: frmtv 

  info = psb_success_
  if (present(filename)) then
    if (filename == '-') then 
      outfile=6
    else
      if (present(iunit)) then 
        outfile=iunit
      else
        outfile=99
      endif
      open(outfile,file=filename, err=901, action='WRITE')
    endif
  else 
    if (present(iunit)) then 
      outfile=iunit
    else
      outfile=6
    endif
  endif

  write(outfile,'(a)') '%%MatrixMarket matrix array real general'
  write(outfile,'(a)') '% '//trim(header)
  write(outfile,'(a)') '% '
  nrow = size(b,1) 
  ncol = size(b,2) 
  write(outfile,*) nrow, ncol

  write(outfile,fmt='(I14,1x)') ((b(i,j), i=1,nrow),j=1,ncol)

  if (outfile /= 6) close(outfile)

  return 
  ! open failed
901 write(psb_err_unit,*) 'mm_vet_write: could not open file ',&
       & outfile,' for output'
  info = -1
  return

end subroutine mm_ivet2_write

subroutine mm_ivet1_write(b, header, info, iunit, filename)   
  use psb_base_mod
  implicit none
  integer(psb_ipk_), intent(in)  :: b(:)
  character(len=*), intent(in) :: header
  integer(psb_ipk_), intent(out)        :: info
  integer(psb_ipk_), optional, intent(in)          :: iunit
  character(len=*), optional, intent(in) :: filename
  integer(psb_ipk_) :: nrow, ncol, i,root, np,  me,  ircode, j, outfile

  character(len=80)                 :: frmtv 

  info = psb_success_
  if (present(filename)) then
    if (filename == '-') then 
      outfile=6
    else
      if (present(iunit)) then 
        outfile=iunit
      else
        outfile=99
      endif
      open(outfile,file=filename, err=901, action='WRITE')
    endif
  else 
    if (present(iunit)) then 
      outfile=iunit
    else
      outfile=6
    endif
  endif

  write(outfile,'(a)') '%%MatrixMarket matrix array real general'
  write(outfile,'(a)') '% '//trim(header)
  write(outfile,'(a)') '% '
  nrow = size(b,1) 
  ncol = 1
  write(outfile,*) nrow,ncol

  write(frmtv,'(a,i0,a)') '(',ncol,'(i14,1x))'

  do i=1,size(b,1) 
    write(outfile,frmtv) b(i)
  end do

  if (outfile /= 6) close(outfile)

  return 
  ! open failed
901 write(psb_err_unit,*) 'mm_vet_write: could not open file ',&
       & outfile,' for output'
  info = -1
  return

end subroutine mm_ivet1_write

