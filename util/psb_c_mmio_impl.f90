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
subroutine mm_cvet_read(b, info, iunit, filename)   
  use psb_base_mod
  implicit none
  complex(psb_spk_), allocatable, intent(out)  :: b(:)
  integer(psb_ipk_), intent(out)        :: info
  integer(psb_ipk_), optional, intent(in)          :: iunit
  character(len=*), optional, intent(in) :: filename
  integer(psb_ipk_) :: nrow, ncol, i,root, np,  me,  ircode, j,infile
  real(psb_spk_)       :: bre, bim 
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

  if ((psb_tolower(type) == 'complex').and.(psb_tolower(sym) == 'general')) then
    allocate(b(nrow),stat = ircode)
    if (ircode /= 0)   goto 993
    do i=1, nrow
      read(infile,fmt=*,end=902) bre,bim
      b(i) = cmplx(bre,bim,kind=psb_spk_)
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
end subroutine mm_cvet_read

subroutine mm_cvet2_read(b, info, iunit, filename)   
  use psb_base_mod
  implicit none
  complex(psb_spk_), allocatable, intent(out)  :: b(:,:)
  integer(psb_ipk_), intent(out)        :: info
  integer(psb_ipk_), optional, intent(in)          :: iunit
  character(len=*), optional, intent(in) :: filename
  integer(psb_ipk_) :: nrow, ncol, i,root, np,  me,  ircode, j,infile
  real(psb_spk_)       :: bre, bim 
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
    do j=1, ncol
      do i=1, nrow
        read(infile,fmt=*,end=902) bre,bim
        b(i,j) = cmplx(bre,bim,kind=psb_spk_)
      end do
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
end subroutine mm_cvet2_read

subroutine mm_cvet2_write(b, header, info, iunit, filename)   
  use psb_base_mod
  implicit none
  complex(psb_spk_), intent(in)  :: b(:,:)
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

  write(outfile,'(a)') '%%MatrixMarket matrix array complex general'
  write(outfile,'(a)') '% '//trim(header)
  write(outfile,'(a)') '% '
  nrow = size(b,1) 
  ncol = size(b,2) 
  write(outfile,*) nrow,ncol

  write(outfile,fmt='(2(es26.18,1x))') ((b(i,j), i=1,nrow),j=1,ncol)

  if (outfile /= 6) close(outfile)

  return 
  ! open failed
901 write(psb_err_unit,*) 'mm_vet_write: could not open file ',&
       & outfile,' for output'
  info = -1
  return

end subroutine mm_cvet2_write

subroutine mm_cvet1_write(b, header, info, iunit, filename)   
  use psb_base_mod
  implicit none
  complex(psb_spk_), intent(in)  :: b(:)
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

  write(outfile,'(a)') '%%MatrixMarket matrix array complex general'
  write(outfile,'(a)') '% '//trim(header)
  write(outfile,'(a)') '% '
  nrow = size(b,1) 
  ncol = 1
  write(outfile,*) nrow,ncol

  write(frmtv,'(a,i0,a)') '(',2*ncol,'(es26.18,1x))'

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

end subroutine mm_cvet1_write

subroutine mm_cvect_read(b, info, iunit, filename)   
  use psb_base_mod
  use psb_mmio_mod, psb_protect_name => mm_cvect_read
  implicit none
  type(psb_c_vect_type), intent(inout)  :: b
  integer(psb_ipk_), intent(out)        :: info
  integer(psb_ipk_), optional, intent(in) :: iunit
  character(len=*), optional, intent(in) :: filename
  !
  complex(psb_spk_), allocatable :: bv(:)

  call mm_array_read(bv, info, iunit, filename)
  call b%bld(bv)
  
end subroutine mm_cvect_read

subroutine mm_cvect_write(b, header, info, iunit, filename)   
  use psb_base_mod
  use psb_mmio_mod, psb_protect_name => mm_cvect_write
  implicit none
  type(psb_c_vect_type), intent(inout)  :: b
  character(len=*), intent(in) :: header
  integer(psb_ipk_), intent(out)        :: info
  integer(psb_ipk_), optional, intent(in)          :: iunit
  character(len=*), optional, intent(in) :: filename
  info = psb_success_
  if (.not.allocated(b%v)) return
  call b%sync()

  call mm_array_write(b%v%v,header,info,iunit,filename)
  
end subroutine mm_cvect_write

subroutine cmm_mat_read(a, info, iunit, filename)   
  use psb_base_mod
  implicit none
  type(psb_cspmat_type), intent(out)  :: a
  integer(psb_ipk_), intent(out)        :: info
  integer(psb_ipk_), optional, intent(in)          :: iunit
  character(len=*), optional, intent(in) :: filename
  character      :: mmheader*15, fmt*15, object*10, type*10, sym*15
  character(1024)      :: line
  integer(psb_ipk_) :: nrow, ncol, nnzero
  integer(psb_ipk_) :: ircode, i,nzr,infile
  type(psb_c_coo_sparse_mat), allocatable :: acoo
  real(psb_spk_) :: are, aim
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

  read(infile,fmt=*,end=902) mmheader, object, fmt, type, sym

  if ( (psb_tolower(object) /= 'matrix').or.(psb_tolower(fmt) /= 'coordinate')) then
    write(psb_err_unit,*) 'READ_MATRIX: input file type not yet supported'
    info=909
    return
  end if

  do 
    read(infile,fmt='(a)') line
    if (line(1:1) /= '%')  exit
  end do
  read(line,fmt=*) nrow,ncol,nnzero

  allocate(acoo, stat=ircode)
  if (ircode /= 0)   goto 993    
  if ((psb_tolower(type) == 'complex').and.(psb_tolower(sym) == 'general')) then
    call acoo%allocate_mnnz(nrow,ncol,nnzero)
    do i=1,nnzero
      read(infile,fmt=*,end=902) acoo%ia(i),acoo%ja(i),are,aim
      acoo%val(i) = cmplx(are,aim,kind=psb_spk_)
    end do
    call acoo%set_nzeros(nnzero)
    call acoo%fix(info)

    call a%mv_from(acoo)

  else if ((psb_tolower(type) == 'complex').and.(psb_tolower(sym) == 'symmetric')) then
    ! we are generally working with non-symmetric matrices, so
    ! we de-symmetrize what we are about to read
    call acoo%allocate_mnnz(nrow,ncol,2*nnzero)
    do i=1,nnzero
      read(infile,fmt=*,end=902,err=905) acoo%ia(i),acoo%ja(i),are,aim
      acoo%val(i) = cmplx(are,aim,kind=psb_spk_)
    end do
    nzr = nnzero
    do i=1,nnzero
      if (acoo%ia(i) /= acoo%ja(i)) then 
        nzr = nzr + 1
        acoo%val(nzr) = acoo%val(i)
        acoo%ia(nzr) = acoo%ja(i)
        acoo%ja(nzr) = acoo%ia(i)
      end if
    end do
    call acoo%set_nzeros(nzr)
    call acoo%fix(info)

    call a%mv_from(acoo)

  else if ((psb_tolower(type) == 'complex').and.(psb_tolower(sym) == 'hermitian')) then
    ! we are generally working with non-symmetric matrices, so
    ! we de-symmetrize what we are about to read
    call acoo%allocate_mnnz(nrow,ncol,2*nnzero)
    do i=1,nnzero
      read(infile,fmt=*,end=902,err=905) acoo%ia(i),acoo%ja(i),are,aim
      acoo%val(i) = cmplx(are,aim,kind=psb_spk_)
    end do
    nzr = nnzero
    do i=1,nnzero
      if (acoo%ia(i) /= acoo%ja(i)) then 
        nzr = nzr + 1
        acoo%val(nzr) = conjg(acoo%val(i))
        acoo%ia(nzr) = acoo%ja(i)
        acoo%ja(nzr) = acoo%ia(i)
      end if
    end do
    call acoo%set_nzeros(nzr)
    call acoo%fix(info)

    call a%mv_from(acoo)

  else
    write(psb_err_unit,*) 'read_matrix: matrix type not yet supported'
    info=904
  end if
  if (infile /= 5) close(infile)
  return 

  ! open failed
901 info=901
  write(psb_err_unit,*) 'read_matrix: could not open file ',filename,' for input'
  return
902 info=902
  write(psb_err_unit,*) 'READ_MATRIX: Unexpected end of file '
  return
905 info=905
  write(psb_err_unit,*) 'READ_MATRIX: Error at line',i
  return
993 info=993
  write(psb_err_unit,*) 'READ_MATRIX: Memory allocation failure'
  return
end subroutine cmm_mat_read


subroutine cmm_mat_write(a,mtitle,info,iunit,filename)
  use psb_base_mod
  implicit none
  type(psb_cspmat_type), intent(in)  :: a
  integer(psb_ipk_), intent(out)        :: info
  character(len=*), intent(in) :: mtitle
  integer(psb_ipk_), optional, intent(in)          :: iunit
  character(len=*), optional, intent(in) :: filename
  integer(psb_ipk_) :: iout


  info = psb_success_

  if (present(filename)) then 
    if (filename == '-') then 
      iout=6
    else
      if (present(iunit)) then 
        iout = iunit
      else
        iout=99
      endif
      open(iout,file=filename, err=901, action='WRITE')
    endif
  else 
    if (present(iunit)) then 
      iout = iunit   
    else
      iout=6
    endif
  endif

  call a%print(iout,head=mtitle)

  if (iout /= 6) close(iout)


  return

901 continue 
  info=901
  write(psb_err_unit,*) 'Error while opening ',filename
  return
end subroutine cmm_mat_write

subroutine lcmm_mat_read(a, info, iunit, filename)   
  use psb_base_mod
  implicit none
  type(psb_lcspmat_type), intent(out)  :: a
  integer(psb_ipk_), intent(out)        :: info
  integer(psb_ipk_), optional, intent(in)          :: iunit
  character(len=*), optional, intent(in) :: filename
  character      :: mmheader*15, fmt*15, object*10, type*10, sym*15
  character(1024)      :: line
  integer(psb_lpk_) :: nrow, ncol, nnzero, i,nzr
  integer(psb_ipk_) :: ircode,infile
  type(psb_lc_coo_sparse_mat), allocatable :: acoo
  real(psb_spk_) :: are, aim
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

  read(infile,fmt=*,end=902) mmheader, object, fmt, type, sym

  if ( (psb_tolower(object) /= 'matrix').or.(psb_tolower(fmt) /= 'coordinate')) then
    write(psb_err_unit,*) 'READ_MATRIX: input file type not yet supported'
    info=909
    return
  end if

  do 
    read(infile,fmt='(a)') line
    if (line(1:1) /= '%')  exit
  end do
  read(line,fmt=*) nrow,ncol,nnzero

  allocate(acoo, stat=ircode)
  if (ircode /= 0)   goto 993    
  if ((psb_tolower(type) == 'complex').and.(psb_tolower(sym) == 'general')) then
    call acoo%allocate_mnnz(nrow,ncol,nnzero)
    do i=1,nnzero
      read(infile,fmt=*,end=902) acoo%ia(i),acoo%ja(i),are,aim
      acoo%val(i) = cmplx(are,aim,kind=psb_spk_)
    end do
    call acoo%set_nzeros(nnzero)
    call acoo%fix(info)

    call a%mv_from(acoo)

  else if ((psb_tolower(type) == 'complex').and.(psb_tolower(sym) == 'symmetric')) then
    ! we are generally working with non-symmetric matrices, so
    ! we de-symmetrize what we are about to read
    call acoo%allocate_mnnz(nrow,ncol,2*nnzero)
    do i=1,nnzero
      read(infile,fmt=*,end=902,err=905) acoo%ia(i),acoo%ja(i),are,aim
      acoo%val(i) = cmplx(are,aim,kind=psb_spk_)
    end do
    nzr = nnzero
    do i=1,nnzero
      if (acoo%ia(i) /= acoo%ja(i)) then 
        nzr = nzr + 1
        acoo%val(nzr) = acoo%val(i)
        acoo%ia(nzr) = acoo%ja(i)
        acoo%ja(nzr) = acoo%ia(i)
      end if
    end do
    call acoo%set_nzeros(nzr)
    call acoo%fix(info)

    call a%mv_from(acoo)

  else if ((psb_tolower(type) == 'complex').and.(psb_tolower(sym) == 'hermitian')) then
    ! we are generally working with non-symmetric matrices, so
    ! we de-symmetrize what we are about to read
    call acoo%allocate_mnnz(nrow,ncol,2*nnzero)
    do i=1,nnzero
      read(infile,fmt=*,end=902,err=905) acoo%ia(i),acoo%ja(i),are,aim
      acoo%val(i) = cmplx(are,aim,kind=psb_spk_)
    end do
    nzr = nnzero
    do i=1,nnzero
      if (acoo%ia(i) /= acoo%ja(i)) then 
        nzr = nzr + 1
        acoo%val(nzr) = conjg(acoo%val(i))
        acoo%ia(nzr) = acoo%ja(i)
        acoo%ja(nzr) = acoo%ia(i)
      end if
    end do
    call acoo%set_nzeros(nzr)
    call acoo%fix(info)

    call a%mv_from(acoo)

  else
    write(psb_err_unit,*) 'read_matrix: matrix type not yet supported'
    info=904
  end if
  if (infile /= 5) close(infile)
  return 

  ! open failed
901 info=901
  write(psb_err_unit,*) 'read_matrix: could not open file ',filename,' for input'
  return
902 info=902
  write(psb_err_unit,*) 'READ_MATRIX: Unexpected end of file '
  return
905 info=905
  write(psb_err_unit,*) 'READ_MATRIX: Error at line',i
  return
993 info=993
  write(psb_err_unit,*) 'READ_MATRIX: Memory allocation failure'
  return
end subroutine lcmm_mat_read


subroutine lcmm_mat_write(a,mtitle,info,iunit,filename)
  use psb_base_mod
  implicit none
  type(psb_lcspmat_type), intent(in)  :: a
  integer(psb_ipk_), intent(out)        :: info
  character(len=*), intent(in) :: mtitle
  integer(psb_ipk_), optional, intent(in)          :: iunit
  character(len=*), optional, intent(in) :: filename
  integer(psb_ipk_) :: iout


  info = psb_success_

  if (present(filename)) then 
    if (filename == '-') then 
      iout=6
    else
      if (present(iunit)) then 
        iout = iunit
      else
        iout=99
      endif
      open(iout,file=filename, err=901, action='WRITE')
    endif
  else 
    if (present(iunit)) then 
      iout = iunit   
    else
      iout=6
    endif
  endif

  call a%print(iout,head=mtitle)

  if (iout /= 6) close(iout)


  return

901 continue 
  info=901
  write(psb_err_unit,*) 'Error while opening ',filename
  return
end subroutine lcmm_mat_write
