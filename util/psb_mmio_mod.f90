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
module psb_mmio_mod
  use psb_base_mod
  public mm_mat_read, mm_mat_write
  interface mm_mat_read
    module procedure dmm_mat_read, zmm_mat_read
  end interface
  interface mm_mat_write
    module procedure dmm_mat_write,zmm_mat_write
  end interface

contains

  subroutine dmm_mat_read(a, iret, iunit, filename)   
    use psb_base_mod
    implicit none
    type(psb_dspmat_type), intent(out)  :: a
    integer, intent(out)        :: iret
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    character      :: mmheader*15, fmt*15, object*10, type*10, sym*15
    character(1024)      :: line
    integer        :: nrow, ncol, nnzero, neltvl, nrhs, nrhsix
    integer        :: ircode, i,iel,nzr,infile, j
    logical, parameter :: debug=.false.

    iret = 0

    if (present(filename)) then
      if (filename=='-') then 
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

    if ( (tolower(object) /= 'matrix').or.(tolower(fmt)/='coordinate')) then
      write(0,*) 'READ_MATRIX: input file type not yet supported'
      iret=909
      return
    end if
    if (debug) write(*,*) mmheader,':', object, ':',fmt,':', type,':', sym

    do 
      read(infile,fmt='(a)') line
      if (line(1:1) /= '%')  exit
    end do
    if (debug) write(*,*) 'Line on input : "',line,'"'
    read(line,fmt=*) nrow,ncol,nnzero
    if (debug) write(*,*) 'Out: ',nrow,ncol,nnzero
    
    if ((tolower(type) == 'real').and.(tolower(sym) == 'general')) then
      call psb_sp_all(nrow,ncol,a,nnzero,ircode)
      a%fida   = 'COO'
      a%descra = 'G'      
      if (ircode /= 0)   goto 993
      do i=1,nnzero
        read(infile,fmt=*,end=902) a%ia1(i),a%ia2(i),a%aspk(i)
      end do
      a%infoa(psb_nnz_) = nnzero
      call psb_ipcoo2csr(a,ircode)

    else if ((tolower(type) == 'real').and.(tolower(sym) == 'symmetric')) then
      ! we are generally working with non-symmetric matrices, so
      ! we de-symmetrize what we are about to read
      call psb_sp_all(nrow,ncol,a,2*nnzero,ircode)
      a%fida   = 'COO'
      a%descra = 'G'      
      if (ircode /= 0)   goto 993
      do i=1,nnzero
        read(infile,fmt=*,end=902) a%ia1(i),a%ia2(i),a%aspk(i)
      end do

      nzr = nnzero
      do i=1,nnzero
        if (a%ia1(i) /= a%ia2(i)) then 
          nzr = nzr + 1
          a%aspk(nzr) = a%aspk(i)
          a%ia1(nzr) = a%ia2(i)
          a%ia2(nzr) = a%ia1(i)
        end if
      end do
      a%infoa(psb_nnz_) = nzr
      call psb_ipcoo2csr(a,ircode)

    else
      write(0,*) 'read_matrix: matrix type not yet supported'
      iret=904
    end if
    if (infile/=5) close(infile)
    return 

    ! open failed
901 iret=901
    write(0,*) 'read_matrix: could not open file ',filename,' for input'
    return
902 iret=902
    write(0,*) 'READ_MATRIX: Unexpected end of file '
    return
993 iret=993
    write(0,*) 'READ_MATRIX: Memory allocation failure'
    return
  end subroutine dmm_mat_read





  subroutine dmm_mat_write(a,mtitle,iret,eiout,filename)
    use psb_base_mod
    implicit none
    type(psb_dspmat_type), intent(in)  :: a
    integer, intent(out)        :: iret
    character(len=*), intent(in) :: mtitle
    integer, optional, intent(in)          :: eiout
    character(len=*), optional, intent(in) :: filename
    integer                     :: iout


    iret = 0

    if (present(filename)) then 
      if (filename=='-') then 
        iout=6
      else
        if (present(eiout)) then 
          iout = eiout
        else
          iout=99
        endif
        open(iout,file=filename, err=901, action='WRITE')
      endif
    else 
      if (present(eiout)) then 
        iout = eiout   
      else
        iout=6
      endif
    endif
    
    call psb_csprt(iout,a,head=mtitle)

    if (iout /= 6) close(iout)


    return

901 continue 
    iret=901
    write(0,*) 'Error while opening ',filename
    return
  end subroutine dmm_mat_write



  subroutine zmm_mat_read(a, iret, iunit, filename)   
    use psb_base_mod
    implicit none
    type(psb_zspmat_type), intent(out)  :: a
    integer, intent(out)        :: iret
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    character      :: mmheader*15, fmt*15, object*10, type*10, sym*15
    character(1024)      :: line
    integer        :: nrow, ncol, nnzero, neltvl, nrhs, nrhsix
    integer        :: ircode, i,iel,nzr,infile,j
    real(kind(1.d0))   :: are, aim
    logical, parameter :: debug=.false.
    

    iret = 0

    if (present(filename)) then
      if (filename=='-') then 
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

    if ( (tolower(object) /= 'matrix').or.(tolower(fmt)/='coordinate')) then
      write(0,*) 'READ_MATRIX: input file type not yet supported'
      iret=909
      return
    end if
    if (debug) write(*,*) mmheader,':', object, ':',fmt,':', type,':', sym

    do 
      read(infile,fmt='(a)') line
      if (line(1:1) /= '%')  exit
    end do
    if (debug) write(*,*) 'Line on input : "',line,'"'
    read(line,fmt=*) nrow,ncol,nnzero
    if (debug) write(*,*) 'Out: ',nrow,ncol,nnzero
    
    if ((tolower(type) == 'complex').and.(tolower(sym) == 'general')) then
      call psb_sp_all(nrow,ncol,a,nnzero,ircode)
      if (ircode /= 0)   goto 993
      a%fida   = 'COO'
      a%descra = 'G'      
      do i=1,nnzero
        read(infile,fmt=*,end=902) a%ia1(i),a%ia2(i),are,aim
        a%aspk(i) = cmplx(are,aim)
      end do
      a%infoa(psb_nnz_) = nnzero
      
      call psb_ipcoo2csr(a,ircode)

    else if ((tolower(type) == 'complex').and.(tolower(sym) == 'symmetric')) then
      ! we are generally working with non-symmetric matrices, so
      ! we de-symmetrize what we are about to read
      call psb_sp_all(nrow,ncol,a,2*nnzero,ircode)
      if (ircode /= 0)   goto 993
      a%fida   = 'COO'
      a%descra = 'G'      
      do i=1,nnzero
        read(infile,fmt=*,end=902) a%ia1(i),a%ia2(i),are,aim
        a%aspk(i) = cmplx(are,aim)
      end do

      nzr = nnzero
      do i=1,nnzero
        if (a%ia1(i) /= a%ia2(i)) then 
          nzr = nzr + 1
          a%aspk(nzr) = a%aspk(i)
          a%ia1(nzr) = a%ia2(i)
          a%ia2(nzr) = a%ia1(i)
        end if
      end do
      a%infoa(psb_nnz_) = nzr
      call psb_ipcoo2csr(a,ircode)

    else if ((tolower(type) == 'complex').and.(tolower(sym) == 'hermitian')) then
      ! we are generally working with non-symmetric matrices, so
      ! we de-symmetrize what we are about to read
      call psb_sp_all(nrow,ncol,a,2*nnzero,ircode)
      if (ircode /= 0)   goto 993
      a%fida   = 'COO'
      a%descra = 'G'      
      do i=1,nnzero
        read(infile,fmt=*,end=902) a%ia1(i),a%ia2(i),are,aim
        a%aspk(i) = cmplx(are,aim)
      end do

      nzr = nnzero
      do i=1,nnzero
        if (a%ia1(i) /= a%ia2(i)) then 
          nzr = nzr + 1
          a%aspk(nzr) = conjg(a%aspk(i))
          a%ia1(nzr)  = a%ia2(i)
          a%ia2(nzr)  = a%ia1(i)
        end if
      end do
      a%infoa(psb_nnz_) = nzr
      call psb_ipcoo2csr(a,ircode)

    else
      write(0,*) 'read_matrix: matrix type not yet supported'
      iret=904
    end if
    if (infile/=5) close(infile)
    return 

    ! open failed
901 iret=901
    write(0,*) 'read_matrix: could not open file ',filename,' for input'
    return
902 iret=902
    write(0,*) 'READ_MATRIX: Unexpected end of file '
    return
993 iret=993
    write(0,*) 'READ_MATRIX: Memory allocation failure'
    return
  end subroutine zmm_mat_read



  subroutine zmm_mat_write(a,mtitle,iret,eiout,filename)
    use psb_base_mod
    implicit none
    type(psb_zspmat_type), intent(in)  :: a
    integer, intent(out)        :: iret
    character(len=*), intent(in) :: mtitle
    integer, optional, intent(in)          :: eiout
    character(len=*), optional, intent(in) :: filename
    integer                     :: iout


    iret = 0

    if (present(filename)) then 
      if (filename=='-') then 
        iout=6
      else
        if (present(eiout)) then 
          iout = eiout
        else
          iout=99
        endif
        open(iout,file=filename, err=901, action='WRITE')
      endif
    else 
      if (present(eiout)) then 
        iout = eiout   
      else
        iout=6
      endif
    endif

    call psb_csprt(iout,a,head=mtitle)

    if (iout /= 6) close(iout)


    return

901 continue 
    iret=901
    write(0,*) 'Error while opening ',filename
    return
  end subroutine zmm_mat_write


end module psb_mmio_mod
