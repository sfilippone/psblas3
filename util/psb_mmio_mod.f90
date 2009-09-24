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
module psb_mmio_mod

  public mm_mat_read, mm_mat_write, mm_vet_read
  interface mm_mat_read
    module procedure smm_mat_read,  dmm_mat_read, cmm_mat_read, zmm_mat_read
  end interface
  interface mm_mat_write
    module procedure smm_mat_write, dmm_mat_write, cmm_mat_write,  zmm_mat_write
  end interface
  interface mm_vet_read
    module procedure mm_svet_read, mm_dvet_read, mm_cvet_read, mm_zvet_read
  end interface

contains

  subroutine mm_svet_read(b, info, iunit, filename)   
    use psb_base_mod
    implicit none
    real(psb_spk_), allocatable, intent(out)  :: b(:,:)
    integer, intent(out)        :: info
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    integer              :: nrow, ncol, i,root, np,  me,  ircode, j,infile
    character            :: mmheader*15, fmt*15, object*10, type*10, sym*15,&
         & line*1024

    info = 0
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

    read(infile,fmt=*, end=902) mmheader, object, fmt, type, sym
    
    if ( (object /= 'matrix').or.(fmt /= 'array')) then
      write(0,*) 'read_rhs: input file type not yet supported'
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

    if (infile/=5) close(infile)

    return 
    ! open failed
901 write(0,*) 'mm_vet_read: could not open file ',&
         & infile,' for input'
    info = -1
    return

902 write(0,*) 'mmv_vet_read: unexpected end of file ',infile,&
         & ' during input'
    info = -2
    return
993 write(0,*) 'mm_vet_read: memory allocation failure'
    info = -3
    return
  end subroutine mm_svet_read


  subroutine mm_dvet_read(b, info, iunit, filename)   
    use psb_base_mod
    implicit none
    real(psb_dpk_), allocatable, intent(out)  :: b(:,:)
    integer, intent(out)        :: info
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    integer              :: nrow, ncol, i,root, np,  me,  ircode, j, infile
    character            :: mmheader*15, fmt*15, object*10, type*10, sym*15,&
         & line*1024

    info = 0
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

    read(infile,fmt=*, end=902) mmheader, object, fmt, type, sym
    
    if ( (object /= 'matrix').or.(fmt /= 'array')) then
      write(0,*) 'read_rhs: input file type not yet supported'
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
    if (infile/=5) close(infile)

    return 
    ! open failed
901 write(0,*) 'mm_vet_read: could not open file ',&
         & infile,' for input'
    info = -1
    return

902 write(0,*) 'mmv_vet_read: unexpected end of file ',infile,&
         & ' during input'
    info = -2
    return
993 write(0,*) 'mm_vet_read: memory allocation failure'
    info = -3
    return
  end subroutine mm_dvet_read


  subroutine mm_cvet_read(b, info, iunit, filename)   
    use psb_base_mod
    implicit none
    complex(psb_spk_), allocatable, intent(out)  :: b(:,:)
    integer, intent(out)        :: info
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    integer              :: nrow, ncol, i,root, np,  me,  ircode, j,infile
    real(psb_spk_)       :: bre, bim 
    character            :: mmheader*15, fmt*15, object*10, type*10, sym*15,&
         & line*1024

    info = 0
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

    read(infile,fmt=*, end=902) mmheader, object, fmt, type, sym
    
    if ( (object /= 'matrix').or.(fmt /= 'array')) then
      write(0,*) 'read_rhs: input file type not yet supported'
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
    if (infile/=5) close(infile)

    return 
    ! open failed
901 write(0,*) 'mm_vet_read: could not open file ',&
         & infile,' for input'
    info = -1
    return

902 write(0,*) 'mmv_vet_read: unexpected end of file ',infile,&
         & ' during input'
    info = -2
    return
993 write(0,*) 'mm_vet_read: memory allocation failure'
    info = -3
    return
  end subroutine mm_cvet_read


  subroutine mm_zvet_read(b, info, iunit, filename)   
    use psb_base_mod
    implicit none
    complex(psb_dpk_), allocatable, intent(out)  :: b(:,:)
    integer, intent(out)        :: info
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    integer              :: nrow, ncol, i,root, np,  me,  ircode, j,infile
    real(psb_dpk_)       :: bre, bim 
    character            :: mmheader*15, fmt*15, object*10, type*10, sym*15,&
         & line*1024

    info = 0
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

    read(infile,fmt=*, end=902) mmheader, object, fmt, type, sym
    
    if ( (object /= 'matrix').or.(fmt /= 'array')) then
      write(0,*) 'read_rhs: input file type not yet supported'
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
          b(i,j) = cmplx(bre,bim,kind=psb_dpk_)
        end do
      end do
      
    end if      ! read right hand sides
    if (infile/=5) close(infile)

    return 
    ! open failed
901 write(0,*) 'mm_vet_read: could not open file ',&
         & infile,' for input'
    info = -1
    return

902 write(0,*) 'mmv_vet_read: unexpected end of file ',infile,&
         & ' during input'
    info = -2
    return
993 write(0,*) 'mm_vet_read: memory allocation failure'
    info = -3
    return
  end subroutine mm_zvet_read


  subroutine smm_mat_read(a, info, iunit, filename)   
    use psb_base_mod
    implicit none
    type(psb_s_sparse_mat), intent(out)  :: a
    integer, intent(out)        :: info
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    character      :: mmheader*15, fmt*15, object*10, type*10, sym*15
    character(1024)      :: line
    integer        :: nrow, ncol, nnzero
    integer        :: ircode, i,nzr,infile
    type(psb_s_coo_sparse_mat), allocatable :: acoo

    info = 0

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

    if ( (psb_tolower(object) /= 'matrix').or.(psb_tolower(fmt)/='coordinate')) then
      write(0,*) 'READ_MATRIX: input file type not yet supported'
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
    if ((psb_tolower(type) == 'real').and.(psb_tolower(sym) == 'general')) then
      call acoo%allocate(nrow,ncol,nnzero)
      do i=1,nnzero
        read(infile,fmt=*,end=902) acoo%ia(i),acoo%ja(i),acoo%val(i)
      end do
      call acoo%set_nzeros(nnzero)
      call acoo%fix(info)

      call a%mv_from(acoo)
      call a%cscnv(ircode,type='csr')

    else if ((psb_tolower(type) == 'real').and.(psb_tolower(sym) == 'symmetric')) then
      ! we are generally working with non-symmetric matrices, so
      ! we de-symmetrize what we are about to read
      call acoo%allocate(nrow,ncol,nnzero)
      do i=1,nnzero
        read(infile,fmt=*,end=902) acoo%ia(i),acoo%ja(i),acoo%val(i)
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
      call a%cscnv(ircode,type='csr')

    else
      write(0,*) 'read_matrix: matrix type not yet supported'
      info=904
    end if
    

    if (infile/=5) close(infile)


    return 

    ! open failed
901 info=901
    write(0,*) 'read_matrix: could not open file ',filename,' for input'
    return
902 info=902
    write(0,*) 'READ_MATRIX: Unexpected end of file '
    return
993 info=993
    write(0,*) 'READ_MATRIX: Memory allocation failure'
    return
  end subroutine smm_mat_read


  subroutine smm_mat_write(a,mtitle,info,iunit,filename)
    use psb_base_mod
    implicit none
    type(psb_s_sparse_mat), intent(in)  :: a
    integer, intent(out)        :: info
    character(len=*), intent(in) :: mtitle
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    integer                     :: iout


    info = 0

    if (present(filename)) then 
      if (filename=='-') then 
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
    write(0,*) 'Error while opening ',filename
    return
  end subroutine smm_mat_write
  


  subroutine dmm_mat_read(a, info, iunit, filename)   
    use psb_base_mod
    implicit none
    type(psb_d_sparse_mat), intent(out)  :: a
    integer, intent(out)        :: info
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    character      :: mmheader*15, fmt*15, object*10, type*10, sym*15
    character(1024)      :: line
    integer        :: nrow, ncol, nnzero
    integer        :: ircode, i,nzr,infile
    type(psb_d_coo_sparse_mat), allocatable :: acoo

    info = 0

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

    if ( (psb_tolower(object) /= 'matrix').or.(psb_tolower(fmt)/='coordinate')) then
      write(0,*) 'READ_MATRIX: input file type not yet supported'
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
    if ((psb_tolower(type) == 'real').and.(psb_tolower(sym) == 'general')) then
      call acoo%allocate(nrow,ncol,nnzero)
      do i=1,nnzero
        read(infile,fmt=*,end=902) acoo%ia(i),acoo%ja(i),acoo%val(i)
      end do
      call acoo%set_nzeros(nnzero)
      call acoo%fix(info)

      call a%mv_from(acoo)
      call a%cscnv(ircode,type='csr')

    else if ((psb_tolower(type) == 'real').and.(psb_tolower(sym) == 'symmetric')) then
      ! we are generally working with non-symmetric matrices, so
      ! we de-symmetrize what we are about to read
      call acoo%allocate(nrow,ncol,nnzero)
      do i=1,nnzero
        read(infile,fmt=*,end=902) acoo%ia(i),acoo%ja(i),acoo%val(i)
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
      call a%cscnv(ircode,type='csr')

    else
      write(0,*) 'read_matrix: matrix type not yet supported'
      info=904
    end if
    

    if (infile/=5) close(infile)


    return 

    ! open failed
901 info=901
    write(0,*) 'read_matrix: could not open file ',filename,' for input'
    return
902 info=902
    write(0,*) 'READ_MATRIX: Unexpected end of file '
    return
993 info=993
    write(0,*) 'READ_MATRIX: Memory allocation failure'
    return
  end subroutine dmm_mat_read


  subroutine dmm_mat_write(a,mtitle,info,iunit,filename)
    use psb_base_mod
    implicit none
    type(psb_d_sparse_mat), intent(in)  :: a
    integer, intent(out)        :: info
    character(len=*), intent(in) :: mtitle
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    integer                     :: iout


    info = 0

    if (present(filename)) then 
      if (filename=='-') then 
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
    write(0,*) 'Error while opening ',filename
    return
  end subroutine dmm_mat_write

  subroutine cmm_mat_read(a, info, iunit, filename)   
    use psb_base_mod
    implicit none
    type(psb_cspmat_type), intent(out)  :: a
    integer, intent(out)        :: info
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    character      :: mmheader*15, fmt*15, object*10, type*10, sym*15
    character(1024)      :: line
    integer        :: nrow, ncol, nnzero
    integer        :: ircode, i,nzr,infile
    real(psb_spk_)   :: are, aim

    info = 0

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

    if ( (psb_tolower(object) /= 'matrix').or.(psb_tolower(fmt)/='coordinate')) then
      write(0,*) 'READ_MATRIX: input file type not yet supported'
      info=909
      return
    end if

    do 
      read(infile,fmt='(a)') line
      if (line(1:1) /= '%')  exit
    end do
    read(line,fmt=*) nrow,ncol,nnzero
    
    if ((psb_tolower(type) == 'complex').and.(psb_tolower(sym) == 'general')) then
      call psb_sp_all(nrow,ncol,a,nnzero,ircode)
      if (ircode /= 0)   goto 993
      a%fida   = 'COO'
      a%descra = 'G'      
      do i=1,nnzero
        read(infile,fmt=*,end=902) a%ia1(i),a%ia2(i),are,aim
        a%aspk(i) = cmplx(are,aim,kind=psb_spk_)
      end do
      a%infoa(psb_nnz_) = nnzero
      
      call psb_spcnv(a,ircode,afmt='csr')

    else if ((psb_tolower(type) == 'complex').and.(psb_tolower(sym) == 'symmetric')) then
      ! we are generally working with non-symmetric matrices, so
      ! we de-symmetrize what we are about to read
      call psb_sp_all(nrow,ncol,a,2*nnzero,ircode)
      if (ircode /= 0)   goto 993
      a%fida   = 'COO'
      a%descra = 'G'      
      do i=1,nnzero
        read(infile,fmt=*,end=902) a%ia1(i),a%ia2(i),are,aim
        a%aspk(i) = cmplx(are,aim,kind=psb_spk_)
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
      call psb_spcnv(a,ircode,afmt='csr')

    else if ((psb_tolower(type) == 'complex').and.(psb_tolower(sym) == 'hermitian')) then
      ! we are generally working with non-symmetric matrices, so
      ! we de-symmetrize what we are about to read
      call psb_sp_all(nrow,ncol,a,2*nnzero,ircode)
      if (ircode /= 0)   goto 993
      a%fida   = 'COO'
      a%descra = 'G'      
      do i=1,nnzero
        read(infile,fmt=*,end=902) a%ia1(i),a%ia2(i),are,aim
        a%aspk(i) = cmplx(are,aim,kind=psb_spk_)
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
      call psb_spcnv(a,ircode,afmt='csr')

    else
      write(0,*) 'read_matrix: matrix type not yet supported'
      info=904
    end if
    if (infile/=5) close(infile)
    return 

    ! open failed
901 info=901
    write(0,*) 'read_matrix: could not open file ',filename,' for input'
    return
902 info=902
    write(0,*) 'READ_MATRIX: Unexpected end of file '
    return
993 info=993
    write(0,*) 'READ_MATRIX: Memory allocation failure'
    return
  end subroutine cmm_mat_read



  subroutine cmm_mat_write(a,mtitle,info,iunit,filename)
    use psb_base_mod
    implicit none
    type(psb_cspmat_type), intent(in)  :: a
    integer, intent(out)        :: info
    character(len=*), intent(in) :: mtitle
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    integer                     :: iout


    info = 0

    if (present(filename)) then 
      if (filename=='-') then 
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

    call psb_csprt(iout,a,head=mtitle)

    if (iout /= 6) close(iout)


    return

901 continue 
    info=901
    write(0,*) 'Error while opening ',filename
    return
  end subroutine cmm_mat_write

  subroutine zmm_mat_read(a, info, iunit, filename)   
    use psb_base_mod
    implicit none
    type(psb_zspmat_type), intent(out)  :: a
    integer, intent(out)        :: info
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    character      :: mmheader*15, fmt*15, object*10, type*10, sym*15
    character(1024)      :: line
    integer        :: nrow, ncol, nnzero
    integer        :: ircode, i,nzr,infile
    real(psb_dpk_)   :: are, aim

    info = 0

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

    if ( (psb_tolower(object) /= 'matrix').or.(psb_tolower(fmt)/='coordinate')) then
      write(0,*) 'READ_MATRIX: input file type not yet supported'
      info=909
      return
    end if

    do 
      read(infile,fmt='(a)') line
      if (line(1:1) /= '%')  exit
    end do
    read(line,fmt=*) nrow,ncol,nnzero
    
    if ((psb_tolower(type) == 'complex').and.(psb_tolower(sym) == 'general')) then
      call psb_sp_all(nrow,ncol,a,nnzero,ircode)
      if (ircode /= 0)   goto 993
      a%fida   = 'COO'
      a%descra = 'G'      
      do i=1,nnzero
        read(infile,fmt=*,end=902) a%ia1(i),a%ia2(i),are,aim
        a%aspk(i) = cmplx(are,aim,kind=psb_dpk_)
      end do
      a%infoa(psb_nnz_) = nnzero
      
      call psb_spcnv(a,ircode,afmt='csr')

    else if ((psb_tolower(type) == 'complex').and.(psb_tolower(sym) == 'symmetric')) then
      ! we are generally working with non-symmetric matrices, so
      ! we de-symmetrize what we are about to read
      call psb_sp_all(nrow,ncol,a,2*nnzero,ircode)
      if (ircode /= 0)   goto 993
      a%fida   = 'COO'
      a%descra = 'G'      
      do i=1,nnzero
        read(infile,fmt=*,end=902) a%ia1(i),a%ia2(i),are,aim
        a%aspk(i) = cmplx(are,aim,kind=psb_dpk_)
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
      call psb_spcnv(a,ircode,afmt='csr')

    else if ((psb_tolower(type) == 'complex').and.(psb_tolower(sym) == 'hermitian')) then
      ! we are generally working with non-symmetric matrices, so
      ! we de-symmetrize what we are about to read
      call psb_sp_all(nrow,ncol,a,2*nnzero,ircode)
      if (ircode /= 0)   goto 993
      a%fida   = 'COO'
      a%descra = 'G'      
      do i=1,nnzero
        read(infile,fmt=*,end=902) a%ia1(i),a%ia2(i),are,aim
        a%aspk(i) = cmplx(are,aim,kind=psb_dpk_)
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
      call psb_spcnv(a,ircode,afmt='csr')

    else
      write(0,*) 'read_matrix: matrix type not yet supported'
      info=904
    end if
    if (infile/=5) close(infile)
    return 

    ! open failed
901 info=901
    write(0,*) 'read_matrix: could not open file ',filename,' for input'
    return
902 info=902
    write(0,*) 'READ_MATRIX: Unexpected end of file '
    return
993 info=993
    write(0,*) 'READ_MATRIX: Memory allocation failure'
    return
  end subroutine zmm_mat_read



  subroutine zmm_mat_write(a,mtitle,info,iunit,filename)
    use psb_base_mod
    implicit none
    type(psb_zspmat_type), intent(in)  :: a
    integer, intent(out)        :: info
    character(len=*), intent(in) :: mtitle
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    integer                     :: iout


    info = 0

    if (present(filename)) then 
      if (filename=='-') then 
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

    call psb_csprt(iout,a,head=mtitle)

    if (iout /= 6) close(iout)


    return

901 continue 
    info=901
    write(0,*) 'Error while opening ',filename
    return
  end subroutine zmm_mat_write


end module psb_mmio_mod
