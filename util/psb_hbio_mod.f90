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
module psb_hbio_mod

  public hb_read, hb_write
  interface hb_read
    module procedure shb_read,  dhb_read, zhb_read !chb_read
  end interface
  interface hb_write
    module procedure shb_write, dhb_write,zhb_write !chb_write
  end interface

contains

  subroutine shb_read(a, iret, iunit, filename,b,mtitle)   
    use psb_base_mod
    implicit none
    type(psb_sspmat_type), intent(out)     :: a
    integer, intent(out)                   :: iret
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    real(psb_spk_), optional, allocatable :: b(:) 
    character(len=72), optional, intent(out) :: mtitle

    character  :: rhstype,type*3,key*8
    character indfmt*16,ptrfmt*16,rhsfmt*20,valfmt*20
    integer        :: indcrd,  ptrcrd, totcrd,&
         & valcrd, rhscrd, nrow, ncol, nnzero, neltvl, nrhs, nrhsix
    integer                     :: ircode, i,nzr,infile, info
    character(len=*), parameter :: fmt10='(a72,a8,/,5i14,/,a3,11x,4i14,/,2a16,2a20)'
    character(len=*), parameter :: fmt11='(a1,13x,2i14)'
    character(len=*), parameter :: fmt111='(1x,a8,1x,i8,1x,a10)'

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

    read (infile,fmt=fmt10) mtitle,key,totcrd,ptrcrd,indcrd,valcrd,rhscrd,&
         & type,nrow,ncol,nnzero,neltvl,ptrfmt,indfmt,valfmt,rhsfmt
    if (rhscrd > 0) read(infile,fmt=fmt11)rhstype,nrhs,nrhsix

    call psb_sp_all(a,nnzero,nrow+1,nnzero,ircode)
    if (ircode /= 0 ) then 
      write(0,*) 'Memory allocation failed'
      goto 993
    end if

    a%m    = nrow
    a%k    = ncol
    a%fida = 'CSR'
    a%descra='G'


    if (psb_tolower(type(1:1)) == 'r') then 
      if (psb_tolower(type(2:2)) == 'u') then 


        read (infile,fmt=ptrfmt) (a%ia2(i),i=1,nrow+1)
        read (infile,fmt=indfmt) (a%ia1(i),i=1,nnzero)
        if (valcrd > 0) read (infile,fmt=valfmt) (a%aspk(i),i=1,nnzero)

        if (present(b)  .and. (rhscrd > 0)) then
          call psb_ensure_size(nrow,b,info)
          if (info == 0) read (infile,fmt=rhsfmt) (b(i),i=1,nrow)
        endif

      else if (psb_tolower(type(2:2)) == 's') then 

        ! we are generally working with non-symmetric matrices, so
        ! we de-symmetrize what we are about to read

        read (infile,fmt=ptrfmt) (a%ia2(i),i=1,nrow+1)
        read (infile,fmt=indfmt) (a%ia1(i),i=1,nnzero)
        if (valcrd > 0) read (infile,fmt=valfmt) (a%aspk(i),i=1,nnzero)

        if (present(b)  .and. (rhscrd > 0)) then
          call psb_ensure_size(nrow,b,info)
          if (info == 0) read (infile,fmt=rhsfmt) (b(i),i=1,nrow)
        endif

        call psb_spcnv(a,ircode,afmt='csr')
        if (ircode /= 0)   goto 993  

        call psb_sp_reall(a,2*nnzero,ircode)
        ! A is now in COO format
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
        if (ircode /= 0)   goto 993

      else
        write(0,*) 'read_matrix: matrix type not yet supported'
        iret=904
      end if
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
    write(0,*) 'DHB_READ: Unexpected end of file '
    return
993 iret=993
    write(0,*) 'DHB_READ: Memory allocation failure'
    return
  end subroutine shb_read


  subroutine shb_write(a,iret,iunit,filename,key,rhs,mtitle)
    use psb_base_mod
    implicit none
    type(psb_sspmat_type), intent(in)  :: a
    integer, intent(out)        :: iret
    character(len=*), optional, intent(in) :: mtitle
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    character(len=*), optional, intent(in) :: key
    real(psb_spk_), optional             :: rhs(:)
    integer                     :: iout

    character(len=*), parameter::  ptrfmt='(10I8)',indfmt='(10I8)'
    integer, parameter :: jptr=10,jind=10
    character(len=*), parameter::  valfmt='(4E20.12)',rhsfmt='(4E20.12)'
    integer, parameter :: jval=4,jrhs=4
    character(len=*), parameter :: fmt10='(a72,a8,/,5i14,/,a3,11x,4i14,/,2a16,2a20)'
    character(len=*), parameter :: fmt11='(a1,13x,2i14)'
    character(len=*), parameter :: fmt111='(1x,a8,1x,i8,1x,a10)'
    
    character(len=72) :: mtitle_
    character(len=8)  :: key_

    character  :: rhstype,type*3

    integer    :: i,indcrd,ptrcrd,rhscrd,totcrd,valcrd,&
         & nrow,ncol,nnzero, neltvl, nrhs, nrhsix

    iret = 0

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
    
    if (present(mtitle)) then 
      mtitle_ = mtitle
    else
      mtitle_ = 'Temporary PSBLAS title '
    endif
    if (present(key)) then 
      key_ = key
    else
      key_ = 'PSBMAT00'
    endif

    if (psb_toupper(a%fida) == 'CSR') then 
      
      nrow   = a%m
      ncol   = a%k
      nnzero = a%ia2(nrow+1)-1
      
      neltvl = 0 

      ptrcrd = (nrow+1)/jptr
      if (mod(nrow+1,jptr) > 0) ptrcrd = ptrcrd + 1
      indcrd = nnzero/jind
      if (mod(nnzero,jind) > 0) indcrd = indcrd + 1
      valcrd = nnzero/jval
      if (mod(nnzero,jval) > 0) valcrd = valcrd + 1
      if (present(rhs)) then 
        if (size(rhs)<nrow) then 
          rhscrd = 0
        else
          rhscrd = nrow/jrhs
          if (mod(nrow,jrhs) > 0) rhscrd = rhscrd + 1
        endif
        nrhs = 1
      else
        rhscrd = 0
        nrhs   = 0
      end if
      totcrd = ptrcrd + indcrd + valcrd + rhscrd
      
      nrhsix  = nrhs*nrow
      rhstype = 'F'
      type    = 'RUA'
      
      write (iout,fmt=fmt10) mtitle_,key_,totcrd,ptrcrd,indcrd,valcrd,rhscrd,&
           &  type,nrow,ncol,nnzero,neltvl,ptrfmt,indfmt,valfmt,rhsfmt
      if (rhscrd > 0) write (iout,fmt=fmt11) rhstype,nrhs,nrhsix
      write (iout,fmt=ptrfmt) (a%ia2(i),i=1,nrow+1)
      write (iout,fmt=indfmt) (a%ia1(i),i=1,nnzero)
      if (valcrd > 0) write (iout,fmt=valfmt) (a%aspk(i),i=1,nnzero)
      if (rhscrd > 0) write (iout,fmt=rhsfmt) (rhs(i),i=1,nrow)


    else

      write(0,*) 'format: ',a%fida,' not yet implemented'

    endif

    if (iout /= 6) close(iout)


    return

901 continue 
    iret=901
    write(0,*) 'Error while opening ',filename
    return
  end subroutine shb_write

  subroutine dhb_read(a, iret, iunit, filename,b,mtitle)   
    use psb_base_mod
    implicit none
    type(psb_dspmat_type), intent(out)     :: a
    integer, intent(out)                   :: iret
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    real(psb_dpk_), optional, allocatable :: b(:) 
    character(len=72), optional, intent(out) :: mtitle

    character  :: rhstype,type*3,key*8
    character indfmt*16,ptrfmt*16,rhsfmt*20,valfmt*20
    integer        :: indcrd,  ptrcrd, totcrd,&
         & valcrd, rhscrd, nrow, ncol, nnzero, neltvl, nrhs, nrhsix
    integer                     :: ircode, i,nzr,infile, info
    character(len=*), parameter :: fmt10='(a72,a8,/,5i14,/,a3,11x,4i14,/,2a16,2a20)'
    character(len=*), parameter :: fmt11='(a1,13x,2i14)'
    character(len=*), parameter :: fmt111='(1x,a8,1x,i8,1x,a10)'

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

    read (infile,fmt=fmt10) mtitle,key,totcrd,ptrcrd,indcrd,valcrd,rhscrd,&
         & type,nrow,ncol,nnzero,neltvl,ptrfmt,indfmt,valfmt,rhsfmt
    if (rhscrd > 0) read(infile,fmt=fmt11)rhstype,nrhs,nrhsix

    call psb_sp_all(a,nnzero,nrow+1,nnzero,ircode)
    if (ircode /= 0 ) then 
      write(0,*) 'Memory allocation failed'
      goto 993
    end if

    a%m    = nrow
    a%k    = ncol
    a%fida = 'CSR'
    a%descra='G'


    if (psb_tolower(type(1:1)) == 'r') then 
      if (psb_tolower(type(2:2)) == 'u') then 


        read (infile,fmt=ptrfmt) (a%ia2(i),i=1,nrow+1)
        read (infile,fmt=indfmt) (a%ia1(i),i=1,nnzero)
        if (valcrd > 0) read (infile,fmt=valfmt) (a%aspk(i),i=1,nnzero)

        if (present(b)  .and. (rhscrd > 0)) then
          call psb_ensure_size(nrow,b,info)
          if (info == 0) read (infile,fmt=rhsfmt) (b(i),i=1,nrow)
        endif

      else if (psb_tolower(type(2:2)) == 's') then 

        ! we are generally working with non-symmetric matrices, so
        ! we de-symmetrize what we are about to read

        read (infile,fmt=ptrfmt) (a%ia2(i),i=1,nrow+1)
        read (infile,fmt=indfmt) (a%ia1(i),i=1,nnzero)
        if (valcrd > 0) read (infile,fmt=valfmt) (a%aspk(i),i=1,nnzero)

        if (present(b)  .and. (rhscrd > 0)) then
          call psb_ensure_size(nrow,b,info)
          if (info == 0) read (infile,fmt=rhsfmt) (b(i),i=1,nrow)
        endif

        call psb_spcnv(a,ircode,afmt='csr')
        if (ircode /= 0)   goto 993  

        call psb_sp_reall(a,2*nnzero,ircode)
        ! A is now in COO format
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
        if (ircode /= 0)   goto 993

      else
        write(0,*) 'read_matrix: matrix type not yet supported'
        iret=904
      end if
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
    write(0,*) 'DHB_READ: Unexpected end of file '
    return
993 iret=993
    write(0,*) 'DHB_READ: Memory allocation failure'
    return
  end subroutine dhb_read




  subroutine dhb_write(a,iret,iunit,filename,key,rhs,mtitle)
    use psb_base_mod
    implicit none
    type(psb_dspmat_type), intent(in)  :: a
    integer, intent(out)        :: iret
    character(len=*), optional, intent(in) :: mtitle
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    character(len=*), optional, intent(in) :: key
    real(psb_dpk_), optional             :: rhs(:)
    integer                     :: iout

    character(len=*), parameter::  ptrfmt='(10I8)',indfmt='(10I8)'
    integer, parameter :: jptr=10,jind=10
    character(len=*), parameter::  valfmt='(4E20.12)',rhsfmt='(4E20.12)'
    integer, parameter :: jval=4,jrhs=4
    character(len=*), parameter :: fmt10='(a72,a8,/,5i14,/,a3,11x,4i14,/,2a16,2a20)'
    character(len=*), parameter :: fmt11='(a1,13x,2i14)'
    character(len=*), parameter :: fmt111='(1x,a8,1x,i8,1x,a10)'
    
    character(len=72) :: mtitle_
    character(len=8)  :: key_

    character  :: rhstype,type*3

    integer    :: i,indcrd,ptrcrd,rhscrd,totcrd,valcrd,&
         & nrow,ncol,nnzero, neltvl, nrhs, nrhsix

    iret = 0

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
    
    if (present(mtitle)) then 
      mtitle_ = mtitle
    else
      mtitle_ = 'Temporary PSBLAS title '
    endif
    if (present(key)) then 
      key_ = key
    else
      key_ = 'PSBMAT00'
    endif

    if (psb_toupper(a%fida) == 'CSR') then 
      
      nrow   = a%m
      ncol   = a%k
      nnzero = a%ia2(nrow+1)-1
      
      neltvl = 0 

      ptrcrd = (nrow+1)/jptr
      if (mod(nrow+1,jptr) > 0) ptrcrd = ptrcrd + 1
      indcrd = nnzero/jind
      if (mod(nnzero,jind) > 0) indcrd = indcrd + 1
      valcrd = nnzero/jval
      if (mod(nnzero,jval) > 0) valcrd = valcrd + 1
      if (present(rhs)) then 
        if (size(rhs)<nrow) then 
          rhscrd = 0
        else
          rhscrd = nrow/jrhs
          if (mod(nrow,jrhs) > 0) rhscrd = rhscrd + 1
        endif
        nrhs = 1
      else
        rhscrd = 0
        nrhs   = 0
      end if
      totcrd = ptrcrd + indcrd + valcrd + rhscrd
      
      nrhsix  = nrhs*nrow
      rhstype = 'F'
      type    = 'RUA'
      
      write (iout,fmt=fmt10) mtitle_,key_,totcrd,ptrcrd,indcrd,valcrd,rhscrd,&
           &  type,nrow,ncol,nnzero,neltvl,ptrfmt,indfmt,valfmt,rhsfmt
      if (rhscrd > 0) write (iout,fmt=fmt11) rhstype,nrhs,nrhsix
      write (iout,fmt=ptrfmt) (a%ia2(i),i=1,nrow+1)
      write (iout,fmt=indfmt) (a%ia1(i),i=1,nnzero)
      if (valcrd > 0) write (iout,fmt=valfmt) (a%aspk(i),i=1,nnzero)
      if (rhscrd > 0) write (iout,fmt=rhsfmt) (rhs(i),i=1,nrow)


    else

      write(0,*) 'format: ',a%fida,' not yet implemented'

    endif

    if (iout /= 6) close(iout)


    return

901 continue 
    iret=901
    write(0,*) 'Error while opening ',filename
    return
  end subroutine dhb_write

!!$  subroutine chb_read(a, iret, iunit, filename,b,mtitle)   
!!$    use psb_base_mod
!!$    implicit none
!!$    type(psb_cspmat_type), intent(out)     :: a
!!$    integer, intent(out)                   :: iret
!!$    integer, optional, intent(in)          :: iunit
!!$    character(len=*), optional, intent(in) :: filename
!!$    real(psb_spk_), optional, allocatable  :: b(:) 
!!$    character(len=72), optional, intent(out) :: mtitle
!!$
!!$    character  :: rhstype,type*3,key*8
!!$    character indfmt*16,ptrfmt*16,rhsfmt*20,valfmt*20
!!$    integer        :: indcrd,  ptrcrd, totcrd,&
!!$         & valcrd, rhscrd, nrow, ncol, nnzero, neltvl, nrhs, nrhsix
!!$    integer                     :: ircode, i,nzr,infile,info
!!$    character(len=*), parameter :: fmt10='(a72,a8,/,5i14,/,a3,11x,4i14,/,2a16,2a20)'
!!$    character(len=*), parameter :: fmt11='(a1,13x,2i14)'
!!$    character(len=*), parameter :: fmt111='(1x,a8,1x,i8,1x,a10)'
!!$
!!$    iret = 0
!!$
!!$    if (present(filename)) then
!!$      if (filename=='-') then 
!!$        infile=5
!!$      else
!!$        if (present(iunit)) then 
!!$          infile=iunit
!!$        else
!!$          infile=99
!!$        endif
!!$        open(infile,file=filename, status='OLD', err=901, action='READ')
!!$      endif
!!$    else 
!!$      if (present(iunit)) then 
!!$        infile=iunit
!!$      else
!!$        infile=5
!!$      endif
!!$    endif
!!$
!!$    read (infile,fmt=fmt10) mtitle,key,totcrd,ptrcrd,indcrd,valcrd,rhscrd,&
!!$         & type,nrow,ncol,nnzero,neltvl,ptrfmt,indfmt,valfmt,rhsfmt
!!$    if (rhscrd > 0) read(infile,fmt=fmt11)rhstype,nrhs,nrhsix
!!$
!!$
!!$    if (psb_tolower(type(1:1)) == 'c') then 
!!$      if (psb_tolower(type(2:2)) == 'u') then 
!!$
!!$        call psb_sp_all(a,nnzero,nrow+1,nnzero,ircode)
!!$        
!!$        if (ircode /= 0 ) then 
!!$          write(0,*) 'Memory allocation failed'
!!$          goto 993
!!$        end if
!!$        
!!$        a%m    = nrow
!!$        a%k    = ncol
!!$        a%fida = 'CSR'
!!$        a%descra='G'
!!$
!!$        read (infile,fmt=ptrfmt) (a%ia2(i),i=1,nrow+1)
!!$        read (infile,fmt=indfmt) (a%ia1(i),i=1,nnzero)
!!$        if (valcrd > 0) read (infile,fmt=valfmt) (a%aspk(i),i=1,nnzero)
!!$
!!$        if (present(b)  .and. (rhscrd > 0)) then
!!$          call psb_ensure_size(nrow,b,info)
!!$          if (info == 0) read (infile,fmt=rhsfmt) (b(i),i=1,nrow)
!!$        endif
!!$
!!$      else if (psb_tolower(type(2:2)) == 's') then 
!!$
!!$        ! we are generally working with non-symmetric matrices, so
!!$        ! we de-symmetrize what we are about to read
!!$        
!!$        call psb_sp_all(nrow,ncol,a,nnzero,ircode)
!!$
!!$        if (ircode /= 0 ) then 
!!$          write(0,*) 'Memory allocation failed'
!!$          goto 993
!!$        end if
!!$        a%fida = 'CSR'
!!$        a%descra='G'
!!$        
!!$
!!$        read (infile,fmt=ptrfmt) (a%ia2(i),i=1,nrow+1)
!!$        read (infile,fmt=indfmt) (a%ia1(i),i=1,nnzero)
!!$        if (valcrd > 0) read (infile,fmt=valfmt) (a%aspk(i),i=1,nnzero)
!!$
!!$        if (present(b)  .and. (rhscrd > 0)) then
!!$          call psb_ensure_size(nrow,b,info)
!!$          if (info == 0) read (infile,fmt=rhsfmt) (b(i),i=1,nrow)
!!$        endif
!!$        
!!$
!!$        call psb_spcnv(a,ircode,afmt='coo')
!!$        if (ircode /= 0) then
!!$          write(0,*) 'ipcsr2coo ',ircode
!!$          goto 993  
!!$        end if
!!$
!!$        call psb_sp_reall(a,2*nnzero,ircode)
!!$      
!!$        if (ircode /= 0 ) then 
!!$          write(0,*) 'Memory allocation failed'
!!$          goto 993
!!$        end if
!!$
!!$        ! A is now in COO format
!!$        nzr = nnzero
!!$        do i=1,nnzero
!!$          if (a%ia1(i) /= a%ia2(i)) then 
!!$            nzr = nzr + 1
!!$            a%aspk(nzr) = a%aspk(i)
!!$            a%ia1(nzr) = a%ia2(i)
!!$            a%ia2(nzr) = a%ia1(i)
!!$          end if
!!$        end do
!!$        a%infoa(psb_nnz_) = nzr
!!$        call psb_spcnv(a,ircode,afmt='csr')
!!$        if (ircode /= 0) then
!!$          write(0,*) 'ipcoo2csr ',ircode
!!$          goto 993  
!!$        end if
!!$
!!$      else if (psb_tolower(type(2:2)) == 'h') then 
!!$
!!$        ! we are generally working with non-symmetric matrices, so
!!$        ! we de-symmetrize what we are about to read
!!$        
!!$
!!$        call psb_sp_all(nrow,ncol,a,2*nnzero,ircode)
!!$        
!!$        if (ircode /= 0 ) then 
!!$          write(0,*) 'Memory allocation failed'
!!$          goto 993
!!$        end if
!!$        a%fida = 'CSR'
!!$        a%descra='G'
!!$        
!!$
!!$        read (infile,fmt=ptrfmt) (a%ia2(i),i=1,nrow+1)
!!$        read (infile,fmt=indfmt) (a%ia1(i),i=1,nnzero)
!!$        if (valcrd > 0) read (infile,fmt=valfmt) (a%aspk(i),i=1,nnzero)
!!$
!!$        if (present(b)  .and. (rhscrd > 0)) then
!!$          call psb_ensure_size(nrow,b,info)
!!$          if (info == 0) read (infile,fmt=rhsfmt) (b(i),i=1,nrow)
!!$        endif
!!$        
!!$        call psb_spcnv(a,ircode,afmt='coo')
!!$        if (ircode /= 0) then
!!$          write(0,*) 'ipcsr2coo ',ircode
!!$          goto 993  
!!$        end if
!!$        call psb_sp_reall(a,2*nnzero,ircode)
!!$      
!!$        if (ircode /= 0 ) then 
!!$          write(0,*) 'Memory allocation failed'
!!$          goto 993
!!$        end if
!!$
!!$        ! A is now in COO format
!!$        nzr = nnzero
!!$        do i=1,nnzero
!!$          if (a%ia1(i) /= a%ia2(i)) then 
!!$            nzr = nzr + 1
!!$            a%aspk(nzr) = conjg(a%aspk(i))
!!$            a%ia1(nzr) = a%ia2(i)
!!$            a%ia2(nzr) = a%ia1(i)
!!$          end if
!!$        end do
!!$        a%infoa(psb_nnz_) = nzr
!!$        call psb_spcnv(a,ircode,afmt='csr')
!!$        if (ircode /= 0) then
!!$          write(0,*) 'ipcoo2csr ',ircode
!!$          goto 993  
!!$        end if
!!$
!!$      else
!!$        write(0,*) 'read_matrix: matrix type not yet supported'
!!$        iret=904
!!$      end if
!!$    else
!!$      write(0,*) 'read_matrix: matrix type not yet supported'
!!$      iret=904
!!$    end if
!!$    if (infile/=5) close(infile)
!!$
!!$    return 
!!$
!!$    ! open failed
!!$901 iret=901
!!$    write(0,*) 'read_matrix: could not open file ',filename,' for input'
!!$    return
!!$902 iret=902
!!$    write(0,*) 'ZHB_READ: Unexpected end of file '
!!$    return
!!$993 iret=993
!!$    write(0,*) 'ZHB_READ: Memory allocation failure'
!!$    return
!!$  end subroutine chb_read
!!$
!!$
!!$
!!$
!!$  subroutine chb_write(a,iret,iunit,filename,key,rhs,mtitle)
!!$    use psb_base_mod
!!$    implicit none
!!$    type(psb_cspmat_type), intent(in)  :: a
!!$    integer, intent(out)        :: iret
!!$    character(len=*), optional, intent(in) :: mtitle
!!$    integer, optional, intent(in)          :: iunit
!!$    character(len=*), optional, intent(in) :: filename
!!$    character(len=*), optional, intent(in) :: key
!!$    complex(psb_spk_), optional             :: rhs(:)
!!$    integer                     :: iout
!!$
!!$    character(len=*), parameter::  ptrfmt='(10I8)',indfmt='(10I8)'
!!$    integer, parameter :: jptr=10,jind=10
!!$    character(len=*), parameter::  valfmt='(4E20.12)',rhsfmt='(4E20.12)'
!!$    integer, parameter :: jval=4,jrhs=4
!!$    character(len=*), parameter :: fmt10='(a72,a8,/,5i14,/,a3,11x,4i14,/,2a16,2a20)'
!!$    character(len=*), parameter :: fmt11='(a1,13x,2i14)'
!!$    character(len=*), parameter :: fmt111='(1x,a8,1x,i8,1x,a10)'
!!$    
!!$    character(len=72) :: mtitle_
!!$    character(len=8)  :: key_
!!$
!!$    character  :: rhstype,type*3
!!$
!!$    integer    :: i,indcrd,ptrcrd,rhscrd,totcrd,valcrd,&
!!$         & nrow,ncol,nnzero, neltvl,nrhs,nrhsix
!!$
!!$    iret = 0
!!$
!!$    if (present(filename)) then 
!!$      if (filename=='-') then 
!!$        iout=6
!!$      else
!!$        if (present(iunit)) then 
!!$          iout = iunit
!!$        else
!!$          iout=99
!!$        endif
!!$        open(iout,file=filename, err=901, action='WRITE')
!!$      endif
!!$    else 
!!$      if (present(iunit)) then 
!!$        iout = iunit   
!!$      else
!!$        iout=6
!!$      endif
!!$    endif
!!$    
!!$    if (present(mtitle)) then 
!!$      mtitle_ = mtitle
!!$    else
!!$      mtitle_ = 'Temporary PSBLAS title '
!!$    endif
!!$    if (present(key)) then 
!!$      key_ = key
!!$    else
!!$      key_ = 'PSBMAT00'
!!$    endif
!!$    if (psb_toupper(a%fida) == 'CSR') then 
!!$      
!!$      nrow   = a%m
!!$      ncol   = a%k
!!$      nnzero = a%ia2(nrow+1)-1
!!$      neltvl = 0 
!!$      ptrcrd = (nrow+1)/jptr
!!$      if (mod(nrow+1,jptr) > 0) ptrcrd = ptrcrd + 1
!!$      indcrd = nnzero/jind
!!$      if (mod(nnzero,jind) > 0) indcrd = indcrd + 1
!!$      valcrd = nnzero/jval
!!$      if (mod(nnzero,jval) > 0) valcrd = valcrd + 1
!!$      if (present(rhs)) then 
!!$        if (size(rhs)<nrow) then 
!!$          rhscrd = 0
!!$        else
!!$          rhscrd = nrow/jrhs
!!$          if (mod(nrow,jrhs) > 0) rhscrd = rhscrd + 1
!!$        endif
!!$        nrhs = 1
!!$      else
!!$        rhscrd = 0
!!$        nrhs = 0
!!$      end if
!!$      totcrd = ptrcrd + indcrd + valcrd + rhscrd
!!$      nrhsix = nrhs * nrow
!!$      rhstype='F'
!!$      type='CUA'
!!$      
!!$      write (iout,fmt=fmt10) mtitle_,key_,totcrd,ptrcrd,indcrd,valcrd,rhscrd,&
!!$           &  type,nrow,ncol,nnzero,neltvl,ptrfmt,indfmt,valfmt,rhsfmt
!!$      if (rhscrd > 0) write (iout,fmt=fmt11) rhstype,nrhs,nrhsix
!!$      write (iout,fmt=ptrfmt) (a%ia2(i),i=1,nrow+1)
!!$      write (iout,fmt=indfmt) (a%ia1(i),i=1,nnzero)
!!$      if (valcrd > 0) write (iout,fmt=valfmt) (a%aspk(i),i=1,nnzero)
!!$      if (rhscrd > 0) write (iout,fmt=rhsfmt) (rhs(i),i=1,nrow)
!!$
!!$
!!$    else
!!$
!!$      write(0,*) 'format: ',a%fida,' not yet implemented'
!!$
!!$    endif
!!$
!!$    if (iout /= 6) close(iout)
!!$
!!$
!!$    return
!!$
!!$901 continue 
!!$    iret=901
!!$    write(0,*) 'Error while opening ',filename
!!$    return
!!$  end subroutine chb_write


  subroutine zhb_read(a, iret, iunit, filename,b,mtitle)   
    use psb_base_mod
    implicit none
    type(psb_zspmat_type), intent(out)     :: a
    integer, intent(out)                   :: iret
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    real(psb_dpk_), optional, allocatable  :: b(:) 
    character(len=72), optional, intent(out) :: mtitle

    character  :: rhstype,type*3,key*8
    character indfmt*16,ptrfmt*16,rhsfmt*20,valfmt*20
    integer        :: indcrd,  ptrcrd, totcrd,&
         & valcrd, rhscrd, nrow, ncol, nnzero, neltvl, nrhs, nrhsix
    integer                     :: ircode, i,nzr,infile,info
    character(len=*), parameter :: fmt10='(a72,a8,/,5i14,/,a3,11x,4i14,/,2a16,2a20)'
    character(len=*), parameter :: fmt11='(a1,13x,2i14)'
    character(len=*), parameter :: fmt111='(1x,a8,1x,i8,1x,a10)'

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

    read (infile,fmt=fmt10) mtitle,key,totcrd,ptrcrd,indcrd,valcrd,rhscrd,&
         & type,nrow,ncol,nnzero,neltvl,ptrfmt,indfmt,valfmt,rhsfmt
    if (rhscrd > 0) read(infile,fmt=fmt11)rhstype,nrhs,nrhsix


    if (psb_tolower(type(1:1)) == 'c') then 
      if (psb_tolower(type(2:2)) == 'u') then 

        call psb_sp_all(a,nnzero,nrow+1,nnzero,ircode)
        
        if (ircode /= 0 ) then 
          write(0,*) 'Memory allocation failed'
          goto 993
        end if
        
        a%m    = nrow
        a%k    = ncol
        a%fida = 'CSR'
        a%descra='G'

        read (infile,fmt=ptrfmt) (a%ia2(i),i=1,nrow+1)
        read (infile,fmt=indfmt) (a%ia1(i),i=1,nnzero)
        if (valcrd > 0) read (infile,fmt=valfmt) (a%aspk(i),i=1,nnzero)

        if (present(b)  .and. (rhscrd > 0)) then
          call psb_ensure_size(nrow,b,info)
          if (info == 0) read (infile,fmt=rhsfmt) (b(i),i=1,nrow)
        endif

      else if (psb_tolower(type(2:2)) == 's') then 

        ! we are generally working with non-symmetric matrices, so
        ! we de-symmetrize what we are about to read
        
        call psb_sp_all(nrow,ncol,a,nnzero,ircode)

        if (ircode /= 0 ) then 
          write(0,*) 'Memory allocation failed'
          goto 993
        end if
        a%fida = 'CSR'
        a%descra='G'
        

        read (infile,fmt=ptrfmt) (a%ia2(i),i=1,nrow+1)
        read (infile,fmt=indfmt) (a%ia1(i),i=1,nnzero)
        if (valcrd > 0) read (infile,fmt=valfmt) (a%aspk(i),i=1,nnzero)

        if (present(b)  .and. (rhscrd > 0)) then
          call psb_ensure_size(nrow,b,info)
          if (info == 0) read (infile,fmt=rhsfmt) (b(i),i=1,nrow)
        endif
        

        call psb_spcnv(a,ircode,afmt='coo')
        if (ircode /= 0) then
          write(0,*) 'ipcsr2coo ',ircode
          goto 993  
        end if

        call psb_sp_reall(a,2*nnzero,ircode)
      
        if (ircode /= 0 ) then 
          write(0,*) 'Memory allocation failed'
          goto 993
        end if

        ! A is now in COO format
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
        if (ircode /= 0) then
          write(0,*) 'ipcoo2csr ',ircode
          goto 993  
        end if

      else if (psb_tolower(type(2:2)) == 'h') then 

        ! we are generally working with non-symmetric matrices, so
        ! we de-symmetrize what we are about to read
        

        call psb_sp_all(nrow,ncol,a,2*nnzero,ircode)
        
        if (ircode /= 0 ) then 
          write(0,*) 'Memory allocation failed'
          goto 993
        end if
        a%fida = 'CSR'
        a%descra='G'
        

        read (infile,fmt=ptrfmt) (a%ia2(i),i=1,nrow+1)
        read (infile,fmt=indfmt) (a%ia1(i),i=1,nnzero)
        if (valcrd > 0) read (infile,fmt=valfmt) (a%aspk(i),i=1,nnzero)

        if (present(b)  .and. (rhscrd > 0)) then
          call psb_ensure_size(nrow,b,info)
          if (info == 0) read (infile,fmt=rhsfmt) (b(i),i=1,nrow)
        endif
        
        call psb_spcnv(a,ircode,afmt='coo')
        if (ircode /= 0) then
          write(0,*) 'ipcsr2coo ',ircode
          goto 993  
        end if
        call psb_sp_reall(a,2*nnzero,ircode)
      
        if (ircode /= 0 ) then 
          write(0,*) 'Memory allocation failed'
          goto 993
        end if

        ! A is now in COO format
        nzr = nnzero
        do i=1,nnzero
          if (a%ia1(i) /= a%ia2(i)) then 
            nzr = nzr + 1
            a%aspk(nzr) = conjg(a%aspk(i))
            a%ia1(nzr) = a%ia2(i)
            a%ia2(nzr) = a%ia1(i)
          end if
        end do
        a%infoa(psb_nnz_) = nzr
        call psb_spcnv(a,ircode,afmt='csr')
        if (ircode /= 0) then
          write(0,*) 'ipcoo2csr ',ircode
          goto 993  
        end if

      else
        write(0,*) 'read_matrix: matrix type not yet supported'
        iret=904
      end if
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
    write(0,*) 'ZHB_READ: Unexpected end of file '
    return
993 iret=993
    write(0,*) 'ZHB_READ: Memory allocation failure'
    return
  end subroutine zhb_read




  subroutine zhb_write(a,iret,iunit,filename,key,rhs,mtitle)
    use psb_base_mod
    implicit none
    type(psb_zspmat_type), intent(in)  :: a
    integer, intent(out)        :: iret
    character(len=*), optional, intent(in) :: mtitle
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    character(len=*), optional, intent(in) :: key
    complex(psb_dpk_), optional             :: rhs(:)
    integer                     :: iout

    character(len=*), parameter::  ptrfmt='(10I8)',indfmt='(10I8)'
    integer, parameter :: jptr=10,jind=10
    character(len=*), parameter::  valfmt='(4E20.12)',rhsfmt='(4E20.12)'
    integer, parameter :: jval=4,jrhs=4
    character(len=*), parameter :: fmt10='(a72,a8,/,5i14,/,a3,11x,4i14,/,2a16,2a20)'
    character(len=*), parameter :: fmt11='(a1,13x,2i14)'
    character(len=*), parameter :: fmt111='(1x,a8,1x,i8,1x,a10)'
    
    character(len=72) :: mtitle_
    character(len=8)  :: key_

    character  :: rhstype,type*3

    integer    :: i,indcrd,ptrcrd,rhscrd,totcrd,valcrd,&
         & nrow,ncol,nnzero, neltvl,nrhs,nrhsix

    iret = 0

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
    
    if (present(mtitle)) then 
      mtitle_ = mtitle
    else
      mtitle_ = 'Temporary PSBLAS title '
    endif
    if (present(key)) then 
      key_ = key
    else
      key_ = 'PSBMAT00'
    endif
    if (psb_toupper(a%fida) == 'CSR') then 
      
      nrow   = a%m
      ncol   = a%k
      nnzero = a%ia2(nrow+1)-1
      neltvl = 0 
      ptrcrd = (nrow+1)/jptr
      if (mod(nrow+1,jptr) > 0) ptrcrd = ptrcrd + 1
      indcrd = nnzero/jind
      if (mod(nnzero,jind) > 0) indcrd = indcrd + 1
      valcrd = nnzero/jval
      if (mod(nnzero,jval) > 0) valcrd = valcrd + 1
      if (present(rhs)) then 
        if (size(rhs)<nrow) then 
          rhscrd = 0
        else
          rhscrd = nrow/jrhs
          if (mod(nrow,jrhs) > 0) rhscrd = rhscrd + 1
        endif
        nrhs = 1
      else
        rhscrd = 0
        nrhs = 0
      end if
      totcrd = ptrcrd + indcrd + valcrd + rhscrd
      nrhsix = nrhs * nrow
      rhstype='F'
      type='CUA'
      
      write (iout,fmt=fmt10) mtitle_,key_,totcrd,ptrcrd,indcrd,valcrd,rhscrd,&
           &  type,nrow,ncol,nnzero,neltvl,ptrfmt,indfmt,valfmt,rhsfmt
      if (rhscrd > 0) write (iout,fmt=fmt11) rhstype,nrhs,nrhsix
      write (iout,fmt=ptrfmt) (a%ia2(i),i=1,nrow+1)
      write (iout,fmt=indfmt) (a%ia1(i),i=1,nnzero)
      if (valcrd > 0) write (iout,fmt=valfmt) (a%aspk(i),i=1,nnzero)
      if (rhscrd > 0) write (iout,fmt=rhsfmt) (rhs(i),i=1,nrow)


    else

      write(0,*) 'format: ',a%fida,' not yet implemented'

    endif

    if (iout /= 6) close(iout)


    return

901 continue 
    iret=901
    write(0,*) 'Error while opening ',filename
    return
  end subroutine zhb_write


end module psb_hbio_mod
