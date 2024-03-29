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
subroutine zhb_read(a, iret, iunit, filename,b,g,x,mtitle)   
  use psb_base_mod
  implicit none
  type(psb_zspmat_type), intent(out)    :: a
  integer(psb_ipk_), intent(out)                   :: iret
  integer(psb_ipk_), optional, intent(in)          :: iunit
  character(len=*), optional, intent(in) :: filename
  complex(psb_dpk_), optional, allocatable, intent(out)  :: b(:,:), g(:,:), x(:,:) 
  character(len=72), optional, intent(out) :: mtitle

  character  :: rhstype*3,type*3,key*8
  character(len=72) :: mtitle_
  character indfmt*16,ptrfmt*16,rhsfmt*20,valfmt*20
  integer(psb_ipk_) :: indcrd,  ptrcrd, totcrd,&
       & valcrd, rhscrd, nrow, ncol, nnzero, neltvl, nrhs, nrhsix
  type(psb_z_csc_sparse_mat) :: acsc
  type(psb_z_coo_sparse_mat) :: acoo
  integer(psb_ipk_) :: ircode, i,nzr,infile, info
  character(len=*), parameter :: fmt10='(a72,a8,/,5i14,/,a3,11x,4i14,/,2a16,2a20)'
  character(len=*), parameter :: fmt11='(a3,11x,2i14)'
  character(len=*), parameter :: fmt111='(1x,a8,1x,i8,1x,a10)'

  iret = 0
  ircode = 0
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

  read (infile,fmt=fmt10) mtitle_,key,totcrd,ptrcrd,indcrd,valcrd,rhscrd,&
       & type,nrow,ncol,nnzero,neltvl,ptrfmt,indfmt,valfmt,rhsfmt
  if (rhscrd > 0) read(infile,fmt=fmt11)rhstype,nrhs,nrhsix

  call acsc%allocate_mnnz(nrow,ncol,nnzero)
  if (ircode /= 0 ) then 
    write(psb_err_unit,*) 'Memory allocation failed'
    goto 993
  end if

  if (present(mtitle)) mtitle=mtitle_


  if (psb_tolower(type(1:1)) == 'c') then 
    if (psb_tolower(type(2:2)) == 'u') then 


      read (infile,fmt=ptrfmt) (acsc%icp(i),i=1,ncol+1)
      read (infile,fmt=indfmt) (acsc%ia(i),i=1,nnzero)
      if (valcrd > 0) read (infile,fmt=valfmt) (acsc%val(i),i=1,nnzero)

      call a%mv_from(acsc)

      if (present(b)) then
        if ((psb_toupper(rhstype(1:1)) == 'F').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,b,info)
          read (infile,fmt=rhsfmt) (b(i,1),i=1,nrow)
        endif
      endif
      if (present(g)) then
        if ((psb_toupper(rhstype(2:2)) == 'G').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,g,info)
          read (infile,fmt=rhsfmt) (g(i,1),i=1,nrow)
        endif
      endif
      if (present(x)) then
        if ((psb_toupper(rhstype(3:3)) == 'X').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,x,info)
          read (infile,fmt=rhsfmt) (x(i,1),i=1,nrow)
        endif
      endif

    else if (psb_tolower(type(2:2)) == 's') then 

      ! we are generally working with non-symmetric matrices, so
      ! we de-symmetrize what we are about to read

      read (infile,fmt=ptrfmt) (acsc%icp(i),i=1,ncol+1)
      read (infile,fmt=indfmt) (acsc%ia(i),i=1,nnzero)
      if (valcrd > 0) read (infile,fmt=valfmt) (acsc%val(i),i=1,nnzero)


      if (present(b)) then
        if ((psb_toupper(rhstype(1:1)) == 'F').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,b,info)
          read (infile,fmt=rhsfmt) (b(i,1),i=1,nrow)
        endif
      endif
      if (present(g)) then
        if ((psb_toupper(rhstype(2:2)) == 'G').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,g,info)
          read (infile,fmt=rhsfmt) (g(i,1),i=1,nrow)
        endif
      endif
      if (present(x)) then
        if ((psb_toupper(rhstype(3:3)) == 'X').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,x,info)
          read (infile,fmt=rhsfmt) (x(i,1),i=1,nrow)
        endif
      endif


      call acoo%mv_from_fmt(acsc,info)
      call acoo%reallocate_nz(2*nnzero)
      ! A is now in COO format
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
      call acoo%fix(ircode)
      if (ircode == 0) call a%mv_from(acoo)
      if (ircode /= 0) goto 993

    else if (psb_tolower(type(2:2)) == 'h') then 

      ! we are generally working with non-symmetric matrices, so
      ! we de-symmetrize what we are about to read

      read (infile,fmt=ptrfmt) (acsc%icp(i),i=1,ncol+1)
      read (infile,fmt=indfmt) (acsc%ia(i),i=1,nnzero)
      if (valcrd > 0) read (infile,fmt=valfmt) (acsc%val(i),i=1,nnzero)


      if (present(b)) then
        if ((psb_toupper(rhstype(1:1)) == 'F').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,b,info)
          read (infile,fmt=rhsfmt) (b(i,1),i=1,nrow)
        endif
      endif
      if (present(g)) then
        if ((psb_toupper(rhstype(2:2)) == 'G').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,g,info)
          read (infile,fmt=rhsfmt) (g(i,1),i=1,nrow)
        endif
      endif
      if (present(x)) then
        if ((psb_toupper(rhstype(3:3)) == 'X').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,x,info)
          read (infile,fmt=rhsfmt) (x(i,1),i=1,nrow)
        endif
      endif


      call acoo%mv_from_fmt(acsc,info)
      call acoo%reallocate_nz(2*nnzero)
      ! A is now in COO format
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
      call acoo%fix(ircode)
      if (ircode == 0) call a%mv_from(acoo)
      if (ircode /= 0) goto 993

    else
      write(psb_err_unit,*) 'read_matrix: matrix type not yet supported'
      iret=904
    end if
  else
    write(psb_err_unit,*) 'read_matrix: matrix type not yet supported'
    iret=904
  end if

  call a%cscnv(ircode,type='csr')
  if (infile /= 5) close(infile)

  return 

  ! open failed
901 iret=901
  write(psb_err_unit,*) 'read_matrix: could not open file ',filename,' for input'
  return
902 iret=902
  write(psb_err_unit,*) 'HB_READ: Unexpected end of file '
  return
993 iret=993
  write(psb_err_unit,*) 'HB_READ: Memory allocation failure'
  return
end subroutine zhb_read

subroutine zhb_write(a,iret,iunit,filename,key,rhs,g,x,mtitle)
  use psb_base_mod
  implicit none
  type(psb_zspmat_type), intent(in), target :: a
  integer(psb_ipk_), intent(out)        :: iret
  character(len=*), optional, intent(in) :: mtitle
  integer(psb_ipk_), optional, intent(in)          :: iunit
  character(len=*), optional, intent(in) :: filename
  character(len=*), optional, intent(in) :: key
  complex(psb_dpk_), optional             :: rhs(:), g(:), x(:)
  integer(psb_ipk_) :: iout

  character(len=*), parameter::  ptrfmt='(10I8)',indfmt='(10I8)'
  integer(psb_ipk_), parameter :: jptr=10,jind=10
  character(len=*), parameter::  valfmt='(4E20.12)',rhsfmt='(4E20.12)'
  integer(psb_ipk_), parameter :: jval=2,jrhs=2
  character(len=*), parameter :: fmt10='(a72,a8,/,5i14,/,a3,11x,4i14,/,2a16,2a20)'
  character(len=*), parameter :: fmt11='(a3,11x,2i14)'
  character(len=*), parameter :: fmt111='(1x,a8,1x,i8,1x,a10)'

  type(psb_z_csc_sparse_mat)  ::  acsc
  character(len=72) :: mtitle_
  character(len=8)  :: key_

  character  :: rhstype*3,type*3

  integer(psb_ipk_) :: i,indcrd,ptrcrd,rhscrd,totcrd,valcrd,&
       & nrow,ncol,nnzero, neltvl, nrhs, nrhsix

  iret = 0

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


  call acsc%cp_from_fmt(a%a, iret)
  if (iret /= 0) return


  nrow   = acsc%get_nrows()
  ncol   = acsc%get_ncols()
  nnzero = acsc%get_nzeros()

  neltvl = 0 

  ptrcrd = (ncol+1)/jptr
  if (mod(ncol+1,jptr) > 0) ptrcrd = ptrcrd + 1
  indcrd = nnzero/jind
  if (mod(nnzero,jind) > 0) indcrd = indcrd + 1
  valcrd = nnzero/jval
  if (mod(nnzero,jval) > 0) valcrd = valcrd + 1
  rhstype = ''
  if (present(rhs)) then 
    if (size(rhs)<nrow) then 
      rhscrd = 0
    else
      rhscrd = nrow/jrhs
      if (mod(nrow,jrhs) > 0) rhscrd = rhscrd + 1
    endif
    nrhs = 1
    rhstype(1:1) = 'F'
  else
    rhscrd = 0
    nrhs   = 0
  end if
  totcrd = ptrcrd + indcrd + valcrd + rhscrd

  nrhsix  = nrhs*nrow

  if (present(g)) then 
    rhstype(2:2) = 'G'
  end if
  if (present(x)) then 
    rhstype(3:3) = 'X'
  end if
  type    = 'CUA'

  write (iout,fmt=fmt10) mtitle_,key_,totcrd,ptrcrd,indcrd,valcrd,rhscrd,&
       &  type,nrow,ncol,nnzero,neltvl,ptrfmt,indfmt,valfmt,rhsfmt
  if (rhscrd > 0) write (iout,fmt=fmt11) rhstype,nrhs,nrhsix
  write (iout,fmt=ptrfmt) (acsc%icp(i),i=1,ncol+1)
  write (iout,fmt=indfmt) (acsc%ia(i),i=1,nnzero)
  if (valcrd > 0) write (iout,fmt=valfmt) (acsc%val(i),i=1,nnzero)
  if (rhscrd > 0) write (iout,fmt=rhsfmt) (rhs(i),i=1,nrow)
  if (present(g).and.(rhscrd>0)) write (iout,fmt=rhsfmt) (g(i),i=1,nrow)
  if (present(x).and.(rhscrd>0)) write (iout,fmt=rhsfmt) (x(i),i=1,nrow)


  if (iout /= 6) close(iout)


  return

901 continue 
  iret=901
  write(psb_err_unit,*) 'Error while opening ',filename
  return
end subroutine zhb_write

subroutine lzhb_read(a, iret, iunit, filename,b,g,x,mtitle)   
  use psb_base_mod
  implicit none
  type(psb_lzspmat_type), intent(out)    :: a
  integer(psb_ipk_), intent(out)                   :: iret
  integer(psb_ipk_), optional, intent(in)          :: iunit
  character(len=*), optional, intent(in) :: filename
  complex(psb_dpk_), optional, allocatable, intent(out)  :: b(:,:), g(:,:), x(:,:) 
  character(len=72), optional, intent(out) :: mtitle

  character  :: rhstype*3,type*3,key*8
  character(len=72) :: mtitle_
  character indfmt*16,ptrfmt*16,rhsfmt*20,valfmt*20
  integer(psb_lpk_) :: indcrd,  ptrcrd, totcrd,&
       & valcrd, rhscrd, nrow, ncol, nnzero, neltvl, nrhs, nrhsix
  type(psb_lz_csc_sparse_mat) :: acsc
  type(psb_lz_coo_sparse_mat) :: acoo
  integer(psb_ipk_) :: ircode, infile, info
  integer(psb_lpk_) :: i,nzr
  character(len=*), parameter :: fmt10='(a72,a8,/,5i14,/,a3,11x,4i14,/,2a16,2a20)'
  character(len=*), parameter :: fmt11='(a3,11x,2i14)'
  character(len=*), parameter :: fmt111='(1x,a8,1x,i8,1x,a10)'

  iret = 0
  ircode = 0
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

  read (infile,fmt=fmt10) mtitle_,key,totcrd,ptrcrd,indcrd,valcrd,rhscrd,&
       & type,nrow,ncol,nnzero,neltvl,ptrfmt,indfmt,valfmt,rhsfmt
  if (rhscrd > 0) read(infile,fmt=fmt11)rhstype,nrhs,nrhsix

  call psb_lz_csc_allocate_mnnz(nrow,ncol,nnzero)
  if (ircode /= 0 ) then 
    write(psb_err_unit,*) 'Memory allocation failed'
    goto 993
  end if

  if (present(mtitle)) mtitle=mtitle_


  if (psb_tolower(type(1:1)) == 'c') then 
    if (psb_tolower(type(2:2)) == 'u') then 


      read (infile,fmt=ptrfmt) (acsc%icp(i),i=1,ncol+1)
      read (infile,fmt=indfmt) (acsc%ia(i),i=1,nnzero)
      if (valcrd > 0) read (infile,fmt=valfmt) (acsc%val(i),i=1,nnzero)

      call a%mv_from(acsc)

      if (present(b)) then
        if ((psb_toupper(rhstype(1:1)) == 'F').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,b,info)
          read (infile,fmt=rhsfmt) (b(i,1),i=1,nrow)
        endif
      endif
      if (present(g)) then
        if ((psb_toupper(rhstype(2:2)) == 'G').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,g,info)
          read (infile,fmt=rhsfmt) (g(i,1),i=1,nrow)
        endif
      endif
      if (present(x)) then
        if ((psb_toupper(rhstype(3:3)) == 'X').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,x,info)
          read (infile,fmt=rhsfmt) (x(i,1),i=1,nrow)
        endif
      endif

    else if (psb_tolower(type(2:2)) == 's') then 

      ! we are generally working with non-symmetric matrices, so
      ! we de-symmetrize what we are about to read

      read (infile,fmt=ptrfmt) (acsc%icp(i),i=1,ncol+1)
      read (infile,fmt=indfmt) (acsc%ia(i),i=1,nnzero)
      if (valcrd > 0) read (infile,fmt=valfmt) (acsc%val(i),i=1,nnzero)


      if (present(b)) then
        if ((psb_toupper(rhstype(1:1)) == 'F').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,b,info)
          read (infile,fmt=rhsfmt) (b(i,1),i=1,nrow)
        endif
      endif
      if (present(g)) then
        if ((psb_toupper(rhstype(2:2)) == 'G').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,g,info)
          read (infile,fmt=rhsfmt) (g(i,1),i=1,nrow)
        endif
      endif
      if (present(x)) then
        if ((psb_toupper(rhstype(3:3)) == 'X').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,x,info)
          read (infile,fmt=rhsfmt) (x(i,1),i=1,nrow)
        endif
      endif


      call acoo%mv_from_fmt(acsc,info)
      call acoo%reallocate_nz(2*nnzero)
      ! A is now in COO format
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
      call acoo%fix(ircode)
      if (ircode == 0) call a%mv_from(acoo)
      if (ircode /= 0) goto 993

    else if (psb_tolower(type(2:2)) == 'h') then 

      ! we are generally working with non-symmetric matrices, so
      ! we de-symmetrize what we are about to read

      read (infile,fmt=ptrfmt) (acsc%icp(i),i=1,ncol+1)
      read (infile,fmt=indfmt) (acsc%ia(i),i=1,nnzero)
      if (valcrd > 0) read (infile,fmt=valfmt) (acsc%val(i),i=1,nnzero)


      if (present(b)) then
        if ((psb_toupper(rhstype(1:1)) == 'F').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,b,info)
          read (infile,fmt=rhsfmt) (b(i,1),i=1,nrow)
        endif
      endif
      if (present(g)) then
        if ((psb_toupper(rhstype(2:2)) == 'G').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,g,info)
          read (infile,fmt=rhsfmt) (g(i,1),i=1,nrow)
        endif
      endif
      if (present(x)) then
        if ((psb_toupper(rhstype(3:3)) == 'X').and.(rhscrd > 0)) then 
          call psb_realloc(nrow,1,x,info)
          read (infile,fmt=rhsfmt) (x(i,1),i=1,nrow)
        endif
      endif


      call acoo%mv_from_fmt(acsc,info)
      call acoo%reallocate_nz(2*nnzero)
      ! A is now in COO format
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
      call acoo%fix(ircode)
      if (ircode == 0) call a%mv_from(acoo)
      if (ircode /= 0) goto 993

    else
      write(psb_err_unit,*) 'read_matrix: matrix type not yet supported'
      iret=904
    end if
  else
    write(psb_err_unit,*) 'read_matrix: matrix type not yet supported'
    iret=904
  end if

  call a%cscnv(ircode,type='csr')
  if (infile /= 5) close(infile)

  return 

  ! open failed
901 iret=901
  write(psb_err_unit,*) 'read_matrix: could not open file ',filename,' for input'
  return
902 iret=902
  write(psb_err_unit,*) 'HB_READ: Unexpected end of file '
  return
993 iret=993
  write(psb_err_unit,*) 'HB_READ: Memory allocation failure'
  return
end subroutine lzhb_read

subroutine lzhb_write(a,iret,iunit,filename,key,rhs,g,x,mtitle)
  use psb_base_mod
  implicit none
  type(psb_lzspmat_type), intent(in), target :: a
  integer(psb_ipk_), intent(out)        :: iret
  character(len=*), optional, intent(in) :: mtitle
  integer(psb_ipk_), optional, intent(in)          :: iunit
  character(len=*), optional, intent(in) :: filename
  character(len=*), optional, intent(in) :: key
  complex(psb_dpk_), optional             :: rhs(:), g(:), x(:)
  integer(psb_ipk_) :: iout

  character(len=*), parameter::  ptrfmt='(10I8)',indfmt='(10I8)'
  integer(psb_ipk_), parameter :: jptr=10,jind=10
  character(len=*), parameter::  valfmt='(4E20.12)',rhsfmt='(4E20.12)'
  integer(psb_ipk_), parameter :: jval=2,jrhs=2
  character(len=*), parameter :: fmt10='(a72,a8,/,5i14,/,a3,11x,4i14,/,2a16,2a20)'
  character(len=*), parameter :: fmt11='(a3,11x,2i14)'
  character(len=*), parameter :: fmt111='(1x,a8,1x,i8,1x,a10)'

  type(psb_lz_csc_sparse_mat)  ::  acsc
  character(len=72) :: mtitle_
  character(len=8)  :: key_

  character  :: rhstype*3,type*3

  integer(psb_ipk_) :: i,indcrd,ptrcrd,rhscrd,totcrd,valcrd,&
       & nrow,ncol,nnzero, neltvl, nrhs, nrhsix

  iret = 0

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


  call acsc%cp_from_fmt(a%a, iret)
  if (iret /= 0) return


  nrow   = acsc%get_nrows()
  ncol   = acsc%get_ncols()
  nnzero = acsc%get_nzeros()

  neltvl = 0 

  ptrcrd = (ncol+1)/jptr
  if (mod(ncol+1,jptr) > 0) ptrcrd = ptrcrd + 1
  indcrd = nnzero/jind
  if (mod(nnzero,jind) > 0) indcrd = indcrd + 1
  valcrd = nnzero/jval
  if (mod(nnzero,jval) > 0) valcrd = valcrd + 1
  rhstype = ''
  if (present(rhs)) then 
    if (size(rhs)<nrow) then 
      rhscrd = 0
    else
      rhscrd = nrow/jrhs
      if (mod(nrow,jrhs) > 0) rhscrd = rhscrd + 1
    endif
    nrhs = 1
    rhstype(1:1) = 'F'
  else
    rhscrd = 0
    nrhs   = 0
  end if
  totcrd = ptrcrd + indcrd + valcrd + rhscrd

  nrhsix  = nrhs*nrow

  if (present(g)) then 
    rhstype(2:2) = 'G'
  end if
  if (present(x)) then 
    rhstype(3:3) = 'X'
  end if
  type    = 'CUA'

  write (iout,fmt=fmt10) mtitle_,key_,totcrd,ptrcrd,indcrd,valcrd,rhscrd,&
       &  type,nrow,ncol,nnzero,neltvl,ptrfmt,indfmt,valfmt,rhsfmt
  if (rhscrd > 0) write (iout,fmt=fmt11) rhstype,nrhs,nrhsix
  write (iout,fmt=ptrfmt) (acsc%icp(i),i=1,ncol+1)
  write (iout,fmt=indfmt) (acsc%ia(i),i=1,nnzero)
  if (valcrd > 0) write (iout,fmt=valfmt) (acsc%val(i),i=1,nnzero)
  if (rhscrd > 0) write (iout,fmt=rhsfmt) (rhs(i),i=1,nrow)
  if (present(g).and.(rhscrd>0)) write (iout,fmt=rhsfmt) (g(i),i=1,nrow)
  if (present(x).and.(rhscrd>0)) write (iout,fmt=rhsfmt) (x(i),i=1,nrow)


  if (iout /= 6) close(iout)


  return

901 continue 
  iret=901
  write(psb_err_unit,*) 'Error while opening ',filename
  return
end subroutine lzhb_write
