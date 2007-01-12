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
!
!  READ_MAT subroutine reads a matrix and its right hand sides,
!  all stored in Matrix Market format file. The B field is optional,.
!
!  Character                            :: filename*20
!     On Entry: name of file to be processed.
!     On Exit : unchanged.
!
!  Type(D_SPMAT)                        :: A
!     On Entry: fresh variable.
!     On Exit : will contain the global sparse matrix as follows:
!        A%AS for coefficient values
!        A%IA1  for column indices
!        A%IA2  for row pointers
!        A%M    for number of global matrix rows
!        A%K    for number of global matrix columns
!
!  Integer                              :: ICTXT
!     On Entry: BLACS context.
!     On Exit : unchanged.
!
!  Real(Kind(1.D0)), Pointer, Optional  :: B(:,:)
!     On Entry: fresh variable.
!     On Exit:  will contain right hand side(s).
!
!  Integer, Optional                    :: inroot
!     On Entry: Index of root processor (default: 0)
!     On Exit : unchanged.
!
module psb_read_mat_mod
  interface read_mat
    module procedure dreadmat, zreadmat
  end interface
  interface read_rhs
    module procedure dread_rhs, zread_rhs
  end interface


contains

  subroutine dreadmat (filename, a, ictxt, inroot)
    use psb_base_mod
    use psb_mmio_mod
    implicit none
    integer                               :: ictxt
    type(psb_dspmat_type)                 :: a
    character(len=*)                      :: filename
    integer, optional                     :: inroot
    integer, parameter          :: infile = 2
    integer                     :: info, root, nprow, npcol, myprow, mypcol

    if (present(inroot)) then
      root = inroot
    else
      root = 0
    end if
    call psb_info(ictxt, myprow, nprow)    
    if (myprow == root) then
      write(*, '("Reading matrix...")')      ! open input file
      call mm_mat_read(a,info,infile,filename)
      if (info /= 0) then 
        write(0,*) 'Error return from MM_MAT_READ ',info
        call psb_abort(ictxt)   ! Unexpected End of File
      endif
    end if 
    return 

  end subroutine dreadmat


  subroutine dread_rhs (filename, b, ictxt, inroot)   
    use psb_base_mod
    implicit none
    integer                               :: ictxt
    character                             :: filename*(*)
    integer, optional                     :: inroot
    integer, parameter   :: infile = 2
    integer              :: nrow, ncol, i,root, nprow, npcol, myprow, mypcol, ircode, j
    character            :: mmheader*15, fmt*15, object*10, type*10, sym*15,&
         & line*1024
    real(kind(1.0d0)), allocatable  :: b(:,:)
    if (present(inroot)) then
      root = inroot
    else
      root = 0
    end if
    call psb_info(ictxt, myprow, nprow)    
    if (myprow == root) then
      write(*, '("Reading rhs...")')      ! open input file
      open(infile,file=filename, status='old', err=901, action="read")
      read(infile,fmt=*, end=902) mmheader, object, fmt, type, sym
      write(0,*)'obj fmt',object, fmt
      if ( (object .ne. 'matrix').or.(fmt.ne.'array')) then
        write(0,*) 'read_rhs: input file type not yet supported'
        call psb_abort(ictxt)
      end if

      do 
         read(infile,fmt='(a)') line
         if (line(1:1) /= '%')  exit
      end do

      read(line,fmt=*)nrow,ncol
      
      if ((tolower(type) == 'real').and.(tolower(sym) == 'general')) then
        allocate(b(nrow,ncol),stat = ircode)
        if (ircode /= 0)   goto 993
        read(infile,fmt=*,end=902) ((b(i,j), i=1,nrow),j=1,ncol)
           
      else
        write(0,*) 'read_rhs: rhs type not yet supported'
        call psb_abort(ictxt)
      end if      ! read right hand sides
      write(*,*) 'end read_rhs'
    end if 
    return 
    ! open failed
901 write(0,*) 'read_rhs: could not open file ',&
         & infile,' for input'
    call psb_abort(ictxt)   ! unexpected end of file
902 write(0,*) 'read_rhs: unexpected end of file ',infile,&
         & ' during input'
    call psb_abort(ictxt)   ! allocation failed
993 write(0,*) 'read_rhs: memory allocation failure'
    call psb_abort(ictxt)
  end subroutine dread_rhs


  subroutine zreadmat (filename, a, ictxt, inroot)
    use psb_base_mod
    use psb_mmio_mod
    implicit none
    integer                               :: ictxt
    type(psb_zspmat_type)                 :: a
    character(len=*)                      :: filename
    integer, optional                     :: inroot
    integer, parameter          :: infile = 2
    integer                     :: info, root, nprow, npcol, myprow, mypcol

    if (present(inroot)) then
      root = inroot
    else
      root = 0
    end if
    call psb_info(ictxt, myprow, nprow)
    if (myprow == root) then
      write(*, '("Reading matrix...")')      ! open input file
      call mm_mat_read(a,info,infile,filename)
      if (info /= 0) then 
        write(0,*) 'Error return from MM_MAT_READ ',info
        call psb_abort(ictxt)   ! Unexpected End of File
      endif
    end if 
    return 

  end subroutine zreadmat


  subroutine zread_rhs (filename, b, ictxt, inroot)   
    use psb_base_mod
    implicit none
    integer                               :: ictxt
    character                             :: filename*(*)
    integer, optional                     :: inroot
    integer, parameter   :: infile = 2
    integer              :: nrow, ncol, i,root, nprow, npcol, myprow, mypcol, ircode, j
    character            :: mmheader*15, fmt*15, object*10, type*10, sym*15,&
         & line*1024
    real(kind(1.d0))     :: bre, bim 
    complex(kind(1.0d0)), allocatable  :: b(:,:)
    if (present(inroot)) then
      root = inroot
    else
      root = 0
    end if
    call psb_info(ictxt, myprow, nprow)
    if (myprow == root) then
      write(*, '("Reading rhs...")')      ! open input file
      open(infile,file=filename, status='old', err=901, action="read")
      read(infile,fmt=*, end=902) mmheader, object, fmt, type, sym
!!$      write(0,*)'obj fmt',object, fmt
      if ( (object .ne. 'matrix').or.(fmt.ne.'array')) then
        write(0,*) 'read_rhs: input file type not yet supported'
        call psb_abort(ictxt)
      end if

      do 
         read(infile,fmt='(a)') line
         if (line(1:1) /= '%')  exit
      end do

      read(line,fmt=*)nrow,ncol
      
      if ((tolower(type) == 'complex').and.(tolower(sym) == 'general')) then
        allocate(b(nrow,ncol),stat = ircode)
        if (ircode /= 0)   goto 993
        do j=1, ncol
          do i=1, nrow
            read(infile,fmt=*,end=902) bre,bim
            b(i,j) = cmplx(bre,bim)
          end do
        end do
      else
        write(0,*) 'read_rhs: rhs type not yet supported'
        call psb_abort(ictxt)
      end if      ! read right hand sides
      write(*,*) 'end read_rhs'
    end if 
    return 
    ! open failed
901 write(0,*) 'read_rhs: could not open file ',&
         & infile,' for input'
    call psb_abort(ictxt)   ! unexpected end of file
902 write(0,*) 'read_rhs: unexpected end of file ',infile,&
         & ' during input'
    call psb_abort(ictxt)   ! allocation failed
993 write(0,*) 'read_rhs: memory allocation failure'
    call psb_abort(ictxt)
  end subroutine zread_rhs



end module psb_read_mat_mod
