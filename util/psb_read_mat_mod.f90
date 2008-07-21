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
!  Real(psb_dpk_), Pointer, Optional  :: B(:,:)
!     On Entry: fresh variable.
!     On Exit:  will contain right hand side(s).
!
!  Integer, Optional                    :: inroot
!     On Entry: Index of root processor (default: 0)
!     On Exit : unchanged.
!
module psb_read_mat_mod
  interface read_mat
    module procedure sreadmat, dreadmat, creadmat, zreadmat 
  end interface
  interface read_rhs
    module procedure sread_rhs, cread_rhs, dread_rhs, zread_rhs
  end interface


contains

  subroutine sreadmat (filename, a, ictxt, inroot)
    use psb_base_mod
    use psb_mmio_mod
    implicit none
    integer                               :: ictxt
    type(psb_sspmat_type)                 :: a
    character(len=*)                      :: filename
    integer, optional                     :: inroot
    integer, parameter          :: infile = 2
    integer                     :: info, root, np,  me

    if (present(inroot)) then
      root = inroot
    else
      root = psb_root_
    end if
    call psb_info(ictxt, me, np)    
    if (me == root) then
      write(*, '("Reading matrix...")')      ! open input file
      call mm_mat_read(a,info,infile,filename)
      if (info /= 0) then 
        write(0,*) 'Error return from MM_MAT_READ ',info
        call psb_abort(ictxt)   ! Unexpected End of File
      endif
    end if 
    return 

  end subroutine sreadmat


  subroutine sread_rhs (filename, b, info,ictxt, inroot)   
    use psb_base_mod
    use psb_mmio_mod
    implicit none
    integer                               :: ictxt
    character                             :: filename*(*)
    real(psb_spk_), allocatable, intent(out)  :: b(:,:)
    integer, intent(out) :: info
    integer, optional                     :: inroot
    integer, parameter   :: infile = 2
    integer              :: nrow, ncol, i,root, np,  me,  ircode, j
    
    if (present(inroot)) then
      root = inroot
    else
      root = psb_root_
    end if
    info = 0 
    call psb_info(ictxt, me, np)    
    if (me == root) then
      write(*, '("Reading rhs...")')      ! open input file
      call mm_vet_read(filename,b,info)
      if (info /= 0) then 
        write(0,*) 'read_rhs: something went wrong.'
        return 
      end if      ! read right hand sides
      write(*,*) 'end read_rhs'
    end if
    return 
    
  end subroutine sread_rhs


  subroutine dreadmat (filename, a, ictxt, inroot)
    use psb_base_mod
    use psb_mmio_mod
    implicit none
    integer                               :: ictxt
    type(psb_dspmat_type)                 :: a
    character(len=*)                      :: filename
    integer, optional                     :: inroot
    integer, parameter          :: infile = 2
    integer                     :: info, root, np,  me

    if (present(inroot)) then
      root = inroot
    else
      root = psb_root_
    end if
    call psb_info(ictxt, me, np)    
    if (me == root) then
      write(*, '("Reading matrix...")')      ! open input file
      call mm_mat_read(a,info,infile,filename)
      if (info /= 0) then 
        write(0,*) 'Error return from MM_MAT_READ ',info
        call psb_abort(ictxt)   ! Unexpected End of File
      endif
    end if 
    return 

  end subroutine dreadmat


  subroutine dread_rhs (filename, b, info,ictxt, inroot)   
    use psb_base_mod
    use psb_mmio_mod
    implicit none
    integer                               :: ictxt
    character                             :: filename*(*)
    real(psb_dpk_), allocatable, intent(out)  :: b(:,:)
    integer, intent(out) :: info
    integer, optional                     :: inroot
    integer, parameter   :: infile = 2
    integer              :: nrow, ncol, i,root, np,  me,  ircode, j
    
    if (present(inroot)) then
      root = inroot
    else
      root = psb_root_
    end if
    info = 0 
    call psb_info(ictxt, me, np)    
    if (me == root) then
      write(*, '("Reading rhs...")')      ! open input file
      call mm_vet_read(filename,b,info)
      if (info /= 0) then 
        write(0,*) 'read_rhs: something went wrong.'
        return 
      end if      ! read right hand sides
      write(*,*) 'end read_rhs'
    end if
    return 
    
  end subroutine dread_rhs

  subroutine creadmat (filename, a, ictxt, inroot)
    use psb_base_mod
    use psb_mmio_mod
    implicit none
    integer                               :: ictxt
    type(psb_cspmat_type)                 :: a
    character(len=*)                      :: filename
    integer, optional                     :: inroot
    integer, parameter          :: infile = 2
    integer                     :: info, root, np,  me

    if (present(inroot)) then
      root = inroot
    else
      root = psb_root_
    end if
    call psb_info(ictxt, me, np)
    if (me == root) then
      write(*, '("Reading matrix...")')      ! open input file
      call mm_mat_read(a,info,infile,filename)
      if (info /= 0) then 
        write(0,*) 'Error return from MM_MAT_READ ',info
        call psb_abort(ictxt)   ! Unexpected End of File
      endif
    end if 
    return 

  end subroutine creadmat


  subroutine cread_rhs (filename, b, info,ictxt, inroot)   
    use psb_base_mod
    use psb_mmio_mod
    implicit none
    integer                               :: ictxt
    character                             :: filename*(*)
    complex(psb_spk_), allocatable, intent(out)  :: b(:,:)
    integer, intent(out) :: info
    integer, optional                     :: inroot
    integer, parameter   :: infile = 2
    integer              :: nrow, ncol, i,root, np,  me,  ircode, j
    
    if (present(inroot)) then
      root = inroot
    else
      root = psb_root_
    end if
    info = 0 
    call psb_info(ictxt, me, np)    
    if (me == root) then
      write(*, '("Reading rhs...")')      ! open input file
      call mm_vet_read(filename,b,info)
      if (info /= 0) then 
        write(0,*) 'read_rhs: something went wrong.'
        return 
      end if      ! read right hand sides
      write(*,*) 'end read_rhs'
    end if
    return 
    
  end subroutine cread_rhs


  subroutine zreadmat (filename, a, ictxt, inroot)
    use psb_base_mod
    use psb_mmio_mod
    implicit none
    integer                               :: ictxt
    type(psb_zspmat_type)                 :: a
    character(len=*)                      :: filename
    integer, optional                     :: inroot
    integer, parameter          :: infile = 2
    integer                     :: info, root, np,  me

    if (present(inroot)) then
      root = inroot
    else
      root = psb_root_
    end if
    call psb_info(ictxt, me, np)
    if (me == root) then
      write(*, '("Reading matrix...")')      ! open input file
      call mm_mat_read(a,info,infile,filename)
      if (info /= 0) then 
        write(0,*) 'Error return from MM_MAT_READ ',info
        call psb_abort(ictxt)   ! Unexpected End of File
      endif
    end if 
    return 

  end subroutine zreadmat


  subroutine zread_rhs (filename, b, info,ictxt, inroot)   
    use psb_base_mod
    use psb_mmio_mod
    implicit none
    integer                               :: ictxt
    character                             :: filename*(*)
    complex(psb_dpk_), allocatable, intent(out)  :: b(:,:)
    integer, intent(out) :: info
    integer, optional                     :: inroot
    integer, parameter   :: infile = 2
    integer              :: nrow, ncol, i,root, np,  me,  ircode, j
    
    if (present(inroot)) then
      root = inroot
    else
      root = psb_root_
    end if
    info = 0 
    call psb_info(ictxt, me, np)    
    if (me == root) then
      write(*, '("Reading rhs...")')      ! open input file
      call mm_vet_read(filename,b,info)
      if (info /= 0) then 
        write(0,*) 'read_rhs: something went wrong.'
        return 
      end if      ! read right hand sides
      write(*,*) 'end read_rhs'
    end if
    return 
    
  end subroutine zread_rhs

end module psb_read_mat_mod
