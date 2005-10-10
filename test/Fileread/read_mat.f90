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
!  Real(Kind(1.D0)), Pointer, Optional  :: indwork(:)
!     On Entry/Exit: Double Precision Work Area.
!
!  Integer, Pointer, Optional           :: iniwork()
!     On Entry/Exit: Integer Work Area.
!
module read_mat
  public readmat
  public read_rhs
  public zreadmat
  public zread_rhs


contains

  subroutine readmat (filename, a, ictxt, inroot)
    use psb_spmat_type
    use mmio
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
    call blacs_gridinfo(ictxt, nprow, npcol, myprow, mypcol)    
    if (myprow == root) then
      write(*, *) 'start read_matrix'      ! open input file
      call mm_mat_read(a,info,infile,filename)
      write(*, *) 'end read_matrix'
      if (info /= 0) then 
        write(0,*) 'Error return from MM_MAT_READ ',info
        call blacs_abort(ictxt, 1)   ! Unexpected End of File
      endif
    end if 
    return 

  end subroutine readmat


  subroutine read_rhs (filename, b, ictxt, inroot,&
       & indwork, iniwork)   
    implicit none
    integer                               :: ictxt
    character                             :: filename*(*)
    integer, optional                     :: inroot
    real(kind(1.0d0)), pointer, optional  :: indwork(:)
    integer, pointer, optional            :: iniwork(:)   ! local variables
    integer, parameter   :: infile = 2
    integer              :: nrow, ncol, i,root, nprow, npcol, myprow, mypcol, ircode, j
    character            :: mmheader*15, fmt*15, object*10, type*10, sym*15,&
         & line*1024
    real(kind(1.0d0)), pointer  :: b(:,:)
    if (present(inroot)) then
      root = inroot
    else
      root = 0
    end if
    call blacs_gridinfo(ictxt, nprow, npcol, myprow, mypcol)    
    if (myprow == root) then
      write(*, *) 'start read_rhs'      ! open input file
      open(infile,file=filename, status='old', err=901, action="read")
      read(infile,fmt=*, end=902) mmheader, object, fmt, type, sym
      write(0,*)'obj fmt',object, fmt
      if ( (object .ne. 'matrix').or.(fmt.ne.'array')) then
        write(0,*) 'read_rhs: input file type not yet supported'
        call blacs_abort(ictxt, 1)
      end if

      do 
         read(infile,fmt='(a)') line
         if (line(1:1) /= '%')  exit
      end do

      read(line,fmt=*)nrow,ncol
      
      call lowerc(type,1,10)
      call lowerc(sym,1,15)
      if ((type == 'real').and.(sym == 'general')) then
        allocate(b(nrow,ncol),stat = ircode)
        if (ircode /= 0)   goto 993
        read(infile,fmt=*,end=902) ((b(i,j), i=1,nrow),j=1,ncol)
           
      else
        write(0,*) 'read_rhs: rhs type not yet supported'
        call blacs_abort(ictxt, 1)
      end if      ! read right hand sides
      write(*,*) 'end read_rhs'
    end if 
    return 
    ! open failed
901 write(0,*) 'read_rhs: could not open file ',&
         & infile,' for input'
    call blacs_abort(ictxt, 1)   ! unexpected end of file
902 write(0,*) 'read_rhs: unexpected end of file ',infile,&
         & ' during input'
    call blacs_abort(ictxt, 1)   ! allocation failed
993 write(0,*) 'read_rhs: memory allocation failure'
    call blacs_abort(ictxt, 1)
  end subroutine read_rhs



end module read_mat
