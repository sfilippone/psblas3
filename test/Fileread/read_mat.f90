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
    use typesp
    use mmio
    implicit none
    integer                               :: ictxt
    type(d_spmat)                         :: a
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

  subroutine zreadmat (filename, a, ictxt, inroot)
    use typesp
    use mmio
    implicit none
    integer                               :: ictxt
    type(z_spmat)                         :: a
    character                             :: filename*(*)
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
      if (info /= 0) then 
        write(0,*) 'Error return from MM_MAT_READ ',info
        call blacs_abort(ictxt, 1)   ! Unexpected End of File
      endif
    end if 
    return 

  end subroutine zreadmat


  SUBROUTINE READ_RHS (FILENAME, B, ICTXT, INROOT,&
       & INDWORK, INIWORK)   
    IMPLICIT NONE
    INTEGER                               :: ICTXT
    CHARACTER                             :: FILENAME*(*)
    INTEGER, OPTIONAL                     :: INROOT
    REAL(KIND(1.0D0)), POINTER, OPTIONAL  :: INDWORK(:)
    INTEGER, POINTER, OPTIONAL            :: INIWORK(:)   ! Local Variables
    INTEGER, PARAMETER   :: INFILE = 2
    INTEGER              :: NROW, NCOL, I,ROOT, NPROW, NPCOL, MYPROW, MYPCOL, IRCODE, J
    CHARACTER            :: MMHEADER*15, FMT*15, OBJECT*10, TYPE*10, SYM*15,&
         & LINE*1024
    REAL(KIND(1.0D0)), POINTER  :: B(:,:)
    IF (PRESENT(INROOT)) THEN
      ROOT = INROOT
    ELSE
      ROOT = 0
    END IF
    CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYPROW, MYPCOL)    
    IF (MYPROW == ROOT) THEN
      WRITE(*, *) 'Start read_rhs'      ! Open Input File
      OPEN(INFILE,FILE=FILENAME, STATUS='OLD', ERR=901, ACTION="READ")
      READ(INFILE,FMT=*, END=902) MMHEADER, OBJECT, FMT, TYPE, SYM
      write(0,*)'obj fmt',object, fmt
      IF ( (OBJECT .NE. 'matrix').OR.(FMT.NE.'array')) THEN
        WRITE(0,*) 'READ_RHS: input file type not yet supported'
        CALL BLACS_ABORT(ICTXT, 1)
      END IF

      do 
         read(infile,fmt='(a)') line
         if (line(1:1) /= '%')  exit
      end do

      READ(LINE,FMT=*)NROW,NCOL
      
      CALL LOWERC(TYPE,1,10)
      CALL LOWERC(SYM,1,15)
      IF ((TYPE == 'real').AND.(SYM == 'general')) THEN
        ALLOCATE(B(NROW,NCOL),STAT = IRCODE)
        IF (IRCODE /= 0)   GOTO 993
        READ(INFILE,FMT=*,END=902) ((B(I,J), I=1,NROW),J=1,NCOL)
           
      ELSE
        WRITE(0,*) 'READ_RHS: RHS type not yet supported'
        CALL BLACS_ABORT(ICTXT, 1)
      END IF      ! Read Right Hand Sides
      WRITE(*,*) 'End READ_RHS'
    END IF 
    RETURN 
    ! Open failed
901 WRITE(0,*) 'READ_RHS: Could not open file ',&
         & INFILE,' for input'
    CALL BLACS_ABORT(ICTXT, 1)   ! Unexpected End of File
902 WRITE(0,*) 'READ_RHS: Unexpected end of file ',INFILE,&
         & ' during input'
    CALL BLACS_ABORT(ICTXT, 1)   ! Allocation Failed
993 WRITE(0,*) 'READ_RHS: Memory allocation failure'
    CALL BLACS_ABORT(ICTXT, 1)
  END SUBROUTINE READ_RHS


  SUBROUTINE ZREAD_RHS(FILENAME, B, ICTXT, INROOT)
    IMPLICIT NONE
    INTEGER                               :: ICTXT
    CHARACTER                             :: FILENAME*(*)
    INTEGER, OPTIONAL                     :: INROOT
    INTEGER, PARAMETER          :: INFILE = 2
    INTEGER                     :: NROW, NCOL, I,ROOT, NPROW, NPCOL, MYPROW, MYPCOL, IRCODE, J
    CHARACTER                   :: MMHEADER*15, FMT*15, OBJECT*10, TYPE*10, SYM*15,&
         & LINE*1024
    COMPLEX(KIND(1.0D0)), POINTER  :: B(:,:)
    IF (PRESENT(INROOT)) THEN
      ROOT = INROOT
    ELSE
      ROOT = 0
    END IF
    CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYPROW, MYPCOL)    
    IF (MYPROW == ROOT) THEN
      WRITE(*, *) 'Start read_rhs'      ! Open Input File
      OPEN(INFILE,FILE=FILENAME, STATUS='OLD', ERR=901, ACTION="READ")
      READ(INFILE,FMT=*, END=902) MMHEADER, OBJECT, FMT, TYPE, SYM
      write(0,*)'obj fmt',object, fmt
      IF ( (OBJECT .NE. 'matrix').OR.(FMT.NE.'array')) THEN
        WRITE(0,*) 'READ_RHS: input file type not yet supported'
        CALL BLACS_ABORT(ICTXT, 1)
      END IF

      do 
         read(infile,fmt='(a)') line
         if (line(1:1) /= '%')  exit
      end do

      READ(LINE,FMT=*)NROW,NCOL

      CALL LOWERC(TYPE,1,10)
      CALL LOWERC(SYM,1,15)
      IF ((TYPE == 'complex').AND.(SYM == 'general')) THEN
        ALLOCATE(B(NROW,NCOL),STAT = IRCODE)
        IF (IRCODE /= 0)   GOTO 993
        READ(INFILE,FMT=*,END=902) ((B(I,J), I=1,NROW),J=1,NCOL)
           
      ELSE
        WRITE(0,*) 'READ_RHS: RHS type not yet supported'
        CALL BLACS_ABORT(ICTXT, 1)
      END IF      ! Read Right Hand Sides
      WRITE(*,*) 'End READ_RHS'
    END IF 
    RETURN 
    ! Open failed
901 WRITE(0,*) 'READ_RHS: Could not open file ',&
         & INFILE,' for input'
    CALL BLACS_ABORT(ICTXT, 1)   ! Unexpected End of File
902 WRITE(0,*) 'READ_RHS: Unexpected end of file ',INFILE,&
         & ' during input'
    CALL BLACS_ABORT(ICTXT, 1)   ! Allocation Failed
993 WRITE(0,*) 'READ_RHS: Memory allocation failure'
    CALL BLACS_ABORT(ICTXT, 1)
  END SUBROUTINE ZREAD_RHS



END MODULE READ_MAT
