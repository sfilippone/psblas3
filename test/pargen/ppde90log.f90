!
! This sample program shows how to build and solve a sparse linear
!
! The program  solves a linear system based on the partial differential
! equation 
!
! 
!
!     the equation generated is:
!   b1 d d (u)  b2  d d (u)    a1 d (u))  a2 d (u)))   
! -   ------   -    ------  +    ----- +  ------     + a3 u = 0
!      dx dx         dy dy         dx        dy        
!
! 
! with  Dirichlet boundary conditions on the unit cube 
!
!    0<=x,y,z<=1
! 
! The equation is discretized with finite differences and uniform stepsize;
! the resulting  discrete  equation is
!
! ( u(x,y,z)(2b1+2b2+a1+a2)+u(x-1,y)(-b1-a1)+u(x,y-1)(-b2-a2)+
!  -u(x+1,y)b1-u(x,y+1)b2)*(1/h**2)
!
! Example taken from: C.T.Kelley
!    Iterative Methods for Linear and Nonlinear Equations
!    SIAM 1995
!
!
! In this sample program the index space of the discretized
! computational domain is first numbered sequentially in a standard way, 
! then the corresponding vector is distributed according to an HPF BLOCK
! distribution directive.
!
! Boundary conditions are set in a very simple way, by adding 
! equations of the form
!
!   u(x,y) = rhs(x,y)
!
Program PDE90
  USE TYPESP
  USE TYPEDESC
  USE TYPEPREC
  USE F90TOOLS
  USE F90METHD
  USE F90PREC
  use mpi
  Implicit none

  interface 
    !.....user passed subroutine.....
    subroutine part_block(glob_index,n,np,pv,nv)
      INTEGER, INTENT(IN)  :: GLOB_INDEX, N, NP
      INTEGER, INTENT(OUT) :: NV
      INTEGER, INTENT(OUT) :: PV(*) 
    end subroutine part_block
  end interface
  ! input parameters
  Character :: CMETHD*10, PREC*10, AFMT*5
  Integer      :: IDIM, IRET

  ! Miscellaneous 
  Integer, Parameter   :: IZERO=0, IONE=1
  Character, PARAMETER :: ORDER='R'
  INTEGER              :: IARGC,CONVERT_DESCR,dim, CHECK_DESCR
  REAL(KIND(1.D0)), PARAMETER :: DZERO = 0.D0, ONE = 1.D0
  REAL(KIND(1.D0)) :: T1, T2, TPREC, TSOLVE, T3, T4 
!!$  EXTERNAL  MPI_WTIME
  integer mpe_log_get_event_number,mpe_Describe_state,mpe_log_event

  ! Sparse Matrix and preconditioner
  TYPE(D_SPMAT) :: A,  L, U, H
  TYPE(PREC_DATA) :: PRE
  ! Descriptor
  TYPE(desc_type)    :: DESC_A, DESC_A_OUT
  ! Dense Matrices
  REAL(KIND(1.d0)), POINTER :: B(:), X(:), D(:), LD(:)
  INTEGER, pointer :: WORK(:)
  ! BLACS parameters
  INTEGER            :: nprow, npcol, icontxt, iam, np, myprow, mypcol

  ! Solver parameters
  INTEGER            :: ITER, ITMAX,IERR,ITRACE, METHD,IPREC, ISTOPC,&
       & IPARM(20), ML
  REAL(KIND(1.D0))   :: ERR, EPS, RPARM(20)

  ! Other variables
  INTEGER            :: I,INFO,iprecb,iprece,islvb,islve
  INTEGER            :: INTERNAL, M,II

  ! Initialize BLACS  
  CALL BLACS_PINFO(IAM, NP)
  CALL BLACS_GET(IZERO, IZERO, ICONTXT)
  iprecb = mpe_log_get_event_number()
  iprece = mpe_log_get_event_number()
  islvb  = mpe_log_get_event_number()
  islve  = mpe_log_get_event_number()
  if (iam==0) then 
    info = mpe_describe_state(iprecb,iprece,"Preconditioner","OrangeRed")
    info = mpe_describe_state(islvb,islve,"Solver","DarkGreen")
  endif


  ! Rectangular Grid,  P x 1

  CALL BLACS_GRIDINIT(ICONTXT, ORDER, NP, IONE)
  CALL BLACS_GRIDINFO(ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL)


  !
  !  Get parameters
  !
  CALL GET_PARMS(ICONTXT,CMETHD,PREC,AFMT,IDIM,ISTOPC,ITMAX,ITRACE,ML)
  
  !
  !  Allocate and fill in the coefficient matrix, RHS and initial guess 
  !

  CALL BLACS_BARRIER(ICONTXT,'All')
  T1 = MPI_WTIME()
  CALL CREATE_MATRIX(IDIM,A,B,X,DESC_A,PART_BLOCK,ICONTXT,AFMT)  
  T2 = MPI_WTIME() - T1

  DIM=SIZE(A%ASPK)

  ALLOCATE(H%ASPK(DIM),H%IA1(DIM),H%IA2(DIM),H%PL(SIZE(A%PL)),&
       & H%PL(SIZE(A%PL)),D(SIZE(A%PL)),&
       & DESC_A_OUT%MATRIX_DATA(SIZE(DESC_A%MATRIX_DATA)),&
       & DESC_A_OUT%HALO_INDEX(SIZE(DESC_A%HALO_INDEX)),&
       & DESC_A_OUT%OVRLAP_INDEX(SIZE(DESC_A%OVRLAP_INDEX)),&
       & DESC_A_OUT%OVRLAP_ELEM(SIZE(DESC_A%OVRLAP_ELEM)),&
       & DESC_A_OUT%LOC_TO_GLOB(SIZE(DESC_A%LOC_TO_GLOB)),&
       & DESC_A_OUT%GLOB_TO_LOC(SIZE(DESC_A%GLOB_TO_LOC)), WORK(1024))
  check_descr=15
!  work(5)=9
!!$  WRITE(0,*)'CALLING VERIFY'
!!$  CALL F90_PSVERIFY(D,A,DESC_A,CHECK_DESCR,CONVERT_DESCR,H,&
!!$       & DESC_A_OUT,WORK)
!!$  WRITE(0,*)'VERIFY DONE',CONVERT_DESCR

  deallocate(work)

  CALL DGAMX2D(ICONTXT,'A',' ',IONE, IONE,T2,IONE,T1,T1,-1,-1,-1)
  IF (IAM.EQ.0) Write(6,*) 'Matrix creation Time : ',T2

  !
  !  Prepare the preconditioner.
  !  
  write(0,*)'PRECONDIZIONATORE=',prec
  SELECT CASE (PREC)     
  CASE ('SCHW6')
     IPREC = 6
  CASE ('SCHW5')
     IPREC = 5
  CASE ('SCHW4')
     IPREC = 4
  CASE ('SCHW3')
     IPREC = 3
  CASE ('ILU')
     IPREC = 2
  CASE ('DIAGSC')
     IPREC = 1
  CASE ('NONE')
     IPREC = 0
  CASE DEFAULT
     WRITE(0,*) 'Unknown preconditioner'
     CALL BLACS_ABORT(ICONTXT,-1)
  END SELECT
  pre%prec=iprec
  CALL BLACS_BARRIER(ICONTXT,'All')
  T1 = MPI_WTIME()
  info = MPE_Log_event( iprecb, 0, "start Precond" )

  CALL PRECONDITIONER(A,PRE,DESC_A,IRET)
!!$  CALL PRECONDITIONER(IPREC,A,L,U,D,DESC_A,IRET)
  info = MPE_Log_event( iprece, 0, "end Precond" )
  TPREC = MPI_WTIME()-T1
  
  CALL DGAMX2D(icontxt,'A',' ',IONE, IONE,TPREC,IONE,t1,t1,-1,-1,-1)
  
  IF (IAM.EQ.0) WRITE(6,*) 'Preconditioner Time : ',TPREC

  IF (IRET.NE.0) THEN
     WRITE(0,*) 'Error on preconditioner',IRET
     CALL BLACS_ABORT(ICONTXT,-1)
     STOP
  END IF
 
  !
  ! Iterative method parameters 
  !
  write(*,*) 'Calling Iterative method', size(b),ml
  CALL BLACS_BARRIER(ICONTXT,'All')
  T1 = MPI_WTIME()  
  EPS   = 1.D-9
  info = MPE_Log_event( islvb, 0, "start Solver" )
  IF (CMETHD.EQ.'BICGSTAB') THEN 
    CALL  F90_BICGSTAB(A,PRE,B,X,EPS,DESC_A,& 
         & ITMAX,ITER,ERR,IERR,ITRACE)     
!!$  ELSE IF (CMETHD.EQ.'BICG') THEN 
!!$    CALL  F90_BICG(A,PRE,B,X,EPS,DESC_A,& 
!!$         & ITMAX,ITER,ERR,IERR,ITRACE)     
!!$  ELSE  IF (CMETHD.EQ.'CGS') THEN 
!!$    CALL  F90_CGS(A,PRE,B,X,EPS,DESC_A,& 
!!$         & ITMAX,ITER,ERR,IERR,ITRACE)     
!!$  ELSE IF (CMETHD.EQ.'BICGSTABL') THEN 
!!$    CALL  F90_BICGSTABL(A,PRE,B,X,EPS,DESC_A,& 
!!$         & ITMAX,ITER,ERR,IERR,ITRACE,ML)     
  ELSE
    write(0,*) 'Unknown method ',cmethd
  end IF
  info = MPE_Log_event( islve, 0, "end Solver" )  
  CALL BLACS_BARRIER(ICONTXT,'All')
  T2 = MPI_WTIME() - T1
  CALL DGAMX2D(ICONTXT,'A',' ',IONE, IONE,T2,IONE,T1,T1,-1,-1,-1)

  IF (IAM.EQ.0) THEN
     WRITE(6,*) 'Time to Solve Matrix : ',T2
     WRITE(6,*) 'Time per iteration : ',T2/ITER
     WRITE(6,*) 'Number of iterations : ',ITER
     WRITE(6,*) 'Error on exit : ',ERR
     WRITE(6,*) 'INFO  on exit : ',IERR
  END IF

  !  
  !  Cleanup storage and exit
  !
  CALL F90_PSDSFREE(B,DESC_A)
  CALL F90_PSDSFREE(X,DESC_A)
!!$  CALL F90_PSDSFREE(D,DESC_A)

  CALL F90_PSSPFREE(A,DESC_A)
!!$  CALL F90_PSSPFREE(L,DESC_A) 
!!$  CALL F90_PSSPFREE(U,DESC_A)          
  CALL F90_PSDSCFREE(DESC_A,info)
  
  CALL BLACS_GRIDEXIT(ICONTXT)
  CALL BLACS_EXIT(0)
  
  STOP
  
CONTAINS
  !
  ! Get iteration parameters from the command line
  !
  SUBROUTINE  GET_PARMS(ICONTXT,CMETHD,PREC,AFMT,IDIM,ISTOPC,ITMAX,ITRACE,ML)
    integer      :: icontxt
    Character    :: CMETHD*10, PREC*10, AFMT*5
    Integer      :: IDIM, IRET, ISTOPC,ITMAX,ITRACE,ML
    Character*40 :: CHARBUF
    INTEGER      :: IARGC, NPROW, NPCOL, MYPROW, MYPCOL
    EXTERNAL     IARGC
    INTEGER      :: INTBUF(10), IP
    
    CALL BLACS_GRIDINFO(ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL)

    IF (MYPROW==0) THEN
       READ(*,*) IP
       IF (IP.GE.3) THEN
          READ(*,*) CMETHD
          READ(*,*) PREC
          READ(*,*) AFMT
         
        ! Convert strings in array
          DO I = 1, LEN(CMETHD)
             INTBUF(I) = IACHAR(CMETHD(I:I))
          END DO
        ! Broadcast parameters to all processors
          CALL IGEBS2D(ICONTXT,'ALL',' ',10,1,INTBUF,10)
        
          DO I = 1, LEN(PREC)
             INTBUF(I) = IACHAR(PREC(I:I))
          END DO
        ! Broadcast parameters to all processors
          CALL IGEBS2D(ICONTXT,'ALL',' ',10,1,INTBUF,10)
        
          DO I = 1, LEN(AFMT)
             INTBUF(I) = IACHAR(AFMT(I:I))
          END DO
        ! Broadcast parameters to all processors
          CALL IGEBS2D(ICONTXT,'ALL',' ',10,1,INTBUF,10)
        
          READ(*,*) IDIM
          IF (IP.GE.4) THEN
             READ(*,*) ISTOPC
          ELSE
             ISTOPC=1        
          ENDIF
          IF (IP.GE.5) THEN
             READ(*,*) ITMAX
          ELSE
             ITMAX=500
          ENDIF
          IF (IP.GE.6) THEN
             READ(*,*) ITRACE
          ELSE
             ITRACE=-1
          ENDIF
          IF (IP.GE.7) THEN
             READ(*,*) ML
          ELSE
             ML=1
          ENDIF
        ! Broadcast parameters to all processors    
          
          INTBUF(1) = IDIM
          INTBUF(2) = ISTOPC
          INTBUF(3) = ITMAX
          INTBUF(4) = ITRACE
          INTBUF(5) = ML
          CALL IGEBS2D(ICONTXT,'ALL',' ',5,1,INTBUF,5)

          WRITE(6,*)'Solving matrix: ELL1'      
          WRITE(6,*)'on  grid',IDIM,'x',IDIM,'x',IDIM
          WRITE(6,*)' with BLOCK data distribution, NP=',Np,&
               & ' Preconditioner=',PREC,&
               & ' Iterative methd=',CMETHD      
       ELSE
        ! Wrong number of parameter, print an error message and exit
          CALL PR_USAGE(0)      
          CALL BLACS_ABORT(ICONTXT,-1)
          STOP 1
       ENDIF
    ELSE
   ! Receive Parameters
       CALL IGEBR2D(ICONTXT,'ALL',' ',10,1,INTBUF,10,0,0)
       DO I = 1, 10
          CMETHD(I:I) = ACHAR(INTBUF(I))
       END DO
       CALL IGEBR2D(ICONTXT,'ALL',' ',10,1,INTBUF,10,0,0)
       DO I = 1, 10
          PREC(I:I) = ACHAR(INTBUF(I))
       END DO
       CALL IGEBR2D(ICONTXT,'ALL',' ',10,1,INTBUF,10,0,0)
       DO I = 1, 5
          AFMT(I:I) = ACHAR(INTBUF(I))
       END DO
       CALL IGEBR2D(ICONTXT,'ALL',' ',5,1,INTBUF,5,0,0)
       IDIM    = INTBUF(1)
       ISTOPC  = INTBUF(2)
       ITMAX   = INTBUF(3)
       ITRACE  = INTBUF(4)
       ML      = INTBUF(5)
    END IF
    RETURN
    
  END SUBROUTINE GET_PARMS
  !
  !  Print an error message 
  !  
  SUBROUTINE PR_USAGE(IOUT)
    INTEGER :: IOUT
    WRITE(IOUT,*)'Incorrect parameter(s) found'
    WRITE(IOUT,*)' Usage:  pde90 methd prec dim &
         &[istop itmax itrace]'  
    WRITE(IOUT,*)' Where:'
    WRITE(IOUT,*)'     methd:    CGSTAB TFQMR CGS' 
    WRITE(IOUT,*)'     prec :    ILU DIAGSC NONE'
    WRITE(IOUT,*)'     dim       number of points along each axis'
    WRITE(IOUT,*)'               the size of the resulting linear '
    WRITE(IOUT,*)'               system is dim**3'
    WRITE(IOUT,*)'     istop     Stopping criterion  1, 2 or 3 [1]  '
    WRITE(IOUT,*)'     itmax     Maximum number of iterations [500] '
    WRITE(IOUT,*)'     itrace    0  (no tracing, default) or '  
    WRITE(IOUT,*)'               >= 0 do tracing every ITRACE'
    WRITE(IOUT,*)'               iterations ' 
  END SUBROUTINE PR_USAGE

!
!  Subroutine to allocate and fill in the coefficient matrix and
!  the RHS. 
!
  SUBROUTINE CREATE_MATRIX(IDIM,A,B,T,DESC_A,PARTS,ICONTXT,AFMT)
    !
    !   Discretize the partial diferential equation
    ! 
    !   b1 dd(u)  b2 dd(u)    b3 dd(u)    a1 d(u)   a2 d(u)  a3 d(u)  
    ! -   ------ -  ------ -  ------ -  -----  -  ------  -  ------ + a4 u 
    !      dxdx     dydy       dzdz        dx       dy         dz   
    !
    !  = 0 
    ! 
    ! boundary condition: Dirichlet
    !    0< x,y,z<1
    !  
    !  u(x,y,z)(2b1+2b2+2b3+a1+a2+a3)+u(x-1,y,z)(-b1-a1)+u(x,y-1,z)(-b2-a2)+
    !  + u(x,y,z-1)(-b3-a3)-u(x+1,y,z)b1-u(x,y+1,z)b2-u(x,y,z+1)b3

    USE TYPESP
    USE TYPEDESC
    USE F90TOOLS
    USE F90METHD
    Implicit None
    INTEGER                  :: IDIM
    integer, parameter       :: nbmax=10
    Real(Kind(1.D0)),Pointer :: B(:),T(:)
    Type (desc_type)  :: DESC_A
    Integer                  :: ICONTXT
    INTERFACE 
      !   .....user passed subroutine.....
      SUBROUTINE PARTS(GLOBAL_INDX,N,NP,PV,NV)
        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: GLOBAL_INDX, N, NP
        INTEGER, INTENT(OUT) :: NV
        INTEGER, INTENT(OUT) :: PV(*) 
      END SUBROUTINE PARTS
    END INTERFACE   ! Local variables
    Type(D_SPMAT)            :: A
    Real(Kind(1.d0))         :: ZT(NBMAX),GLOB_X,GLOB_Y,GLOB_Z
    Integer                  :: M,N,NNZ,GLOB_ROW,J
    Type (D_SPMAT)           :: ROW_MAT
    Integer                  :: X,Y,Z,COUNTER,IA,I,INDX_OWNER
    INTEGER                  :: NPROW,NPCOL,MYPROW,MYPCOL
    Integer                  :: ELEMENT
    INTEGER                  :: INFO, NV, INV
    INTEGER, ALLOCATABLE     :: PRV(:)
    INTEGER, pointer         :: ierrv(:)
    Real(Kind(1.d0)), pointer ::  DWORK(:)
    INTEGER,POINTER        ::  IWORK(:)
    character              :: afmt*5
    ! deltah dimension of each grid cell
    ! deltat discretization time
    Real(Kind(1.D0))         :: DELTAH
    Real(Kind(1.d0)),Parameter   :: RHS=0.d0,ONE=1.d0,ZERO=0.d0
    Real(Kind(1.d0))   :: MPI_WTIME, T1, T2, T3, TINS
    Real(Kind(1.d0))   :: a1, a2, a3, a4, b1, b2, b3 
    external            mpi_wtime,a1, a2, a3, a4, b1, b2, b3
    integer            :: nb, ir1, ir2, ipr
    logical            :: own
    ! common area


    CALL BLACS_GRIDINFO(ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL)

    DELTAH = 1.D0/(IDIM-1)

    ! Initialize array descriptor and sparse matrix storage. Provide an
    ! estimate of the number of non zeroes 
    CALL SETERR(2)
    allocate(ierrv(6))

    ierrv(:) = 0
    M   = IDIM*IDIM*IDIM
    N   = M
    NNZ = ((N*9)/(NPROW*NPCOL))
    write(*,*) 'Size: n ',n
    Call F90_PSDSCALL(N,N,PARTS,ICONTXT,IERRV,DESC_A)
    write(*,*) 'Allocating A : nnz',nnz
    Call F90_PSSPALL(A,IERRV,DESC_A,NNZ=NNZ)
    ! Define  RHS from boundary conditions; also build initial guess 
    write(*,*) 'Allocating B'
    Call F90_PSDSALL(N,B,IERRV,DESC_A)
    write(*,*) 'Allocating T'
    Call F90_PSDSALL(N,T,IERRV,DESC_A)

    ! We build an auxiliary matrix consisting of one row at a
    ! time; just a small matrix. Might be extended to generate 
    ! a bunch of rows per call. 
    ! 
    ROW_MAT%DESCRA(1:1) = 'G'
    ROW_MAT%FIDA        = 'CSR'
    write(*,*) 'Allocating ROW_MAT',20*nbmax
    ALLOCATE(ROW_MAT%ASPK(20*nbmax),ROW_MAT%IA1(20*nbmax),&
         &ROW_MAT%IA2(20*nbmax),PRV(NPROW),stat=info)
    if (info.ne.0 ) then 
      write(*,*) 'Memory allocation error'
      call blacs_abort(icontxt,-1)
    endif

    TINS = 0.D0
    CALL BLACS_BARRIER(ICONTXT,'ALL')
    T1 = MPI_WTIME()

    ! Loop over rows belonging to current process in a BLOCK
    ! distribution.

    ROW_MAT%IA2(1)=1    
    DO GLOB_ROW = 1, N
      CALL PARTS(GLOB_ROW,N,NPROW,PRV,NV)
      DO INV = 1, NV
        INDX_OWNER = PRV(INV)
        IF (INDX_OWNER == MYPROW) THEN
          ! Local matrix pointer 
          ELEMENT=1
          ! Compute gridpoint Coordinates
          IF (MOD(GLOB_ROW,(IDIM*IDIM)).EQ.0) THEN
            X = GLOB_ROW/(IDIM*IDIM)
          ELSE
            X = GLOB_ROW/(IDIM*IDIM)+1
          ENDIF
          IF (MOD((GLOB_ROW-(X-1)*IDIM*IDIM),IDIM).EQ.0) THEN
            Y = (GLOB_ROW-(X-1)*IDIM*IDIM)/IDIM
          ELSE
            Y = (GLOB_ROW-(X-1)*IDIM*IDIM)/IDIM+1
          ENDIF
          Z = GLOB_ROW-(X-1)*IDIM*IDIM-(Y-1)*IDIM
          ! GLOB_X, GLOB_Y, GLOB_X coordinates
          GLOB_X=X*DELTAH
          GLOB_Y=Y*DELTAH
          GLOB_Z=Z*DELTAH

          ! Check on boundary points 
          IF (X.EQ.1) THEN
            ROW_MAT%ASPK(ELEMENT)=ONE
            ROW_MAT%IA1(ELEMENT)=(X-1)*IDIM*IDIM+(Y-1)*IDIM+(Z)
            ELEMENT=ELEMENT+1
          ELSE IF (Y.EQ.1) THEN
            ROW_MAT%ASPK(ELEMENT)=ONE
            ROW_MAT%IA1(ELEMENT)=(X-1)*IDIM*IDIM+(Y-1)*IDIM+(Z)
            ELEMENT=ELEMENT+1
          ELSE IF (Z.EQ.1) THEN
            ROW_MAT%ASPK(ELEMENT)=ONE
            ROW_MAT%IA1(ELEMENT)=(X-1)*IDIM*IDIM+(Y-1)*IDIM+(Z)
            ELEMENT=ELEMENT+1
          ELSE IF (X.EQ.IDIM) THEN
            ROW_MAT%ASPK(ELEMENT)=ONE
            ROW_MAT%IA1(ELEMENT)=(X-1)*IDIM*IDIM+(Y-1)*IDIM+(Z)
            ELEMENT=ELEMENT+1
          ELSE IF (Y.EQ.IDIM) THEN
            ROW_MAT%ASPK(ELEMENT)=ONE
            ROW_MAT%IA1(ELEMENT)=(X-1)*IDIM*IDIM+(Y-1)*IDIM+(Z)
            ELEMENT=ELEMENT+1
          ELSE IF (Z.EQ.IDIM) THEN
            ROW_MAT%ASPK(ELEMENT)=ONE
            ROW_MAT%IA1(ELEMENT)=(X-1)*IDIM*IDIM+(Y-1)*IDIM+(Z)
            ELEMENT=ELEMENT+1
          ELSE
            ! Internal point: build discretization
            !   
            !  Term depending on   (x-1,y,z)
            !
            ROW_MAT%ASPK(ELEMENT)=-B1(GLOB_X,GLOB_Y,GLOB_Z)&
                 & -A1(GLOB_X,GLOB_Y,GLOB_Z)
            ROW_MAT%ASPK(ELEMENT) = ROW_MAT%ASPK(ELEMENT)/(DELTAH*&
                 & DELTAH)
            ROW_MAT%IA1(ELEMENT)=(X-2)*IDIM*IDIM+(Y-1)*IDIM+(Z)
            ELEMENT=ELEMENT+1
            !  Term depending on     (x,y-1,z)
            ROW_MAT%ASPK(ELEMENT)=-B2(GLOB_X,GLOB_Y,GLOB_Z)&
                 & -A2(GLOB_X,GLOB_Y,GLOB_Z)
            ROW_MAT%ASPK(ELEMENT) = ROW_MAT%ASPK(ELEMENT)/(DELTAH*&
                 & DELTAH)
            ROW_MAT%IA1(ELEMENT)=(X-1)*IDIM*IDIM+(Y-2)*IDIM+(Z)
            ELEMENT=ELEMENT+1
            !  Term depending on     (x,y,z-1)
            ROW_MAT%ASPK(ELEMENT)=-B3(GLOB_X,GLOB_Y,GLOB_Z)&
                 & -A3(GLOB_X,GLOB_Y,GLOB_Z)
            ROW_MAT%ASPK(ELEMENT) = ROW_MAT%ASPK(ELEMENT)/(DELTAH*&
                 & DELTAH)
            ROW_MAT%IA1(ELEMENT)=(X-1)*IDIM*IDIM+(Y-1)*IDIM+(Z-1)
            ELEMENT=ELEMENT+1
            !  Term depending on     (x,y,z)
            ROW_MAT%ASPK(ELEMENT)=2*B1(GLOB_X,GLOB_Y,GLOB_Z)&
                 & +2*B2(GLOB_X,GLOB_Y,GLOB_Z)&
                 & +2*B3(GLOB_X,GLOB_Y,GLOB_Z)&
                 & +A1(GLOB_X,GLOB_Y,GLOB_Z)&
                 & +A2(GLOB_X,GLOB_Y,GLOB_Z)&
                 & +A3(GLOB_X,GLOB_Y,GLOB_Z)
            ROW_MAT%ASPK(ELEMENT) = ROW_MAT%ASPK(ELEMENT)/(DELTAH*&
                 & DELTAH)
            ROW_MAT%IA1(ELEMENT)=(X-1)*IDIM*IDIM+(Y-1)*IDIM+(Z)
            ELEMENT=ELEMENT+1                  
            !  Term depending on     (x,y,z+1)
            ROW_MAT%ASPK(ELEMENT)=-B1(GLOB_X,GLOB_Y,GLOB_Z)
            ROW_MAT%ASPK(ELEMENT) = ROW_MAT%ASPK(ELEMENT)/(DELTAH*&
                 & DELTAH)
            ROW_MAT%IA1(ELEMENT)=(X-1)*IDIM*IDIM+(Y-1)*IDIM+(Z+1)
            ELEMENT=ELEMENT+1
            !  Term depending on     (x,y+1,z)
            ROW_MAT%ASPK(ELEMENT)=-B2(GLOB_X,GLOB_Y,GLOB_Z)
            ROW_MAT%ASPK(ELEMENT) = ROW_MAT%ASPK(ELEMENT)/(DELTAH*&
                 & DELTAH)
            ROW_MAT%IA1(ELEMENT)=(X-1)*IDIM*IDIM+(Y)*IDIM+(Z)
            ELEMENT=ELEMENT+1
            !  Term depending on     (x+1,y,z)
            ROW_MAT%ASPK(ELEMENT)=-B3(GLOB_X,GLOB_Y,GLOB_Z)
            ROW_MAT%ASPK(ELEMENT) = ROW_MAT%ASPK(ELEMENT)/(DELTAH*&
                 & DELTAH)
            ROW_MAT%IA1(ELEMENT)=(X)*IDIM*IDIM+(Y-1)*IDIM+(Z)
            ELEMENT=ELEMENT+1
          ENDIF
          ROW_MAT%M=1
          ROW_MAT%K=N
          ROW_MAT%IA2(2)=ELEMENT       
          ! IA== GLOBAL ROW INDEX
          IA=GLOB_ROW
!!$             IA=(X-1)*IDIM*IDIM+(Y-1)*IDIM+(Z)
!!$        write(0,*) 'Inserting row ',ia,' On proc',myprow
          T3 = MPI_WTIME()
          CALL F90_PSSPINS(A,IA,1,ROW_MAT,IERRV,DESC_A)       
          if (ierrv(1).ne.0) then 
            write(0,*) 'On row ',ia,' IERRV:',ierrv(:)
          endif
          TINS = TINS + (MPI_WTIME()-T3)
          ! Build RHS  
          IF (X==1) THEN     
            GLOB_Y=(Y-IDIM/2)*DELTAH
            GLOB_Z=(Z-IDIM/2)*DELTAH        
            ZT(1) = EXP(-GLOB_Y**2-GLOB_Z**2)
          ELSE IF ((Y==1).OR.(Y==IDIM).OR.(Z==1).OR.(Z==IDIM)) THEN 
            GLOB_X=3*(X-1)*DELTAH
            GLOB_Y=(Y-IDIM/2)*DELTAH
            GLOB_Z=(Z-IDIM/2)*DELTAH        
            ZT(1) = EXP(-GLOB_Y**2-GLOB_Z**2)*EXP(-GLOB_X)
          ELSE
            ZT(1) = 0.D0
          ENDIF
          CALL F90_PSDSINS(1,B,IA,ZT(1:1),IERRV,DESC_A)
          ZT(1)=0.D0
          CALL F90_PSDSINS(1,T,IA,ZT(1:1),IERRV,DESC_A)
        END IF
      END DO
    END DO

    CALL BLACS_BARRIER(ICONTXT,'ALL')    
    T2 = MPI_WTIME()

    WRITE(*,*) '     pspins  time',TINS
    WRITE(*,*) '   Insert time',(T2-T1)

    DEALLOCATE(ROW_MAT%ASPK,ROW_MAT%IA1,ROW_MAT%IA2)

    write(*,*) 'Calling SPASB'
    CALL BLACS_BARRIER(ICONTXT,'ALL')
    T1 = MPI_WTIME()

    CALL F90_PSSPASB(A,IERRV,DESC_A,AFMT=AFMT,DUP=2) 

    CALL BLACS_BARRIER(ICONTXT,'ALL')
    T2 = MPI_WTIME()

    WRITE(0,*) '   Assembly  time',(T2-T1),' ',a%fida(1:4)

    CALL F90_PSDSASB(B,IERRV,DESC_A)
    CALL F90_PSDSASB(T,IERRV,DESC_A)
    IF (MYPROW.EQ.0) THEN
      WRITE(0,*) '   End CREATE_MATRIX'
    ENDIF
    RETURN

  END SUBROUTINE CREATE_MATRIX
END PROGRAM PDE90
!
! Functions parametrizing the differential equation 
!  
FUNCTION A1(X,Y,Z)
  REAL(KIND(1.D0)) :: A1
  REAL(KIND(1.D0)) :: X,Y,Z
  A1=1.D0
END FUNCTION A1
FUNCTION A2(X,Y,Z)
  REAL(KIND(1.D0)) ::  A2
  REAL(KIND(1.D0)) :: X,Y,Z
  A2=2.D1*Y
END FUNCTION A2
FUNCTION A3(X,Y,Z)
  REAL(KIND(1.D0)) ::  A3
  REAL(KIND(1.D0)) :: X,Y,Z      
  A3=1.D0
END FUNCTION A3
FUNCTION A4(X,Y,Z)
  REAL(KIND(1.D0)) ::  A4
  REAL(KIND(1.D0)) :: X,Y,Z      
  A4=1.D0
END FUNCTION A4
FUNCTION B1(X,Y,Z)
  REAL(KIND(1.D0)) ::  B1   
  REAL(KIND(1.D0)) :: X,Y,Z
  B1=1.D0
END FUNCTION B1
FUNCTION B2(X,Y,Z)
  REAL(KIND(1.D0)) ::  B2
  REAL(KIND(1.D0)) :: X,Y,Z
  B2=1.D0
END FUNCTION B2
FUNCTION B3(X,Y,Z)
  REAL(KIND(1.D0)) ::  B3
  REAL(KIND(1.D0)) :: X,Y,Z
  B3=1.D0
END FUNCTION B3


