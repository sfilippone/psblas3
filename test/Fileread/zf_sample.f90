PROGRAM ZF_SAMPLE
  USE TYPESP
  USE TYPEDESC
  USE F90TOOLS
  USE F90PSBLAS
  USE F90METHD
  USE MAT_DIST
  USE READ_MAT
  USE PARTGRAPH
  USE GETP
  IMPLICIT NONE

  ! Input parameters
  CHARACTER*20 :: CMETHD, PREC, MTRX_FILE, RHS_FILE
  CHARACTER*80 :: CHARBUF

  DOUBLE PRECISION DDOT
  EXTERNAL DDOT
  INTERFACE 
    !   .....user passed subroutine.....
    SUBROUTINE PART_BLOCK(GLOBAL_INDX,N,NP,PV,NV)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: GLOBAL_INDX, N, NP
      INTEGER, INTENT(OUT) :: NV
      INTEGER, INTENT(OUT) :: PV(*) 
    END SUBROUTINE PART_BLOCK
  END INTERFACE   ! Local variables

  INTEGER, PARAMETER    :: IZERO=0, IONE=1
  CHARACTER, PARAMETER  :: ORDER='R'
  COMPLEX(KIND(1.D0)), POINTER,SAVE :: B_COL(:), X_COL(:), R_COL(:), &
       & B_COL_GLOB(:), X_COL_GLOB(:), R_COL_GLOB(:), B_GLOB(:,:), &
       &Z(:), Q(:),Z1(:)
  INTEGER              :: IARGC
  Real(Kind(1.d0)), Parameter :: Dzero = 0.d0, One = 1.d0
  Real(Kind(1.d0)) :: MPI_WTIME, T1, T2, TPREC, R_AMAX, B_AMAX,bb(1,1),lambda,scale,resmx,resmxp
  integer :: nrhs, nrow, nx1, nx2, n_row
  External IARGC, MPI_WTIME
  integer bsze,overlap
  common/part/bsze,overlap

  ! Sparse Matrices
  TYPE(Z_SPMAT) :: A, AUX_A, L, U
!!$  TYPE(D_PRECN) :: APRC
  ! Dense Matrices
  COMPLEX(KIND(1.D0)), POINTER ::  AUX_B(:,:) , AUX1(:), AUX2(:), VDIAG(:),&
       & AUX_G(:,:), AUX_X(:,:)

  ! Communications data structure
  TYPE(desc_type)    :: DESC_A

  ! BLACS parameters
  INTEGER            :: NPROW, NPCOL, ICTXT, IAM, NP, MYPROW, MYPCOL
  logical            :: amroot

  ! Solver paramters
  INTEGER            :: ITER, ITMAX, IERR, ITRACE, IRCODE, IPART,&
       & IPREC, METHD, ISTOPC, ML
  integer, pointer :: ierrv(:)
  REAL(KIND(1.D0))   :: ERR, EPS
  integer   iparm(20)
  real(kind(1.d0)) rparm(20)

  ! Other variables
  INTEGER            :: I,INFO,J
  INTEGER            :: INTERNAL, M,II,NNZERO

  ! common area
  INTEGER M_PROBLEM, NPROC

  allocate(ierrv(6))
  ! Initialize BLACS
  CALL BLACS_PINFO(IAM, NP)
  CALL BLACS_GET(IZERO, IZERO, ICTXT)

  ! Rectangular Grid,  Np x 1

  CALL BLACS_GRIDINIT(ICTXT, ORDER, NP, IONE)
  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYPROW, MYPCOL)
  amroot = (myprow==0).and.(mypcol==0)
  !
  !  Get parameters
  !
  CALL GET_PARMS(ICTXT,MTRX_FILE,RHS_FILE,CMETHD,PREC,&
       & IPART,ISTOPC,ITMAX,ITRACE,ML,IPREC,EPS)
  CALL BLACS_BARRIER(ICTXT,'A')
  T1 = MPI_WTIME()  
  ! Read the input matrix to be processed and (possibly) the RHS 
  NRHS = 1
  IF (amroot) THEN
    NULLIFY(AUX_B)
    CALL ZREADMAT(MTRX_FILE, AUX_A, ICTXT)
    M_PROBLEM = AUX_A%M
    CALL IGEBS2D(ICTXT,'A',' ',1,1,M_PROBLEM,1)

    IF(RHS_FILE.NE.'NONE') THEN
       !  Reading an RHS
       CALL ZREAD_RHS(RHS_FILE,AUX_B,ICTXT)
    END IF

    IF (ASSOCIATED(AUX_B).and.SIZE(AUX_B,1)==M_PROBLEM) THEN
      ! If any RHS were present, broadcast the first one
      write(0,*) 'Ok, got an RHS ',aux_b(m_problem,1)
      B_COL_GLOB =>AUX_B(:,1)
    ELSE
      write(0,*) 'Inventing an RHS '
      ALLOCATE(AUX_B(M_PROBLEM,1), STAT=IRCODE)
      IF (IRCODE /= 0) THEN
        WRITE(0,*) 'Memory allocation failure in ZF_SAMPLE'
        CALL BLACS_ABORT(ICTXT,-1)
        STOP
      ENDIF
      write(0,*) 'Check on AUX_B ',size(aux_b,1),size(aux_b,2),m_problem
      B_COL_GLOB => AUX_B(:,1)
      
      DO I=1, M_PROBLEM
        B_COL_GLOB(I) = CMPLX(I*2.0/M_PROBLEM,0)
      ENDDO
    ENDIF
    CALL ZGEBS2D(ICTXT,'A',' ',M_PROBLEM,1,B_COL_GLOB,M_PROBLEM) 

  ELSE
    CALL IGEBR2D(ICTXT,'A',' ',1,1,M_PROBLEM,1,0,0)
    WRITE(0,*) 'Receiving AUX_B'
    ALLOCATE(AUX_B(M_PROBLEM,1), STAT=IRCODE)
    IF (IRCODE /= 0) THEN
      WRITE(0,*) 'Memory allocation failure in ZF_SAMPLE'
      CALL BLACS_ABORT(ICTXT,-1)
      STOP
    ENDIF
    write(0,*) 'Check on AUX_B ',size(aux_b,1),size(aux_b,2),m_problem
    B_COL_GLOB =>AUX_B(:,1)
    CALL ZGEBR2D(ICTXT,'A',' ',M_PROBLEM,1,B_COL_GLOB,M_PROBLEM,0,0) 
  END IF
  NPROC = NPROW 

  ! Switch over different partition types
  IF (IPART.EQ.0) THEN 
    CALL BLACS_BARRIER(ICTXT,'A')
    WRITE(6,*) 'Partition type: BLOCK'
    CALL ZMATDIST(AUX_A, A, PART_BLOCK, ICTXT, &
         & DESC_A,B_COL_GLOB,B_COL)
  ELSE IF (IPART.EQ.2) THEN 
    WRITE(6,*) amroot,' Partition type: GRAPH'
    IF (amroot) THEN 
      CALL BUILD_GRPPART(AUX_A%M,AUX_A%FIDA,AUX_A%IA1,AUX_A%IA2,NP)
    ENDIF
    call blacs_barrier(ictxt,'All')
    CALL DISTR_GRPPART(0,0,ICTXT)
    CALL ZMATDIST(AUX_A, A, PART_GRAPH, ICTXT, &
         & DESC_A,B_COL_GLOB,B_COL)
  ELSE 
    WRITE(6,*) 'Partition type: BLOCK'
   CALL ZMATDIST(AUX_A, A, PART_BLOCK, ICTXT, &
        & DESC_A,B_COL_GLOB,B_COL)
  END IF
  
  write(*,*) amroot,' Done matdist'
  CALL F90_PSDSALL(M_PROBLEM,X_COL,IERRV,DESC_A)
  X_COL(:) =0.0
  CALL F90_PSDSASB(X_COL,IERRV,DESC_A)
  CALL F90_PSDSALL(M_PROBLEM,R_COL,IERRV,DESC_A)
  R_COL(:) =0.0
  CALL F90_PSDSASB(R_COL,IERRV,DESC_A)
  T2 = MPI_WTIME() - T1
  
  CALL DGAMX2D(ICTXT, 'A', ' ', IONE, IONE, T2, IONE,&
       & T1, T1, -1, -1, -1)
  
  IF (amroot) THEN
     WRITE(6,*) 'Time to Read and Partition Matrix : ',T2
  END IF
  
  !
  !  Prepare the preconditioning matrix. Note the availability
  !  of optional parameters
  !

  IF (amroot) WRITE(6,*) 'Preconditioner : "',PREC(1:6),'"  ',iprec
  T1 = MPI_WTIME()
  CALL PRECONDITIONER(IPREC,A,L,U,VDIAG,DESC_A,INFO)
  
  TPREC = MPI_WTIME()-T1
  
  
  CALL DGAMX2D(ICTXT,'A',' ',IONE, IONE,TPREC,IONE,T1,T1,-1,-1,-1)
  
  WRITE(0,*) 'Preconditioner Time : ',TPREC,'  ',&
       &prec,iprec
  IF (INFO /= 0) THEN
    WRITE(0,*) 'Error in preconditioner :',INFO
    CALL BLACS_ABORT(ICTXT,-1)
    STOP
  END IF
  IPARM = 0
  RPARM = 0.D0
  write(0,*) amroot,'Starting method'
  CALL BLACS_BARRIER(ICTXT,'All')
  T1 = MPI_WTIME()
  IF (CMETHD.EQ.'BICGSTAB') Then
     CALL  F90_BICGSTAB(A,IPREC,L,U,VDIAG,B_COL,X_COL,EPS,DESC_A,& 
       & ITMAX,ITER,ERR,IERR,ITRACE)     
!!$  ELSE IF (CMETHD.EQ.'BICG') Then
!!$    CALL  F90_BICG(A,IPREC,L,U,VDIAG,B_COL,X_COL,EPS,DESC_A,& 
!!$       & ITMAX,ITER,ERR,IERR,ITRACE)     
!!$  ELSE IF (CMETHD.EQ.'CGS') Then
!!$    CALL  F90_CGS(A,IPREC,L,U,VDIAG,B_COL,X_COL,EPS,DESC_A,& 
!!$       & ITMAX,ITER,ERR,IERR,ITRACE)     
!!$  ELSE IF (CMETHD.EQ.'BICGSTABL') Then
!!$    CALL  F90_BICGSTABL(A,IPREC,L,U,VDIAG,B_COL,X_COL,EPS,DESC_A,& 
!!$       & ITMAX,ITER,ERR,IERR,ITRACE,ML)     
  ENDIF
  CALL BLACS_BARRIER(ICTXT,'All')
  T2 = MPI_WTIME() - T1
  CALL DGAMX2D(ICTXT,'A',' ',IONE, IONE,T2,IONE,T1,T1,-1,-1,-1)
  call f90_psaxpby((1.d0,0.d0),b_col,(0.d0,0.d0),r_col,desc_A)
  call f90_psspmm((-1.d0,0.d0),a,x_col,(1.d0,0.d0),r_col,desc_a)
  call f90_amax(resmx,r_col,desc_a)
  where (b_col/= 0)
    r_col = r_col/b_col
  end where
  call f90_amax(resmxp,r_col,desc_a)

!!$  ITER=IPARM(5)
!!$  ERR = RPARM(2)
  if (amroot) then
     write(6,*) 'methd iprec istopc   : ',iprec, istopc
     write(6,*) 'Number of iterations : ',iter
     write(6,*) 'Time to Solve Matrix : ',T2
     WRITE(6,*) 'Time per iteration   : ',T2/(ITER)
     WRITE(6,*) 'Error on exit        : ',ERR
  END IF
  
  allocate(x_col_glob(m_problem),r_col_glob(m_problem),stat=ierr)
  if (ierr.ne.0) then 
    write(0,*) 'Allocation error: no data collection'
  else
    call f90_psdgatherm(x_col_glob,x_col,desc_a,iroot=0)
    call f90_psdgatherm(r_col_glob,r_col,desc_a,iroot=0)
    if (amroot) then
      write(0,*) 'Saving X on file'
      write(20,*) 'Matrix: ',mtrx_file
      write(20,*) 'Computed solution on ',NPROW,' processors.'
      write(20,*) 'Iterations to convergence: ',iter
      write(20,*) 'Error indicator (infinity norm) on exit:', &
           & ' ||r||/(||A||||x||+||b||) = ',err
      write(20,*) 'Max residual = ',resmx, resmxp
      do i=1,m_problem
        write(20,998) i,x_col_glob(i),r_col_glob(i),b_col_glob(i)
      enddo
    end if
  end if
998 format(i8,4(2x,g20.14))
993 format(i6,4(1x,e12.6))

  
!!$  !
!!$  ! Raleygh quotients for first eigenvalue
!!$  ! 
!!$  CALL F90_PSDSall(M_problem,Q,ierrv,DESC_A)
!!$  CALL F90_PSDSall(M_problem,Z,ierrv,DESC_A)
!!$  CALL F90_PSDSall(M_problem,Z1,ierrv,DESC_A)
!!$  CALL F90_PSDSasb(Q,ierrv,DESC_A)
!!$  CALL F90_PSDSasb(Z,ierrv,DESC_A)
!!$  CALL F90_PSDSasb(Z1,ierrv,DESC_A)
!!$  scale = f90_psnrm2(x_col,desc_a)
!!$  scale = one/scale
!!$  call f90_psaxpby(scale,x_col,dzero,q,desc_A)
!!$  call f90_psspmm(one,a,q,dzero,z,desc_a)
!!$  do i=1, itmax   
!!$    scale  = f90_psnrm2(z,desc_a)
!!$    scale = one/scale
!!$    call f90_psaxpby(one,z,dzero,z1,desc_a)
!!$    call f90_psaxpby(scale,z,dzero,q,desc_a)
!!$    call f90_psspmm(one,a,q,dzero,z,desc_a)
!!$    lambda = f90_psdot(q,z,desc_A)
!!$    scale  = f90_psnrm2(z,desc_A)
!!$    if (amroot) write(0,*) 'Lambda: ',i,lambda, scale
!!$  enddo
!!$  call f90_psaxpby(-one,z,one,z1,desc_a)
!!$  scale  = f90_psnrm2(z1,desc_A)
!!$  if (amroot) write(0,*) 'Final check: ',i,lambda, scale
!!$  do i=1, desc_a%matrix_data(n_row_)
!!$    scale=z(i)/q(i)
!!$    write(0,*) 'Vector check: ',i,lambda, scale,abs(scale-lambda)
!!$  enddo

  CALL F90_PSDSFREE(B_COL, DESC_A)
  CALL F90_PSDSFREE(X_COL, DESC_A)
  CALL F90_PSSPFREE(A, DESC_A)
  IF (IPREC.GE.2) THEN
     CALL F90_PSSPFREE(L, DESC_A)
     CALL F90_PSSPFREE(U, DESC_A)
  END IF
  CALL F90_PSDSCFREE(DESC_A,info)
  CALL BLACS_GRIDEXIT(ICTXT)
  CALL BLACS_EXIT(0)
  
END PROGRAM ZF_SAMPLE
  




