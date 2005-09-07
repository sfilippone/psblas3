PROGRAM DF_SAMPLE
  USE F90SPARSE
  USE MAT_DIST
  USE READ_MAT
  USE PARTGRAPH
  USE GETP
  use mpi
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
  REAL(KIND(1.D0)), POINTER, SAVE :: B_COL(:), X_COL(:), R_COL(:), &
       & B_COL_GLOB(:), X_COL_GLOB(:), R_COL_GLOB(:), B_GLOB(:,:), &
       &Z(:), Q(:),Z1(:)  
  INTEGER              :: IARGC, CHECK_DESCR, CONVERT_DESCR
  Real(Kind(1.d0)), Parameter :: Dzero = 0.d0, One = 1.d0
  Real(Kind(1.d0)) :: T1, T2, TPREC, R_AMAX, B_AMAX,bb(1,1),&
       &lambda,scale,resmx,resmxp
  integer :: nrhs, nrow, nx1, nx2, n_row, dim,iread
  logical :: amroot
  integer mpe_log_get_event_number,mpe_Describe_state,mpe_log_event

  External IARGC
  integer bsze,overlap
  common/part/bsze,overlap
  INTEGER, POINTER :: WORK(:)
  ! Sparse Matrices
  TYPE(D_SPMAT) :: A, AUX_A, H
  TYPE(D_PREC) :: PRE
!!$  TYPE(D_PRECN) :: APRC
  ! Dense Matrices
  REAL(KIND(1.D0)), POINTER ::  AUX_B(:,:) , AUX1(:), AUX2(:), VDIAG(:), &
       &  AUX_G(:,:), AUX_X(:,:), D(:)

  ! Communications data structure
  TYPE(desc_type)    :: DESC_A, DESC_A_OUT

  ! BLACS parameters
  INTEGER            :: NPROW, NPCOL, ICTXT, IAM, NP, MYPROW, MYPCOL

  ! Solver paramters
  INTEGER            :: ITER, ITMAX, IERR, ITRACE, IRCODE, IPART,&
       & METHD, ISTOPC, ML
  integer, pointer :: ierrv(:)
  REAL(KIND(1.D0))   :: ERR, EPS
  integer   iparm(20)
  real(kind(1.d0)) rparm(20)

  ! Other variables
  INTEGER            :: I,INFO,J, iprecb,iprece,islvb,islve
  INTEGER            :: INTERNAL, M,II,NNZERO

  ! common area
  INTEGER M_PROBLEM, NPROC
  

  allocate(ierrv(6))
  ! Initialize BLACS
  CALL BLACS_PINFO(IAM, NP)
  CALL BLACS_GET(IZERO, IZERO, ICTXT)
  iprecb = mpe_log_get_event_number()
  iprece = mpe_log_get_event_number()
  islvb  = mpe_log_get_event_number()
  islve  = mpe_log_get_event_number()
  if (iam==0) then 
    info = mpe_describe_state(iprecb,iprece,"Preconditioner","OrangeRed")
    info = mpe_describe_state(islvb,islve,"Solver","DarkGreen")
  endif

  ! Rectangular Grid,  Np x 1

  CALL BLACS_GRIDINIT(ICTXT, ORDER, NP, IONE)
  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYPROW, MYPCOL)
  AMROOT = (MYPROW==0).AND.(MYPCOL==0)

  !
  !  Get parameters
  !
  CALL GET_PARMS(ICTXT,MTRX_FILE,RHS_FILE,CMETHD,PREC,&
       & IPART,ISTOPC,ITMAX,ITRACE,ML,PRE%PREC,EPS)
  
  CALL BLACS_BARRIER(ICTXT,'A')
  T1 = MPI_WTIME()  
  ! Read the input matrix to be processed and (possibly) the RHS 
  NRHS = 1
  NPROC = NPROW 

  IF (AMROOT) THEN
    NULLIFY(AUX_B)
    CALL READMAT(MTRX_FILE, AUX_A, ICTXT)

    M_PROBLEM = AUX_A%M
    CALL IGEBS2D(ICTXT,'A',' ',1,1,M_PROBLEM,1)

    IF(RHS_FILE /= 'NONE') THEN
       !  Reading an RHS
       CALL READ_RHS(RHS_FILE,AUX_B,ICTXT)
    END IF

    IF (ASSOCIATED(AUX_B).and.SIZE(AUX_B,1)==M_PROBLEM) THEN
      ! If any RHS were present, broadcast the first one
      write(0,*) 'Ok, got an RHS ',aux_b(m_problem,1)
      B_COL_GLOB =>AUX_B(:,1)
    ELSE
      write(0,*) 'Inventing an RHS '
      ALLOCATE(AUX_B(M_PROBLEM,1), STAT=IRCODE)
      IF (IRCODE /= 0) THEN
        WRITE(0,*) 'Memory allocation failure in DF_SAMPLE'
        CALL BLACS_ABORT(ICTXT,-1)
        STOP
      ENDIF
      B_COL_GLOB =>AUX_B(:,1)
      DO I=1, M_PROBLEM
        B_COL_GLOB(I) = REAL(I)*2.0/REAL(M_PROBLEM)
      ENDDO      
    ENDIF
    CALL DGEBS2D(ICTXT,'A',' ',M_PROBLEM,1,B_COL_GLOB,M_PROBLEM) 
  ELSE
    CALL IGEBR2D(ICTXT,'A',' ',1,1,M_PROBLEM,1,0,0)
    WRITE(0,*) 'Receiving AUX_B'
    ALLOCATE(AUX_B(M_PROBLEM,1), STAT=IRCODE)
    IF (IRCODE /= 0) THEN
      WRITE(0,*) 'Memory allocation failure in DF_SAMPLE'
      CALL BLACS_ABORT(ICTXT,-1)
      STOP
    ENDIF
    B_COL_GLOB =>AUX_B(:,1)
    CALL DGEBR2D(ICTXT,'A',' ',M_PROBLEM,1,B_COL_GLOB,M_PROBLEM,0,0) 
  END IF

  ! Switch over different partition types
  IF (IPART.EQ.0) THEN 
    CALL BLACS_BARRIER(ICTXT,'A')
    WRITE(6,*) 'Partition type: BLOCK'
    CALL MATDIST(AUX_A, A, PART_BLOCK, ICTXT, &
         & DESC_A,B_COL_GLOB,B_COL)
  ELSE IF (IPART.EQ.2) THEN 
    WRITE(0,*) 'Partition type: GRAPH'
    IF (AMROOT) THEN 
!!$      WRITE(0,*) 'Call BUILD',size(aux_a%ia1),size(aux_a%ia2),np
      CALL BUILD_GRPPART(AUX_A%M,AUX_A%FIDA,AUX_A%IA1,AUX_A%IA2,NP)
    ENDIF
    CALL BLACS_BARRIER(ICTXT,'A')
!!$    WRITE(0,*) myprow,'Done BUILD_GRPPART'
    CALL DISTR_GRPPART(0,0,ICTXT)
!!$    WRITE(0,*) myprow,'Done DISTR_GRPPART'
    CALL MATDIST(AUX_A, A, PART_GRAPH, ICTXT, &
         & DESC_A,B_COL_GLOB,B_COL)
  ELSE 
    WRITE(6,*) 'Partition type: BLOCK'
    CALL MATDIST(AUX_A, A, PART_BLOCK, ICTXT, &
         & DESC_A,B_COL_GLOB,B_COL)
  END IF
  
  CALL F90_PSDSALL(M_PROBLEM,X_COL,IERRV,DESC_A)
  X_COL(:) =0.0
  CALL F90_PSDSASB(X_COL,IERRV,DESC_A)
  CALL F90_PSDSALL(M_PROBLEM,R_COL,IERRV,DESC_A)
  R_COL(:) =0.0
  CALL F90_PSDSASB(R_COL,IERRV,DESC_A)
  T2 = MPI_WTIME() - T1
  
  
  DIM=SIZE(A%ASPK)

!!$  ALLOCATE(H%ASPK(DIM),H%IA1(DIM),H%IA2(DIM),H%PL(SIZE(A%PL)),&
!!$       & H%PL(SIZE(A%PL)),D(SIZE(A%PL)),&
!!$       & DESC_A_OUT%MATRIX_DATA(SIZE(DESC_A%MATRIX_DATA)),&
!!$       & DESC_A_OUT%HALO_INDEX(SIZE(DESC_A%HALO_INDEX)),&
!!$       & DESC_A_OUT%OVRLAP_INDEX(SIZE(DESC_A%OVRLAP_INDEX)),&
!!$       & DESC_A_OUT%OVRLAP_ELEM(SIZE(DESC_A%OVRLAP_ELEM)),&
!!$       & DESC_A_OUT%LOC_TO_GLOB(SIZE(DESC_A%LOC_TO_GLOB)),&
!!$       & DESC_A_OUT%GLOB_TO_LOC(SIZE(DESC_A%GLOB_TO_LOC)), WORK(dim))
!!$  check_descr=15
!  work(5)=9
!!$  WRITE(0,*)'CALLING VERIFY'
!!$  CALL F90_PSVERIFY(D,A,DESC_A,CHECK_DESCR,CONVERT_DESCR,H,&
!!$       & DESC_A_OUT,WORK)
!!$  WRITE(0,*)'VERIFY DONE',CONVERT_DESCR

!  deallocate(work)


  CALL DGAMX2D(ICTXT, 'A', ' ', IONE, IONE, T2, IONE,&
       & T1, T1, -1, -1, -1)
  
  IF (AMROOT) THEN
     WRITE(6,*) 'Time to Read and Partition Matrix : ',T2
  END IF
  
  !
  !  Prepare the preconditioning matrix. Note the availability
  !  of optional parameters
  !

  IF (AMROOT) WRITE(6,*) 'Preconditioner : "',PREC(1:6),'"  ',PRE%PREC


!!$  do i=1,a%m
!!$     do j=a%ia2(i),a%ia2(i+1)-1
!!$        write(0,*)'a ',i,a%ia1(j),a%aspk(j)
!!$     end do
!!$  end do
!!$
!!$  write(0,*)'halo_index',desc_a%halo_index(:)
!!$  write(0,*)'ovrlap_index',desc_a%ovrlap_index(:)
!!$  write(0,*)'ovrlap_elem',desc_a%ovrlap_elem(:)

  info = MPE_Log_event( iprecb, 0, "start Precond" )
  T1 = MPI_WTIME()

  CALL PRECONDITIONER(A,PRE,DESC_A,INFO)!,'F')
  TPREC = MPI_WTIME()-T1
  info = MPE_Log_event( iprece, 0, "end Precond" )
  
  
  CALL DGAMX2D(ICTXT,'A',' ',IONE, IONE,TPREC,IONE,T1,T1,-1,-1,-1)
  
  WRITE(0,*) 'Preconditioner Time :',TPREC,' ',&
       &prec,pre%prec
  IF (INFO /= 0) THEN
    WRITE(0,*) 'Error in preconditioner :',INFO
    CALL BLACS_ABORT(ICTXT,-1)
    STOP
  END IF
  
  IPARM = 0
  RPARM = 0.D0   
  CALL BLACS_BARRIER(ICTXT,'All')
  info = MPE_Log_event( islvb, 0, "start Solver" )
  T1 = MPI_WTIME()
  IF (CMETHD.EQ.'BICGSTAB') Then
    CALL  F90_BICGSTAB(A,PRE,B_COL,X_COL,EPS,DESC_A,& 
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
  T2 = MPI_WTIME() - T1
  info = MPE_Log_event( islve, 0, "end Solver" )  
  CALL DGAMX2D(ICTXT,'A',' ',IONE, IONE,T2,IONE,T1,T1,-1,-1,-1)
  call f90_psaxpby(1.d0,b_col,0.d0,r_col,desc_A)
  call f90_psspmm(-1.d0,a,x_col,1.d0,r_col,desc_a)
  call f90_amax(resmx,r_col,desc_a)
  where (b_col /= 0.d0)
    r_col = r_col/b_col
  end where
  call f90_amax(resmxp,r_col,desc_a)

!!$  ITER=IPARM(5)
!!$  ERR = RPARM(2)
  if (amroot) then
     write(6,*) 'methd iprec istopc   : ',pre%prec, istopc
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
  CALL F90_PSPRECFREE(PRE,info)
  CALL F90_PSDSCFREE(DESC_A,info)
  CALL BLACS_GRIDEXIT(ICTXT)
  CALL BLACS_EXIT(0)
  
END PROGRAM DF_SAMPLE
  




