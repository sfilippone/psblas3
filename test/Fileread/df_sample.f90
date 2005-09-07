PROGRAM DF_SAMPLE
  USE F90SPARSE
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
  INTERFACE 
    !   .....user passed subroutine.....
    SUBROUTINE PART_BLK2(GLOBAL_INDX,N,NP,PV,NV)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: GLOBAL_INDX, N, NP
      INTEGER, INTENT(OUT) :: NV
      INTEGER, INTENT(OUT) :: PV(*) 
    END SUBROUTINE PART_BLK2
  END INTERFACE   ! Local variables


  INTEGER, PARAMETER    :: IZERO=0, IONE=1
  CHARACTER, PARAMETER  :: ORDER='R'
  REAL(KIND(1.D0)), POINTER, SAVE :: B_COL(:), X_COL(:), R_COL(:), &
       & B_COL_GLOB(:), X_COL_GLOB(:), R_COL_GLOB(:), B_GLOB(:,:), &
       &Z(:), Q(:),Z1(:)  
  INTEGER              :: IARGC, CHECK_DESCR, CONVERT_DESCR
  Real(Kind(1.d0)), Parameter :: Dzero = 0.d0, One = 1.d0
  Real(Kind(1.d0)) :: MPI_WTIME, T1, T2, TPREC, R_AMAX, B_AMAX,bb(1,1),&
       &lambda,scale,resmx,resmxp
  integer :: nrhs, nrow, nx1, nx2, n_row, dim,iread
  logical :: amroot
  External IARGC, MPI_WTIME
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
       & METHD, ISTOPC, ML, iprec, novr
  integer, pointer :: ierrv(:)
  character(len=5)   :: afmt
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
  AMROOT = (MYPROW==0).AND.(MYPCOL==0)

  !
  !  Get parameters
  !
  CALL GET_PARMS(ICTXT,MTRX_FILE,RHS_FILE,CMETHD,PREC,&
       & IPART,AFMT,ISTOPC,ITMAX,ITRACE,novr,iprec,EPS)

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
!!$        B_COL_GLOB(I) = REAL(I)*2.0/REAL(M_PROBLEM)        
        B_COL_GLOB(I) = 1.D0
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
         & DESC_A,B_COL_GLOB,B_COL,fmt=afmt)
  ELSE  IF (IPART.EQ.1) THEN 
    CALL BLACS_BARRIER(ICTXT,'A')
    WRITE(6,*) 'Partition type: BLK2'
    CALL MATDIST(AUX_A, A, PART_BLK2, ICTXT, &
         & DESC_A,B_COL_GLOB,B_COL,fmt=afmt)
  ELSE IF (IPART.EQ.2) THEN 
    WRITE(0,*) 'Partition type: GRAPH'
    IF (AMROOT) THEN 
!!$      WRITE(0,*) 'Call BUILD',size(aux_a%ia1),size(aux_a%ia2),np
      WRITE(0,*) 'Build type: GRAPH ',aux_a%fida,&
           &aux_a%m
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
         & DESC_A,B_COL_GLOB,B_COL,fmt=afmt)
  END IF
  
  CALL F90_PSDSALL(M_PROBLEM,X_COL,IERRV,DESC_A)
  X_COL(:) =0.0
  CALL F90_PSDSASB(X_COL,IERRV,DESC_A)
  CALL F90_PSDSALL(M_PROBLEM,R_COL,IERRV,DESC_A)
  R_COL(:) =0.0
  CALL F90_PSDSASB(R_COL,IERRV,DESC_A)
  T2 = MPI_WTIME() - T1
  
  
  CALL DGAMX2D(ICTXT, 'A', ' ', IONE, IONE, T2, IONE,&
       & T1, T1, -1, -1, -1)
  
  IF (AMROOT) THEN
     WRITE(6,*) 'Time to Read and Partition Matrix : ',T2
  END IF
  
  !
  !  Prepare the preconditioning matrix. Note the availability
  !  of optional parameters
  !

  IF (AMROOT) WRITE(6,*) 'Preconditioner : "',PREC(1:6),'"  ',iprec
  
  ! Zero initial guess.
  select case(iprec)
  case(noprec_)
    call psb_precset(pre,'noprec')
  case(diagsc_)             
    call psb_precset(pre,'diagsc')
  case(ilu_)             
    call psb_precset(pre,'ilu')
  case(asm_)             
    call psb_precset(pre,'asm',iv=(/novr,halo_,sum_/))
  case(ash_)             
    call psb_precset(pre,'asm',iv=(/novr,nohalo_,sum_/))
  case(ras_)             
    call psb_precset(pre,'asm',iv=(/novr,halo_,none_/))
  case(rash_)             
    call psb_precset(pre,'asm',iv=(/novr,nohalo_,none_/))
  end select


  T1 = MPI_WTIME()

  CALL psb_precbld(A,PRE,DESC_A,INFO)!,'F')
  TPREC = MPI_WTIME()-T1
  
  
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
  T1 = MPI_WTIME()
  IF (CMETHD.EQ.'BICGSTAB') Then
    CALL  F90_BICGSTAB(A,PRE,B_COL,X_COL,EPS,DESC_A,& 
       & ITMAX,ITER,ERR,IERR,ITRACE,istop=istopc)     
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
  call f90_psaxpby(1.d0,b_col,0.d0,r_col,desc_A)
  call f90_psspmm(-1.d0,a,x_col,1.d0,r_col,desc_a)
  call f90_nrm2(resmx,r_col,desc_a)
!!$  where (b_col /= 0.d0)
!!$    r_col = r_col/b_col
!!$  end where
  call f90_amax(resmxp,r_col,desc_a)

!!$  ITER=IPARM(5)
!!$  ERR = RPARM(2)
  if (amroot) then 
    call prec_descr(6,pre)
    call csprt(60+myprow,a)
!!$    write(6,*) 'Number of iterations : ',iter
!!$    WRITE(6,*) 'Error on exit        : ',ERR
    write(6,*) 'Matrix: ',mtrx_file
    write(6,*) 'Computed solution on ',NPROW,' processors.'
    write(6,*) 'Iterations to convergence: ',iter
    write(6,*) 'Error indicator on exit:',err
    write(6,*) 'Time to Buil Prec.   : ',TPREC
    write(6,*) 'Time to Solve Matrix : ',T2
    WRITE(6,*) 'Time per iteration   : ',T2/(ITER)
    write(6,*) 'Total Time           : ',T2+TPREC
    write(6,*) 'Residual norm 2   = ',resmx
    write(6,*) 'Residual norm inf = ',resmxp
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

  
  CALL F90_PSDSFREE(B_COL, DESC_A)
  CALL F90_PSDSFREE(X_COL, DESC_A)
  CALL F90_PSSPFREE(A, DESC_A)
  CALL psb_precfree(PRE,info)
  CALL F90_PSDSCFREE(DESC_A,info)
  CALL BLACS_GRIDEXIT(ICTXT)
  CALL BLACS_EXIT(0)
  
END PROGRAM DF_SAMPLE
  




