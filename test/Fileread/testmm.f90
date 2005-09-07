PROGRAM TESTMM
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


  INTEGER, PARAMETER    :: IZERO=0, IONE=1
  CHARACTER, PARAMETER  :: ORDER='R'
  REAL(KIND(1.D0)), POINTER, SAVE :: B_COL(:), X_COL(:), R_COL(:), &
       & B_COL_GLOB(:), X_COL_GLOB(:), R_COL_GLOB(:), B_GLOB(:,:), &
       &Z(:), Q(:),Z1(:), XM(:,:), YM1(:,:), YMM(:,:)
  INTEGER              :: IARGC, CHECK_DESCR, CONVERT_DESCR
  Real(Kind(1.d0)), Parameter :: Dzero = 0.d0, One = 1.d0
  Real(Kind(1.d0)) :: MPI_WTIME, T1, T2, TPREC, R_AMAX, B_AMAX,bb(1,1),&
       &lambda,scale,resmx,resmxp, tlpm1, tlpmm, tt, tnc1
  integer :: nrhs, nrow, nx1, nx2, n_row, dim,iread, itry
  integer, parameter :: ntry=16
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
       & METHD, ISTOPC, ML, NCOLS, nc
  integer, pointer :: ierrv(:)
  character(len=5) :: afmt
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
       & IPART,AFMT,NCOLS,ITMAX,ITRACE,PRE%N_OVR,PRE%PREC,EPS)

  CALL BLACS_BARRIER(ICTXT,'A')
  T1 = MPI_WTIME()  
  ! Read the input matrix to be processed and (possibly) the RHS 
  NRHS = 1
  NPROC = NPROW 

  IF (AMROOT) THEN
    NULLIFY(AUX_B)
    CALL READMAT(MTRX_FILE, AUX_A, ICTXT)
    WRITE(0,*) 'From readmat:  ',aux_a%fida,aux_a%m,':',&
         &aux_a%ia2(aux_a%m+1)-1,':',aux_a%ia1(1:10)
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
        WRITE(0,*) 'Memory allocation failure in TESTMM'
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
      WRITE(0,*) 'Memory allocation failure in TESTMM'
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
         & DESC_A,B_COL_GLOB,B_COL,FMT=AFMT)
  ELSE IF (IPART.EQ.2) THEN 
    IF (AMROOT) THEN 
!!$      WRITE(0,*) 'Call BUILD',size(aux_a%ia1),size(aux_a%ia2),np
      WRITE(0,*) 'Build type: GRAPH ',aux_a%fida,&
           &aux_a%m
      CALL BUILD_GRPPART(AUX_A%M,AUX_A%FIDA,AUX_A%IA1,AUX_A%IA2,NP)
    ENDIF

    CALL DISTR_GRPPART(0,0,ICTXT)

    CALL MATDIST(AUX_A, A, PART_GRAPH, ICTXT, &
         & DESC_A,B_COL_GLOB,B_COL,FMT=AFMT)
  ELSE 
    WRITE(6,*) 'Partition type: BLOCK'
    CALL MATDIST(AUX_A, A, PART_BLOCK, ICTXT, &
         & DESC_A,B_COL_GLOB,B_COL,FMT=AFMT)
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

  T1 = MPI_WTIME()

  CALL PRECONDITIONER(A,PRE,DESC_A,INFO)!,'F')
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
    write(6,*) 'methd iprec          : ',pre%prec
    write(6,*) 'Number of iterations : ',iter
    write(6,*) 'Time to Solve Matrix : ',t2
    write(6,*) 'Time per iteration   : ',t2/(iter)
    write(6,*) 'Error on exit        : ',err
  end if


  do nc=1, ncols
    call f90_psdsall(m_problem,nc,xm,ierrv,desc_a)
    call f90_psdsall(m_problem,nc,ym1,ierrv,desc_a)
    call f90_psdsall(m_problem,nc,ymm,ierrv,desc_a)
    ym1(:,:) = 0.d0
    ymm(:,:) = 0.d0
    do j=1,nc
      xm(:,j) = j
    end do
    call f90_psdsasb(xm,ierrv,desc_a)  
    call f90_psdsasb(ym1,ierrv,desc_a)  
    call f90_psdsasb(ymm,ierrv,desc_a)  

    tlpm1 = 1.d200
    do itry=1,ntry
      call blacs_barrier(ictxt,'All')
      T1 = MPI_WTIME()
      do i=1, nc
        call f90_psspmm(1.d0,a,xm(:,i),1.d0,ym1(:,i),desc_a)      
      enddo
      t2 = mpi_wtime()-t1
      call dgamx2d(ictxt,'a',' ',ione, ione,t2,ione,t1,t1,-1,-1,-1)
      tlpm1 = min(tlpm1,t2)
!!$      write(0,*) 'Timing for loop ',nc,itry,t2
    enddo

    tlpmm = 1.d200
    do itry=1,ntry
      call blacs_barrier(ictxt,'All')
      T1 = MPI_WTIME()
      call f90_psspmm(1.d0,a,xm,1.d0,ymm,desc_a)      
      t2 = mpi_wtime()-t1
      call dgamx2d(ictxt,'a',' ',ione, ione,t2,ione,t1,t1,-1,-1,-1)
      tlpmm = min(tlpmm,t2)
!!$      write(0,*) 'Timing for mm   ',nc,itry,t2
    enddo

!!$    ymm = ymm - ym1
    if (nc == 1) tnc1 = tlpm1 
    if (amroot) then
!!$      write(6,*) 'Size         : ',ncols,size(xm,2),size(ym1,2)
!!$      write(6,*) 'Loop         : ',tlpm1   
!!$      write(6,*) 'Single call  : ',tlpmm
      write(6,997) nc, tlpm1, tlpmm, tlpm1/(nc*tnc1),tlpmm/(nc*tnc1)
997 format(i8,4(2x,g16.10))
    end if

!!$    write(6,*) 'maxdiff      : ',maxval(ymm)

    call f90_psdsfree(xm,desc_a)
    call f90_psdsfree(ymm,desc_a)
    call f90_psdsfree(ym1,desc_a)
  end do

  if (.false.) then   
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
  endif

  CALL F90_PSDSFREE(B_COL, DESC_A)
  CALL F90_PSDSFREE(X_COL, DESC_A)
  CALL F90_PSSPFREE(A, DESC_A)
  CALL F90_PSPRECFREE(PRE,info)
  CALL F90_PSDSCFREE(DESC_A,info)
  CALL BLACS_GRIDEXIT(ICTXT)
  CALL BLACS_EXIT(0)


END PROGRAM TESTMM
  




