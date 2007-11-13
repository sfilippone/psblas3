C   SUBROUTINE ZCSMM(TRANS,M,N,K,ALPHA,PL,FIDA,DESCRA,A,IA1,IA2,
C                      INFOA,PR,B,LDB,BETA,C,LDC,WORK,LWORK,IERROR)
C     Purpose
C     =======
C
C     Computing matrix-matrix product
C                 C <-- ALPHA PL A   PR B + BETA C    or
C                 C <-- ALPHA PL At  PR B + BETA C    or
C                 C <-- ALPHA PL Atc PR B + BETA C 
C
C     Parameters
C     ==========
C
C     TRANS    - CHARACTER*1
C             On entry TRANS specifies if the routine operates with matrix A
C             or with the transpose of A as follows:
C                TRANS = 'N'         ->  use matrix A
C                TRANS = 'T'         ->  use A' (transpose of matrix A)
C                TRANS = 'C'         ->  use conjugate transpose of A
C             Unchanged on exit.
C
C     M        - INTEGER
C             On entry: number of rows of matrix A (A') and
C                       number of rows of matrix C
C             Unchanged on exit.
C
C     N        - INTEGER
C             On entry: number of columns of matrix B
C             and number of columns of matrix C.
C             Unchanged on exit.
C
C     K        - INTEGER
C             On entry: number of columns of matrix A (A') and
C                       number of rows of matrix B
C             Unchanged on exit.
C
C     ALPHA    - COMPLEX*16
C             On entry: multiplicative constant.
C             Unchanged on exit.
C
C     PL       - INTEGER array of dimension (M)
C             On entry PL specifies the row permutation of matrix A
C             (PL(1) == 0 if no permutation).
C             Unchanged on exit.
C
C     FIDA     - CHARACTER*5
C             On entry FIDA defines the format of the input sparse matrix.
C             Unchanged on exit.
C
C     DESCRA   - CHARACTER*1 array of DIMENSION (9)
C             On entry DESCRA describes the characteristics of the input
C             sparse matrix.
C             Unchanged on exit.
C
C     A        - COMPLEX*16 array of DIMENSION (*)
C             On entry A specifies the values of the input sparse
C             matrix.
C             Unchanged on exit.
C
C     IA1      - INTEGER array of dimension (*)
C             On entry IA1 holds integer information on input sparse
C             matrix.  Actual information will depend on data format used.
C             Unchanged on exit.
C
C     IA2      - INTEGER array of dimension (*)
C             On entry IA2 holds integer information on input sparse
C             matrix.  Actual information will depend on data format used.
C             Unchanged on exit.
C
C     INFOA     - INTEGER array of length 10.
C             On entry can hold auxiliary information on input matrices
C             formats or environment of subsequent calls.
C             Might be changed on exit.
C
C     PR       - INTEGER array of dimension (K)
C             On entry PR specifies the column permutation of matrix A
C             (PR(1) == 0 if no permutation).
C             Unchanged on exit.
C
C     B        - COMPLEX*16 matrix of dimension (LDB,*)
C             On entry: dense matrix.
C             Unchanged on exit.
C
C     LDB      - INTEGER
C             On entry: leading dimension of B
C             Unchanged on exit.
C
C     BETA     - COMPLEX*16
C             On entry: multiplicative constant.
C             Unchanged on exit.
C
C     C        - COMPLEX*16 matrix of dimension (LDC,*)
C             On entry: dense matrix.
C             On exit is updated with the matrix-matrix product.
C
C     LDC      - INTEGER
C             On entry: leading dimension of C
C             Unchanged on exit.
C
C     WORK     - COMPLEX*16 array of dimension (LWORK)
C             On entry: work area.
C             On exit INT(WORK(1)) contains the minimum value
C             for LWORK satisfying ZCSMM memory requirements.
C
C     LWORK    - INTEGER
C             On entry LWORK specifies the dimension of WORK
C             LWORK should be set as follows:
C                LWORK = (LWORK for DxxxMM) + Pr*K*N + Pl*M*N
C             where Pr íPlù = 1 if right íleftù permutation has to
C             be performed, 0 otherwise.
C             Unchanged on exit.
C
C     IERROR   - INTEGER
C             On exit IERROR contains the value of error flag as follows:
C             IERROR = 0   no error
C             IERROR > 0   warning
C             IERROR < 0   fatal error
C
C     Local Variables
C     ===============
C
C     LWORKM   - INTEGER
C             Minimum work area dimension for ZCSMM
C
C     LWORKB   - INTEGER
C             Work area dimension for matrix B in subroutine DLPUPD
C
C     LWORKC   - INTEGER
C             Work area dimension for matrix C in subroutine DLPUPD
C
C     LWORKS   - INTEGER
C             Work area dimension for subroutine DSWMM
C
C     P        - INTEGER
C             Pointer to work area
C
C     LP       - LOGICAL
C             LP is true if left permutation is required
C
C     RP       - LOGICAL
C             RP is true if right permutation is required
C
C     Notes
C     =====
C       Some tests have shown that it is more efficient to divide the
C     sparse matrix-dense matrix multiplication step and the dense
C     matrix permutation step, and it is more efficient to put
C     together the left permutation and update (C <- xxx + BETA C)
C     steps. So, the sequence of operations is:
C                      Right permutation             ZLPUPD
C                      Matrix-Matrix product         ZSWMM
C                      Left permutation and update   ZLPUPD
C       In order to avoid useless memory transfer, the above scheme is
C     simplified according to whether right and left permutation have to
C     be performed. If left permutation is not required, the update step
C     is performed in the sparse matrix-dense matrix multiplication kernel.
C
C     It is not possible to call this subroutine with LWORK=0 to get     #
C     the minimal value for LWORK. This functionality needs a better     #
C     connection with DxxxMM                                             #
C
C
      SUBROUTINE ZCSMM(TRANS,M,N,K,ALPHA,PL,FIDA,DESCRA,A,IA1,IA2,   
     &                 INFOA,PR,B,LDB,BETA,C,LDC,WORK,LWORK,IERROR)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER    M,N,K,LDB,LDC,LWORK, IERROR
      CHARACTER  TRANS
      COMPLEX*16 ALPHA,BETA
C     .. Array Arguments ..
      INTEGER    IA1(*),IA2(*),INFOA(*),PL(*),PR(*)
      CHARACTER  DESCRA*11, FIDA*5
      COMPLEX*16 A(*),B(LDB,*),C(LDC,*),WORK(*)
C     .. Local Scalars ..
      INTEGER           LWORKM,  LWORKB, LWORKC, LWORKS, P, ERR_ACT
      LOGICAL           LP, RP
C     .. Local Array..
      INTEGER    INT_VAL(5)
      CHARACTER*30 NAME,  STRINGS(2)
C     .. Parameters ..
      COMPLEX*16 ZERO
      INTEGER    IONE
      PARAMETER  (ZERO=(0.D0, 0.D0), IONE=1)
C     .. External Subroutines ..
      EXTERNAL          ZSWMM, ZLPUPD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, IDINT
C     .. Executable Statements ..
C
C     Check for argument errors
C
      NAME = 'ZCSMM\0'
      IERROR=0

      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)
      
      IF     (M.LT.0) THEN
         IERROR = 10
         INT_VAL(1) = 2
         INT_VAL(2) = M
      ELSE IF (K.LT.0) THEN
         IERROR = 10
         INT_VAL(1) = 4
         INT_VAL(2) = K
      ELSE IF (N.LT.0) THEN
         IERROR = 10
         INT_VAL(1) = 3
         INT_VAL(2) = N
      ELSE IF (TRANS.NE.'T' .AND. TRANS.NE.'N' .AND. TRANS.NE.'C') THEN
         IERROR = 40
         INT_VAL(1) = 1
         STRINGS(1) = TRANS//'\0'
      ELSE IF (LDB.LT.K) THEN
         IERROR = 50
         INT_VAL(1) = 15
         INT_VAL(2) = 4
         INT_VAL(3) = LDB
         INT_VAL(4) = K
      ELSE IF (LDC.LT.M) THEN
         IERROR = 50
         INT_VAL(1) = 18
         INT_VAL(2) = 2
         INT_VAL(3) = LDC
         INT_VAL(4) = M
      ENDIF

C
C     Error handling
C
      IF(IERROR.NE.0) THEN
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      ENDIF

C
C     Inizializations
C
      LP = PL(1).NE.0
      RP = PR(1).NE.0
      LWORKB = K*N
      LWORKC = M*N
      LWORKM = 0
      IF (RP) LWORKM = LWORKB
      IF (LP) LWORKM = LWORKM + LWORKC
      IF (LWORK.LT.LWORKM) THEN
         IERROR = 60
         INT_VAL(1) = 20
         INT_VAL(2) = LWORKM
         INT_VAL(3) = LWORK
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      ENDIF
      LWORKS = LWORK - LWORKM

C
C     Check for M, N, K
C
      IF(M.LE.0 .OR. N.LE.0 .OR. K.LE.0) THEN
         GOTO 9999
      ENDIF
C
C     Switching on PR and PL
C
      IF     (LP .AND. RP) THEN
C
C        Both right and left permutation required
C
         P=LWORKB+1
         CALL ZLPUPD(K,N,PR,B,LDB,ZERO,WORK,K)
         CALL ZSWMM(TRANS,M,N,K,ALPHA,FIDA,DESCRA,A,IA1,IA2,INFOA,  
     &           WORK,K,ZERO,WORK(P),M,WORK(P+LWORKC),LWORKS,IERROR)
         LWORKS = IDINT(DBLE(WORK(P+LWORKC)))
         IF(IERROR .NE. 0) THEN
            IF (IERROR.EQ.3010) THEN
               STRINGS(1) = FIDA//'\0'
               CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
               GOTO 9999
            ENDIF
         ENDIF
         CALL ZLPUPD(M,N,PL,WORK(P),M,BETA,C,LDC)
      ELSE IF(.NOT.LP .AND. RP) THEN
C
C        Only right permutation required
C
         P=LWORKB+1
         CALL ZLPUPD(K,N,PR,B,LDB,ZERO,WORK,K)
         CALL ZSWMM(TRANS,M,N,K,ALPHA,FIDA,DESCRA,A,IA1,IA2,INFOA,  
     &           WORK,K,BETA,C,LDC,WORK(P),LWORKS,IERROR)
         LWORKS = IDINT(DBLE(WORK(P)))
         IF(IERROR .NE. 0) THEN
            IF (IERROR.EQ.3010) THEN
               STRINGS(1) = FIDA//'\0'
               CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
               GOTO 9999
            ENDIF
         ENDIF
      ELSE IF(.NOT.RP .AND. LP) THEN
C
C        Only left permutation required
C
         P=LWORKC+1
         CALL ZSWMM(TRANS,M,N,K,ALPHA,FIDA,DESCRA,A,IA1,IA2,INFOA,  
     &           B,LDB,ZERO,WORK,M,WORK(P),LWORKS,IERROR)
         LWORKS = IDINT(DBLE(WORK(P)))
         IF(IERROR .NE. 0) THEN
            IF (IERROR.EQ.3010) THEN
               STRINGS(1) = FIDA//'\0'
               CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
               GOTO 9999
            ENDIF
         ENDIF

         CALL ZLPUPD(M,N,PL,WORK,M,BETA,C,LDC)
      ELSE IF(.NOT.RP .AND. .NOT.LP) THEN
C
C        No permutations required
C
         CALL ZSWMM(TRANS,M,N,K,ALPHA,FIDA,DESCRA,A,IA1,IA2,INFOA,
     &           B,LDB,BETA,C,LDC,WORK,LWORKS,IERROR)
         LWORKS = IDINT(DBLE(WORK(1)))
         IF(IERROR .NE. 0) THEN
            IF (IERROR.EQ.3010) THEN
               STRINGS(1) = FIDA//'\0'
               CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)

               GOTO 9999
            ENDIF
         ENDIF
      ENDIF
C
C     Return minimum workarea dimension
C
      LWORKM = LWORKM + LWORKS
      WORK(1) = DBLE(LWORKM)
      CALL FCPSB_ERRACTIONRESTORE(ERR_ACT)
      RETURN

 9999 CONTINUE
      CALL FCPSB_ERRACTIONRESTORE(ERR_ACT)

      IF ( ERR_ACT .NE. 0 ) THEN 
         CALL FCPSB_SERROR()
         RETURN
      ENDIF


      RETURN
      END


