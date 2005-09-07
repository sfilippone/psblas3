C     SUBROUTINE DSWSM(TRANS,M,N,ALPHA,UNITD,D,FIDT,DESCRT,T,IT1,IT2,
C                      INFOT,B,LDB,BETA,C,LDC,WORK,LWORK,IERROR)
C
C     Purpose
C     =======
C
C     Solving triangular systems of equations with multiple right-hand sides
C                 C <-- ALPHA D T-1 B + BETA C   or
C                 C <-- ALPHA D T-t B + BETA C   or
C                 C <-- ALPHA T-1 D B + BETA C   or
C                 C <-- ALPHA T-t D B + BETA C
C     Actual computing performed by sparse Toolkit kernels.
C     This routine selects the proper kernel for each
C     data structure.
C
C     Parameters
C     ==========
C
C     TRANS    - CHARACTER*1
C             On entry TRANS specifies whether the routine operates with
C             matrix T or with the transpose of T as follows:
C                TRANS = 'N'         ->  use matrix T
C                TRANS = 'T' or 'C'  ->  use T' (transpose of matrix T)
C             Unchanged on exit.
C
C     M        - INTEGER
C             On entry: number of rows and columns of matrix T
C             and number of rows of matrices B and C.
C             Unchanged on exit.
C
C     N        - INTEGER
C             On entry: number of columns of matrices B and C
C             (number of right-hand sides).
C             Unchanged on exit.
C
C     ALPHA    - DOUBLE PRECISION
C             On entry: multiplicative constant.
C             Unchanged on exit.
C
C     UNITD    - CHARACTER*1
C             On entry UNITD specifies whether the diagonal matrix is unit
C             or whether row or column scaling has to be performed, as follows:
C                UNITD = 'U'         ->  unit matrix (no scaling)
C                UNITD = 'L'         ->  scale on the left (row scaling)
C                UNITD = 'R'         ->  scale on the right (column scaling)
C                UNITD = 'B'         ->  scale on the right and on the left
C                                             with D^1/2
C             Unchanged on exit.
C
C     D        - DOUBLE PRECISION array of dimension (M)
C             On entry D specifies the main diagonal of the matrix used
C             for scaling.
C             Unchanged on exit.
C
C     FIDT     - CHARACTER*5
C             On entry FIDT defines the format of the input sparse matrix.
C             Unchanged on exit.
C
C     DESCRT   - CHARACTER*1 array of DIMENSION (9)
C             On entry DESCRT describes the characteristics of the input
C             sparse matrix.
C             Unchanged on exit.
C
C
C     T        - DOUBLE PRECISION array of DIMENSION (*)
C             On entry T specifies the values of the input sparse
C             matrix.
C             Unchanged on exit.
C
C     IT1      - INTEGER array of dimension (*)
C             On entry IT1 holds integer information on input sparse
C             matrix.  Actual information will depend on data format used.
C             Unchanged on exit.
C
C     IT2      - INTEGER array of dimension (*)
C             On entry IT2 holds integer information on input sparse
C             matrix.  Actual information will depend on data format used.
C             Unchanged on exit.
C
C     INFOT     - INTEGER array of dimension (10)
C             On entry can hold auxiliary information on input matrices
C             formats or environment of subsequent calls.
C             Might be changed on exit.
C
C     B        - DOUBLE PRECISION array of dimension (LDB,*)
C             On entry: matrix of right-hand sides
C             Unchanged on exit.
C
C     LDB      - INTEGER
C             On entry: leading dimension of B.
C             Unchanged on exit.
C
C     BETA     - DOUBLE PRECISION
C             On entry: multiplicative constant.
C             Unchanged on exit.
C
C     C        - DOUBLE PRECISION array of dimension (LDC,*)
C             On exit: solutions of triangular systems
C
C     LDC      - INTEGER
C             On entry: leading dimension of C.
C             Unchanged on exit.
C
C     WORK     - DOUBLE PRECISION array of dimension (LWORK)
C             On entry: work area.
C             On exit INT(WORK(1)) contains the minimum value
C             for LWORK satisfying DSWSM memory requirements.
C
C     LWORK    - INTEGER
C             On entry LWORK specifies the dimension of WORK
C             Unchanged on exit.
C
C     IERROR   - INTEGER
C             On exit IERROR contains the value of error flag as follows:
C             IERROR = 0   no error
C             IERROR > 0   warning
C             IERROR < 0   fatal error
C
C     Note
C     ====
C     All checks on argument are performed in the calling routine.
C
C
      SUBROUTINE DSWSM(TRANS,M,N,ALPHA,UNITD,D,FIDT,DESCRT,T,IT1,IT2, 
     &                 INFOT,B,LDB,BETA,C,LDC,WORK,LWORK,IERROR)
C     .. Scalar Arguments ..
      INTEGER           M, N, LDB, LDC, LWORK, IERROR
      CHARACTER         UNITD, TRANS
      DOUBLE PRECISION  ALPHA, BETA
C     .. Array Arguments ..
      INTEGER           IT1(*), IT2(*), INFOT(*)
      CHARACTER         DESCRT*11, FIDT*5
      DOUBLE PRECISION  T(*), B(LDB,*), C(LDC,*), D(*), WORK(*)
C     .. Local Scalars ..
      INTEGER           ONE
C     .. Parameters ..
      PARAMETER        (ONE=1)
C     .. External Subroutines ..
      EXTERNAL         DCSRSM, DCOPY
      LOGICAL           DEBUG
      PARAMETER         (DEBUG=.FALSE.)

      CHARACTER*20      NAME

C     .. Executable Statements ..

      NAME = 'DSWSM\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

C
C     Check for identity matrix
C
      IF(DESCRT(1:1).EQ.'D' .AND. DESCRT(3:3).EQ.'U') THEN
         CALL DCOPY(M,B,ONE,C,ONE)
         GOTO 9998
      ENDIF

      if (debug) write(*,*) 'DSWSM ',m,n,ierror,' ',unitd
C
C     Switching on FIDT: proper sparse BLAS routine is selected
C     according to data structure
C
      IF (FIDT(1:3).EQ.'CSR') THEN
C
C        T, IT1, IT2 --->  AR,   JA,   IA
C                         VAL, INDX, PNTR
C        INFOT(*) not used
C
         CALL  DCSRSM(TRANS,M,N,UNITD,D,ALPHA,DESCRT,T,IT1,          
     &                IT2,B,LDB,BETA,C,LDC,WORK,LWORK)
      ELSE IF (FIDT(1:3).EQ.'JAD') THEN
         
         CALL  DJADSM(TRANS,M,N,D,UNITD,0,ALPHA,DESCRT,T,IT1,IT2,
     +    0,B,LDB,BETA,C,LDC,WORK)
         
      ELSE IF (FIDT(1:3).EQ.'COO') THEN
        
        CALL  DCOOSM(TRANS,M,N,UNITD,D,ALPHA,DESCRT,T,IT1,IT2,INFOT,
     +     B,LDB,BETA,C,LDC,WORK,LWORK,IERROR)
        
      ELSE
C
C     This data structure not yet considered
C
         IERROR = 3010
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999

      END IF

 9998 CONTINUE
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
