C
C     Purpose
C     =======
C
C     Performing column permutation of a sparse matrix.
C
C     Parameters
C     ==========
C
C     TRANS    - CHARACTER*1
C             On entry TRANS specifies whether the routine will use
C             matrix P or the transpose of P for the permutation as follows:
C                TRANS = 'N'         ->  permute with matrix P
C                TRANS = 'T' or 'C'  ->  permute the transpose of P
C             Unchanged on exit.
C
C     M        - INTEGER
C             On entry: number of rows of matrix A.
C             Unchanged on exit.
C
C     N        - INTEGER
C             On entry: number of columns of matrix A.
C             Unchanged on exit.
C
C     DESCRA   - CHARACTER*5 array of DIMENSION (10)
C             On entry DESCRA defines the format of the input sparse matrix.
C             Unchanged on exit.
C
C     IA1      - INTEGER array of dimension (*)
C             On entry IA1 holds integer information on input sparse
C             matrix.  Actual information will depend on data format used.
C             On exit contain integer information on permuted matrix.
C
C     IA2      - INTEGER array of dimension (*)
C             On entry IA2 holds integer information on input sparse
C             matrix.  Actual information will depend on data format used.
C             On exit contain integer information on permuted matrix.
C
C     INFOA     - INTEGER array of dimension (10)
C             On entry can hold auxiliary information on input matrices
C             formats or environment of subsequent calls.
C             Might be changed on exit.
C
C     P        - INTEGER array of dimension (M)
C             On entry P specifies the column permutation of matrix A
C             (P(1) == 0 if no permutation).
C             Unchanged on exit.
C
C     WORK     - DOUBLE PRECISION array of dimension (LWORK)
C             On entry: work area.
C             On exit INT(WORK(1)) contains the minimum value
C             for LWORK satisfying DSPRP memory requirements.
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
C     WORK     - DOUBLE PRECISION array of dimension (LWORK)
C             Work area.
C
C     LWORK    - INTEGER
C             On entry LWORK specifies the dimension of WORK.
C             LWORK must be greater than zero.
C             On exit LWORK is the maximum between the initial value and
C             the minimum value satisfying DSPRP memory requirements.
C
C     IERROR    - INTEGER
C             On exit IERROR contains the value of error flag as follows:
C             IERROR = 0      no error
C             IERROR = 4      error on dimension of vector WORK
C             IERROR = 32     unknown flag TRANS
C             IERROR = 64     LWORK  <=  0
C             IERROR = 128    this data structure not yet considered
C                                                                        
C     Notes                                                              
C     =====                                                              
C     It is not possible to call this subroutine with LWORK=0 to get     
C     the minimal value for LWORK. This functionality needs a better     
C     connection with DxxxMM                                             
C
C
      SUBROUTINE DCSRP1(TRANS,M,N,DESCRA,JA,IA,
     +   P,WORK,IWORK,LWORK,IERROR)
      IMPLICIT NONE
      INCLUDE 'sparker.fh'
C     .. Scalar Arguments ..
      INTEGER          LWORK,M, N, IERROR
      CHARACTER        TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION WORK(*)
      INTEGER          JA(*), IA(*), P(*), IWORK(*)
      CHARACTER        DESCRA*11
C     .. Local Scalars ..
      INTEGER          I, J
c     .. Local Arrays ..
      CHARACTER*20       NAME
      INTEGER            INT_VAL(5)
C     .. Intrinsic Functions ..
      INTRINSIC        DBLE
C
C     .. Executable Statements ..
C
      NAME = 'DCSRP1\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

      IF(TRANS.EQ.'N') THEN
         DO 30 I=1,M
            DO 10 J=IA(I),IA(I+1)-1
               JA(J) = P(JA(J))
 10         CONTINUE
 30      CONTINUE
         WORK(1) = 0.D0
      ELSE IF(TRANS.EQ.'T') THEN
C
C        LWORK refers here to INTEGER IWORK (alias for WORK)
C
         IF(LWORK.LT.M) THEN
          IERROR = 60
          INT_VAL(1) = 18
          INT_VAL(2) = NNZ+2
          INT_VAL(3) = LAUX
          CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
          GOTO 9999
         ENDIF
C
C        Transpose permutation matrix
C
         DO 20 I=1,N
            IWORK(P(I)) = I
 20      CONTINUE
C
C        Permute columns
C
         DO 50 I=1,M
            DO 40 J=IA(I),IA(I+1)-1
               JA(J) = IWORK(JA(J))
 40         CONTINUE
 50      CONTINUE
C
C        WORK(1) refers here to a value for a DOUBLE PRECISION WORK
C
         WORK(1) = DBLE((M+1)/DBLEINT_)
      ENDIF

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
