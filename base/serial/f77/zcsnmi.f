C     DOUBLE PRECISION FUNCTION ZCSNMI(TRANS,M,N,FIDA,DESCRA,A,IA1,IA2,      &
C    &                 INFOA,IERROR)
C     Purpose
C     =======
C
C     Computing matrix infinity norm
C                 nrmi <-- ||A||infty      or
C                 nrmi <-- ||At||infty
C
C     Parameters
C     ==========
C
C     TRANS    - CHARACTER*1
C             On entry TRANS specifies if the routine operates with matrix A
C             or with the transpose of A as follows:
C                TRANS = 'N'         ->  use matrix A
C                TRANS = 'T' or 'C'  ->  use A' (transpose of matrix A)
C             Unchanged on exit.
C
C     M        - INTEGER
C             On entry: number of rows of matrix A (A')
C             Unchanged on exit.
C
C     N        - INTEGER
C             On entry: number of columns of matrix A.
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
C     A        - DOUBLE PRECISION array of DIMENSION (*)
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
C     IERROR   - INTEGER
C             On exit IERROR contains the value of error flag as follows:
C             IERROR = 0   no error
C             IERROR > 0   warning
C             IERROR < 0   fatal error
C
C     Notes
C     =====
C
      DOUBLE PRECISION FUNCTION ZCSNMI(TRANS,M,N,FIDA,DESCRA,A,IA1,IA2,
     &                 INFOA,IERROR)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           M,N, IERROR
      CHARACTER         TRANS
C     .. Array Arguments ..
      INTEGER           IA1(*),IA2(*),INFOA(*)
      CHARACTER         DESCRA*11, FIDA*5
      COMPLEX*16        A(*)
C     .. Local Array..
      INTEGER           INT_VAL(5), ERR_ACT
      CHARACTER*30      NAME, STRINGS(2)
C     .. External Subroutines ..
      DOUBLE PRECISION  ZCRNRMI, ZCOONRMI
      EXTERNAL          ZCRNRMI, ZCOONRMI
C     .. Executable Statements ..
C
C     Check for argument errors
C
      IERROR = 0
      NAME = 'ZCSNMI\0'
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

      IF     (M.LT.0) THEN
         IERROR = 10
         INT_VAL(1) = 2
         INT_VAL(2) = M
      ELSE IF (N.LT.0) THEN
         IERROR = 10
         INT_VAL(1) = 3
         INT_VAL(2) = N
      ELSE IF (TRANS.NE.'T' .AND. TRANS.NE.'N' .AND. TRANS.NE.'C') THEN
         IERROR = 40
         INT_VAL(1) = 1
         STRINGS(1) = TRANS//'\0'
      ENDIF

C
C     Error handling
C
      IF(IERROR.NE.0) THEN
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      ENDIF

C
C     Check for M, N, K
C
      IF(M.LE.0 .OR. N.LE.0) THEN
         ZCSNMI = 0.D0
         GOTO 9999
      ENDIF

C     ... Compute infinity norm for matrix A ...
      IF (FIDA(1:3).EQ.'CSR') THEN
         ZCSNMI = ZCRNRMI(TRANS,M,N,DESCRA,A,IA1,IA2,
     +      INFOA,IERROR)
c$$$      ELSE IF (FIDA(1:3).EQ.'JAD') THEN
c$$$         ZCSNMI = ZJDNRMI(TRANS,M,N,DESCRA,A,IA1,IA2,
c$$$     +      INFOA,IERROR)
      ELSE IF (FIDA(1:3).EQ.'COO') THEN
         ZCSNMI = ZCOONRMI(TRANS,M,N,DESCRA,A,IA1,IA2,
     +      INFOA,IERROR)
      ELSE
C
C     This data structure not yet considered
C
         IERROR = 3010
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
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
