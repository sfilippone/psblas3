C     Subroutine DLPUPD(M,N,PERM,B,LDB,BETA,C,LDC,IERROR)
C     Purpose
C     =======
C
C     Computing   C  <--  PERM B + BETA C
C     where PERM is a permutation matrix.
C
C     Parameters
C     ==========
C
C     M        - INTEGER
C             On entry M specifies the number of rows of matrices
C             B and C.
C             Unchanged on exit.
C
C     N        - INTEGER
C             On entry N specifies the number of columns of matrices
C             B and C.
C             Unchanged on exit.
C
C     PERM     - INTEGER array of dimension (N)
C             On entry PERM specifies the values of a permutation matrix.
C             Unchanged on exit.
C
C     B        - DOUBLE PRECISION matrix of dimension (LDB,*)
C             On entry: dense matrix.
C             Unchanged on exit.
C
C     LDB      - INTEGER
C             On entry LDB holds the value of the leading dimension of B
C             Unchanged on exit.
C
C     BETA     - DOUBLE PRECISION
C             On entry: multiplicative constant.
C             Unchanged on exit.
C
C     C        - DOUBLE PRECISION matrix of dimension (LDC,*)
C             On entry: dense matrix.
C             On exit is updated as shown above.
C
C     LDC      - INTEGER
C             On entry LDC holds the value of the leading dimension of C
C             Unchanged on exit.
C
C     Note
C     ====
C     All checks on argument are performed in the calling routines.
C
C
      SUBROUTINE DLPUPD(M,N,PERM,B,LDB,BETA,C,LDC)
C     .. Scalar Arguments ..
      INTEGER           M, N, LDB, LDC
      DOUBLE PRECISION  BETA
C     .. Array Arguments ..
      INTEGER           PERM(*)
      DOUBLE PRECISION  B(LDB,*), C(LDC,*)
C     .. Local Scalars ..
      INTEGER           I,J
C
C     .. Executable Statements ..
C
C      Switching on BETA
C
      IF     (BETA.NE.0.D0) THEN
C
C        Performing left permutation and update
C
         DO 40 J = 1, N
            DO 30 I = 1, M
               C(I,J) = B(PERM(I),J) + BETA*C(I,J)
   30       CONTINUE
   40    CONTINUE
      ELSE IF(BETA.EQ.0.D0) THEN
C
C        Performing right or left permutation
C
         DO 160 J = 1, N
            DO 150 I = 1, M
               C(I,J) = B(PERM(I),J)
  150       CONTINUE
  160    CONTINUE
      ENDIF
      RETURN
      END
