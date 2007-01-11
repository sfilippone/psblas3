c
c       What if a wrong DESCRA is passed?
c
c
*
*
      SUBROUTINE ZCOOMM(TRANSA,M,K,N,ALPHA,DESCRA,AR,
     *   IA,JA,INFOA,B,LDB,BETA,C,LDC,WORK,LWORK)
C
C
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA, BETA
      INTEGER           K, LDB, LDC, M, N, LWORK
      CHARACTER         TRANSA
C     .. Array Arguments ..
      COMPLEX*16        AR(*), B(LDB,*), C(LDC,*),  WORK(*)
      INTEGER           IA(*), JA(*),INFOA(*)
      CHARACTER         DESCRA*11
C     .. Local Scalars ..
      INTEGER           I, J
      CHARACTER         DIAG, TRANS


C     .. External Subroutines ..
      EXTERNAL          ZCOOMV
C     .. Executable Statements ..
C
C
C
      IF (DESCRA(1:1).EQ.'G') TRANS = TRANSA
      IF (DESCRA(1:1).EQ.'S' .AND. DESCRA(2:2).EQ.'U') TRANS = 'U'
      IF (DESCRA(1:1).EQ.'S' .AND. DESCRA(2:2).EQ.'L') TRANS = 'L'
c
c     Does DSRMV manage this case too? <????????????????????????????
c
      IF (DESCRA(1:1).EQ.'D') THEN
         IF (DESCRA(3:3).EQ.'U') THEN
            DO 40 I = 1, K
               DO 20 J = 1, M
                  C(J,I) = BETA*C(J,I) + ALPHA*B(J,I)
   20          CONTINUE
   40       CONTINUE
         ELSE
            DO 80 I = 1, K
               DO 60 J = 1, M
                  C(J,I) = BETA*C(J,I) + ALPHA*AR(J)*B(J,I)
   60          CONTINUE
   80       CONTINUE
         END IF
         RETURN
      END IF
C
      IF (DESCRA(3:3).EQ.'N') DIAG = 'N'
      IF (DESCRA(3:3).EQ.'U') DIAG = 'U'
C
C     C = A*B  OR C=A'*B
C
C

      DO 100 I = 1, K
         CALL ZCOOMV(TRANS,DIAG,M,N,ALPHA,AR,IA,JA,INFOA,
     +     B(1,I),BETA,C(1,I),WORK)
 100  CONTINUE
      RETURN
      END
