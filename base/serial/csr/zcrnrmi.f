      DOUBLE PRECISION FUNCTION ZCRNRMI(TRANS,M,N,DESCRA,A,IA1,IA2,      
     &   INFOA,IERROR)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           M,N, IERROR
      CHARACTER         TRANS
C     .. Array Arguments ..
      INTEGER    IA1(*),IA2(*),INFOA(*)
      CHARACTER  DESCRA*11
      COMPLEX*16 A(*)
C     .. Local scalars ..
      INTEGER I, J
      DOUBLE PRECISION NRMI, SUM

      NRMI = 0.D0
      DO I = 1, M
         SUM = 0.D0
         DO J = IA2(I), IA2(I+1)-1
C
C            .. definition coerent abs
C
C            SUM = SUM + ABS(A(J))
C
C            .. essl_way abs
C
             SUM = SUM + ABS(DBLE(A(J))) + ABS(AIMAG(A(J)))
         ENDDO         
         NRMI = MAX(NRMI, SUM)
      ENDDO

      ZCRNRMI = NRMI
      END
