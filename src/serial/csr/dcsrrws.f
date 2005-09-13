      SUBROUTINE DCSRRWS(TRANS,M,N,DESCRA,A,IA1,IA2,
     &  INFOA,ROWSUM,IERROR)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           M,N, IERROR
      CHARACTER         TRANS
C     .. Array Arguments ..
      INTEGER           IA1(*),IA2(*),INFOA(*)
      CHARACTER         DESCRA*11
      DOUBLE PRECISION  A(*), ROWSUM(*)
C     .. Local scalars ..
      INTEGER I, J

      IF (TRANS.EQ.'N') THEN
        DO I = 1, M
          ROWSUM(I) = 0.0D0
          DO J = IA2(I), IA2(I + 1) - 1
            ROWSUM(I) = ROWSUM(I) + ABS(A(J))
          ENDDO
        ENDDO
      ELSE IF ((TRANS.EQ.'T').OR.(TRANS.EQ.'C')) THEN
        DO J = 1, N
          ROWSUM(J) = 0.0D0
        ENDDO
        DO I = 1, M
          DO J = IA2(I), IA2(I + 1) - 1
            ROWSUM(IA1(J)) = ROWSUM(IA1(J)) + ABS(A(J))
          ENDDO
        ENDDO          
      ENDIF
      END



