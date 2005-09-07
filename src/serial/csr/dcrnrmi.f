C     ... Compute Infinity norm for sparse matrix in CSR Format ...
      DOUBLE PRECISION FUNCTION DCRNRMI(TRANS,M,N,DESCRA,A,IA1,IA2,
     +  INFOA,IERROR)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           M,N, IERROR
      CHARACTER         TRANS
C     .. Array Arguments ..
      INTEGER           IA1(*),IA2(*),INFOA(*)
      CHARACTER         DESCRA*11
      DOUBLE PRECISION  A(*)
C     .. Local scalars ..
      INTEGER I, J
      DOUBLE PRECISION NRMI, SUM

      IERROR=0
      NRMI = 0.0
      DO I = 1, M
        SUM = 0.0
        DO J = IA2(I), IA2(I+1)-1
          SUM = SUM + ABS(A(J))
        ENDDO

        IF (SUM.GT.NRMI) THEN
          NRMI = SUM
        ENDIF
      ENDDO

      DCRNRMI = NRMI
      END
