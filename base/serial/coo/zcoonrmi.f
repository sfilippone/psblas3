C     ... Compute Infinity norm for sparse matrix in CSR Format ...
      DOUBLE PRECISION FUNCTION ZCOONRMI(TRANS,M,N,DESCRA,A,IA1,IA2,
     +   INFOA,IERROR)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           M,N, IERROR
      CHARACTER         TRANS
C     .. Array Arguments ..
      INTEGER           IA1(*),IA2(*),INFOA(*)
      CHARACTER         DESCRA*11
      COMPLEX*16        A(*)
C     .. Local scalars ..
      INTEGER I, J, K, NNZ
      DOUBLE PRECISION NRMI, SUM

      NRMI = 0.0
      NNZ  = INFOA(1)
      I    = 1
      J    = I
      DO WHILE (I.LE.NNZ)

        DO WHILE ((IA1(J).EQ.IA1(I)).AND.
     +     (J.LE.NNZ))
          J = J+1
        ENDDO
        
        SUM = 0.0
        DO K = I, J-1
          SUM = SUM + ABS(DBLE(A(K))) + ABS(AIMAG(A(K)))
        ENDDO        
        IF (SUM.GT.NRMI) THEN
          NRMI = SUM
        ENDIF
        I = J 
      ENDDO

      ZCOONRMI = NRMI
      END
