      SUBROUTINE ZCOORWS(TRANS,M,N,DESCRA,A,IA1,IA2,
     &  INFOA,ROWSUM,IERROR)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           M,N, IERROR
      CHARACTER         TRANS
C     .. Array Arguments ..
      INTEGER           IA1(*),IA2(*),INFOA(*)
      CHARACTER         DESCRA*11
      COMPLEX*16        A(*), ROWSUM(*)
C     .. Local scalars ..
      INTEGER I, J, NNZ, K
      DOUBLE PRECISION  SUM
      logical psb_lsame
      external psb_lsame

      NNZ = INFOA(1)
      IF (psb_lsame(TRANS,'N')) THEN
        DO I=1, M 
          ROWSUM(I) = 0.0D0
        ENDDO
        I    = 1
        J    = I
        DO WHILE (I.LE.NNZ)
          
          DO WHILE ((IA1(J).EQ.IA1(I)).AND.
     +      (J.LE.NNZ))
            J = J+1
          ENDDO
          
          SUM = 0.0
          DO K = I, J-1
            SUM = SUM + ABS(REAL(A(J))) + ABS(AIMAG(A(J)))
          ENDDO        
          ROWSUM(IA1(I)) = ROWSUM(IA1(I)) + SUM
          I = J 
        ENDDO
        
      ELSE IF (psb_lsame(TRANS,'T').OR.psb_lsame(TRANS,'C')) THEN
        DO J = 1, N
          ROWSUM(J) = 0.0D0
        ENDDO
        DO I = 1, NNZ
          ROWSUM(IA2(I)) = ROWSUM(IA2(I)) +
     +      ABS(REAL(A(J))) + ABS(AIMAG(A(J)))
        ENDDO
      ELSE
        ierror = -1
      ENDIF
      RETURN
      END



