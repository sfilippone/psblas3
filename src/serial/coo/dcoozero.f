      SUBROUTINE DCOOZERO(M,N,DESCRA,A,IA1,IA2,
     +   INFOA,IA,JA,MZ,NZ,IERROR)
C
C     This subroutione performs the operation:
C
C     A(IA : IA + MZ - 1, JA : JA + NZ - 1) = 0
C
C     This isn't accomplished by removing elements 
C     from sparse matrix representation, but assigning them
C     the zero value.
C     Columns are supposed to be ordered 
C     into the same row. This subroutine will
C     not work properly otherwise.
C
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           M,N,IA,JA,MZ,NZ,IERROR
C     .. Array Arguments ..
      INTEGER           IA1(*),IA2(*),INFOA(*)
      CHARACTER         DESCRA*11
      DOUBLE PRECISION  A(*)
C     .. Local scalars ..
      INTEGER I, J, JBEGIN, JEND, AUX, NNZ
      DOUBLE PRECISION 

      IERROR=0
      IF (((JA + NZ - 1) .GT. N) .OR.
     +    ((IA + MZ - 1) .GT. M) .OR.
     +    (IA .LT. 1) .OR. (JA .LT. 1)) THEN
         IERROR = 1
         GOTO 9999
      ENDIF   
      NNZ = INFOA(1)
      I   = 1
      DO WHILE ((IA1(I).LT.IA).AND.(I.LE.NNZ))
        I = I + 1
      ENDDO
      DO WHILE ((IA1(I).LE.(IA+MZ-1)).AND.(I.LE.NNZ))
        IF ((JA.LE.IA2(I)).AND.(IA2(I).LE.(JA+NZ-1))) THEN
          A(I) = 0.0D0
        ENDIF
        I = I + 1 
      ENDDO

      RETURN
      END
