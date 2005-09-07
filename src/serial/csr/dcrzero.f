      SUBROUTINE DCRZERO(M,N,DESCRA,A,IA1,IA2,      
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
C     not work properly if this hypotesis is not 
C     verified.
C
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           M,N,IA,JA,MZ,NZ,IERROR
C     .. Array Arguments ..
      INTEGER           IA1(*),IA2(*),INFOA(*)
      CHARACTER         DESCRA*11
      DOUBLE PRECISION  A(*)
C     .. Local scalars ..
      INTEGER I, J, JBEGIN, JEND, AUX
      DOUBLE PRECISION 

      IERROR=0
      IF (((JA + NZ - 1) .GT. N) .OR.
     +    ((IA + MZ - 1) .GT. M) .OR.
     +    (IA .LT. 1) .OR. (JA .LT. 1)) THEN
         IERROR = 1
         GOTO 9999
      ENDIF      
      DO I = IA, IA + M - 1
C        Scan current line until found first element
         DO JBEGIN = IA2(I), IA2(I + 1) - 1
            IF (IA1(JBEGIN) .GE. JA)   EXIT
         ENDDO
C        If reached last column end not yet
C        encountered proper column, skip this row
         IF ((JBEGIN .EQ. IA2(I + 1) - 1) .AND.
     +       (IA1(JBEGIN) .LT. JA))   CYCLE
C        Now I'm sure there's at least one element
C        to process: scan until found last element
         AUX = JA + N - 1
         DO JEND = JBEGIN, IA2(I + 1) - 1
            IF (IA1(JEND) .GE. AUX)   EXIT
         ENDDO
         IF (IA1(JEND) .GT. AUX)   JEND = JEND - 1      
         DO J = JBEGIN, JEND
            A(J) = 0.0D0
         ENDDO
      ENDDO
 9999 RETURN
      END
