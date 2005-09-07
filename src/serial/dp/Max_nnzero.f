      INTEGER FUNCTION MAX_NNZERO(M, IA2)
      IMPLICIT NONE

      INTEGER M
      INTEGER IA2(*)

      INTEGER I,MAX_NZ

      MAX_NZ = 0

      DO I = 1, M
         MAX_NZ = MAX(MAX_NZ,IA2(I+1)- IA2(I))
      ENDDO

      MAX_NNZERO = MAX_NZ
      END
         
