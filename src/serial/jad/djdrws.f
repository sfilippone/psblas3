C     ... Compute infinity norm for sparse matrix in CSR Format ...
      SUBROUTINE DJDRWS(TRANS,M,N,DESCRA,A,JA,IA,
     +   INFOA,ROWSUM,IERROR)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           M,N, IERROR
      CHARACTER         TRANS
C     .. Array Arguments ..
      INTEGER           JA(*),IA(*),INFOA(*)
      CHARACTER         DESCRA*11
      DOUBLE PRECISION  A(*), ROWSUM(*)
C     .. Local scalars ..
      INTEGER PNG, PIA, PJA
C     .. External routines ..
      DOUBLE PRECISION DJADNR
      EXTERNAL DJADNR

      IERROR = 0
      PNG = IA(1)
      PIA = IA(2)
      PJA = IA(3)

      IF (DESCRA(1:1).EQ.'G') THEN
        CALL DJADRWS(TRANS,M,N,IA(PNG),
     +     A,JA,IA(PJA),IA(PIA),
     +     INFOA,ROWSUM,IERROR)
      ELSE 
        IERROR = 3011        
      ENDIF
      END
