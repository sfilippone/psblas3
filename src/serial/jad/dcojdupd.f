      SUBROUTINE DCOJDUPD(M, N, DESCRA, A, IA1,
     +  IA2, INFOA, IA, JA, DESCRH, H, IH1, IH2,
     +  INFOH, IH, JH, FLAG, GLOB_TO_LOC,
     +  IWORK, LIWORK, IERROR)
C
C     .. Matrix A to be updated is required to be stored with
C     .. column indices belonging to the same row ordered.
C     .. Block H to be inserted don't need to be stored in such a way.
C
C     Flag = 0: put elements to 0.0D0;
C     Flag = 1: replace elements with new value;
C     Flag = 2: sum block value to elements;
C
      IMPLICIT NONE
      include 'sparker.fh'
C     .. Scalar Arguments ..
      INTEGER           IA, JA, IH, JH, M, N,
     +  IERROR, FLAG, LIWORK 
C     .. Array Arguments ..
      INTEGER           IA1(*),IA2(*),IH1(*),IH2(*),
     +  INFOA(*),INFOH(*),IWORK(*),
     +  GLOB_TO_LOC(*)
      CHARACTER         DESCRA*11,DESCRH*11
      DOUBLE PRECISION  A(*),H(*)
C     .. Local scalars ..
      INTEGER           J, NNZ, IP1, NNZI
C     .. Local arrays ..
      IERROR = 0
      IF (IBITS(INFOA(UPD_),2,1).EQ.1) THEN 
C
C     Smart update capability
C       
        IP1 = INFOA(UPD_PNT_)
        NNZ = IA1(IP1+NNZ_)
        NNZI = INFOH(1) 
        DO J = 1, NNZI
          NNZ = NNZ + 1 
          A(NNZ) = H(J)
        ENDDO
        IA1(IP1+NNZ_) = NNZ
      ELSE 
        IERROR = 2
      ENDIF
 9999 CONTINUE 
      RETURN
      END


