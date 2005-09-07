      SUBROUTINE DCRCRUPD(M, N, DESCRA, A, IA1,
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
      INTEGER           I, J, XBLCK, XMATR,
     +  NRC, IPH, JPH, JPA, LPA, IRET, LNK, NNZ, IP1
C     .. Local arrays ..
      IERROR = 0
c$$$      write(0,*) 'dcrcrupd ',infoa(upd_),ibits(infoa(upd_),2,1)
      IF (IBITS(INFOA(UPD_),2,1).EQ.1) THEN 
C
C     Smart update capability
C       
        IP1 = INFOA(UPD_PNT_)
        NNZ = IA2(IP1+NNZ_)
        DO I = 1, M
          XBLCK = IH + I - 1
          DO J = IH2(XBLCK),IH2(XBLCK+1) - 1
            NNZ = NNZ + 1 
            A(NNZ) = H(J)
          ENDDO
        ENDDO
        IA2(IP1+NNZ_) = NNZ
      ELSE 
        IF (FLAG.EQ.0) THEN 
          DO I = 1, M
            XBLCK = IH + I - 1
            XMATR = IA + I - 1
            NRC   = IH2(XBLCK+1) - IH2(XBLCK)
            IPH   = IH2(XBLCK)
            IF (LIWORK.LT.2*NRC+2) THEN
              IERROR = 10
              RETURN
            ENDIF
            DO J = 1, NRC
              IWORK(J) =  GLOB_TO_LOC(JA - JH + IH1(IPH+J-1))
            ENDDO 
            CALL MRGSRT(NRC,IWORK(1),IWORK(NRC+1),IRET)
            
            JPA      = IA2(XMATR)
            LPA      = IA2(XMATR+1)
            LNK = IWORK(NRC+1)
            DO J = 1, NRC 
              JPH = IWORK(LNK)           
              DO WHILE ((IA1(JPA).NE.JPH).AND.(JPA.LT.LPA)) 
                JPA = JPA + 1
              ENDDO
              IF (IA1(JPA).EQ.JPH) THEN 
                A(JPA) = 0.0D0
                LNK = IWORK(NRC+1+LNK)
              ELSE
                IERROR = 1
                GOTO 9999
              ENDIF            
            enddo
          ENDDO
        ELSE IF (FLAG.EQ.1) THEN 
          DO I = 1, M
            XBLCK = IH + I - 1
            XMATR = IA + I - 1
            NRC   = IH2(XBLCK+1) - IH2(XBLCK)
            IPH   = IH2(XBLCK)
            IF (LIWORK.LT.2*NRC+2) THEN
              IERROR = 10
              RETURN
            ENDIF
            DO J = 1, NRC
              IWORK(J) =  GLOB_TO_LOC(JA - JH + IH1(IPH+J-1))
            ENDDO 
            CALL MRGSRT(NRC,IWORK(1),IWORK(NRC+1),IRET)
            
            JPA      = IA2(XMATR)
            LPA      = IA2(XMATR+1)
            LNK = IWORK(NRC+1)
            DO J = 1, NRC 
              JPH = IWORK(LNK)           
              DO WHILE((IA1(JPA).NE.JPH).AND.(JPA.LT.LPA)) 
                JPA = JPA + 1
              ENDDO
              IF (IA1(JPA).EQ.JPH) THEN 
                A(JPA) = H(IPH+LNK-1)
                LNK = IWORK(NRC+1+LNK)
              ELSE     
                IERROR = 1
                GOTO 9999
              ENDIF
            ENDDO
          enddo
        ELSE IF (FLAG.EQ.2) THEN 
          DO I = 1, M
            XBLCK = IH + I - 1
            XMATR = IA + I - 1
            NRC   = IH2(XBLCK+1) - IH2(XBLCK)
            IPH   = IH2(XBLCK)
            IF (LIWORK.LT.2*NRC+2) THEN
              IERROR = 10
              RETURN
            ENDIF
            DO J = 1, NRC
              IWORK(J) =  GLOB_TO_LOC(JA - JH + IH1(IPH+J-1))
            ENDDO 
            CALL MRGSRT(NRC,IWORK(1),IWORK(NRC+1),IRET)
            
            JPA      = IA2(XMATR)
            LPA      = IA2(XMATR+1)
            LNK = IWORK(NRC+1)
            DO J = 1, NRC 
              JPH = IWORK(LNK)           
              DO WHILE((IA1(JPA).NE.JPH).AND.(JPA.LT.LPA)) 
                JPA = JPA + 1
              ENDDO
              IF (IA1(JPA).EQ.JPH) THEN 
                A(JPA) = A(JPA) + H(IPH+LNK-1)
                LNK = IWORK(NRC+1+LNK)
              ELSE     
                IERROR = 1
                GOTO 9999
              ENDIF
            ENDDO
          enddo
        ELSE
          IERROR = 2
        ENDIF
      ENDIF
 9999 CONTINUE 
      RETURN
      END




