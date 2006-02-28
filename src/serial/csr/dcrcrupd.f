C
C             Parallel Sparse BLAS  v2.0
C   (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
C                      Alfredo Buttari        University of Rome Tor Vergata
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions
C are met:
C   1. Redistributions of source code must retain the above copyright
C      notice, this list of conditions and the following disclaimer.
C   2. Redistributions in binary form must reproduce the above copyright
C      notice, this list of conditions, and the following disclaimer in the
C      documentation and/or other materials provided with the distribution.
C   3. The name of the PSBLAS group or the names of its contributors may
C      not be used to endorse or promote products derived from this
C      software without specific written permission.
C
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
C TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
C PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
C BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
C CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
C SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
C INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
C CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
C ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
C POSSIBILITY OF SUCH DAMAGE.
C
C 
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
      include 'psb_const.fh'
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
      IF (IBITS(INFOA(PSB_UPD_),2,1).EQ.1) THEN 
C
C     Smart update capability
C       
        IP1 = INFOA(PSB_UPD_PNT_)
        NNZ = IA2(IP1+PSB_NNZ_)
        DO I = 1, M
          XBLCK = IH + I - 1
          DO J = IH2(XBLCK),IH2(XBLCK+1) - 1
            NNZ = NNZ + 1 
            A(NNZ) = H(J)
          ENDDO
        ENDDO
        IA2(IP1+PSB_NNZ_) = NNZ
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




