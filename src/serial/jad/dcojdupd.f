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
      INTEGER           J, NNZ, IP1, NNZI
C     .. Local arrays ..
      IERROR = 0
      IF (IBITS(INFOA(PSB_UPD_),2,1).EQ.1) THEN 
C
C     Smart update capability
C       
        IP1 = INFOA(PSB_UPD_PNT_)
        NNZ = IA1(IP1+PSB_NNZ_)
        NNZI = INFOH(1) 
        DO J = 1, NNZI
          NNZ = NNZ + 1 
          A(NNZ) = H(J)
        ENDDO
        IA1(IP1+PSB_NNZ_) = NNZ
      ELSE 
        IERROR = 2
      ENDIF
 9999 CONTINUE 
      RETURN
      END


