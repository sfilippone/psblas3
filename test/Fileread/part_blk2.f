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
C
C User defined function corresponding to an HPF  BLOCK partition 
C
      SUBROUTINE PART_BLK2(IDX,N,NP,PV,NV)
      
      IMPLICIT NONE
      
      INTEGER  IDX, N, NP
      INTEGER  NV
      INTEGER  PV(*)
      DOUBLE PRECISION   DDIFF
      INTEGER  IB1, IB2, IP, NB, NB1, NNB1
       
      NV = 1
      NB   = N/NP
      NB1  = NB+1 
      NNB1 = MOD(N,NP) 
      IF (IDX .LE. (NNB1*NB1))  THEN 
        PV(1) = (IDX - 1) / NB1
      ELSE 
        IF (NB > 0) THEN  
          IP    = ( (IDX-NNB1*NB1) - 1)/NB
          PV(1) = NNB1 + IP
        ELSE
          write(0,*) 'Impossible ??? '
          PV(1) = NNB1
        ENDIF
      ENDIF
            
      RETURN
      END 
      

      SUBROUTINE BLD_PARTBLK2(N,NP,IVG)

      INTEGER  N, IVG(*),NP
      INTEGER  IB1, IB2, IP, NB, NB1, NNB1, I

      NB   = N/NP
      NB1  = NB+1 
      NNB1 = MOD(N,NP) 
      DO I=1,N
         IF (I .LE. (NNB1*NB1))  THEN 
            IVG(I) = (I - 1) / NB1
         ELSE 
            IF (NB > 0) THEN  
               IP    = ( (I-NNB1*NB1) - 1)/NB
               IVG(I) = NNB1 + IP
            ELSE
               write(0,*) 'Impossible ??? '
               IVG(I) = NNB1
            ENDIF
         ENDIF
      ENDDO
      
      END
