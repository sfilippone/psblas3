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
      SUBROUTINE PART_BLOCK(GLOBAL_INDX,N,NP,PV,NV)
      
      IMPLICIT NONE
      
      INTEGER  GLOBAL_INDX, N, NP
      INTEGER  NV
      INTEGER  PV(*)
      INTEGER  DIM_BLOCK
      DOUBLE PRECISION   DDIFF
      INTEGER  IB1, IB2, IPV
      
      double precision PC
      PARAMETER   (PC=0.0D0)

      DIM_BLOCK = (N + NP - 1)/NP
      NV = 1  
      PV(NV) = (GLOBAL_INDX - 1) / DIM_BLOCK
      
      IPV = PV(1)
      IB1 = IPV * DIM_BLOCK + 1
      IB2 = (IPV+1) * DIM_BLOCK
      
      DDIFF = DBLE(ABS(GLOBAL_INDX-IB1))/DBLE(DIM_BLOCK)
      IF (DDIFF .lt. PC/2) THEN
C
C     Overlap at the beginning of a block, with the previous proc
C         
         IF (IPV.gt.0) THEN 
            NV     = NV + 1
            PV(NV) = IPV - 1
         ENDIF
      ENDIF

      DDIFF = DBLE(ABS(GLOBAL_INDX-IB2))/DBLE(DIM_BLOCK)
      IF (DDIFF .lt. PC/2) THEN
C
C     Overlap at the end of a block, with the next proc
C         
         IF (IPV.lt.(NP-1)) THEN 
            NV     = NV + 1
            PV(NV) = IPV + 1
         ENDIF
      ENDIF
      
      RETURN
      END 
      


      SUBROUTINE BLD_PARTBLOCK(N,NP,IVG)
      
      INTEGER N,NP,IVG(*)

      INTEGER  DIM_BLOCK,I


      DIM_BLOCK = (N + NP - 1)/NP
      DO I=1,N
         IVG(I) = (I - 1) / DIM_BLOCK
      ENDDO

      END


