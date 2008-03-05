C
C             Parallel Sparse BLAS  version 2.2
C   (C) Copyright 2006/2007/2008
C                      Salvatore Filippone    University of Rome Tor Vergata
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
      SUBROUTINE DCOORWS(TRANS,M,N,DESCRA,A,IA1,IA2,
     &   INFOA,ROWSUM,IERROR)
      use psb_const_mod
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           M,N, IERROR
      CHARACTER         TRANS
C     .. Array Arguments ..
      INTEGER           IA1(*),IA2(*),INFOA(*)
      CHARACTER         DESCRA*11
      real(psb_dpk_)  A(*), ROWSUM(*)
C     .. Local scalars ..
      INTEGER I, J, NNZ, K
      real(psb_dpk_)  SUM
      logical psb_lsame
      external psb_lsame

      NNZ = INFOA(1)
      IF (psb_lsame(TRANS,'N')) THEN
        DO I=1, M 
          ROWSUM(I) = 0.0D0
        ENDDO
        I    = 1
        J    = I
        DO WHILE (I.LE.NNZ)
          
          DO WHILE ((IA1(J).EQ.IA1(I)).AND.
     +       (J.LE.NNZ))
            J = J+1
          ENDDO
          
          SUM = 0.0
          DO K = I, J-1
            SUM = SUM + ABS(A(K))
          ENDDO        
          ROWSUM(IA1(I)) = ROWSUM(IA1(I)) + SUM
          I = J 
        ENDDO
        
      ELSE IF (psb_lsame(TRANS,'T').OR.psb_lsame(TRANS,'C')) THEN
        DO J = 1, N
          ROWSUM(J) = 0.0D0
        ENDDO
        DO I = 1, NNZ
          ROWSUM(IA2(I)) = ROWSUM(IA2(I)) + ABS(A(I))
        ENDDO
      ELSE
        ierror = -1
      ENDIF
      RETURN
      END



