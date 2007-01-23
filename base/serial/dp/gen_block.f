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
      SUBROUTINE GEN_BLOCK(M,NG,IA,AUX)
      use psb_const_mod
      use psb_spmat_type
      IMPLICIT NONE

      INTEGER M, NG
      INTEGER IA(3,*), AUX(*)

      INTEGER BLOCK, I, N_ROWS

      N_ROWS = IA(1,2) - IA(1,1)
      I = 2
      BLOCK = 2
      AUX(1) = 1
      
      DO WHILE(.TRUE.)
        IF (N_ROWS.GT.PSB_MAXJDROWS_) THEN
          AUX(BLOCK) = AUX(BLOCK-1)+PSB_MAXJDROWS_
          N_ROWS = N_ROWS-PSB_MAXJDROWS_
          BLOCK = BLOCK+1
        ELSE IF (N_ROWS.GT.0) THEN
          AUX(BLOCK) = AUX(BLOCK-1)+N_ROWS
          N_ROWS = 0
          BLOCK = BLOCK+1
        ELSE IF (I.LE.NG) THEN
          N_ROWS = IA(1,I+1) - IA(1,I)
          I = I+1
        ELSE
          GOTO 998
        ENDIF
      ENDDO
 998  CONTINUE 

C     ... Copy AUX in IA(1,*)

      NG = BLOCK - 2
      DO I = 1, NG+1
        IA(1,I) = AUX(I)
      ENDDO

      END
