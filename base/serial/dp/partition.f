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
      SUBROUTINE PARTITION(M, WORK, IA, N_BLOCK)
      use psb_const_mod
      use psb_spmat_type
      IMPLICIT NONE


C     ...Scalar arguments...

      INTEGER M, N_BLOCK

C     ...Array arguments...
      
      INTEGER IA(3,*), WORK(*)

C     ...Local scalars...

      INTEGER I, NNZ_ROW, N_ROWS, N_ROWS_EQ, BLOCK

      I = 1      
      N_ROWS = 0
      BLOCK = 2

      WORK(M+1) = 1

      IA(1,1) = 1

      DO WHILE(.TRUE.) 
        IF (N_ROWS.GT.PSB_MAXJDROWS_) THEN
          IA(1,BLOCK) = IA(1,BLOCK-1)+PSB_MAXJDROWS_
          N_ROWS = N_ROWS-PSB_MAXJDROWS_
          BLOCK = BLOCK+1
        ELSE IF (N_ROWS.GE.PSB_MINJDROWS_) THEN
          IA(1,BLOCK) = IA(1,BLOCK-1)+N_ROWS
          N_ROWS = 0
          BLOCK = BLOCK+1
        ELSE IF (I.LE.M) THEN
          N_ROWS_EQ = 0
          NNZ_ROW = -WORK(I)
          DO WHILE (NNZ_ROW.EQ.-WORK(I))
            N_ROWS_EQ = N_ROWS_EQ+1
            I=I+1
          ENDDO
          N_ROWS = N_ROWS + N_ROWS_EQ
        ELSE IF (N_ROWS.NE.0) THEN ! (I.GT.M)
          IA(1,BLOCK) = IA(1,BLOCK-1)+N_ROWS
          BLOCK = BLOCK+1
          GOTO 998
        ELSE
          GOTO 998
        ENDIF
      ENDDO
 998  CONTINUE

      N_BLOCK = BLOCK - 2
      
      if (ia(1,n_block+1)-1 .ne. m) then 
        write(0,*) 'PARTITION: Something wrong',m,
     +    n_block,ia(1,n_block+1),ia(1,n_block)
      endif
      END

