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
      SUBROUTINE CHECK_DIM(M, N, IA, NG, IA2, 
     +   NZ, LARN, LIAN1, LIAN2, IERRV)

      use psb_const_mod
      IMPLICIT NONE

C
C     .. Scalar Arguments ..
      INTEGER M,N,NG,LARN,LIAN1,LIAN2, NZ

C     .. Array Arguments ..
      INTEGER IA(3,*), IA2(*), IERRV(*)

C     Local scalars
      INTEGER NNZ, BLOCK, DIM_BLOCK, LIMIT
      INTEGER MAX_NNZERO, MAX_NZ
      
      EXTERNAL MAX_NNZERO

      MAX_NZ = MAX_NNZERO(M,IA2)
      
      NNZ = NZ
      
c$$$      LIMIT = INT(DIM_BLOCK*PSB_PERCENT_)
      
      DO BLOCK = 1, NG
         DIM_BLOCK = IA(1,BLOCK+1)-IA(1,BLOCK)
         LIMIT = INT(DIM_BLOCK*PSB_PERCENT_)

         NNZ = NNZ+(DIM_BLOCK-LIMIT)*MAX_NZ
      END DO

      IERRV(1)=0
      IERRV(2) = NNZ
      IERRV(3) = NNZ
      IERRV(4) = 6+3*(NG+1)+M+MAX_NZ*NG+1
      IF (6+3*(NG+1)+M+MAX_NZ*NG+1.GT.LIAN2) THEN
         IERRV(1) = 30
c$$$         write(0,*) 'check_dim: error 1',
c$$$     +     6+3*(NG+1)+M+MAX_NZ*NG+1,LIAN2
      ENDIF
      
      IF (NNZ.GT.LIAN1) THEN
c$$$        write(0,*) 'check_dim: error 2',nnz,lian1
         IERRV(1) = 31
      ENDIF
      
      IF (NNZ.GT.LARN) THEN
c$$$        write(0,*) 'check_dim: error 3',nnz,larn
         IERRV(1) = 32
      ENDIF

      RETURN
      END



