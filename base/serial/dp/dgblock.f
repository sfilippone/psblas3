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
      SUBROUTINE DGBLOCK(M,IA2,IPERM,IA,N_BLOCK,WORK,LWORK)
      IMPLICIT NONE

C     ...Scalar arguments...

      INTEGER M, N_BLOCK, LWORK

C     ...Array arguments...
      
      INTEGER IA2(*), IPERM(*), IA(3,*), WORK(*)

C     ...Local scalars...

      INTEGER I, SWAP, KK, LP, IRET
C     Compute number of nnzero elements per row

      IPERM(1) = 0

      DO I = 1, M
        WORK(I) = - IA2(I+1) + IA2(I)
      ENDDO

C     Sorting Array work
C ........................

      CALL MSORT_UP(M,WORK,WORK(M+1),IRET)
      IF (IRET.EQ.0) THEN
C     Construct IPERM Vector
        LP = WORK(M+1)
        
        DO I = 1, M
          IPERM(LP) = I
          LP = WORK(M+1+LP)
        ENDDO

        LP = WORK(M+1)                           
        KK = 1                              
        DO WHILE (.NOT.((LP.EQ.0).OR.(KK.GT.M)))
          DO WHILE (LP.LT.KK)
            LP = WORK(M+1+LP)
          ENDDO
C        Swap values of array work         
          SWAP = WORK(KK)                      
          WORK(KK) = WORK(LP)                     
          WORK(LP) = SWAP                      

C        Swap values of index array work(m+1)
          SWAP = WORK(M+1+LP)                     
          WORK(M+1+LP) = WORK(M+1+KK)                     
          WORK(M+1+KK) = LP

          LP    = SWAP                     
          KK = KK+1                         
        ENDDO

      ENDIF
C     Partitioning Matrix in blocks of rows
      CALL PARTITION(M, WORK, IA, N_BLOCK)

      END



