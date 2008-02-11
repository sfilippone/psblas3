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
      SUBROUTINE ZCSRRWS(TRANS,M,N,DESCRA,A,IA1,IA2,
     &  INFOA,ROWSUM,IERROR)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           M,N, IERROR
      CHARACTER         TRANS
C     .. Array Arguments ..
      INTEGER           IA1(*),IA2(*),INFOA(*)
      CHARACTER         DESCRA*11
      COMPLEX*16        A(*), ROWSUM(*)
C     .. Local scalars ..
      INTEGER I, J
      DOUBLE PRECISION  SUM

      IF (TRANS.EQ.'N') THEN
        DO I = 1, M
          ROWSUM(I) = (0.0D0, 0.0D0)
          DO J = IA2(I), IA2(I + 1) - 1
            ROWSUM(I) = ROWSUM(I) +
     +        ABS(REAL(A(J))) + ABS(AIMAG(A(J)))
          ENDDO
        ENDDO
      ELSE IF ((TRANS.EQ.'T').OR.(TRANS.EQ.'C')) THEN
        DO J = 1, N
          ROWSUM(J) = (0.0D0, 0.0D0)
        ENDDO
        DO I = 1, M
          DO J = IA2(I), IA2(I + 1) - 1
            ROWSUM(IA1(J)) = ROWSUM(IA1(J)) +
     +        ABS(REAL(A(J))) + ABS(AIMAG(A(J)))
          ENDDO
        ENDDO          
      ENDIF
      END



