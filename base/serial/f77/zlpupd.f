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
C     Subroutine ZLPUPD(M,N,PERM,B,LDB,BETA,C,LDC,IERROR)
C     Purpose
C     =======
C
C     Computing   C  <--  PERM B + BETA C
C     where PERM is a permutation matrix.
C
C     Parameters
C     ==========
C
C     M        - INTEGER
C             On entry M specifies the number of rows of matrices
C             B and C.
C             Unchanged on exit.
C
C     N        - INTEGER
C             On entry N specifies the number of columns of matrices
C             B and C.
C             Unchanged on exit.
C
C     PERM     - INTEGER array of dimension (N)
C             On entry PERM specifies the values of a permutation matrix.
C             Unchanged on exit.
C
C     B        - DOUBLE PRECISION matrix of dimension (LDB,*)
C             On entry: dense matrix.
C             Unchanged on exit.
C
C     LDB      - INTEGER
C             On entry LDB holds the value of the leading dimension of B
C             Unchanged on exit.
C
C     BETA     - DOUBLE PRECISION
C             On entry: multiplicative constant.
C             Unchanged on exit.
C
C     C        - DOUBLE PRECISION matrix of dimension (LDC,*)
C             On entry: dense matrix.
C             On exit is updated as shown above.
C
C     LDC      - INTEGER
C             On entry LDC holds the value of the leading dimension of C
C             Unchanged on exit.
C
C     Note
C     ====
C     All checks on argument are performed in the calling routines.
C
C
      SUBROUTINE ZLPUPD(M,N,PERM,B,LDB,BETA,C,LDC)
C     .. Scalar Arguments ..
      INTEGER           M, N, LDB, LDC
      complex(kind(1.d0)) BETA
C     .. Array Arguments ..
      INTEGER           PERM(*)
      complex(kind(1.d0)) B(LDB,*), C(LDC,*)
C     .. Local Scalars ..
      INTEGER           I,J
C
C     .. Executable Statements ..
C
C      Switching on BETA
C
      IF     (BETA.NE.(0.D0,0.d0)) THEN
C
C        Performing left permutation and update
C
         DO 40 J = 1, N
            DO 30 I = 1, M
               C(I,J) = B(PERM(I),J) + BETA*C(I,J)
   30       CONTINUE
   40    CONTINUE
      ELSE IF(BETA.EQ.(0.D0,0.d0)) THEN
C
C        Performing right or left permutation
C
         DO 160 J = 1, N
            DO 150 I = 1, M
               C(I,J) = B(PERM(I),J)
  150       CONTINUE
  160    CONTINUE
      ENDIF
      RETURN
      END
