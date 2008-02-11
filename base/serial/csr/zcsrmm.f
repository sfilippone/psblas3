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
C     SUBROUTINE ZCSRMM(TRANSA,M,K,N,ALPHA,DESCRA,AR,
C    *                  JA,IA,B,LDB,BETA,C,LDC,WORK,LWORK)
C
C     Purpose
C     =======
C
C     Computing   C <-- ALPHA A   B + BETA C    or
C                 C <-- ALPHA At  B + BETA C    or
C                 C <-- ALPHA Atc B + BETA C
C     Called by ZSWMM
C     This routine calls kernel for CSR data structure.
C
C     Parameters
C     ==========
C
C     TRANSA   - CHARACTER*1
C             On entry TRANS specifies if the routine operates with matrix A
C             or with the transpose of A as follows:
C                TRANS = 'N'         ->  use matrix A
C                TRANS = 'T'         ->  use A' (transpose of matrix A)
C                TRANS = 'C'         ->  use conjugate transpose of A
C             Unchanged on exit.
C      
C             N.B.: M, K for C matrix
C                   M, N for A matrix
C                   N, K for B matrix
C             In the calling subroutine, ZSWMM, it was:
C                   M, N for C matrix
C                   M, K for A matrix
C                   K, N for B matrix
C             Check the parameters order!
C 
C     M        - INTEGER
C             On entry: number of rows of matrix A (A') and
C                       number of rows of matrix C
C             Unchanged on exit.
C
C     K        - INTEGER
C             On entry: number of columns of matrix B
C             and number of columns of matrix C.
C             Unchanged on exit.
C
C     N        - INTEGER
C             On entry: number of columns of matrix A (A') and
C                       number of rows of matrix B
C             Unchanged on exit.
C
C     ALPHA    - COMPLEX*16
C             On entry: multiplicative constant.
C             Unchanged on exit.
C
C     DESCRA   - CHARACTER*1 array of DIMENSION (9)
C             On entry DESCRA describes the characteristics of the input
C             sparse matrix.
C             Unchanged on exit.
C
C     AR       - COMPLEX*16 array of DIMENSION (*)
C             On entry AR specifies the values of the input sparse
C             matrix.
C             Unchanged on exit.
C
C     JA       - INTEGER array of dimension (*)
C             On entry JA holds integer information on input sparse
C             matrix.  Actual information will depend on data format used.
C             Unchanged on exit.
C
C     IA       - INTEGER array of dimension (*)
C             On entry IA holds integer information on input sparse
C             matrix.  Actual information will depend on data format used.
C             Unchanged on exit.
C
C     INFOA     - INTEGER array of length 10.
C             On entry can hold auxiliary information on input matrices
C             formats or environment of subsequent calls.
C             Might be changed on exit.
C
C     B        - COMPLEX*16 matrix of dimension (LDB,*)
C             On entry: dense matrix.
C             Unchanged on exit.
C
C     LDB      - INTEGER
C             On entry: leading dimension of B
C             Unchanged on exit.
C
C     BETA     - COMPLEX*16
C             On entry: multiplicative constant.
C             Unchanged on exit.
C
C     C        - COMPLEX*16 matrix of dimension (LDC,*)
C             On entry: dense matrix.
C             On exit is updated with the matrix-matrix product.
C
C     LDC      - INTEGER
C             On entry: leading dimension of C
C             Unchanged on exit.
C
C     WORK     - COMPLEX*16 array of dimension (LWORK)
C             On entry: work area.
C             On exit INT(WORK(1)) contains the minimum value
C             for LWORK satisfying ZSWMM memory requirements.
C
C     LWORK    - INTEGER
C             On entry LWORK specifies the dimension of WORK
C             Unchanged on exit.
C
      SUBROUTINE ZCSRMM(TRANSA,M,K,N,ALPHA,DESCRA,AR,
     *                  JA,IA,B,LDB,BETA,C,LDC,WORK,LWORK)
C
C
C     .. Scalar Arguments ..
      COMPLEX*16 ALPHA, BETA
      INTEGER    K, LDB, LDC, M, N, LWORK
      CHARACTER  TRANSA
C     .. Array Arguments ..
      COMPLEX*16 AR(*), B(LDB,*), C(LDC,*),  WORK(*)
      INTEGER    IA(*), JA(*)
      CHARACTER  DESCRA*11
C     .. Local Scalars ..
      INTEGER           I, J
      CHARACTER         DIAG, TRANS

C     .. External Subroutines ..
      EXTERNAL          ZSRMV
C     .. Executable Statements ..
C
C
C      IF (DESCRA(1).EQ.'G') TRANS = TRANSA
C
C     .. Why to loose TRANSA for H, T, A matrices?
C
      TRANS = TRANSA
C        
      IF ((DESCRA(1:1).EQ.'S').AND.(DESCRA(2:2).EQ.'U')) THEN
         IF (TRANSA.EQ.'C') THEN
            TRANS = 'V'
         ELSE
            TRANS = 'U'
         ENDIF
      ENDIF
      IF ((DESCRA(1:1).EQ.'S').AND.(DESCRA(2:2).EQ.'L')) THEN
         IF (TRANSA.EQ.'C') THEN
            TRANS = 'M'
         ELSE
            TRANS = 'L'
         ENDIF
      ENDIF
      
C     .. Diagonal matrix
      IF (DESCRA(1:1).EQ.'D') THEN
C        .. Diagonal matrix with unitary values
         IF (DESCRA(3:3).EQ.'U') THEN
            DO 40 I = 1, K
               DO 20 J = 1, M
                  C(J,I) = BETA*C(J,I) + ALPHA*B(J,I)
   20          CONTINUE
   40       CONTINUE
            RETURN
C        .. Diagonal matrix to be conjugated 
         ELSE IF (TRANSA.EQ.'C') THEN
            DO 80 I = 1, K
               DO 60 J = 1, M
                  C(J,I) = BETA*C(J,I) + ALPHA * 
     +                     CONJG(AR(J)) * B(J,I)
   60          CONTINUE
   80       CONTINUE
            RETURN
C        .. Generic diagonal matrix
         ELSE
            DO 91 I = 1, K
               DO 90 J = 1, M
                  C(J,I) = BETA*C(J,I) + ALPHA * 
     +                     AR(J) * B(J,I)
   90          CONTINUE
   91       CONTINUE
            RETURN
         ENDIF
      END IF
C
      IF (DESCRA(3:3).EQ.'N') DIAG = 'N'
      IF (DESCRA(3:3).EQ.'U') DIAG = 'U'
C
C     C = A*B  or  C=A'*B  or  C=conjug(A')*B
C
C     TRANS =
C        'N': Compute using A as it is.
C        'T': Compute using the transpose of A.
C        'C': Compute using conjugate transpose of A.
C        'U': A is symmetric, stored upper, compute as it is;
C             (transposition makes no sense)
C        'V': A is symmetric, stored upper, to be conjugated.
C        'L': A is symmetric, stored lower, compute as it is.
C        'M': A is symmetric, stored lower, to be conjugated.
C     DIAG =
C        'U': Unitary diagonal.
C        'N': Generic diagonal.

C
      DO 100 I = 1, K
         CALL ZSRMV(TRANS,DIAG,M,N,ALPHA,AR,JA,IA,B(1,I),
     +      BETA,C(1,I),WORK)
  100    CONTINUE
      RETURN
      END

