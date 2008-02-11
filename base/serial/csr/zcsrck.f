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
C
C     Purpose
C     =======
C
C     Performing checks on sparse matrix.
C
C     Parameters
C     ==========
C
C     TRANS    - CHARACTER*1
C             On entry TRANS specifies whether the routine will use
C             matrix P or the transpose of P for the permutation as follows:
C                TRANS = 'N'         ->  permute with matrix P
C                TRANS = 'T' or 'C'  ->  permute the transpose of P
C             Unchanged on exit.
C
C     M        - INTEGER
C             On entry: number of rows of matrix A.
C             Unchanged on exit.
C
C     N        - INTEGER
C             On entry: number of columns of matrix A.
C             Unchanged on exit.
C
C     DESCRA   - CHARACTER*5 array of DIMENSION (10)
C             On entry DESCRA defines the format of the input sparse matrix.
C             Unchanged on exit.
C
C     A        - DOUBLE PRECISION array of DIMENSION (*)
C             On entry A specifies the values of the input sparse
C             matrix.
C             Unchanged on exit.
C
C     IA1      - INTEGER array of dimension (*)
C             On entry IA1 holds integer information on input sparse
C             matrix.  Actual information will depend on data format used.
C             Unchanged on exit.
C
C     IA2      - INTEGER array of dimension (*)
C             On entry IA2 holds integer information on input sparse
C             matrix.  Actual information will depend on data format used.
C             Unchanged on exit.
C
C     WORK     - DOUBLE PRECISION array of dimension (LWORK)
C             On entry: work area.
C             On exit INT(WORK(1)) contains the minimum value
C             for LWORK satisfying DSPRP memory requirements.
C
C     LWORK    - INTEGER
C             On entry LWORK specifies the dimension of WORK
C             Unchanged on exit.
C
C     IERROR   - INTEGER
C             On exit specify the error code.
C             IERROR = 0 no errors
C             IERROR > 0 error in integrity check


      SUBROUTINE ZCSRCK(TRANS,M,N,DESCRA,A,IA1,IA2,                       
     +   WORK,LWORK,IERROR)
      IMPLICIT NONE                                                     
C     .. Scalar Arguments ..
      INTEGER          LWORK,M, N, IERROR
      CHARACTER        TRANS
C     .. Array Arguments ..
      complex(kind(1.d0)) A(*), WORK(*)
      INTEGER          IA1(*), IA2(*)
      CHARACTER        DESCRA*11
C     .. Local Scalars ..
      INTEGER          I, J, nrow, nind
C
C     .. Executable Statements ..
C

C
C      Check #1: Character descriptor have valid values
C
      IERROR = 0

      IF ((DESCRA(1:1).NE.'G').AND.(DESCRA(1:1).NE.'S').AND.
     &   (DESCRA(1:1).NE.'H').AND.(DESCRA(1:1).NE.'T').AND.
     &   (DESCRA(1:1).NE.'A').AND.(DESCRA(1:1).NE.'D'))  THEN
         IERROR = 11
         GOTO 9999
      END IF
      IF ((DESCRA(2:2).NE.'U').AND.(DESCRA(2:2).NE.'L')) THEN
         IERROR = 12
         GOTO 9999
      END IF
      IF ((DESCRA(3:3).NE.'U').AND.(DESCRA(3:3).NE.'N')) THEN
         IERROR = 13
         GOTO 9999
      END IF
C
C      Check #2: Pointers have non decreasing order
C
      IF (IA2(1).LE.0) THEN
         IERROR = 14
         GOTO 9999
      ENDIF
      
      NROW = 0
      DO 10 I = 1, M
         IF (IA2(I) .GT. IA2(I+1)) THEN
            NROW = NROW + 1
         END IF
 10   CONTINUE
      IF (NROW .GT. 0) THEN
         IERROR = 15
         GOTO 9999
      END IF
C
C      Check #3: Indices are within problem dimension
C
      NIND = 0
      DO 20 I = 1, M
         DO 30 J = IA2(I), IA2(I+1) - 1
            IF ((IA1(J).LT.0) .OR. (IA1(J).GT.N)) THEN
               NIND = NIND + 1
            END IF
 30      CONTINUE
 20   CONTINUE
      IF (NIND .GT. 0) THEN
         IERROR = 16
         GOTO 9999
      END IF
 9999 CONTINUE
      RETURN
      END
