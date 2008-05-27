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
C     SUBROUTINE  DCSRS(TRANS,M,N,FIDA,DESCRA,A,IA1,IA2,      &
C    &                 INFOA,ROWSUM,IERROR)
C     Purpose
C     =======
C
C     Computing sum of absolute values for rows of distributed matrix
C                 ROWSUM(IX) = ASUM(A(IX, 1..N))
C                 IX = 1..M
C
C     Parameters
C     ==========
C
C     TRANS    - CHARACTER*1
C             On entry TRANS specifies if the routine operates with matrix A
C             or with the transpose of A as follows:
C                TRANS = 'N'         ->  use matrix A
C                TRANS = 'T' or 'C'  ->  use A' (transpose of matrix A)
C             Unchanged on exit.
C
C     M        - INTEGER
C             On entry: number of rows of matrix A (A') and
C                       number of rows of matrix C
C             Unchanged on exit.
C
C     N        - INTEGER
C             On entry: number of columns of matrix B
C             and number of columns of matrix C.
C             Unchanged on exit.
C
C     FIDA     - CHARACTER*5
C             On entry FIDA defines the format of the input sparse matrix.
C             Unchanged on exit.
C
C     DESCRA   - CHARACTER*1 array of DIMENSION (9)
C             On entry DESCRA describes the characteristics of the input
C             sparse matrix.
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
C     INFOA     - INTEGER array of length 10.
C             On entry can hold auxiliary information on input matrices
C             formats or environment of subsequent calls.
C             Might be changed on exit.
C
C     IERROR   - INTEGER
C             On exit IERROR contains the value of error flag as follows:
C             IERROR = 0   no error
C             IERROR > 0   warning
C             IERROR < 0   fatal error
C
C     ROWSUM   - DOUBLE PRECISION array of dimension (*)
C             On exit this vector contains the sum of absolute values 
C             of elements of a row (AMAX of row array).
C             It is required that it has dimension:
C                ROWSUM(M) if the subroutine in called with the 'N' option
C                ROWSUM(N) in other cases ('T' or 'C' options).
C
      SUBROUTINE  ZCSRWS(TRANS,M,N,FIDA,DESCRA,A,IA1,IA2,
     &  INFOA,ROWSUM,IERROR)
      use psb_const_mod
      use psb_string_mod
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           M,N, IERROR
      CHARACTER         TRANS
C     .. Array Arguments ..
      INTEGER           IA1(*),IA2(*),INFOA(*)
      CHARACTER         DESCRA*11, FIDA*5
      complex(psb_dpk_)        A(*), ROWSUM(*)
C     .. Local Array..
      INTEGER           INT_VAL(5), ERR_ACT
      CHARACTER*30      NAME,STRINGS(2)
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, IDINT
C     .. Executable Statements ..
C
C     Check for argument errors
C
      
      IERROR = 0
      NAME = 'ZCSRWS\0'
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

      IF     (M.LT.0) THEN
        IERROR = 10
        INT_VAL(1) = 2
        INT_VAL(2) = M
      ELSE IF (N.LT.0) THEN
        IERROR = 10
        INT_VAL(1) = 3
        INT_VAL(2) = N
      ELSE IF (psb_toupper(TRANS).NE.'T' .AND.
     +    psb_toupper(TRANS).NE.'N' .AND.
     +    psb_toupper(TRANS).NE.'C') THEN
        IERROR = 40
        INT_VAL(1) = 1
        STRINGS(1) = TRANS//'\0'
      ENDIF

C
C     Error handling
C
      IF(IERROR.NE.0) THEN
        CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
        GOTO 9999
      ENDIF

      IF(M.LE.0 .OR. N.LE.0) THEN
        GOTO 9999
      ENDIF

      IF (psb_toupper(FIDA(1:3)).EQ.'CSR') THEN
        CALL ZCSRRWS(TRANS,M,N,DESCRA,A,IA1,IA2,
     +    INFOA,ROWSUM,IERROR)
      ELSE IF (psb_toupper(FIDA(1:3)).EQ.'COO') THEN
        CALL ZCOORWS(TRANS,M,N,DESCRA,A,IA1,IA2,
     +    INFOA,ROWSUM,IERROR)
c$$$      ELSE IF (psb_toupper(FIDA(1:3)).EQ.'JAD') THEN
c$$$        CALL DJDRWS(TRANS,M,N,DESCRA,A,IA1,IA2,
c$$$     +     INFOA,ROWSUM,IERROR)
      ELSE
C
C     This data structure not yet considered
C
        IERROR = 3010
        strings(1) = fida//'\0'
      ENDIF

      IF(IERROR.NE.0) THEN
        CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
        GOTO 9999 
      ENDIF
      CALL FCPSB_ERRACTIONRESTORE(ERR_ACT)
      RETURN

 9999 CONTINUE
      CALL FCPSB_ERRACTIONRESTORE(ERR_ACT)

      IF ( ERR_ACT .NE. 0 ) THEN 
         CALL FCPSB_SERROR()
         RETURN
      ENDIF
      
      RETURN
      END


