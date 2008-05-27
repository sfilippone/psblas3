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
C     SUBROUTINE CSWSM(TRANS,M,N,ALPHA,UNITD,D,FIDT,DESCRT,T,IT1,IT2,
C                      INFOT,B,LDB,BETA,C,LDC,WORK,LWORK,IERROR)
C
C     Purpose
C     =======
C
C     Solving triangular systems of equations with multiple right-hand sides
C                 C <-- ALPHA D T-1 B + BETA C   or
C                 C <-- ALPHA D T-t B + BETA C   or
C                 C <-- ALPHA T-1 D B + BETA C   or
C                 C <-- ALPHA T-t D B + BETA C
C     Actual computing performed by sparse Toolkit kernels.
C     This routine selects the proper kernel for each
C     data structure.
C
C     Parameters
C     ==========
C
C     TRANS    - CHARACTER*1
C             On entry TRANS specifies whether the routine operates with
C             matrix T or with the transpose of T as follows:
C                TRANS = 'N'         ->  use matrix T
C                TRANS = 'T'         ->  use T' (transpose of matrix T)
C                TRANS = 'C'         ->  use conjugate transpose of T
C             Unchanged on exit.
C
C     M        - INTEGER
C             On entry: number of rows and columns of matrix T
C             and number of rows of matrices B and C.
C             Unchanged on exit.
C
C     N        - INTEGER
C             On entry: number of columns of matrices B and C
C             (number of right-hand sides).
C             Unchanged on exit.
C
C     ALPHA    - COMPLEX*16
C             On entry: multiplicative constant.
C             Unchanged on exit.
C
C     UNITD    - CHARACTER*1
C             On entry UNITD specifies whether the diagonal matrix is unit
C             or whether row or column scaling has to be performed, as follows:
C                UNITD = 'U'         ->  unit matrix (no scaling)
C                UNITD = 'L'         ->  scale on the left (row scaling)
C                UNITD = 'R'         ->  scale on the right (column scaling)
C                UNITD = 'B'         ->  scale on the right and on the left
C                                             with D^1/2
C             Unchanged on exit.
C
C     D        - COMPLEX*16 array of dimension (M)
C             On entry D specifies the main diagonal of the matrix used
C             for scaling.
C             Unchanged on exit.
C
C     FIDT     - CHARACTER*5
C             On entry FIDT defines the format of the input sparse matrix.
C             Unchanged on exit.
C
C     DESCRT   - CHARACTER*1 array of DIMENSION (9)
C             On entry DESCRT describes the characteristics of the input
C             sparse matrix.
C             Unchanged on exit.
C
C
C     T        - COMPLEX*16 array of DIMENSION (*)
C             On entry T specifies the values of the input sparse
C             matrix.
C             Unchanged on exit.
C
C     IT1      - INTEGER array of dimension (*)
C             On entry IT1 holds integer information on input sparse
C             matrix.  Actual information will depend on data format used.
C             Unchanged on exit.
C
C     IT2      - INTEGER array of dimension (*)
C             On entry IT2 holds integer information on input sparse
C             matrix.  Actual information will depend on data format used.
C             Unchanged on exit.
C
C     INFOT     - INTEGER array of dimension (10)
C             On entry can hold auxiliary information on input matrices
C             formats or environment of subsequent calls.
C             Might be changed on exit.
C
C     B        - COMPLEX*16 array of dimension (LDB,*)
C             On entry: matrix of right-hand sides
C             Unchanged on exit.
C
C     LDB      - INTEGER
C             On entry: leading dimension of B.
C             Unchanged on exit.
C
C     BETA     - COMPLEX*16
C             On entry: multiplicative constant.
C             Unchanged on exit.
C
C     C        - COMPLEX*16 array of dimension (LDC,*)
C             On exit: solutions of triangular systems
C
C     LDC      - INTEGER
C             On entry: leading dimension of C.
C             Unchanged on exit.
C
C     WORK     - COMPLEX*16 array of dimension (LWORK)
C             On entry: work area.
C             On exit INT(WORK(1)) contains the minimum value
C             for LWORK satisfying DSWSM memory requirements.
C
C     LWORK    - INTEGER
C             On entry LWORK specifies the dimension of WORK
C             Unchanged on exit.
C
C     IERROR   - INTEGER
C             On exit IERROR contains the value of error flag as follows:
C             IERROR = 0   no error
C             IERROR > 0   warning
C             IERROR < 0   fatal error
C
C     Note
C     ====
C     All checks on argument are performed in the calling routine.
C
C
      SUBROUTINE CSWSM(TRANS,M,N,ALPHA,UNITD,D,FIDT,DESCRT,T,IT1,IT2,
     &                 INFOT,B,LDB,BETA,C,LDC,WORK,LWORK,IERROR)
      use psb_const_mod
      use psb_string_mod
C     .. Scalar Arguments ..
      INTEGER    M, N, LDB, LDC, LWORK, IERROR
      CHARACTER  UNITD, TRANS
      complex(psb_spk_) ALPHA, BETA
C     .. Array Arguments ..
      INTEGER    IT1(*), IT2(*), INFOT(*)
      CHARACTER  DESCRT*11, FIDT*5
      complex(psb_spk_) T(*), B(LDB,*), C(LDC,*), D(*), WORK(*)
      integer    ERR_ACT
C     .. Local Scalars ..
      INTEGER    ONE
C     .. Parameters ..
      PARAMETER  (ONE=1)
C     .. External Subroutines ..
      EXTERNAL   CCSRSM, CCOPY
      CHARACTER*20      NAME

C     .. Executable Statements ..
      NAME   = 'CSWSM\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)
C
C     Check for identity matrix
C
      IF(psb_toupper(DESCRT(1:1)).EQ.'D' .AND.
     +  psb_toupper(DESCRT(3:3)).EQ.'U') THEN
         CALL CCOPY(M,B,ONE,C,ONE)
         GOTO 9999
      ENDIF
C
C     Switching on FIDT: proper sparse BLAS routine is selected
C     according to data structure
C
      IF (psb_toupper(FIDT(1:3)).EQ.'CSR') THEN
C
C        T, IT1, IT2 --->  AR,   JA,   IA
C                         VAL, INDX, PNTR
C        INFOT(*) not used
C
         CALL  CCSRSM(TRANS,M,N,UNITD,D,ALPHA,DESCRT,T,IT1,      
     &                IT2,B,LDB,BETA,C,LDC,WORK,LWORK)
      ELSE IF ((psb_toupper(FIDT(1:3)).EQ.'JAD')
     +     .AND.(.FALSE.)) THEN
C
C        .. JAD format not yet supported
C         
C         CALL  CJADSM(TRANS,M,N,D,UNITD,0,ALPHA,DESCRT,T,IT1,IT2,0,
C     &      B,LDB,BETA,C,LDC,WORK)
C         
      ELSE
C
C     This data structure not yet considered
C
         IERROR = 3010
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      END IF
      CALL FCPSB_ERRACTIONRESTORE(ERR_ACT)
      RETURN

 9999 CONTINUE
      CALL FCPSB_ERRACTIONRESTORE(ERR_ACT)

      IF ( ERR_ACT .NE. 0 ) THEN 
         CALL FCPSB_SERROR()
      ENDIF

      RETURN
      END


