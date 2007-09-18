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
C     SUBROUTINE ZCRCR(TRANS,M,N,UNITD,D,DESCRA,A,IA1,IA2,INFOA,IP1,
C                      DESCRN,AN,IAN1,IAN2,INFON,IP2,LAN,LIAN1,LIAN2,
C                      WORK,LWORK,IERROR)
C
C     Purpose: CSR to CSR format conversion
C     =======
C
C     Parameter:
C     =========
C
C     TRANS    - CHARACTER*1
C             On entry TRANS specifies whether the routine will convert
C             matrix A or the transpose of A as follows:
C                TRANS = 'N'         ->  convert matrix A
C                TRANS = 'T' or 'C'  ->  convert A' (the transpose of A)
C             Unchanged on exit.
C
C     M        - INTEGER
C             On entry: number of rows of matrix A (A')
C             and number of rows of matrix H
C             Unchanged on exit.
C
C     N        - INTEGER
C             On entry: number of columns of matrix A (A')
C             and number of columns of matrix H
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
C     D        - DOUBLE PRECISION array of dimension (M)
C             On entry D specifies the main diagonal of the matrix used
C             for scaling.
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
C     INFOA    - INTEGER array of dimension (10)
C             On entry can hold auxiliary information on input matrices
C             formats or environment of subsequent calls.
C             Might be changed on exit.
C
C     IP1       - INTEGER array of dimension (M)
C             On exit IP1 specifies the row permutation of matrix AN
C             (IP1(1) == 0 if no permutation).
C
C     DESCRN   - CHARACTER*1 array of DIMENSION (9)
C             On exit DESCRN describes the characteristics of the input
C             sparse matrix.
C             Unchanged on exit.
C
C     AN       - DOUBLE PRECISION array of DIMENSION (LAN)
C             On exit AN specifies the values of the output sparse
C             matrix. If LAN=0, INT(AN(1)) is the minimum value for LAN
C             satisfying DSPDP memory requirements.
C
C     IAN1     - INTEGER array of dimension (LIAN1)
C             On exit IAN1 holds integer information on output sparse
C             matrix.  Actual information will depend on data format used.
C             If LIAN1=0, INT(IAN1(1)) is the minimum value for LIAN1
C             satisfying DSPDP memory requirements.
C
C     IAN2     - INTEGER array of dimension (LIAN2)
C             On exit IAN2 holds integer information on output sparse
C             matrix.  Actual information will depend on data format used.
C             If LIAN2=0, INT(IAN2(1)) is the minimum value for LIAN2
C             satisfying DSPDP memory requirements.
C
C     INFON    - INTEGER array of dimension (10)
C             On exit can hold auxiliary information on output matrices
C             formats or environment of subsequent calls.
C
C     IP2      - INTEGER array of dimension (M)
C             On exit IP2 specifies the column permutation of matrix AN
C             (IP2(1) == 0 if no permutation).
C
C     LAN      - INTEGER
C             On entry LAN specifies the dimension of AN
C             LAN must satisfy memory required from the new data structure.
C             Unchanged on exit.
C
C     LIAN1    - INTEGER
C             On entry LH1 specifies the dimension of IAN1
C             LH1 must satisfy memory required from the new data structure.
C             Unchanged on exit.
C
C     LIAN2    - INTEGER
C             On entry LIAN2 specifies the dimension of IAN2
C             LIAN2 must satisfy memory required from the new data structure.
C             Unchanged on exit.
C
C     WORK     - DOUBLE PRECISION array of dimension (LWORK)
C             On entry: work area.
C             On exit INT(WORK(1)) contains the minimum value
C             for LWORK satisfying DSPDP memory requirements.
C
C     LWORK    - INTEGER
C             On entry LWORK specifies the dimension of WORK
C             LWORK must satisfy memory necessary for the data conversion.
C             Unchanged on exit.
C
C     IERROR   - INTEGER
C             On exit IERROR contains the value of error flag as follows:
C             IERROR = 0   no error
C             IERROR > 0   error
C
C
      SUBROUTINE ZCRCR(TRANS,M,N,UNITD,D,DESCRA,A,IA1,IA2,INFOA,IP1,
     *  DESCRN,AN,IAN1,IAN2,INFON,IP2,LAN,LIAN1,LIAN2,
     *  WORK,LWORK,IERROR)
      use psb_string_mod
      IMPLICIT NONE                                                      
C
C     .. Scalar Arguments ..
      INTEGER          M, N, LAN, LIAN1, LIAN2, LWORK, IERROR
      CHARACTER        TRANS, UNITD
C     .. Array Arguments ..
      complex(kind(1.d0))  A(*), AN(*), D(*), WORK(LWORK)
      INTEGER          IA1(*), IA2(*), IAN1(*), IAN2(*), IP1(*), IP2(*),
     *  INFOA(*), INFON(*)
      CHARACTER        DESCRA*11, DESCRN*11
C     .. Local Scalars ..
      INTEGER          I, J, ERR_ACT
      LOGICAL          EXIT
c     .. Local Arrays ..
      character  idescra*11
      CHARACTER*20       NAME
      INTEGER            INT_VAL(5)

C     .. Intrinsic Functions ..
      INTRINSIC        DBLE, DSQRT

C     .. Executable Statements ..
C
      EXIT=.FALSE.
      NAME = 'DCOCO\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)
C
C       Check for argument errors
C
      idescra=toupper(descra)
      IF(((IDESCRA(1:1) .EQ. 'S' .OR. IDESCRA(1:1) .EQ. 'H' .OR.
     &  IDESCRA(1:1) .EQ. 'A') .AND. (toupper(UNITD) .NE. 'B'))  .OR.
     &  (.NOT.((IDESCRA(3:3).EQ.'N').OR.(IDESCRA(3:3).EQ.'L').OR.
     +  (IDESCRA(3:3).EQ.'U'))) .OR.
     +  toupper(TRANS).NE.'N') THEN
        IERROR = 20
      ENDIF
      IF(LAN.LT.(IA2(M+1)-1)) THEN
        IF     (LAN.LE.0) THEN
          EXIT=.TRUE.
          AN(1) = DBLE(IA2(M+1)-1)
        ELSE
          IERROR = 21
        ENDIF
      ENDIF
      IF(LIAN1.LT.(IA2(M+1)-1)) THEN
        IF     (LAN.LE.0) THEN
          EXIT=.TRUE.
          IAN1(1) = IA2(M+1)-1
        ELSE
          IERROR = 22
        ENDIF
      ENDIF
      IF(LIAN2.LT.(M+1)) THEN
        IF     (LAN.LE.0) THEN
          EXIT=.TRUE.
          IAN2(1) = M+1
        ELSE
          IERROR = 23
        ENDIF
      ENDIF
      IF ((IDESCRA(1:1) .EQ. 'S' .OR. IDESCRA(1:1) .EQ. 'H' .OR.
     &  IDESCRA(1:1) .EQ. 'A') .AND. (toupper(UNITD) .EQ. 'B')) THEN
        IF (LWORK.LT.M) THEN
          IF     (LWORK.LE.0) THEN
            EXIT=.TRUE.
          ELSE
            IERROR = 25
          ENDIF
          WORK(1) = DBLE(M)
        ENDIF
      ELSE
        IF (LWORK.LT.0) THEN
          WORK(1) = 0.D0
        ENDIF
      ENDIF
C
C     Error handling
C
      IF(IERROR.NE.0) THEN
        CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
        GOTO 9999
      END IF

      IF (EXIT)    goto 9998
C
C     Set DESCRN, IP1, IP2
C
      DESCRN(1:3) = IDESCRA(1:3)

      IP1(1)=0
      IP2(1)=0
C
C     Compute output matrix
C
      DO 20 I = 1, M+1
        IAN2(I) = IA2(I)
 20   CONTINUE
      IF ((IDESCRA(1:1) .EQ. 'S' .OR. IDESCRA(1:1) .EQ. 'H' .OR.
     &  IDESCRA(1:1) .EQ. 'A') .AND. (toupper(UNITD) .EQ. 'B')) THEN
        DO 30 I = 1, M
          WORK(I) = DBLE(DSQRT(ABS(D(I))))
 30     CONTINUE
        DO 40 I = 1, M
          DO 50 J = IA2(I), IA2(I+1)-1
            AN(J)   = WORK(I) * A(J) * WORK(IA1(J))
            IAN1(J) = IA1(J)
 50       CONTINUE
 40     CONTINUE
      ELSE IF (toupper(UNITD) .EQ. 'L') THEN
        DO 60 I = 1, M
          DO 70 J = IA2(I), IA2(I+1)-1
            AN(J)   = D(I) * A(J)
            IAN1(J) = IA1(J)
 70       CONTINUE
 60     CONTINUE
      ELSE IF (toupper(UNITD) .EQ. 'R') THEN
        DO 80 I = 1, M
          DO 90 J = IA2(I), IA2(I+1)-1
            AN(J)   = A(J) * D(IA1(J))
            IAN1(J) = IA1(J)
 90       CONTINUE
 80     CONTINUE
      ELSE IF (toupper(UNITD) .EQ. 'U') THEN
        DO 100 J = 1, IA2(M+1)-1
          AN(J)   = A(J)
          IAN1(J) = IA1(J)
 100    CONTINUE
      ENDIF

 9998 CONTINUE
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


