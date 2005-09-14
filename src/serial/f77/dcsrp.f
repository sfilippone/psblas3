C     SUBROUTINE DCSRP(TRANS,M,N,FIDA,DESCRA,IA1,IA2,INFOA,
C                      P,WORK,LWORK,IERROR)
C
C     Purpose
C     =======
C
C     Performing column permutation of a sparse matrix.
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
C     IA1      - INTEGER array of dimension (*)
C             On entry IA1 holds integer information on input sparse
C             matrix.  Actual information will depend on data format used.
C             On exit contain integer information on permuted matrix.
C
C     IA2      - INTEGER array of dimension (*)
C             On entry IA2 holds integer information on input sparse
C             matrix.  Actual information will depend on data format used.
C             On exit contain integer information on permuted matrix.
C
C     INFOA     - INTEGER array of dimension (10)
C             On entry can hold auxiliary information on input matrices
C             formats or environment of subsequent calls.
C             Might be changed on exit.
C
C     P        - INTEGER array of dimension (M)
C             On entry P specifies the column permutation of matrix A
C             (P(1) == 0 if no permutation).
C             Unchanged on exit.
C
C     WORK     - DOUBLE PRECISION array of dimension (LWORK)
C             On entry: work area.
C             On exit INT(WORK(1)) contains the minimum value
C             for LWORK satisfying DCSRP memory requirements.
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
C     WORK     - DOUBLE PRECISION array of dimension (LWORK)
C             Work area.
C
C     LWORK    - INTEGER
C             On entry LWORK specifies the dimension of WORK.
C             LWORK must be greater than zero.
C             On exit LWORK is the maximum between the initial value and
C             the minimum value satisfying DCSRP memory requirements.
C
C     IERROR    - INTEGER
C             On exit IERROR contains the value of error flag as follows:
C             IERROR = 0      no error
C             IERROR = 4      error on dimension of vector WORK
C             IERROR = 32     unknown flag TRANS
C             IERROR = 64     LWORK  <=  0
C             IERROR = 128    this data structure not yet considered
C                                                                        
C     Notes                                                              
C     =====                                                              
C     It is not possible to call this subroutine with LWORK=0 to get     
C     the minimal value for LWORK. This functionality needs a better     
C     connection with DxxxMM                                             
C
C
      SUBROUTINE DCSRP(TRANS,M,N,FIDA,DESCRA,IA1,IA2,INFOA,
     +  P,WORK,LWORK,IERROR)
      IMPLICIT NONE                                                      
C     .. Scalar Arguments ..
      INTEGER          LWORK, M, N, IERROR
      CHARACTER        TRANS
C     .. Local Scalars..
      INTEGER          ERR_ACT
C     .. Array Arguments ..
      DOUBLE PRECISION WORK(LWORK)
      INTEGER          IA1(*), IA2(*), INFOA(*), P(*), INT_VAL(5)
      CHARACTER        DESCRA*11, FIDA*5
C     .. External Subroutines ..
      EXTERNAL          DCSRRP
      logical          debug
      parameter        (debug=.false.)

      CHARACTER*20      NAME
C
C     .. Executable Statements ..
C
C
C     Check on M, N, TRANS
C
      NAME = 'DCSRP\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

      IF     (M.LE.0) THEN
         IERROR = 1
      ELSE IF(N.LE.0) THEN
         IERROR = 3
      ELSE IF(TRANS.NE.'N' .AND. TRANS.NE.'T' .AND. TRANS.NE.'C') THEN
         IERROR = -1
      ENDIF
C
C     Error handling
C
      IF(IERROR.LT.0) THEN
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      ENDIF
C
C     Check for M, N, P
C
      IF(M.LE.0 .OR. N.LE.0 .OR. P(1).EQ.0) THEN
         GOTO 9999
      ENDIF
C
C     Switching on FIDA
C
c$$$      write(0,*) 'DCSRP FORMAT: ',fida
      IF  (FIDA(1:3).EQ.'CSR') THEN
C
C        Permuting CSR structure
C
         CALL DCSRRP(TRANS,M,N,DESCRA,IA1,IA2,P,WORK,LWORK)
      ELSE IF (FIDA(1:3).EQ.'JAD') THEN
        if (debug) write(0,*) 'Calling djadrp',m,p(1),lwork
         CALL DJADRP(TRANS,M,N,DESCRA,IA1,IA2,P,WORK,LWORK)
      ELSE
C
C        This data structure not yet considered
C
         IERROR  = 4
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
