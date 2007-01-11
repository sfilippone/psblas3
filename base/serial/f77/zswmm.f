C     SUBROUTINE ZSWMM(TRANS,M,N,K,ALPHA,FIDA,DESCRA,A,IA1,IA2,
C                      INFOA,B,LDB,BETA,C,LDC,WORK,LWORK,IERROR)
C     Purpose
C     =======
C
C     Computing   C <-- ALPHA A   B + BETA C    or
C                 C <-- ALPHA At  B + BETA C    or
C                 C <-- ALPHA Atc B + BETA C
C     Called by ZCSMM
C     Actual computing performed by sparse Toolkit kernels.
C     This routine selects the proper kernel for each
C     data structure.
C
C     Parameters
C     ==========
C
C     TRANS    - CHARACTER*1
C             On entry TRANS specifies if the routine operates with matrix A
C             or with the transpose of A as follows:
C                TRANS = 'N'         ->  use matrix A
C                TRANS = 'T'         ->  use A' (transpose of matrix A)
C                TRANS = 'C'         ->  use conjugate transpose of A
C             Unchanged on exit.
C
C             N.B.: M, N for C matrix
C                   M, K for A matrix
C                   K, N for B matrix
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
C     K        - INTEGER
C             On entry: number of columns of matrix A (A') and
C                       number of rows of matrix B
C             Unchanged on exit.
C
C     ALPHA    - COMPLEX*16
C             On entry: multiplicative constant.
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
C     A        - COMPLEX*16 array of DIMENSION (*)
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
      SUBROUTINE ZSWMM(TRANS,M,N,K,ALPHA,FIDA,DESCRA,A,IA1,IA2,
     &                 INFOA,B,LDB,BETA,C,LDC,WORK,LWORK,IERROR)
C     .. Scalar Arguments ..
      INTEGER    M,N,K,LDB,LDC,LWORK,IERROR
      CHARACTER  TRANS
      COMPLEX*16 ALPHA,BETA
C     .. Array Arguments ..
      INTEGER    IA1(*),IA2(*),INFOA(*)
      CHARACTER  DESCRA*11, FIDA*5
      COMPLEX*16 A(*),B(LDB,*),C(LDC,*),WORK(*)
C     .. Local Array ..
      INTEGER    INT_VAL(5), ERR_ACT
      COMPLEX*16 Z_VAL(5)
      CHARACTER*30 NAME,  STRINGS(2)
C     .. External Subroutines ..
      EXTERNAL      ZCSRMM
C     +              ZJADMM
C     .. Executable Statements ..
C
C     Switching on FIDA: proper sparse BLAS routine is selected
C     according to data structure
C
      NAME = 'ZSWMM\0'
      IERROR = 0 
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)
      
      IF (FIDA(1:3).EQ.'CSR') THEN
C
C        A, IA1, IA2 --->  AR,   JA,   IA
C                         VAL, INDX, PNTR
C        INFOA(*) not used
        CALL  ZCSRMM(TRANS,M,N,K,ALPHA,DESCRA,A,IA1,                   
     &    IA2,B,LDB,BETA,C,LDC,WORK,LWORK)
C
C     JAD format not jet implemented
C
      ELSE IF (FIDA(1:3).EQ.'COO') THEN
C      
C        INFOA(*) not used    
C
        CALL  ZCOOMM(TRANS,M,N,K,ALPHA,DESCRA,A,IA1,
     +    IA2,INFOA,B,LDB,BETA,C,LDC,WORK)                         
      ELSE IF ((FIDA(1:3).EQ.'JAD').AND.(.FALSE.)) THEN
C      
C        INFOA(*) not used    
C
C         CALL  ZJADMM(TRANS,M,N,K,ALPHA,DESCRA,A,IA1,  
C     &                IA2,B,LDB,BETA,C,LDC,WORK)     
      ELSE
C
C     This data structure not yet considered
C
        IERROR = 3010
        STRINGS(1) = FIDA(1:4)
        CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
        GOTO 9999 
        
      END IF
      CALL FCPSB_ERRACTIONRESTORE(ERR_ACT)

      RETURN
      
 9999 CONTINUE
      CALL FCPSB_ERRACTIONRESTORE(ERR_ACT)

      IF ( ERR_ACT .NE. 0 ) THEN 
         CALL FCPSB_SERROR()
         RETURN
      ENDIF

      END

