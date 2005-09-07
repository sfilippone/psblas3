      SUBROUTINE DCOOSM(TRANST,M,N,UNITD,D,ALPHA,DESCRA,A,IA,JA,INFOA,
     *                  B,LDB,BETA,C,LDC,WORK,LWORK,IERROR)
      IMPLICIT NONE
      LOGICAL           DEBUG
      PARAMETER         (DEBUG=.FALSE.)
      DOUBLE PRECISION  ALPHA, BETA
      INTEGER           LDB, LDC, LWORK, M, N, IERROR
      CHARACTER         UNITD, TRANST
      DOUBLE PRECISION  A(*), B(LDB,*), C(LDC,*), D(*), WORK(*)
      INTEGER           IA(*), JA(*), INFOA(*), INT_VAL(5)
      CHARACTER         DESCRA*11
      INTEGER           I, K, ERR_ACT
      CHARACTER         DIAG, UPLO
      EXTERNAL          XERBLA
      INTRINSIC         DBLE, IDINT
      CHARACTER*20      NAME

      NAME = 'DCOOSM\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

      IF((ALPHA.NE.1.D0) .OR. (BETA.NE.0.D0))then
         IERROR=5
         CALL PSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      ENDIF
      if (debug) write(*,*) 'DCOOSM ',m
      if (debug) write(*,*) 'DCOOSM ',m,ierror
      
      UPLO = '?'
      IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'U') UPLO = 'U'
      IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'L') UPLO = 'L'
      IF (UPLO.EQ.'?') THEN
         IERROR=5
         CALL PSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      END IF
      IF (DESCRA(3:3).EQ.'N') DIAG = 'N'
      IF (DESCRA(3:3).EQ.'U') DIAG = 'U'
      IF(UNITD.EQ.'B') THEN
         IERROR=5
         CALL PSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      ENDIF
      IF (UNITD.EQ.'R') THEN
         DO 40 I = 1, N
            DO 20 K = 1, M
               B(K,I) = B(K,I)*D(K)
   20       CONTINUE
   40    CONTINUE
      END IF

      DO 60 I = 1, N
         CALL DCOOSV(UPLO,TRANST,DIAG,M,A,IA,JA,INFOA,
     +     B(1,I),C(1,I),IERROR)
   60 CONTINUE
      IF(IERROR.NE.0) THEN
         INT_VAL(1)=IERROR
         CALL FCPSB_ERRPUSH(4012,NAME,INT_VAL)
         GOTO 9999
      END IF

      IF (UNITD.EQ.'L') THEN
         DO 45 I = 1, N
            DO 25 K = 1, M
               C(K,I) = C(K,I)*D(K)
   25       CONTINUE
   45    CONTINUE
      END IF
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
