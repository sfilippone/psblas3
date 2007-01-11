      SUBROUTINE ZCOOSM(TRANST,M,N,UNITD,D,ALPHA,DESCRA,A,IA,JA,INFOA,
     *                  B,LDB,BETA,C,LDC,WORK,LWORK,IERROR)
      LOGICAL           DEBUG
      PARAMETER         (DEBUG=.FALSE.)
      COMPLEX*16        ALPHA, BETA
      INTEGER           LDB, LDC, LWORK, M, N, IERROR
      CHARACTER         UNITD, TRANST
      COMPLEX*16        A(*), B(LDB,*), C(LDC,*), D(*), WORK(*)
      INTEGER           IA(*), JA(*), INFOA(*)
      CHARACTER         DESCRA*11
      INTEGER           I, K
      CHARACTER         DIAG, UPLO
      EXTERNAL          XERBLA

      IF((ALPHA.NE.1.D0) .OR. (BETA.NE.0.D0))then
         call xerbla('DCSSM ',9)
         RETURN
      ENDIF
      if (debug) write(*,*) 'ZCOOSM ',m
      if (debug) write(*,*) 'ZCOOSM ',m,ierror
      
      UPLO = '?'
      IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'U') UPLO = 'U'
      IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'L') UPLO = 'L'
      IF (UPLO.EQ.'?') THEN
         CALL XERBLA('ZCSSM ',10)
         RETURN
      END IF
      IF (DESCRA(3:3).EQ.'N') DIAG = 'N'
      IF (DESCRA(3:3).EQ.'U') DIAG = 'U'
      IF(UNITD.EQ.'B') THEN
         CALL XERBLA('ZCSSM ',11)
         RETURN
      ENDIF
      IF (UNITD.EQ.'R') THEN
         DO 40 I = 1, N
            DO 20 K = 1, M
               B(K,I) = B(K,I)*D(K)
   20       CONTINUE
   40    CONTINUE
      END IF

      DO 60 I = 1, N
         CALL ZCOOSV(UPLO,TRANST,DIAG,M,A,IA,JA,INFOA,
     +     B(1,I),C(1,I),IERROR)
   60 CONTINUE
      IF (UNITD.EQ.'L') THEN
         DO 45 I = 1, N
            DO 25 K = 1, M
               C(K,I) = C(K,I)*D(K)
   25       CONTINUE
   45    CONTINUE
      END IF
      RETURN
      END
