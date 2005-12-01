      SUBROUTINE DCSRSM(TRANST,M,N,UNITD,D,ALPHA,DESCRA,A,JA,IA,
     *                  B,LDB,BETA,C,LDC,WORK,LWORK,IERROR)
      DOUBLE PRECISION  ALPHA, BETA
      INTEGER           LDB, LDC, LWORK, M, N, IERROR
      CHARACTER         UNITD, TRANST
      DOUBLE PRECISION  A(*), B(LDB,*), C(LDC,*), D(*), WORK(*)
      INTEGER           IA(*), JA(*)
      CHARACTER         DESCRA*11
      INTEGER           I, K
      CHARACTER         DIAG, UPLO
      LOGICAL DEBUG
      PARAMETER (DEBUG=.FALSE.)
C     .. Local Arrays ..
      CHARACTER*20       NAME
      INTEGER            INT_VAL(5)

      NAME = 'DCSRSM\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)
      int_Val(1)=0
      UPLO = '?'
      IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'U') UPLO = 'U'
      IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'L') UPLO = 'L'
      IF (UPLO.EQ.'?') THEN
         int_val(1) = 7
         IERROR=5
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      END IF
      IF (DESCRA(3:3).EQ.'N') DIAG = 'N'
      IF (DESCRA(3:3).EQ.'U') DIAG = 'U'
      IF(UNITD.EQ.'B') THEN
         IERROR=5
         int_val(1) = 4
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      ENDIF
      IF (UNITD.EQ.'R') THEN
         DO 40 I = 1, N
            DO 20 K = 1, M
               B(K,I) = B(K,I)*D(K)
   20       CONTINUE
   40    CONTINUE
      END IF
      if ((alpha.ne.1.d0) .or.(beta .ne.0.0d0)) then 
        if (lwork .lt. m) then           
          int_val(1) = 17
          IERROR=5
          CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
          GOTO 9999
        END IF
        DO I = 1, N
          CALL DCSRSV(UPLO,TRANST,DIAG,M,A,JA,IA,B(1,I),work)
          do k=1,m
            c(k,i) = beta*c(k,i) + alpha*work(k)
          enddo
        enddo        

      else
        DO 60 I = 1, N
          CALL DCSRSV(UPLO,TRANST,DIAG,M,A,JA,IA,B(1,I),C(1,I))
 60     CONTINUE
      endif
      IF(IERROR.NE.0) THEN
         INT_VAL(1)=IERROR
         CALL FCPSB_ERRPUSH(4012,NAME,INT_VAL)
         GOTO 9999
      END IF

      if (debug) then 
        write(0,*) 'Check from DCSRSM'
        do k=1,m
          write(0,*) k, b(k,1),c(k,1)
        enddo
      endif

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
