      SUBROUTINE DJADSM(TRANST,M,N,VDIAG,TDIAG,PERMQ,ALPHA,DESCRA,
     +   AR,JA,IA,PERMP,B,LDB,BETA,C,LDC,WORK)
C
C
C     .. Scalar Arguments ..
      INTEGER           LDB, LDC, M, N
      CHARACTER         TDIAG, TRANST
      DOUBLE PRECISION  ALPHA, BETA
C     .. Array Arguments ..
      DOUBLE PRECISION  AR(*), B(LDB,*), C(LDC,*), VDIAG(*), WORK(*)
      INTEGER           IA(*), JA(*), PERMP(*), PERMQ(*)
      CHARACTER         DESCRA*11
C     .. Local Scalars ..
      INTEGER           PIA, PJA, PNG
      INTEGER           I, K, ERR_ACT
      CHARACTER         UPLO,UNITD
      logical debug
      parameter (debug=.false.)
      CHARACTER*20      NAME
      INTEGER           INT_VAL(5)
C     .. Executable Statements ..
C
      NAME = 'DJADSM\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

      IF((ALPHA.NE.1.D0) .OR. (BETA.NE.0.D0))then
         IERROR=5
         CALL PSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      ENDIF
      UPLO = '?'
      IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'U') UPLO = 'U'
      IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'L') UPLO = 'L'
C
      IF (UPLO.EQ.'?') THEN
         IERROR=5
         CALL PSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      END IF

      IF (DESCRA(3:3).NE.'U') THEN
         IERROR=5
         CALL PSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      END IF
      UNITD=DESCRA(3:3)
C
C        B = INV(A)*B  OR B=INV(A')*B
C
      if (debug) write(0,*) 'DJADSM : ',m,n,' ',tdiag

      IF (TDIAG.EQ.'R') THEN
        if (debug) write(0,*) 'DJADSM : Right Scale',m,n
        DO  I = 1, N
          DO  K = 1, M
            B(K,I) = B(K,I)*VDIAG(K)
          ENDDO
        ENDDO
      END IF
      
      PNG = IA(1)
      PIA = IA(2)
      PJA = IA(3)

      DO I = 1, N
         CALL DJADSV(UNITD,M,IA(PNG),
     +      AR,JA,IA(PIA),IA(PJA),B(1,I),C(1,I),IERROR)
      ENDDO
      IF(IERROR.NE.0) THEN
         INT_VAL(1)=IERROR
         CALL FCPSB_ERRPUSH(4012,NAME,INT_VAL)
         GOTO 9999
      END IF


      if (debug) then 
        write(0,*) 'Check from DJADSM'
        do k=1,m
          write(0,*) k, b(k,1),c(k,1)
        enddo
      endif

      IF (TDIAG.EQ.'L') THEN
         DO I = 1, N
            DO K = 1, M
               C(K,I) = C(K,I)*VDIAG(K)
            ENDDO
         ENDDO
      END IF
c      write(*,*) 'exit djadsm'
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
