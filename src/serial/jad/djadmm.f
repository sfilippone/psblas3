      SUBROUTINE DJADMM(TRANSA,M,K,N,ALPHA,DESCRA,AR,
     *   JA,IA,B,LDB,BETA,C,LDC,WORK,IERROR)

      IMPLICIT NONE
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, BETA
      INTEGER           K, LDB, LDC, M, N,IERROR
      CHARACTER         TRANSA
C     .. Array Arguments ..
      DOUBLE PRECISION  AR(*), B(LDB,*), C(LDC,*),  WORK(*)
      INTEGER           IA(*), JA(*)
      CHARACTER         DESCRA*11
C     .. Local Scalars ..
      integer           PIA, PJA, PNG
      INTEGER           I, J, KB,NB, ERR_ACT
      CHARACTER         DIAG, TRANS
c     .. Local Arrays ..
      CHARACTER*20       NAME
      INTEGER            INT_VAL(5)
C     .. Executable Statements ..
C

      NAME = 'DJADMM\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

      IF (DESCRA(1:1).EQ.'G')  TRANS = TRANSA
      IF (DESCRA(1:1).EQ.'S' .AND. DESCRA(2:2).EQ.'U') TRANS = 'U'
      IF (DESCRA(1:1).EQ.'S' .AND. DESCRA(2:2).EQ.'L') TRANS = 'L'
      IF (DESCRA(1:1).EQ.'D') THEN
        IF (DESCRA(3:3).EQ.'U') THEN
          DO 40 I = 1, K
            DO 20 J = 1, M
              C(J,I) = BETA*C(J,I) + ALPHA*B(J,I)
 20         CONTINUE
 40       CONTINUE
        ELSE
          DO 80 I = 1, K
            DO 60 J = 1, M
              C(J,I) = BETA*C(J,I) + ALPHA*AR(J)*B(J,I)
 60         CONTINUE
 80       CONTINUE
        END IF
        RETURN
      END IF
      IF (TRANS.EQ.'T'.or.TRANS.EQ.'C') THEN
        IERROR = 3015
        RETURN
      ENDIF
C
      IF (DESCRA(3:3).EQ.'N') DIAG = 'N'
      IF (DESCRA(3:3).EQ.'U') DIAG = 'U'
      PNG = IA(1)
      PIA = IA(2)
      PJA = IA(3)
*      write(6,*)'M N K',m,n,k

      NB = 4
      KB = MOD(K,NB) 
      
      SELECT CASE(KB) 
      CASE (1) 
        IF (K==1) THEN 
          CALL DJADMV(DIAG,M,N,ALPHA,IA(PNG),
     +      AR,JA,IA(PJA),IA(PIA),B(1,1),BETA,C(1,1),IERROR)
        ELSE
          CALL DJADMV2(DIAG,M,N,ALPHA,IA(PNG),
     +      AR,JA,IA(PJA),IA(PIA),B(1,1),LDB,BETA,C(1,1),LDC,IERROR)
          CALL DJADMV3(DIAG,M,N,ALPHA,IA(PNG),
     +      AR,JA,IA(PJA),IA(PIA),B(1,3),LDB,BETA,C(1,3),LDC,IERROR)
          KB = KB + 4
        ENDIF
      CASE(2)
        CALL DJADMV2(DIAG,M,N,ALPHA,IA(PNG),
     +    AR,JA,IA(PJA),IA(PIA),B(1,1),LDB,BETA,C(1,1),LDC,IERROR)
      CASE(3)
        CALL DJADMV3(DIAG,M,N,ALPHA,IA(PNG),
     +    AR,JA,IA(PJA),IA(PIA),B(1,1),LDB,BETA,C(1,1),LDC,IERROR)
        
      END SELECT

      IF(IERROR.NE.0) THEN
         INT_VAL(1)=IERROR
         CALL FCPSB_ERRPUSH(4012,NAME,INT_VAL)
         GOTO 9999
      END IF

      DO I=KB+1,K,NB
        CALL DJADMV4(DIAG,M,N,ALPHA,IA(PNG),
     +    AR,JA,IA(PJA),IA(PIA),B(1,I),LDB,BETA,C(1,I),LDC,IERROR)
      END DO
      
      IF(IERROR.NE.0) THEN
         INT_VAL(1)=IERROR
         CALL FCPSB_ERRPUSH(4012,NAME,INT_VAL)
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

      RETURN
      END
