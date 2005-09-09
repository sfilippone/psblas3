      SUBROUTINE DCRCO(TRANS,M,N,UNITD,D,DESCRA,AR,IA1,IA2,INFO,
     *   IP1,DESCRN,ARN,IAN1,IAN2,INFON,IP2,LARN,LIAN1,
     *   LIAN2,AUX,LAUX,IERROR)

      IMPLICIT NONE
      INCLUDE  'psb_const.fh'

C
C     .. Scalar Arguments ..
      INTEGER            LARN, LAUX, LIAN1, LIAN2, M, N, IERROR
      CHARACTER          TRANS,UNITD
C     .. Array Arguments ..
      DOUBLE PRECISION   AR(*), ARN(*), D(*), AUX(LAUX)
      INTEGER            IA1(*), IA2(*), INFO(*), IAN1(*), IAN2(*),
     *   INFON(*), IP1(*), IP2(*)
      CHARACTER          DESCRA*11, DESCRN*11
C     .. Local Scalars ..
      INTEGER            NNZ, K, ROW, J
      INTEGER            ELEM, ERR_ACT
      LOGICAL            SCALE
      INTEGER MAX_NNZERO
c     .. Local Arrays ..
      CHARACTER*20       NAME
      INTEGER            INT_VAL(5)

C     .. External Subroutines ..
      EXTERNAL           MAX_NNZERO
C     .. Executable Statements ..
C

      NAME = 'DCRCO\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

      IF (TRANS.EQ.'N') THEN
         SCALE  = (UNITD.EQ.'L') ! meaningless
         IP1(1) = 0
         IP2(1) = 0
         NNZ = IA2(M+1)-1
         IF (LARN.LT.NNZ) THEN
          IERROR = 60
          INT_VAL(1) = 18
          INT_VAL(2) = NNZ
          INT_VAL(3) = LARN
         ELSE IF (LIAN1.LT.NNZ) THEN
          IERROR = 60
          INT_VAL(1) = 19
          INT_VAL(2) = NNZ
          INT_VAL(3) = LIAN1
         ELSE IF (LIAN2.LT.NNZ) THEN
          IERROR = 60
          INT_VAL(1) = 20
          INT_VAL(2) = NNZ
          INT_VAL(3) = LIAN2
         ENDIF
         
         IF(IERROR.NE.0) THEN
            CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
            GOTO 9999
         END IF
         
         IF (DESCRA(1:1).EQ.'G') THEN
C        ... Construct COO Representation...
            ELEM = 1

            DO ROW = 1, M
               DO J = IA2(ROW), IA2(ROW+1)-1
                  IAN1(ELEM) = ROW
                  IAN2(ELEM) = IA1(J)
                  ARN(ELEM) = AR(J)
                  ELEM = ELEM + 1
               ENDDO
            ENDDO
            INFON(1) = IA2(M+1)-1
         ELSE IF (DESCRA(1:1).EQ.'S' .AND. DESCRA(2:2).EQ.'U') THEN

            DO 20 K = 1, M
               IP2(K) = K
 20         CONTINUE

c$$$            CALL DVSSG(M,IA1,IA2,IP2,IAN2(PNG),IP1,IP2,AUX(IWLEN),
c$$$     *                 AUX(IWORK1))
c$$$            CALL DVSMR(M,AR,IA1,IA2,IAN2(PNG),AUX(IWLEN),IP1,IP2,
c$$$     *                 IAN2(PIA),IAN2(PJA),IAN1,ARN,AUX(IWORK1),
c$$$     *                 AUX(IWORK2),NJA,IER,SCALE)
C
         ELSE IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'U') THEN

c$$$            CALL DVTFG('U',M,IA1,IA2,IAN2(PNG),IP1,IP2,AUX(IWLEN),
c$$$c    *                 AUX(IWORK1),AUX(IWORK2),IAN1(M+1))
c$$$     *                 AUX(IWORK1),IAN1(1),IAN1(M+5))
c$$$            CALL DVTMR(M,AR,IA1,IA2,ISTROW,IAN2(PNG),AUX(IWLEN),IP1,IP2,
c$$$     *                 IAN2(PIA),IAN2(PJA),IAN1,ARN,NJA,IER,SCALE)
C

         ELSE IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'L') THEN

c$$$            CALL DVTFG('L',M,IA1,IA2,IAN2(PNG),IP1,IP2,AUX(IWLEN),
c$$$c     *                 AUX(IWORK1),AUX(IWORK2),IAN1(M+1))
c$$$     *                 AUX(IWORK1),IAN1(1),IAN1(M+5))
c$$$            CALL DVTMR(M,AR,IA1,IA2,ISTROW,IAN2(PNG),AUX(IWLEN),IP1,IP2,
c$$$     *                 IAN2(PIA),IAN2(PJA),IAN1,ARN,NJA,IER,SCALE)

         END IF
C
      ELSE IF (TRANS.NE.'N') THEN
C
C           TO DO
C     
         IERROR = 3021
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

      RETURN
      END
