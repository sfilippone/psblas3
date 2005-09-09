      SUBROUTINE DJDCO(TRANS,M,N,DESCRA,AR,IA1,IA2,IPERM,INFO,
     *     IP1,DESCRN,ARN,IA1N,IA2N,INFON,IP2,LARN,LIA1N,
     *     LIA2N,AUX,LAUX,IERROR)      
      IMPLICIT NONE
      INCLUDE  'psb_const.fh'      
C     
C     .. Scalar Arguments ..
      INTEGER            LARN, LAUX, LIA1N, LIA2N, M, N, IERROR
      CHARACTER          TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION   AR(*), ARN(*)
      INTEGER            AUX(0:LAUX-1),IPERM(*)
      INTEGER            IA1(*), IA2(*), INFO(*), IA1N(*), 
     *     IA2N(*), INFON(*), IP1(*), IP2(*)
      CHARACTER          DESCRA*11, DESCRN*11
C     .. Local Scalars .. 
      INTEGER            PIA, PJA, PNG, ERR_ACT
      logical     debug
      parameter   (debug=.false.)
c     .. Local Arrays ..
      CHARACTER*20       NAME
      INTEGER            INT_VAL(5)

      NAME = 'DJDCO\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)
            
      PNG = IA2(1)
      PIA = IA2(2)
      PJA = IA2(3)

       if(debug) write(*,*) 'On entry to DJDCO: NNZ LAUX ',
     +     info(1),laux,larn,lia1n,lia2n
      
      CALL DJDCOX(TRANS,M,N,DESCRA,AR,IA2(PIA),IA2(PJA),
     *  IA1,IA2(PNG),IPERM, INFO, IP1,DESCRN,ARN,IA1N,IA2N,INFON,
     *  IP2,LARN,LIA1N, LIA2N,AUX,LAUX,IERROR)
        IF(IERROR.NE.0) THEN
           IERROR=4011
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
      
