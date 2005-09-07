C     Covert matrix from COO format to COO Format
C
      SUBROUTINE DCOCO(TRANS,M,N,UNITD,D,DESCRA,AR,IA1,IA2,INFO,
     *  P1,DESCRN,ARN,IA1N,IA2N,INFON,P2,LARN,LIA1N,
     *  LIA2N,AUX,LAUX,IERROR)

      IMPLICIT NONE
      INCLUDE  'sparker.fh'

C     .. Scalar Arguments ..
      INTEGER            LARN, LAUX, LIA1N, LIA2N, 
     +  M, N, IERROR
      CHARACTER          TRANS,UNITD
C     .. Array Arguments ..
      DOUBLE PRECISION   AR(*), ARN(*), D(*)
      INTEGER            AUX(0:LAUX-1)
      INTEGER            IA1(*), IA2(*), INFO(*), IA1N(*), IA2N(*),
     *  INFON(*), P1(*), P2(*)
      CHARACTER          DESCRA*11, DESCRN*11
C     .. Local Scalars ..
      INTEGER            IPX, IP1, IP2, CHECK_FLAG
      INTEGER            NNZ, K, I, J, NZL, IRET
      INTEGER            ELEM_IN, ELEM_OUT
      LOGICAL            SCALE
      INTEGER MAX_NNZERO
      logical     debug
      parameter   (debug=.false.)
c     .. Local Arrays ..
      CHARACTER*20       NAME
      INTEGER            INT_VAL(5)
C
C     ...Common variables...
C     This flag describe the action to do
      
C     .. External Subroutines ..
      EXTERNAL           MAX_NNZERO
C     .. Executable Statements ..
C

      NAME = 'DCOCO\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

      CHECK_FLAG=IBITS(info(upd_),1,2)
      IF (TRANS.EQ.'N') THEN
        SCALE  = (UNITD.EQ.'L') ! meaningless
        P1(1) = 0
        P2(1) = 0

        NNZ = INFO(nnz_)
        if (debug) then 
          write(*,*) 'On entry to DCOCO: NNZ LAUX ',
     +      nnz,laux,larn,lia1n,lia2n
        endif
        IF (LAUX.LT.NNZ+2) THEN
          IERROR = 60
          INT_VAL(1) = 22
          INT_VAL(2) = NNZ+2
          INT_VAL(3) = LAUX
        ELSE IF (LARN.LT.NNZ) THEN
          IERROR = 60
          INT_VAL(1) = 18
          INT_VAL(2) = NNZ+2
          INT_VAL(3) = LAUX
        ELSE IF (LIA1N.LT.NNZ) THEN
          IERROR = 60
          INT_VAL(1) = 19
          INT_VAL(2) = NNZ+2
          INT_VAL(3) = LAUX
        ELSE IF (LIA2N.LT.M+1) THEN
          IERROR = 60
          INT_VAL(1) = 20
          INT_VAL(2) = NNZ+2
          INT_VAL(3) = LAUX
        ENDIF
        
C
C     Error handling
C
        IF(IERROR.NE.0) THEN
           CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
           GOTO 9999
        END IF

        IF (DESCRA(1:1).EQ.'G') THEN
C
C     Sort COO data structure
C     
          if (debug) write(*,*)'First sort',nnz
          do k=1, nnz
            arn(k)  = ar(k)
            ia1n(k) = ia1(k)
            ia2n(k) = ia2(k)
          enddo
          
          if (debug) write(*,*)'Second sort'            

          if ((lia2n.ge.(2*nnz+ireg_flgs+1))
     +      .and.(laux.ge.2*(2+nnz))) then 
C     
C     Prepare for smart regeneration
c     
            ipx = nnz+3            
            do i=1, nnz
              aux(ipx+i-1) = i
            enddo
            ip1              = nnz+2
            infon(upd_pnt_)  = ip1
            ip2              = ip1+ireg_flgs
            ia2n(ip1+ip2_)   = ip2
            ia2n(ip1+iflag_) = check_flag
            ia2n(ip1+nnzt_)  = nnz
            ia2n(ip1+nnz_)   = 0
            ia2n(ip1+ichk_)  = nnz+check_flag
            if (debug) write(0,*) 'Build check :',ia2n(ip1+nnzt_) 
            
C     .... Order with key IA1N ...
            CALL MRGSRT(NNZ,IA1N,AUX,IRET)
            IF (IRET.EQ.0) CALL REORDVN3(NNZ,ARN,IA1N,IA2N,AUX(IPX),AUX)
C     .... Order with key IA2N ...
            
            I    = 1
            J    = I
            DO WHILE (I.LE.NNZ)
              DO WHILE ((IA1N(J).EQ.IA1N(I)).AND.
     +          (J.LE.NNZ))
                J = J+1
              ENDDO
              NZL = J - I
              CALL MRGSRT(NZL,IA2N(I),AUX,IRET)
              IF (IRET.EQ.0) CALL REORDVN3(NZL,ARN(I),IA1N(I),IA2N(I),
     +          AUX(IPX+I-1),AUX)
              I = J
            ENDDO
            
            ia2n(ip2+aux(ipx+1-1)-1) = 1

C     ... Construct final COO  Representation...
            ELEM_OUT = 1
C     ... Insert remaining element ...
            DO ELEM_IN  = 2, NNZ
              IF ((IA1N(ELEM_IN).EQ.IA1N(ELEM_OUT)).AND.
     +          (IA2N(ELEM_IN).EQ.IA2N(ELEM_OUT))) THEN 
                IF (CHECK_FLAG.EQ.1) THEN
C     ... Error, there are duplicated elements ...
                  IERROR = 130
                  CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                  GOTO 9999
                ELSE IF (CHECK_FLAG.EQ.2) THEN
C     ... Insert only the first duplicated element ...
                  ia2n(ip2+aux(ipx+elem_in-1)-1) = elem_out
                ELSE IF (CHECK_FLAG.EQ.3) THEN
C     ... Sum the duplicated element ...
                  ARN(ELEM_OUT) = ARN(ELEM_OUT) + ARN(ELEM_IN)
                  ia2n(ip2+aux(ipx+elem_in-1)-1) = elem_out
                END IF
              ELSE
                ELEM_OUT = ELEM_OUT + 1
                ARN(ELEM_OUT)  = ARN(ELEM_IN)
                ia2n(ip2+aux(ipx+elem_in-1)-1) = elem_out
                IA1N(ELEM_OUT) = IA1N(ELEM_IN)
                IA2N(ELEM_OUT) = IA2N(ELEM_IN)
              ENDIF
            ENDDO
            
          ELSE
            
C     .... Order with key IA1N ...
            CALL MRGSRT(NNZ,IA1N,AUX,IRET)
            IF (IRET.EQ.0) CALL REORDVN(NNZ,ARN,IA1N,IA2N,AUX)
C     .... Order with key IA2N ...
            
            I    = 1
            J    = I
            DO WHILE (I.LE.NNZ)
              DO WHILE ((IA1N(J).EQ.IA1N(I)).AND.
     +          (J.LE.NNZ))
                J = J+1
              ENDDO
              NZL = J - I
              CALL MRGSRT(NZL,IA2N(I),AUX,IRET)
              IF (IRET.EQ.0) CALL REORDVN(NZL,ARN(I),IA1N(I),IA2N(I),
     +          AUX)
              I = J
            ENDDO
C     ... Construct final COO  Representation...
            ELEM_OUT = 1
C     ... Insert remaining element ...
            DO ELEM_IN  = 2, NNZ
              IF ((IA1N(ELEM_IN).EQ.IA1N(ELEM_OUT)).AND.
     +          (IA2N(ELEM_IN).EQ.IA2N(ELEM_OUT))) THEN 
                IF (CHECK_FLAG.EQ.1) THEN
C     ... Error, there are duplicated elements ...
                  IERROR = 130
                  CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                  GOTO 9999
                ELSE IF (CHECK_FLAG.EQ.2) THEN
C     ... Insert only the first duplicated element ...
                ELSE IF (CHECK_FLAG.EQ.3) THEN
C     ... Sum the duplicated element ...
                  ARN(ELEM_OUT) = ARN(ELEM_OUT) + ARN(ELEM_IN)
                END IF
              ELSE
                ELEM_OUT = ELEM_OUT + 1
                ARN(ELEM_OUT)  = ARN(ELEM_IN)
                IA1N(ELEM_OUT) = IA1N(ELEM_IN)
                IA2N(ELEM_OUT) = IA2N(ELEM_IN)
              ENDIF
            ENDDO
          ENDIF
          INFON(nnz_)  = ELEM_OUT
          infon(srtd_) = isrtdcoo
          
          if (debug) write(*,*)'Done Rebuild COO',infon(1)
          
        ELSE IF (DESCRA(1:1).EQ.'S' .AND. DESCRA(2:2).EQ.'U') THEN

          DO 20 K = 1, M
            P2(K) = K
 20       CONTINUE

        ELSE IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'U') THEN

        ELSE IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'L') THEN
           
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
