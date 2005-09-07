

C     Covert matrix from JAD format to COO Format
C     
      SUBROUTINE DJDCOX(TRANS,M,N,DESCRA,AR,IA,JA,KA,NG,IPERM,INFO,
     *     IP1,DESCRN,ARN,IA1N,IA2N,INFON,IP2,LARN,LIA1N,
     *     LIA2N,AUX,LAUX,IERROR)

      IMPLICIT NONE
      INCLUDE  'sparker.fh'

C     
C     .. Scalar Arguments ..
      INTEGER            NG, LARN, LAUX, LIA1N, LIA2N, M, N, IERROR
      CHARACTER          TRANS,UNITD
C     .. Array Arguments ..
      DOUBLE PRECISION   AR(*), ARN(*) 
      INTEGER            AUX(0:LAUX/2-1),IPERM(*)
      INTEGER            IA(3,*), JA(*), KA(*), INFO(*), IA1N(*), 
     *     IA2N(*), INFON(*), IP1(*), IP2(*)
      CHARACTER          DESCRA*11, DESCRN*11
C     .. Local Scalars .. 
      INTEGER            IPX, IPG, NNZ, K, ROW, 
     *     I, J, NZL, IRET
      LOGICAL            SCALE
      logical     debug
      parameter   (debug=.false.)
c     .. Local Arrays ..
      CHARACTER*20       NAME
      INTEGER            INT_VAL(5)

C     
C     .. Executable Statements ..
C     
      NAME = 'DJDCOX\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

      IF (TRANS.EQ.'N') THEN
C     SCALE  = (UNITD.EQ.'L') ! meaningless
         IP1(1) = 0
         IP2(1) = 0

         IF (IPERM(1).NE.0) THEN
            DO I = 1, M
               AUX(IPERM(I)) = I
            ENDDO
         ENDIF
         
         NNZ = JA(IA(2,NG+1)-1 +1)-1
         
         if (debug) then 
            write(0,*) 'On entry to DJDCOX: NNZ LAUX ',
     +           nnz,laux,larn,lia1n,lia2n
         endif
         IF (LAUX.LT.NNZ+2) THEN
            IERROR = 60
            INT_VAL(1) = 23
            INT_VAL(2) = NNZ+2
            INT_VAL(3) = LAUX
         ELSE IF (LARN.LT.NNZ) THEN
            IERROR = 60
            INT_VAL(1) = 19
            INT_VAL(2) = NNZ+2
            INT_VAL(3) = LAUX
         ELSE IF (LIA1N.LT.NNZ) THEN
            IERROR = 60
            INT_VAL(1) = 20
            INT_VAL(2) = NNZ+2
            INT_VAL(3) = LAUX
         ELSE IF (LIA2N.LT.NNZ) THEN
            IERROR = 60
            INT_VAL(1) = 21
            INT_VAL(2) = NNZ+2
            INT_VAL(3) = LAUX
         ENDIF
         
         IF (DESCRA(1:1).EQ.'G') THEN
            
            DO 200 IPG = 1, NG                                                  
               DO 50 K = IA(2,IPG), IA(3,IPG)-1                                
                  IPX = IA(1,IPG)                                              
                  DO 40 I = JA(K), JA(K+1) - 1                                 
                     ARN(I)  = AR(I) 
                     IA1N(I) = AUX(IPX)
                     IA2N(I) = KA(I)                 
                     IPX = IPX + 1                                    
 40               CONTINUE                                            
 50            CONTINUE                                                            
               
               IPX = IA(1,IPG)                                        
               DO 70 K = IA(3,IPG), IA(2,IPG+1)-1                     
                  DO 60 I = JA(K), JA(K+1) - 1                        
                     ARN(I)  = AR(I) 
                     IA1N(I) = AUX(IPX)
                     IA2N(I) = KA(I)                 
 60               CONTINUE                                            
                  IPX = IPX + 1                                       
 70            CONTINUE                                               
 200        CONTINUE                                     

            
            
C     .... Order with key IA1N....
            CALL MRGSRT(NNZ,IA1N,AUX,IRET)
            IF (IRET.EQ.0) CALL REORDVN(NNZ,ARN,IA1N,IA2N,AUX)           
            
C     .... Order with key IA2N ...
            I    = 1
            J    = I
            DO WHILE (I.LE.NNZ)
               DO WHILE ((IA1N(J).EQ.IA1N(I)).AND.
     +              (J.LE.NNZ))
                  J = J+1
               ENDDO
               NZL = J - I
               CALL MRGSRT(NZL,IA2N(I),AUX,IRET)
               IF (IRET.EQ.0) CALL REORDVN(NZL,ARN(I),IA1N(I),IA2N(I),
     +              AUX)
               I = J
            ENDDO
            INFON(1)=nnz

         ELSE IF (DESCRA(1:1).EQ.'S' .AND. DESCRA(2:2).EQ.'U') THEN

            DO 20 K = 1, M
               IP2(K) = K
 20         CONTINUE

         ELSE IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'U') THEN

         ELSE IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'L') THEN

         END IF
C     
      ELSE IF (TRANS.NE.'N') THEN 
C     
C     TO DO
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

      







