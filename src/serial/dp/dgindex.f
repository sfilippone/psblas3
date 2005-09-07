      SUBROUTINE DGINDEX(M,N_BLOCKS,A,IA1,IA2,ARN,KA,IA,JA, INFON,
     +     LARN,LKA,LJA,IPERM,WORK, LWORK, SIZE_REQ, IERROR)

      IMPLICIT NONE
      INCLUDE  'sparker.fh'

C     ... Scalar arguments ...
      INTEGER          M, LWORK,N_BLOCKS,LARN,LKA,LJA,
     +     SIZE_REQ,IERROR

C     ... Array arguments ...

      DOUBLE PRECISION A(*), ARN(*)
      INTEGER          IA1(*), IA2(*), KA(*), 
     +     IA(3,*), IPERM(*), JA(*), WORK(*),INFON(*)
      
C     .... Local scalars ...
      INTEGER          I, J, BLOCK, ROW, COL, POINT_AR, POINT_JA, IP1,
     +     IP2, IPX, NNZ, DIM_BLOCK, LIMIT, IPW,COUNT, IPC,CHECK_FLAG
      LOGICAL          CSR
c     .. Local Arrays ..
      CHARACTER*20       NAME
      INTEGER            INT_VAL(5)

      NAME = 'DGINDEX\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

      POINT_AR = 1
      POINT_JA = 0
      CHECK_FLAG=IBITS(INFON(UPD_),1,2)

      IF ((LARN.LT.POINT_AR).OR.(LKA.LT.POINT_AR)) THEN
         IERROR = 60
         INT_VAL(1) = 11
         INT_VAL(2) = POINT_AR
         INT_VAL(3) = LARN
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      ENDIF

      NNZ = IA2(M + 1) - 1
      COUNT            = 0

C     .... Invert Permutation Matrix...
      IF (IPERM(1).NE.0) THEN
         DO I = 1, M
            WORK(IPERM(I)) = I
         ENDDO
      ENDIF

      IF ( (LKA .GE.( SIZE_REQ ))  
     +     .AND. (LWORK .GE. (M + NNZ+2))) THEN
C     
C     Prepare for smart regeneration
C     

         IPW              = M + 2
         IP1              = (LKA-IREG_FLGS-2)/2
         IP2              = IP1+IREG_FLGS
         IPC              = IP2 + NNZ + 1
         KA(IP1 + IPC_)   = IPC
         KA(IP1+IP2_)     = IP2
         INFON(UPD_PNT_)  = IP1
         KA(IP1+IFLAG_)   = CHECK_FLAG
         KA(IP1+NNZT_)    = NNZ
         KA(IP1+NNZ_)     = 0
         KA(IP1+ICHK_)    = NNZ+CHECK_FLAG
         I                = M+2
         IPX              = IA2(I+IP2_)

C     Invert permutation for smart regeneration

         DO I = 1, NNZ
            WORK(IPW + IA2(IPX + I -1) - 1) = I
         ENDDO

C     Construct JAD matrix...

         DO BLOCK = 1, N_BLOCKS
            COL = 1
            DIM_BLOCK = IA(1,BLOCK+1)-IA(1,BLOCK)
c$$$  write(0,*) 'DGINDEX: BLOCK LOOP ',block,n_blocks,dim_block
            if (dim_block .gt. maxjdrows) then 
               write(0,*) 'Wrong value for dim_block',block,
     +              IA(1,BLOCK+1),IA(1,BLOCK)
               return
            endif
            LIMIT = INT(DIM_BLOCK*PERCENT)
            POINT_JA = POINT_JA+1
            IF (LJA.LT.POINT_JA) THEN
               IERROR = 60
               INT_VAL(1) = 13
               INT_VAL(2) = POINT_JA
               INT_VAL(3) = LJA
               CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
               GOTO 9999
            ENDIF
            
            IA(2,BLOCK) = POINT_JA
            JA(POINT_JA) = POINT_AR
            CSR = .FALSE.
            
            IF (DIM_BLOCK.NE.0) THEN
C     ... If current block is not empty ...
C     ... For each Column belonging to Block ...
               DO WHILE(.TRUE.) 
C     ... For each row belonging to the block BLOCK ...
                  DO I = IA(1,BLOCK), IA(1,BLOCK+1)-1
                     IF (IPERM(1).EQ.0) THEN
                        ROW = I
                     ELSE
                        ROW = WORK(I)
                     ENDIF
                     
C     ... If the current row is too short ...
                     IF (IA2(ROW)+COL-1.GE.IA2(ROW+1)) THEN
C     ... Switch to CSR representation ...
                        IF (I.LE.IA(1,BLOCK)+LIMIT) THEN
                           CSR=.TRUE.
                           POINT_AR = POINT_AR - I + IA(1,BLOCK)
                           GOTO 998
                        ELSE

                           COUNT = COUNT + 1
                           ARN(POINT_AR) = 0.D0                           
                           KA(POINT_AR) = KA(POINT_AR-1)                           
                           IF(POINT_AR.LT.IP1) THEN
                              KA(IPC + COUNT -1) = POINT_AR
                           ENDIF
C     
C     The following statement assumes that we never get here with POINT_AR=1
C     
                           POINT_AR = POINT_AR+1
                           IF ((LARN.LT.POINT_AR).OR.(LKA.LT.POINT_AR))
     +                          THEN
                              IERROR = 60
                              INT_VAL(1) = 11
                              INT_VAL(2) = POINT_AR
                              INT_VAL(3) = LARN
                              CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                              GOTO 9999
                           ENDIF
                        ENDIF
                     ELSE
                        ARN(POINT_AR) = A(IA2(ROW)+COL-1)
                        KA(POINT_AR) = IA1(IA2(ROW)+COL-1)
                        IF(POINT_AR.LT.IP1) THEN
                           KA(IP2 +  WORK(IPW  + 
     +                          IA2(ROW) +COL -1-1)-1) = POINT_AR
                        ENDIF
                        
                        POINT_AR = POINT_AR+1
                        IF ((LARN.LT.POINT_AR).OR.(LKA.LT.POINT_AR)) 
     +                       THEN
                              IERROR = 60
                              INT_VAL(1) = 11
                              INT_VAL(2) = POINT_AR
                              INT_VAL(3) = LARN
                              CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                              GOTO 9999
                        ENDIF
                     ENDIF
                  ENDDO
                  
                  IF (CSR) GOTO 998
                  
                  IF (LJA.LT.POINT_JA+COL) THEN
                     IERROR = 60
                     INT_VAL(1) = 13
                     INT_VAL(2) = POINT_JA
                     INT_VAL(3) = LJA
                     CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                     GOTO 9999
                  ENDIF
                  IF(IERROR.EQ.0) THEN
                     JA(POINT_JA+COL) = POINT_AR
                  ENDIF
                  COL = COL+1
               ENDDO
               
            ENDIF
 998        CONTINUE 
            
            POINT_JA = POINT_JA+COL-1
            
            IF (LJA.LT.POINT_JA) THEN
               IERROR = 60
               INT_VAL(1) = 13
               INT_VAL(2) = POINT_JA
               INT_VAL(3) = LJA
               CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
               GOTO 9999
            ENDIF
            
            IA(3,BLOCK) = POINT_JA
            
C     ... Start CSR Format ...
            
C     ... For each row belonging to the block BLOCK ...
            DO I = IA(1,BLOCK), IA(1,BLOCK+1)-1
               IF (IPERM(1).EQ.0) THEN
                  ROW = I
               ELSE
                  ROW = WORK(I)
               ENDIF
               
C     ... For each nnzero elements belonging to current row ...
               DO J = IA2(ROW)+COL-1, IA2(ROW+1)-1
                  
                  ARN(POINT_AR) = A(J)
                  KA (POINT_AR) = IA1(J)
                  IF (POINT_AR.LT.IP1) THEN
                     KA(IP2 +  WORK(IPW  + J-1)-1) = POINT_AR
                  ENDIF

                  POINT_AR = POINT_AR+1
                  IF ((LARN.LT.POINT_AR).OR.(LKA.LT.POINT_AR)) THEN
                     IERROR = 60
                     INT_VAL(1) = 11
                     INT_VAL(2) = POINT_AR
                     INT_VAL(3) = LARN
                     CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                     GOTO 9999
                  ENDIF
               ENDDO
               
               POINT_JA = POINT_JA+1
               IF (LJA.LT.POINT_JA) THEN
                  IERROR = 60
                  INT_VAL(1) = 13
                  INT_VAL(2) = POINT_JA
                  INT_VAL(3) = LJA
                  CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                  GOTO 9999
               ENDIF
               JA(POINT_JA) = POINT_AR
            ENDDO
         ENDDO
         
         
      ELSE
c$$$  c      write(*,*)'inizio a ciclare sui blocchi'
         DO BLOCK = 1, N_BLOCKS
            COL = 1
            DIM_BLOCK = IA(1,BLOCK+1)-IA(1,BLOCK)
            LIMIT = INT(DIM_BLOCK*PERCENT)
            POINT_JA = POINT_JA+1
            IF (LJA.LT.POINT_JA) THEN
               IERROR = 60
               INT_VAL(1) = 13
               INT_VAL(2) = POINT_JA
               INT_VAL(3) = LJA
               CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
               GOTO 9999
            ENDIF
            
            IA(2,BLOCK) = POINT_JA
            JA(POINT_JA) = POINT_AR
            CSR = .FALSE.
            
            IF (DIM_BLOCK.NE.0) THEN
C     ... If current block is not empty ...
C     ... For each Column belonging to Block ...
               DO WHILE(.TRUE.) 
C     ... For each row belonging to the block BLOCK ...
                  DO I = IA(1,BLOCK), IA(1,BLOCK+1)-1
                     IF (IPERM(1).EQ.0) THEN
                        ROW = I
                     ELSE
                        ROW = WORK(I)
                     ENDIF
                     
C     ... If the current row is too short ...
                     IF (IA2(ROW)+COL-1.GE.IA2(ROW+1)) THEN
C     ... Switch to CSR representation ...
                        IF (I.LE.IA(1,BLOCK)+LIMIT) THEN
                           CSR=.TRUE.
                           POINT_AR = POINT_AR - I + IA(1,BLOCK)
                           GOTO 999
                        ELSE
                           COUNT= COUNT+1
                           ARN(POINT_AR) = 0.D0
                           KA (POINT_AR) = KA(POINT_AR-1)
C     
C     The following statement assumes that we never get here with POINT_AR=1
C     
                           POINT_AR = POINT_AR+1
                           IF ((LARN.LT.POINT_AR).OR.(LKA.LT.POINT_AR))
     +                          THEN
                              IERROR = 60
                              INT_VAL(1) = 11
                              INT_VAL(2) = POINT_AR
                              INT_VAL(3) = LARN
                              CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                              GOTO 9999
                           ENDIF
                        ENDIF
                     ELSE
                        ARN(POINT_AR) = A(IA2(ROW)+COL-1)
                        KA (POINT_AR) = IA1(IA2(ROW)+COL-1)
                        POINT_AR = POINT_AR+1
                        IF ((LARN.LT.POINT_AR).OR.(LKA.LT.POINT_AR)) 
     +                       THEN
                           IERROR = 60
                           INT_VAL(1) = 11
                           INT_VAL(2) = POINT_AR
                           INT_VAL(3) = LARN
                           CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                           GOTO 9999
                        ENDIF
                     ENDIF
                  ENDDO
                  
                  IF (CSR) GOTO 999
                  
                  IF (LJA.LT.POINT_JA+COL) THEN
                     IERROR = 60
                     INT_VAL(1) = 13
                     INT_VAL(2) = POINT_JA
                     INT_VAL(3) = LJA
                     CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                     GOTO 9999
                  ENDIF
                  
                  JA(POINT_JA+COL) = POINT_AR
                  COL = COL+1
               ENDDO
               
            ENDIF
 999        CONTINUE 
            
            POINT_JA = POINT_JA+COL-1
            
            IF (LJA.LT.POINT_JA) THEN
               IERROR = 60
               INT_VAL(1) = 13
               INT_VAL(2) = POINT_JA
               INT_VAL(3) = LJA
               CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
               GOTO 9999
            ENDIF
            
            IA(3,BLOCK) = POINT_JA
            
C     ... Start CSR Format ...
            
C     ... For each row belonging to the block BLOCK ...
            DO I = IA(1,BLOCK), IA(1,BLOCK+1)-1
               IF (IPERM(1).EQ.0) THEN
                  ROW = I
               ELSE
                  ROW = WORK(I)
               ENDIF
               
C     ... For each nnzero elements belonging to current row ...
               DO J = IA2(ROW)+COL-1, IA2(ROW+1)-1
                  ARN(POINT_AR) = A(J)
                  KA (POINT_AR) = IA1(J)
                  POINT_AR = POINT_AR+1
                  IF ((LARN.LT.POINT_AR).OR.(LKA.LT.POINT_AR)) THEN
                     IERROR = 60
                     INT_VAL(1) = 11
                     INT_VAL(2) = POINT_AR
                     INT_VAL(3) = LARN
                     CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                     GOTO 9999
                  ENDIF
               ENDDO
               
               POINT_JA = POINT_JA+1
               IF (LJA.LT.POINT_JA) THEN
                  IERROR = 60
                  INT_VAL(1) = 13
                  INT_VAL(2) = POINT_JA
                  INT_VAL(3) = LJA
                  CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                  GOTO 9999
               ENDIF
               JA(POINT_JA) = POINT_AR
            ENDDO
         ENDDO
         
         
         
         
      ENDIF
      IA(2,N_BLOCKS+1) = POINT_JA
      KA(IP1 + ZERO_) = COUNT

      IF(POINT_AR.GE.IP1) THEN
         SIZE_REQ=NNZ+COUNT
      ELSE
         SIZE_REQ=0
      ENDIF
      infon(1)=point_ar-1

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

