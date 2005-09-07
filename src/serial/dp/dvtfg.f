C                                                
C  Compute level numbers for the triangular matrix. 
C                                                                               
      SUBROUTINE DVTFG (UPLO,M,JA,IA,NG,IPA,IPAT,KLEN,IWORK1,IWORK2,            
     *                  IWORK3)                                                 
C     .. Scalar Arguments ..                                                    
      INTEGER           M, NG                                                   
      CHARACTER         UPLO                                                    
C     .. Array Arguments ..                                                     
      INTEGER           IA(*), IPA(*), IPAT(*), IWORK3(*), JA(*),               
     *                  KLEN(*), IWORK2(*), IWORK1(*)                           
C     .. Local Scalars ..                                                       
      INTEGER           I, J, L, L0, L1, LEV, NP, iret
C     .. Intrinsic Functions ..                                                 
      INTRINSIC         MAX                                                     
      logical debug
      parameter (debug=.false.)
C     .. Executable Statements ..                                               
C                                                                               
      NG = 0                                                                    
C                                                                               
C     CHECK ON THE NUMBERS OF THE ELEMENTS OF THE MATRIX                        
C                                                                               
      IF ((IA(M+1)-1).EQ.0) THEN                                                
C                                                                               
C        THE MATRIX HASN'T ELEMENTS                                             
C        THE OUTPUT PERMUTATIONS ARE POSED TO THE IDENTITY MATRIX               
C                                                                               
         DO 20 I = 1, M                                                         
            IPA(I) = I                                                          
            IPAT(I) = I                                                         
            KLEN(I) = 0                                                         
   20    CONTINUE                                                               
      ELSE                                                                      
C                                                                               
C          COMPUTE LEVEL NUMBER FOR EACH ROWS                                   
C                                                                               
C          IWORK1:   AUXILIARY VECTOR WHICH CONTAINS                            
C                    LEVEL NUMBERS OF EACH ROW
C                                                                               
         DO 40 I = 1, M                                                         
            IWORK1(I) = 0                                                       
            IWORK3(I) = 0                                                    
   40    CONTINUE                                                               
         IF (UPLO.EQ.'L') THEN                                                  
C                                                                               
C           LOWER TRIANGULAR SPARSE MATRIX                                      
C                                                                               
            DO 80 I = 1, M                                                      
               IWORK1(I) = 1                                                    
               DO 60 J = IA(I), IA(I+1) - 1                        
                  IWORK1(I) = MAX(IWORK1(I),IWORK1(JA(J))+1)                    
   60          CONTINUE                                                         
   80       CONTINUE                                                            
         ELSE IF (UPLO.EQ.'U') THEN                                             
C                                                                               
C           UPPER TRIANGULAR SPARSE MATRIX                                      
C                                                                               
            DO 120 I = M, 1, -1                                                 
               IWORK1(I) = 1                                                    
               DO 100 J = IA(I), IA(I+1) - 1                                    
                  IWORK1(I) = MAX(IWORK1(I),IWORK1(JA(J))+1)                    
  100          CONTINUE                                                         
  120       CONTINUE                                                            
         END IF                                                                 
C                                                                               
C          COUNT NUMBER OF ROWS IN EACH EQUIVALENCE GROUPS                      
C                                                                               
C          NOTE: GROUP = SET OF ROWS HAVING THE SAME LEVEL NUMBER               
C                                                                               
C          IWORK3: AUXILIARY VECTOR WHICH CONTAINS                              
C                  THE NUMBER OF ROWS FOR EACH GROUPS                           
C                                                                               
         DO 140 I = 1, M                                                        
            IWORK3(IWORK1(I)) = IWORK3(IWORK1(I)) + 1                           
  140    CONTINUE                                                               
C                                                                               
C        SET UP IWORK2:                                                         
C        IWORK2(I) POINTS TO THE BEGINNING OF I-TH GROUP                        
C                                                                               
         IWORK2(1) = 1                                                          
         DO 160 I = 2, M + 1                                                    
            IWORK2(I) = IWORK3(I-1) + IWORK2(I-1)                               
            IF (IWORK2(I).EQ.M+1) THEN                                          
               NG = I - 1                                                       
               GO TO 180                                                        
            END IF                                                              
  160    CONTINUE                                                               
  180    CONTINUE                                                               
C                                                                               
C        NG : TOTAL NUMBER OF LEVELS    
C        IWORK3: VECTOR CONTAINING THE NUMBER OF THE ROWS SORTED BY             
C                EQUIVALENCE GROUPS.                                            
C                                                                               
         DO 200 I = 1, M                                                        
            IWORK3(IWORK2(IWORK1(I))) = I                                       
            IWORK2(IWORK1(I)) = IWORK2(IWORK1(I)) + 1                           
  200    CONTINUE                                                               
C                                                                               
C        REGENERATE IWORK2: POINTER INTO IWORK3 FOR EQUIVALENCE GROUPS          
C                                                                               
         DO 220 I = NG + 1, 2, -1                                               
            IWORK2(I) = IWORK2(I-1)                                             
  220    CONTINUE                                                               
         IWORK2(1) = 1                                                          
C                                                                               
C        IWORK1: ROWS LENGTH IN NEW ORDERING                                    
C                                                                               
         DO 240 L = 1, M                                                        
            IWORK1(L) = IA(IWORK3(L)) - IA(IWORK3(L)+1)                         
  240    CONTINUE                                                               
C                                                                               
C        SORT ROWS BY DECREASING NUMBER OF NONZERO ELEMENTS.                    
C                                                                               
C        IPA:  VECTOR OF NUMBER OF ROWS SORTED BY EQUIVALENCE                   
C              GROUPS AND BY ROW LENGTH                                         
C        IPAT: AUXILIARY VECTOR NEED TO THE SORTER ROUTINES                     
C                                                                               
         L1 = IWORK2(2) - IWORK2(1)                                             
         DO 260 L = 1, L1                                                       
            IPA(L) = IWORK3(L)                                            
  260    CONTINUE                                                               
         if (debug) write(0,*) 'DVTFG: Group ',1,':',(ipa(l),l=1,l1)
         DO 360 LEV = 2, NG                                                     
C                                                                               
C           LOOP ON GROUPS                                                      
C           L1: LENGTH OF CURRENT GROUP                                         
C           L0: POINTER TO IPA TO THE FIRST LOCATIONS RESERVED                  
C               FOR CURRENT GROUP                                               
C                                                                               
            L1 = IWORK2(LEV+1) - IWORK2(LEV)                                    
            L0 = IWORK2(LEV) - 1                                                
            CALL MRGSRT(L1,IWORK1(IWORK2(LEV)),IPAT,IRET)                        
            IF (IRET.EQ.0) THEN
              NP = IPAT(1)                                                        
              DO 280 L = 1, L1                                                    
                IPA(L0+L) = IWORK3(L0+NP)                                        
                NP = IPAT(1+NP)                                                  
 280          CONTINUE                                              
            ELSE
C                                                                               
C           VECTOR ALREADY SORTED. NO CHANGE IS NEED                            
C                                                                               
              DO 320 L = 1, L1                                                    
                IPA(L0+L) = IWORK3(L0+L)
 320          CONTINUE                                                            
            ENDIF
              if (debug) write(0,*) 'DVTFG: Group ',lev,
     +          ':',(ipa(l0+l),l=1,l1)
  360    CONTINUE                                                               
C                                                                               
C        IPAT = IPA-1                                                           
C                                                                               
         DO 380 I = 1, M                                                        
            IPAT(IPA(I)) = I                                                    
  380    CONTINUE                                                               
C        DO 400 I = 1, NG                                                       
C           KLEN(I) = IWORK2(I+1) - IWORK2(I)                                   
         DO 400 I = 1, NG+1                                                     
            KLEN(I) = IWORK2(I)                                                 
  400    CONTINUE                                                               
      END IF                                                                    
      RETURN                                                                    
      END                                                                       
