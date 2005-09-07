***********************************************************************
*                                                                     *
*    REFERENCES   =                                                   *
*                                                                     *
*          [1]  D. E. Knuth                                           *
*               The art of computer programming Vol. 1                *
*               Sect. 1.3.3                                           *
*               Addison-Wesley                                        *
*                                                                     *
*                                                                     *
*    FUNCTION: Checks whether a vector really is a permutation        *
*              works by exploiting the cycle structure of the         *
*              permutation.                                           *
*                                                                     *
*                                                                     *
***********************************************************************
      LOGICAL FUNCTION ISAPERM(N,IP)               
C     .. Scalar Arguments ..                                                    
      INTEGER N                                                                 
C     ..                                                                        
C     .. Array Arguments ..                                                     
      INTEGER IP(N)
C     ..                                                                        
C     .. Local Scalars ..                                                       
      INTEGER I,J,M
C     ..          
 
      ISAPERM = .TRUE.
C
C   Sanity check first 
C     
      DO I=1, N 
        IF ((IP(I).LT.1).OR.(IP(I).GT.N)) THEN
          ISAPERM = .FALSE.
          RETURN
        ENDIF
      ENDDO
      
C
C Now work through the cycles, by marking each successive item as negative.
C No cycle should intersect with any other, hence the .GE.1 check. 
C
      DO M = 1, N    
        I = IP(M) 
        IF (I.LT.0) THEN      
          IP(M) = -I          
        ELSE IF (I.NE.M) THEN 
          J     = IP(I)               
          IP(I) = -J          
          I     = J
          DO WHILE ((J.GE.1).AND.(J.NE.M))
            J     = IP(I)               
            IP(I) = -J 
            I     = J               
          ENDDO   
          IP(M) = IABS(IP(M))
          IF (J.NE.M) THEN 
            ISAPERM = .FALSE.
            DO I=1, N
              IP(I) = IABS(IP(I))
            ENDDO
            GOTO 9999
          ENDIF
        END IF             
      ENDDO
 9999 CONTINUE 

      RETURN                                                                    
      END                                                                       
