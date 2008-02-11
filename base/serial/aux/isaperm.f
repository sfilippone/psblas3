C
C             Parallel Sparse BLAS  version 2.2
C   (C) Copyright 2006/2007/2008
C                      Salvatore Filippone    University of Rome Tor Vergata
C                      Alfredo Buttari        University of Rome Tor Vergata
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions
C are met:
C   1. Redistributions of source code must retain the above copyright
C      notice, this list of conditions and the following disclaimer.
C   2. Redistributions in binary form must reproduce the above copyright
C      notice, this list of conditions, and the following disclaimer in the
C      documentation and/or other materials provided with the distribution.
C   3. The name of the PSBLAS group or the names of its contributors may
C      not be used to endorse or promote products derived from this
C      software without specific written permission.
C
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
C TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
C PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
C BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
C CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
C SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
C INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
C CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
C ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
C POSSIBILITY OF SUCH DAMAGE.
C
C 
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
      implicit none

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
