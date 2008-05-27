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
*                PROCEDURAL LOGIC SECTION                             *         
*      SUBROUTINE DJADMV (DIAG,NROW,NCOL,ALPHA,NG,A,KA,JA,IA,X,BETA,Y) *         
*      DOUBLE PRECISION  ZERO                                         *         
*      PARAMETER         (ZERO=0.0D0)                                 *         
*      DOUBLE PRECISION  ACC                                          *         
*      INTEGER           I, J, K, IPX, IPG                            *         
*      LOGICAL           UNI                                          *         
*C     .. Executable Statements ..                                    *         
*C                                                                    *         
*C                                                                    *         
*      IF (DIAG.EQ.'U') THEN                                          *         
*            DO  10 I = 1, M                                          *         
*                Y(I) = BETA*Y(I) + ALPHA*X(I)                        *         
*   10       CONTINUE                                                 *         
*      ELSE                                                           *         
*            DO 20 I = 1, M                                           *         
*               Y(I) = BETA*Y(I)                                      *         
*   20       CONTINUE                                                 *         
*      END IF                                                         *         
*                                                                     *         
*      IF (ALPHA.EQ.ZERO) THEN                                        *         
*         RETURN                                                      *         
*      END IF                                                         *         
*C                                                                    *         
*C                                                                    *         
*C        DO 200 IPG = 1, NG                                          *         
*            DO 50 K = IA(2,IPG), IA(3,IPG)-1                         *         
*               IPX = IA(1,IPG)                                       *         
*               DO 40 I = JA(K), JA(K+1) - 1                          *         
*                  Y(IPX) = Y(IPX) + ALPHA*A(I)*X(KA(I))              *         
*                  IPX = IPX + 1                                      *         
*   40          CONTINUE                                              *         
*   50       CONTINUE                                                 *         
*C                                                                    *         
*C                                                                    *         
*            IPX = IA(1,IPG)                                          *         
*            DO 70 K = IA(3,IPG), IA(2,IPG+1)-1                       *         
*               DO 60 I = JA(K), JA(K+1) - 1                          *         
*                  Y(IPX) = Y(IPX) + ALPHA*A(I)*X(KA(I))              *         
*   60          CONTINUE                                              *         
*               IPX = IPX + 1                                         *         
*   70       CONTINUE                                                 *         
*  200    CONTINUE                                                    *         
*C                                                                    *         
*      RETURN                                                         *         
*C                                                                    *         
*C                                                                    *         
*      END                                                            *         
*                                                                     *         
*                                                                     *         
***********************************************************************       
      SUBROUTINE SJADMV (DIAG,NROW,NCOL,ALPHA,NG,A,KA,JA,IA,X,
     +     BETA,Y,IERROR)
      use psb_const_mod
      IMPLICIT NONE
      real(psb_spk_)  A(*),X(*),Y(*),ALPHA,BETA
      INTEGER           IA(3,*),KA(*),JA(*),NCOL,NROW,NG,IERROR
      CHARACTER         DIAG                                                   
      INTEGER           I, K, IPX, IPG, I0, IN
      INTEGER           NPG
      real(psb_spk_)  Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, 
     +  Y8, Y9, Y10, Y11, Y12, Y13, Y14, Y15
c      .. Executable Statements ..                                              
c                                                                               
c                                                              
      IERROR=0                 
      IF (DIAG.EQ.'U') THEN
        IF (BETA.EQ.SZERO) THEN
          DO I = 1, NROW                                                 
            Y(I) = ALPHA*X(I)                                  
          ENDDO
        ELSE
          DO  10 I = 1, NROW                                                 
            Y(I) = BETA*Y(I) + ALPHA*X(I)                                  
 10       CONTINUE  
        ENDIF
      ELSE                                                                     
        IF (BETA.EQ.SZERO) THEN
          DO I = 1, NROW                                                  
            Y(I) = 0.D0
          ENDDO
        ELSE
          DO 20 I = 1, NROW                                                  
            Y(I) = BETA*Y(I)                                                
 20       CONTINUE                                                           
        END IF
      ENDIF

      IF (ALPHA.EQ.SZERO) THEN                                                  
        RETURN                                                                
      END IF                                                                   
c                                                                               
c                                                                               
      DO 200 IPG = 1, NG                                                    
        K   = IA(2,IPG)
        NPG = JA(K+1)-JA(K)

        IF (NPG.EQ.4) THEN
          IPX = IA(1,IPG)                                                 
          Y0  = SZERO
          Y1  = SZERO
          Y2  = SZERO
          Y3  = SZERO
          K   = IA(2,IPG)
          I0  = JA(K)
          K   = IA(3,IPG)-1
          IN  = JA(K)      
          DO I = I0, IN+3, 4
            Y0 = Y0 + A(I+0)*X(KA(I+0))  
            Y1 = Y1 + A(I+1)*X(KA(I+1))  
            Y2 = Y2 + A(I+2)*X(KA(I+2))  
            Y3 = Y3 + A(I+3)*X(KA(I+3))              
          ENDDO
          Y(IPX+0)  = Y(IPX+0) + ALPHA*Y0 
          Y(IPX+1)  = Y(IPX+1) + ALPHA*Y1 
          Y(IPX+2)  = Y(IPX+2) + ALPHA*Y2  
          Y(IPX+3)  = Y(IPX+3) + ALPHA*Y3 
          
        ELSE IF (NPG.EQ.5) THEN
          
          IPX = IA(1,IPG)                                                 
          Y0  = SZERO
          Y1  = SZERO
          Y2  = SZERO
          Y3  = SZERO
          Y4  = SZERO
          K   = IA(2,IPG)
          I0  = JA(K)
          K   = IA(3,IPG)-1
          IN  = JA(K)
          DO I = I0, IN+4, 5
            Y0 = Y0 + A(I+0)*X(KA(I+0))  
            Y1 = Y1 + A(I+1)*X(KA(I+1))  
            Y2 = Y2 + A(I+2)*X(KA(I+2))  
            Y3 = Y3 + A(I+3)*X(KA(I+3))  
            Y4 = Y4 + A(I+4)*X(KA(I+4))  
          ENDDO
          Y(IPX+0)  = Y(IPX+0) + ALPHA*Y0 
          Y(IPX+1)  = Y(IPX+1) + ALPHA*Y1 
          Y(IPX+2)  = Y(IPX+2) + ALPHA*Y2  
          Y(IPX+3)  = Y(IPX+3) + ALPHA*Y3 
          Y(IPX+4)  = Y(IPX+4) + ALPHA*Y4 

        ELSE IF (NPG.EQ.6) THEN

          IPX = IA(1,IPG)                                                 
          Y0  = SZERO
          Y1  = SZERO
          Y2  = SZERO
          Y3  = SZERO
          Y4  = SZERO
          Y5  = SZERO
          K   = IA(2,IPG)
          I0  = JA(K)
          K   = IA(3,IPG)-1
          IN  = JA(K)
          DO I = I0, IN+5, 6
            Y0 = Y0 + A(I+0)*X(KA(I+0))  
            Y1 = Y1 + A(I+1)*X(KA(I+1))  
            Y2 = Y2 + A(I+2)*X(KA(I+2))  
            Y3 = Y3 + A(I+3)*X(KA(I+3))  
            Y4 = Y4 + A(I+4)*X(KA(I+4))  
            Y5 = Y5 + A(I+5)*X(KA(I+5))  
          ENDDO
          Y(IPX+0)  = Y(IPX+0) + ALPHA*Y0 
          Y(IPX+1)  = Y(IPX+1) + ALPHA*Y1 
          Y(IPX+2)  = Y(IPX+2) + ALPHA*Y2  
          Y(IPX+3)  = Y(IPX+3) + ALPHA*Y3 
          Y(IPX+4)  = Y(IPX+4) + ALPHA*Y4 
          Y(IPX+5)  = Y(IPX+5) + ALPHA*Y5 

        ELSE IF (NPG.EQ.7) THEN 

          IPX = IA(1,IPG)                                                 
          Y0  = SZERO
          Y1  = SZERO
          Y2  = SZERO
          Y3  = SZERO
          Y4  = SZERO
          Y5  = SZERO
          Y6  = SZERO
          K   = IA(2,IPG)
          I0  = JA(K)
          K   = IA(3,IPG)-1
          IN  = JA(K)
          DO I = I0, IN+6, 7
            Y0 = Y0 + A(I+0)*X(KA(I+0))  
            Y1 = Y1 + A(I+1)*X(KA(I+1))  
            Y2 = Y2 + A(I+2)*X(KA(I+2))  
            Y3 = Y3 + A(I+3)*X(KA(I+3))  
            Y4 = Y4 + A(I+4)*X(KA(I+4))  
            Y5 = Y5 + A(I+5)*X(KA(I+5))  
            Y6 = Y6 + A(I+6)*X(KA(I+6))  
          ENDDO
          Y(IPX+0)  = Y(IPX+0) + ALPHA*Y0 
          Y(IPX+1)  = Y(IPX+1) + ALPHA*Y1 
          Y(IPX+2)  = Y(IPX+2) + ALPHA*Y2  
          Y(IPX+3)  = Y(IPX+3) + ALPHA*Y3 
          Y(IPX+4)  = Y(IPX+4) + ALPHA*Y4 
          Y(IPX+5)  = Y(IPX+5) + ALPHA*Y5 
          Y(IPX+6)  = Y(IPX+6) + ALPHA*Y6  

        ELSE IF (NPG.EQ.8) THEN

          IPX = IA(1,IPG)                                                 
          Y0  = SZERO
          Y1  = SZERO
          Y2  = SZERO
          Y3  = SZERO
          Y4  = SZERO
          Y5  = SZERO
          Y6  = SZERO
          Y7  = SZERO
          K   = IA(2,IPG)
          I0  = JA(K)
          K   = IA(3,IPG)-1
          IN  = JA(K)
          DO I = I0, IN+7, 8
            Y0 = Y0 + A(I+0)*X(KA(I+0))  
            Y1 = Y1 + A(I+1)*X(KA(I+1))  
            Y2 = Y2 + A(I+2)*X(KA(I+2))  
            Y3 = Y3 + A(I+3)*X(KA(I+3))  
            Y4 = Y4 + A(I+4)*X(KA(I+4))  
            Y5 = Y5 + A(I+5)*X(KA(I+5))  
            Y6 = Y6 + A(I+6)*X(KA(I+6))  
            Y7 = Y7 + A(I+7)*X(KA(I+7))  
          ENDDO           
          Y(IPX+0)  = Y(IPX+0) + ALPHA*Y0 
          Y(IPX+1)  = Y(IPX+1) + ALPHA*Y1 
          Y(IPX+2)  = Y(IPX+2) + ALPHA*Y2  
          Y(IPX+3)  = Y(IPX+3) + ALPHA*Y3 
          Y(IPX+4)  = Y(IPX+4) + ALPHA*Y4 
          Y(IPX+5)  = Y(IPX+5) + ALPHA*Y5 
          Y(IPX+6)  = Y(IPX+6) + ALPHA*Y6  
          Y(IPX+7)  = Y(IPX+7) + ALPHA*Y7 

        ELSE IF (NPG.EQ.16) THEN
          
          IPX = IA(1,IPG)                                                 
          Y0  = SZERO
          Y1  = SZERO
          Y2  = SZERO
          Y3  = SZERO
          Y4  = SZERO
          Y5  = SZERO
          Y6  = SZERO
          Y7  = SZERO
          Y8  = SZERO
          Y9  = SZERO
          Y10  = SZERO
          Y11  = SZERO
          Y12  = SZERO
          Y13  = SZERO
          Y14  = SZERO
          Y15  = SZERO
          K   = IA(2,IPG)
          I0  = JA(K)
          K   = IA(3,IPG)-1
          IN  = JA(K)
          DO I = I0, IN+15, 16
            Y0 = Y0 + A(I+0)*X(KA(I+0))  
            Y1 = Y1 + A(I+1)*X(KA(I+1))  
            Y2 = Y2 + A(I+2)*X(KA(I+2))  
            Y3 = Y3 + A(I+3)*X(KA(I+3))  
            Y4 = Y4 + A(I+4)*X(KA(I+4))  
            Y5 = Y5 + A(I+5)*X(KA(I+5))  
            Y6 = Y6 + A(I+6)*X(KA(I+6))  
            Y7 = Y7 + A(I+7)*X(KA(I+7))  
            Y8 = Y8 + A(I+8)*X(KA(I+8))  
            Y9 = Y9 + A(I+9)*X(KA(I+9))  
            Y10 = Y10 + A(I+10)*X(KA(I+10))  
            Y11 = Y11 + A(I+11)*X(KA(I+11))  
            Y12 = Y12 + A(I+12)*X(KA(I+12))  
            Y13 = Y13 + A(I+13)*X(KA(I+13))  
            Y14 = Y14 + A(I+14)*X(KA(I+14))  
            Y15 = Y15 + A(I+15)*X(KA(I+15))  
          ENDDO
          Y(IPX+0)  = Y(IPX+0) + ALPHA*Y0 
          Y(IPX+1)  = Y(IPX+1) + ALPHA*Y1 
          Y(IPX+2)  = Y(IPX+2) + ALPHA*Y2  
          Y(IPX+3)  = Y(IPX+3) + ALPHA*Y3 
          Y(IPX+4)  = Y(IPX+4) + ALPHA*Y4 
          Y(IPX+5)  = Y(IPX+5) + ALPHA*Y5 
          Y(IPX+6)  = Y(IPX+6) + ALPHA*Y6  
          Y(IPX+7)  = Y(IPX+7) + ALPHA*Y7 
          Y(IPX+8)  = Y(IPX+8) + ALPHA*Y8 
          Y(IPX+9)  = Y(IPX+9) + ALPHA*Y9 
          Y(IPX+10)  = Y(IPX+10) + ALPHA*Y10  
          Y(IPX+11)  = Y(IPX+11) + ALPHA*Y11 
          Y(IPX+12)  = Y(IPX+12) + ALPHA*Y12 
          Y(IPX+13)  = Y(IPX+13) + ALPHA*Y13 
          Y(IPX+14)  = Y(IPX+14) + ALPHA*Y14  
          Y(IPX+15)  = Y(IPX+15) + ALPHA*Y15

        ELSE

          DO  K = IA(2,IPG), IA(3,IPG)-1                                   
            IPX = IA(1,IPG)                                                 
            DO  I = JA(K), JA(K+1) - 1                                    
              Y(IPX) = Y(IPX) + ALPHA*A(I)*X(KA(I))                        
              IPX = IPX + 1                                                
            ENDDO
          ENDDO
        END IF

c        CSR Product

        IPX = IA(1,IPG)                            
        DO 70 K = IA(3,IPG), IA(2,IPG+1)-1         
          DO 60 I = JA(K), JA(K+1) - 1
            Y(IPX) = Y(IPX) + ALPHA*A(I)*X(KA(I))
 60       CONTINUE                                
          IPX = IPX + 1                           
 70     CONTINUE                                   
 200  CONTINUE                                                              
c                                                                               
      RETURN                                                                   
      END                                                                      

