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
      SUBROUTINE DJADMV2(DIAG,NROW,NCOL,ALPHA,NG,A,KA,JA,IA,
     +  X,LDX,BETA,Y,LDY, IERROR)
      IMPLICIT NONE
      INTEGER           IA(3,*),KA(*),JA(*),NCOL,NROW,NG,LDX,LDY,IERROR
      DOUBLE PRECISION  A(*),X(LDX,*),Y(LDY,*),ALPHA,BETA
      CHARACTER         DIAG                                                   
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)                                           
      INTEGER           I, J, K, IPX, IPG, I0, IN
      INTEGER           NPG
      integer           nb
      parameter         (nb=2)
      DOUBLE PRECISION  Y0(NB), Y1(NB), Y2(NB), Y3(NB), Y4(NB),
     +  Y5(NB), Y6(NB), Y7(NB), Y8(NB), Y9(NB), Y10(NB), Y11(NB),
     +  Y12(NB), Y13(NB), Y14(NB), Y15(NB)
c      .. Executable Statements ..                                              
c                                                                               
c                                                                               
c$$$      write(0,*) 'djadmv2:',diag,alpha,beta,nb
      IERROR=0
      IF (DIAG.EQ.'U') THEN
        IF (BETA.EQ.ZERO) THEN
          DO I = 1, NROW                                                 
            Y(I,1:NB) = ALPHA*X(I,1:NB)                                  
          ENDDO
        ELSE
          DO  10 I = 1, NROW                                                 
            Y(I,1:NB) = BETA*Y(I,1:NB) + ALPHA*X(I,1:NB)
 10       CONTINUE  
        ENDIF
      ELSE                                                                     
        IF (BETA.EQ.ZERO) THEN
          DO I = 1, NROW                                                  
            Y(I,1:NB) = 0.D0
          ENDDO
        ELSE
          DO 20 I = 1, NROW                                                  
            Y(I,1:NB) = BETA*Y(I,1:NB)                                                
 20       CONTINUE                                                           
        END IF
      ENDIF

      IF (ALPHA.EQ.ZERO) THEN                                                  
        RETURN                                                                
      END IF                                                                   
c                                 
c$$$         write(0,*) 'djadmv2:',diag,alpha,beta

      do 200 ipg = 1, ng                                                    
        k   = ia(2,ipg)
        npg = ja(k+1)-ja(k)

c$$$         write(0,*) 'djadmv2:',npg

        if (npg.eq.4) then
          ipx = ia(1,ipg)                                                 
          y0(1:nb)  = zero
          y1(1:nb)  = zero
          y2(1:nb)  = zero
          y3(1:nb)  = zero
          k   = ia(2,ipg)
          i0  = ja(k)
          k   = ia(3,ipg)-1
          in  = ja(k)      
          do i = i0, in+3, 4
            y0(1:nb) = y0(1:nb) + a(i+0)*x(ka(i+0),1:nb)  
            y1(1:nb) = y1(1:nb) + a(i+1)*x(ka(i+1),1:nb)  
            y2(1:nb) = y2(1:nb) + a(i+2)*x(ka(i+2),1:nb)  
            y3(1:nb) = y3(1:nb) + a(i+3)*x(ka(i+3),1:nb)              
          enddo
          y(ipx+0,1:nb)  = y(ipx+0,1:nb) + alpha*y0(1:nb) 
          y(ipx+1,1:nb)  = y(ipx+1,1:nb) + alpha*y1(1:nb) 
          y(ipx+2,1:nb)  = y(ipx+2,1:nb) + alpha*y2(1:nb)  
          y(ipx+3,1:nb)  = y(ipx+3,1:nb) + alpha*y3(1:nb) 
          
        else if (npg.eq.5) then
          
          ipx = ia(1,ipg)                                                 
          y0(1:nb)  = zero
          y1(1:nb)  = zero
          y2(1:nb)  = zero
          y3(1:nb)  = zero
          y4(1:nb)  = zero
          k   = ia(2,ipg)
          i0  = ja(k)
          k   = ia(3,ipg)-1
          in  = ja(k)
          do i = i0, in+4, 5
            y0(1:nb) = y0(1:nb) + a(i+0)*x(ka(i+0),1:nb)  
            y1(1:nb) = y1(1:nb) + a(i+1)*x(ka(i+1),1:nb)  
            y2(1:nb) = y2(1:nb) + a(i+2)*x(ka(i+2),1:nb)  
            y3(1:nb) = y3(1:nb) + a(i+3)*x(ka(i+3),1:nb)  
            y4(1:nb) = y4(1:nb) + a(i+4)*x(ka(i+4),1:nb)  
          enddo
          y(ipx+0,1:nb)  = y(ipx+0,1:nb) + alpha*y0(1:nb) 
          y(ipx+1,1:nb)  = y(ipx+1,1:nb) + alpha*y1(1:nb) 
          y(ipx+2,1:nb)  = y(ipx+2,1:nb) + alpha*y2(1:nb)  
          y(ipx+3,1:nb)  = y(ipx+3,1:nb) + alpha*y3(1:nb) 
          y(ipx+4,1:nb)  = y(ipx+4,1:nb) + alpha*y4(1:nb) 

        else if (npg.eq.6) then

          ipx = ia(1,ipg)                                                 
          y0(1:nb)  = zero
          y1(1:nb)  = zero
          y2(1:nb)  = zero
          y3(1:nb)  = zero
          y4(1:nb)  = zero
          y5(1:nb)  = zero
          k   = ia(2,ipg)
          i0  = ja(k)
          k   = ia(3,ipg)-1
          in  = ja(k)
          do i = i0, in+5, 6
            y0(1:nb) = y0(1:nb) + a(i+0)*x(ka(i+0),1:nb)  
            y1(1:nb) = y1(1:nb) + a(i+1)*x(ka(i+1),1:nb)  
            y2(1:nb) = y2(1:nb) + a(i+2)*x(ka(i+2),1:nb)  
            y3(1:nb) = y3(1:nb) + a(i+3)*x(ka(i+3),1:nb)  
            y4(1:nb) = y4(1:nb) + a(i+4)*x(ka(i+4),1:nb)  
            y5(1:nb) = y5(1:nb) + a(i+5)*x(ka(i+5),1:nb)  
          enddo
          y(ipx+0,1:nb)  = y(ipx+0,1:nb) + alpha*y0(1:nb) 
          y(ipx+1,1:nb)  = y(ipx+1,1:nb) + alpha*y1(1:nb) 
          y(ipx+2,1:nb)  = y(ipx+2,1:nb) + alpha*y2(1:nb)  
          y(ipx+3,1:nb)  = y(ipx+3,1:nb) + alpha*y3(1:nb) 
          y(ipx+4,1:nb)  = y(ipx+4,1:nb) + alpha*y4(1:nb) 
          y(ipx+5,1:nb)  = y(ipx+5,1:nb) + alpha*y5(1:nb) 

        else if (npg.eq.7) then 

          ipx = ia(1,ipg)                                                 
          y0(1:nb)  = zero
          y1(1:nb)  = zero
          y2(1:nb)  = zero
          y3(1:nb)  = zero
          y4(1:nb)  = zero
          y5(1:nb)  = zero
          y6(1:nb)  = zero
          k   = ia(2,ipg)
          i0  = ja(k)
          k   = ia(3,ipg)-1
          in  = ja(k)
          do i = i0, in+6, 7
            y0(1:nb) = y0(1:nb) + a(i+0)*x(ka(i+0),1:nb)  
            y1(1:nb) = y1(1:nb) + a(i+1)*x(ka(i+1),1:nb)  
            y2(1:nb) = y2(1:nb) + a(i+2)*x(ka(i+2),1:nb)  
            y3(1:nb) = y3(1:nb) + a(i+3)*x(ka(i+3),1:nb)  
            y4(1:nb) = y4(1:nb) + a(i+4)*x(ka(i+4),1:nb)  
            y5(1:nb) = y5(1:nb) + a(i+5)*x(ka(i+5),1:nb)  
            y6(1:nb) = y6(1:nb) + a(i+6)*x(ka(i+6),1:nb)  
          enddo
          y(ipx+0,1:nb)  = y(ipx+0,1:nb) + alpha*y0(1:nb) 
          y(ipx+1,1:nb)  = y(ipx+1,1:nb) + alpha*y1(1:nb) 
          y(ipx+2,1:nb)  = y(ipx+2,1:nb) + alpha*y2(1:nb)  
          y(ipx+3,1:nb)  = y(ipx+3,1:nb) + alpha*y3(1:nb) 
          y(ipx+4,1:nb)  = y(ipx+4,1:nb) + alpha*y4(1:nb) 
          y(ipx+5,1:nb)  = y(ipx+5,1:nb) + alpha*y5(1:nb) 
          y(ipx+6,1:nb)  = y(ipx+6,1:nb) + alpha*y6(1:nb)  

        else if (npg.eq.8) then

          ipx = ia(1,ipg)                                                 
          y0(1:nb)  = zero
          y1(1:nb)  = zero
          y2(1:nb)  = zero
          y3(1:nb)  = zero
          y4(1:nb)  = zero
          y5(1:nb)  = zero
          y6(1:nb)  = zero
          y7(1:nb)  = zero
          k   = ia(2,ipg)
          i0  = ja(k)
          k   = ia(3,ipg)-1
          in  = ja(k)
          do i = i0, in+7, 8
            y0(1:nb) = y0(1:nb) + a(i+0)*x(ka(i+0),1:nb)  
            y1(1:nb) = y1(1:nb) + a(i+1)*x(ka(i+1),1:nb)  
            y2(1:nb) = y2(1:nb) + a(i+2)*x(ka(i+2),1:nb)  
            y3(1:nb) = y3(1:nb) + a(i+3)*x(ka(i+3),1:nb)  
            y4(1:nb) = y4(1:nb) + a(i+4)*x(ka(i+4),1:nb)  
            y5(1:nb) = y5(1:nb) + a(i+5)*x(ka(i+5),1:nb)  
            y6(1:nb) = y6(1:nb) + a(i+6)*x(ka(i+6),1:nb)  
            y7(1:nb) = y7(1:nb) + a(i+7)*x(ka(i+7),1:nb)  
          enddo           
          y(ipx+0,1:nb)  = y(ipx+0,1:nb) + alpha*y0(1:nb) 
          y(ipx+1,1:nb)  = y(ipx+1,1:nb) + alpha*y1(1:nb) 
          y(ipx+2,1:nb)  = y(ipx+2,1:nb) + alpha*y2(1:nb)  
          y(ipx+3,1:nb)  = y(ipx+3,1:nb) + alpha*y3(1:nb) 
          y(ipx+4,1:nb)  = y(ipx+4,1:nb) + alpha*y4(1:nb) 
          y(ipx+5,1:nb)  = y(ipx+5,1:nb) + alpha*y5(1:nb) 
          y(ipx+6,1:nb)  = y(ipx+6,1:nb) + alpha*y6(1:nb)  
          y(ipx+7,1:nb)  = y(ipx+7,1:nb) + alpha*y7(1:nb) 

        else if (npg.eq.16) then
          
          ipx = ia(1,ipg)                                                 
          y0(1:nb)  = zero
          y1(1:nb)  = zero
          y2(1:nb)  = zero
          y3(1:nb)  = zero
          y4(1:nb)  = zero
          y5(1:nb)  = zero
          y6(1:nb)  = zero
          y7(1:nb)  = zero
          y8(1:nb)  = zero
          y9(1:nb)  = zero
          y10(1:nb)  = zero
          y11(1:nb)  = zero
          y12(1:nb)  = zero
          y13(1:nb)  = zero
          y14(1:nb)  = zero
          y15(1:nb)  = zero
          k   = ia(2,ipg)
          i0  = ja(k)
          k   = ia(3,ipg)-1
          in  = ja(k)
          do i = i0, in+15, 16
            y0(1:nb) = y0(1:nb) + a(i+0)*x(ka(i+0),1:nb)  
            y1(1:nb) = y1(1:nb) + a(i+1)*x(ka(i+1),1:nb)  
            y2(1:nb) = y2(1:nb) + a(i+2)*x(ka(i+2),1:nb)  
            y3(1:nb) = y3(1:nb) + a(i+3)*x(ka(i+3),1:nb)  
            y4(1:nb) = y4(1:nb) + a(i+4)*x(ka(i+4),1:nb)  
            y5(1:nb) = y5(1:nb) + a(i+5)*x(ka(i+5),1:nb)  
            y6(1:nb) = y6(1:nb) + a(i+6)*x(ka(i+6),1:nb)  
            y7(1:nb) = y7(1:nb) + a(i+7)*x(ka(i+7),1:nb)  
            y8(1:nb) = y8(1:nb) + a(i+8)*x(ka(i+8),1:nb)  
            y9(1:nb) = y9(1:nb) + a(i+9)*x(ka(i+9),1:nb)  
            y10(1:nb) = y10(1:nb) + a(i+10)*x(ka(i+10),1:nb)  
            y11(1:nb) = y11(1:nb) + a(i+11)*x(ka(i+11),1:nb)  
            y12(1:nb) = y12(1:nb) + a(i+12)*x(ka(i+12),1:nb)  
            y13(1:nb) = y13(1:nb) + a(i+13)*x(ka(i+13),1:nb)  
            y14(1:nb) = y14(1:nb) + a(i+14)*x(ka(i+14),1:nb)  
            y15(1:nb) = y15(1:nb) + a(i+15)*x(ka(i+15),1:nb)  
          enddo
          y(ipx+0,1:nb)  = y(ipx+0,1:nb) + alpha*y0(1:nb)  
          y(ipx+1,1:nb)  = y(ipx+1,1:nb) + alpha*y1(1:nb)   
          y(ipx+2,1:nb)  = y(ipx+2,1:nb) + alpha*y2(1:nb)
          y(ipx+3,1:nb)  = y(ipx+3,1:nb) + alpha*y3(1:nb)   
          y(ipx+4,1:nb)  = y(ipx+4,1:nb) + alpha*y4(1:nb)   
          y(ipx+5,1:nb)  = y(ipx+5,1:nb) + alpha*y5(1:nb)   
          y(ipx+6,1:nb)  = y(ipx+6,1:nb) + alpha*y6(1:nb)   
          y(ipx+7,1:nb)  = y(ipx+7,1:nb) + alpha*y7(1:nb)   
          y(ipx+8,1:nb)  = y(ipx+8,1:nb) + alpha*y8(1:nb)   
          y(ipx+9,1:nb)  = y(ipx+9,1:nb) + alpha*y9(1:nb)   
          y(ipx+10,1:nb)  = y(ipx+10,1:nb) + alpha*y10(1:nb)  
          y(ipx+11,1:nb)  = y(ipx+11,1:nb) + alpha*y11(1:nb) 
          y(ipx+12,1:nb)  = y(ipx+12,1:nb) + alpha*y12(1:nb)   
          y(ipx+13,1:nb)  = y(ipx+13,1:nb) + alpha*y13(1:nb)   
          y(ipx+14,1:nb)  = y(ipx+14,1:nb) + alpha*y14(1:nb)     
          y(ipx+15,1:nb)  = y(ipx+15,1:nb) + alpha*y15(1:nb)   

        else

          do  k = ia(2,ipg), ia(3,ipg)-1                                   
            ipx = ia(1,ipg)                                                 
            do  i = ja(k), ja(k+1) - 1                                    
              y(ipx,1:nb) = y(ipx,1:nb) + alpha*a(i)*x(ka(i),1:nb)                     
              ipx = ipx + 1                                                
            enddo
          enddo
        end if

c        csr product

        ipx = ia(1,ipg)                            
        do 70 k = ia(3,ipg), ia(2,ipg+1)-1         
          do 60 i = ja(k), ja(k+1) - 1
            y(ipx,1:nb) = y(ipx,1:nb) + alpha*a(i)*x(ka(i),1:nb)
 60       continue                                
          ipx = ipx + 1                           
 70     continue                                   
 200  continue                                                              
c                                                                               
      return                                                                   
      end                                                                      

