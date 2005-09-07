C
C User defined function corresponding to an HPF  BLOCK partition 
C
      SUBROUTINE PART_BLK2(IDX,N,NP,PV,NV)
      
      IMPLICIT NONE
      
      INTEGER  IDX, N, NP
      INTEGER  NV
      INTEGER  PV(*)
      DOUBLE PRECISION   DDIFF
      INTEGER  IB1, IB2, IP, NB, NB1, NNB1
       
      NV = 1
      NB   = N/NP
      NB1  = NB+1 
      NNB1 = MOD(N,NP) 
      IF (IDX .LE. (NNB1*NB1))  THEN 
        PV(1) = (IDX - 1) / NB1
      ELSE 
        IF (NB > 0) THEN  
          IP    = ( (IDX-NNB1*NB1) - 1)/NB
          PV(1) = NNB1 + IP
        ELSE
          write(0,*) 'Impossible ??? '
          PV(1) = NNB1
        ENDIF
      ENDIF
            
      RETURN
      END 
      

      SUBROUTINE BLD_PARTBLK2(N,NP,IVG)

      INTEGER  N, IVG(*),NP
      INTEGER  IB1, IB2, IP, NB, NB1, NNB1, I

      NB   = N/NP
      NB1  = NB+1 
      NNB1 = MOD(N,NP) 
      DO I=1,N
         IF (I .LE. (NNB1*NB1))  THEN 
            IVG(I) = (I - 1) / NB1
         ELSE 
            IF (NB > 0) THEN  
               IP    = ( (I-NNB1*NB1) - 1)/NB
               IVG(I) = NNB1 + IP
            ELSE
               write(0,*) 'Impossible ??? '
               IVG(I) = NNB1
            ENDIF
         ENDIF
      ENDDO
      
      END
