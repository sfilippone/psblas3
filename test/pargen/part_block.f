C
C User defined function corresponding to an HPF  BLOCK partition 
C
      SUBROUTINE PART_BLOCK(GLOBAL_INDX,N,NP,PV,NV)
      
      IMPLICIT NONE
      
      INTEGER  GLOBAL_INDX, N, NP
      INTEGER  NV
      INTEGER  PV(*)
      INTEGER  DIM_BLOCK
      DOUBLE PRECISION   DDIFF
      INTEGER  IB1, IB2, IPV
      
      double precision PC
      PARAMETER   (PC=0.0D0)

      DIM_BLOCK = (N + NP - 1)/NP
      NV = 1  
      PV(NV) = (GLOBAL_INDX - 1) / DIM_BLOCK
      
      IPV = PV(1)
      IB1 = IPV * DIM_BLOCK + 1
      IB2 = (IPV+1) * DIM_BLOCK
      
      DDIFF = DBLE(ABS(GLOBAL_INDX-IB1))/DBLE(DIM_BLOCK)
      IF (DDIFF .LT. PC/2) THEN
C
C     Overlap at the beginning of a block, with the previous proc
C         
         IF (IPV.GT.0) THEN 
            NV     = NV + 1
            PV(NV) = IPV - 1
         ENDIF
      ENDIF

      DDIFF = DBLE(ABS(GLOBAL_INDX-IB2))/DBLE(DIM_BLOCK)
      IF (DDIFF .LT. PC/2) THEN
C
C     Overlap at the end of a block, with the next proc
C         
         IF (IPV.LT.(NP-1)) THEN 
            NV     = NV + 1
            PV(NV) = IPV + 1
         ENDIF
      ENDIF
      
      RETURN
      END 
      
