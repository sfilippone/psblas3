      SUBROUTINE GEN_BLOCK(M,NG,IA,AUX)
      IMPLICIT NONE

      INCLUDE 'psb_const.fh'
      INTEGER M, NG
      INTEGER IA(3,*), AUX(*)

      INTEGER BLOCK, I, N_ROWS

      N_ROWS = IA(1,2) - IA(1,1)
      I = 2
      BLOCK = 2
      AUX(1) = 1
      
      DO WHILE(.TRUE.)
        IF (N_ROWS.GT.PSB_MAXJDROWS_) THEN
          AUX(BLOCK) = AUX(BLOCK-1)+PSB_MAXJDROWS_
          N_ROWS = N_ROWS-PSB_MAXJDROWS_
          BLOCK = BLOCK+1
        ELSE IF (N_ROWS.GT.0) THEN
          AUX(BLOCK) = AUX(BLOCK-1)+N_ROWS
          N_ROWS = 0
          BLOCK = BLOCK+1
        ELSE IF (I.LE.NG) THEN
          N_ROWS = IA(1,I+1) - IA(1,I)
          I = I+1
        ELSE
          GOTO 998
        ENDIF
      ENDDO
 998  CONTINUE 

C     ... Copy AUX in IA(1,*)

      NG = BLOCK - 2
      DO I = 1, NG+1
        IA(1,I) = AUX(I)
      ENDDO

      END
