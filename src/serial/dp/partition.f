      SUBROUTINE PARTITION(M, WORK, IA, N_BLOCK)
      IMPLICIT NONE

      INCLUDE 'sparker.fh'

C     ...Scalar arguments...

      INTEGER M, N_BLOCK

C     ...Array arguments...
      
      INTEGER IA(3,*), WORK(*)

C     ...Local scalars...

      INTEGER I, NNZ_ROW, N_ROWS, N_ROWS_EQ, BLOCK

      I = 1      
      N_ROWS = 0
      BLOCK = 2

      WORK(M+1) = 1

      IA(1,1) = 1

      DO WHILE(.TRUE.) 
        IF (N_ROWS.GT.MAXJDROWS) THEN
          IA(1,BLOCK) = IA(1,BLOCK-1)+MAXJDROWS
          N_ROWS = N_ROWS-MAXJDROWS
          BLOCK = BLOCK+1
        ELSE IF (N_ROWS.GE.MINJDROWS) THEN
          IA(1,BLOCK) = IA(1,BLOCK-1)+N_ROWS
          N_ROWS = 0
          BLOCK = BLOCK+1
        ELSE IF (I.LE.M) THEN
          N_ROWS_EQ = 0
          NNZ_ROW = -WORK(I)
          DO WHILE (NNZ_ROW.EQ.-WORK(I))
            N_ROWS_EQ = N_ROWS_EQ+1
            I=I+1
          ENDDO
          N_ROWS = N_ROWS + N_ROWS_EQ
        ELSE IF (N_ROWS.NE.0) THEN ! (I.GT.M)
          IA(1,BLOCK) = IA(1,BLOCK-1)+N_ROWS
          BLOCK = BLOCK+1
          GOTO 998
        ELSE
          GOTO 998
        ENDIF
      ENDDO
 998  CONTINUE

      N_BLOCK = BLOCK - 2
      
      if (ia(1,n_block+1)-1 .ne. m) then 
        write(0,*) 'PARTITION: Something wrong',m,
     +    n_block,ia(1,n_block+1),ia(1,n_block)
      endif
      END

