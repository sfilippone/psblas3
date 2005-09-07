      SUBROUTINE DGBLOCK(M,IA2,IPERM,IA,N_BLOCK,WORK,LWORK)
      IMPLICIT NONE

C     ...Scalar arguments...

      INTEGER M, N_BLOCK, LWORK

C     ...Array arguments...
      
      INTEGER IA2(*), IPERM(*), IA(3,*), WORK(*)

C     ...Local scalars...

      INTEGER I, SWAP, KK, LP, IRET
C     Compute number of nnzero elements per row

      IPERM(1) = 0

      DO I = 1, M
        WORK(I) = - IA2(I+1) + IA2(I)
      ENDDO

C     Sorting Array work
C ........................

      CALL MRGSRT(M,WORK,WORK(M+1),IRET)
      IF (IRET.EQ.0) THEN
C     Construct IPERM Vector
        LP = WORK(M+1)
        
        DO I = 1, M
          IPERM(LP) = I
          LP = WORK(M+1+LP)
        ENDDO

        LP = WORK(M+1)                           
        KK = 1                              
        DO WHILE (.NOT.((LP.EQ.0).OR.(KK.GT.M)))
          DO WHILE (LP.LT.KK)
            LP = WORK(M+1+LP)
          ENDDO
C        Swap values of array work         
          SWAP = WORK(KK)                      
          WORK(KK) = WORK(LP)                     
          WORK(LP) = SWAP                      

C        Swap values of index array work(m+1)
          SWAP = WORK(M+1+LP)                     
          WORK(M+1+LP) = WORK(M+1+KK)                     
          WORK(M+1+KK) = LP

          LP    = SWAP                     
          KK = KK+1                         
        ENDDO

      ENDIF
C     Partitioning Matrix in blocks of rows
      CALL PARTITION(M, WORK, IA, N_BLOCK)

      END



