      SUBROUTINE CHECK_DIM(M, N, IA, NG, IA2, 
     +   NZ, LARN, LIAN1, LIAN2, IERRV)

      IMPLICIT NONE
      INCLUDE  'psb_const.fh'

C
C     .. Scalar Arguments ..
      INTEGER M,N,NG,LARN,LIAN1,LIAN2, NZ

C     .. Array Arguments ..
      INTEGER IA(3,*), IA2(*), IERRV(*)

C     Local scalars
      INTEGER NNZ, BLOCK, DIM_BLOCK, LIMIT
      INTEGER MAX_NNZERO, MAX_NZ
      
      EXTERNAL MAX_NNZERO

      MAX_NZ = MAX_NNZERO(M,IA2)
      
      NNZ = NZ
      
      LIMIT = INT(DIM_BLOCK*PSB_PERCENT_)
      
      DO BLOCK = 1, NG
         DIM_BLOCK = IA(1,BLOCK+1)-IA(1,BLOCK)
         LIMIT = INT(DIM_BLOCK*PSB_PERCENT_)

         NNZ = NNZ+(DIM_BLOCK-LIMIT)*MAX_NZ
      END DO

      IERRV(1)=0
      IERRV(2) = NNZ
      IERRV(3) = NNZ
      IERRV(4) = 6+3*(NG+1)+M+MAX_NZ*NG+1
      IF (6+3*(NG+1)+M+MAX_NZ*NG+1.GT.LIAN2) THEN
         IERRV(1) = 30
c$$$         write(0,*) 'check_dim: error 1',
c$$$     +     6+3*(NG+1)+M+MAX_NZ*NG+1,LIAN2
      ENDIF
      
      IF (NNZ.GT.LIAN1) THEN
c$$$        write(0,*) 'check_dim: error 2',nnz,lian1
         IERRV(1) = 31
      ENDIF
      
      IF (NNZ.GT.LARN) THEN
c$$$        write(0,*) 'check_dim: error 3',nnz,larn
         IERRV(1) = 32
      ENDIF

      RETURN
      END



