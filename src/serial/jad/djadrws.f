C     ... Compute infinity norma for sparse matrix in CSR Format ...
      SUBROUTINE  DJADRWS(TRANS,M,N,NG,A,KA,JA,IA,
     +   INFOA,ROWSUM,IERROR)
      IMPLICIT NONE
      INCLUDE  'psb_const.fh'
C     .. Scalar Arguments ..
      INTEGER           M,N, IERROR, NG
      CHARACTER         TRANS
C     .. Array Arguments ..
      INTEGER           KA(*),JA(*),IA(3,*),INFOA(*)
      DOUBLE PRECISION  A(*), rowsum(*)
C     ... Local Scalars ..
      DOUBLE PRECISION NRMI
      INTEGER I, IR, K, IPG, NPG, IPX

      NRMI = 0.0
      IR   = 0
      DO IPG = 1, NG
         K = IA(2,IPG)
         NPG = JA(K+1)- JA(K)

C        ...  ...
         DO I = 1, NPG
            ROWSUM(IR+I) = 0.0
         ENDDO

         DO K = IA(2,IPG), IA(3,IPG)-1
            IPX = 1                                                 
            DO  I = JA(K), JA(K+1) - 1                                    
               ROWSUM(IR+IPX) = ROWSUM(IR+IPX) + ABS(A(I))
               IPX = IPX + 1                                                
            ENDDO
         ENDDO

C       ... CSR Representation ...
         
         IPX = 1
         DO K = IA(3,IPG), IA(2,IPG+1)-1         
            DO I = JA(K), JA(K+1) - 1
               ROWSUM(IR+IPX) = ROWSUM(IR+IPX) + ABS(A(I))
            ENDDO
            IPX = IPX + 1                           
         ENDDO
         
         IR = IR + NPG
      ENDDO

      END
