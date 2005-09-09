C     ... Compute infinity norma for sparse matrix in CSR Format ...
      DOUBLE PRECISION FUNCTION DJADNR(TRANS,M,N,NG,A,KA,JA,IA,
     +  INFOA,IERROR)
      IMPLICIT NONE
      INCLUDE  'psb_const.fh'
C     .. Scalar Arguments ..
      INTEGER           M,N, IERROR, NG
      CHARACTER         TRANS
C     .. Array Arguments ..
      INTEGER           KA(*),JA(*),IA(3,*),INFOA(*)
      DOUBLE PRECISION  A(*)
C     ... Local Array ..
      DOUBLE PRECISION NRMI_BLOCK(PSB_MAXJDROWS_)
C     ... Local Scalars ..
      DOUBLE PRECISION NRMI
      INTEGER I, K, IPG, NPG, IPX

      IERROR=0
      NRMI = 0.0

      DO IPG = 1, NG
        K = IA(2,IPG)
        NPG = JA(K+1)- JA(K)

C        ... Initialize NRMI_BLOCK ...
        DO I = 1, NPG
          NRMI_BLOCK(I) = 0.0
        ENDDO

        DO K = IA(2,IPG), IA(3,IPG)-1
          IPX = 1                                                 
          DO  I = JA(K), JA(K+1) - 1                                    
            NRMI_BLOCK(IPX) = NRMI_BLOCK(IPX) + ABS(A(I))
            IPX = IPX + 1                                                
          ENDDO
        ENDDO

C       ... CSR Representation ...
        
        IPX = 1
        DO K = IA(3,IPG), IA(2,IPG+1)-1         
          DO I = JA(K), JA(K+1) - 1
            NRMI_BLOCK(IPX) = NRMI_BLOCK(IPX) + ABS(A(I))
          ENDDO
          IPX = IPX + 1                           
        ENDDO
        
C        ... Compute Max in Block ...
        DO I = 1, NPG
          IF (NRMI_BLOCK(I).GT.NRMI) THEN
            NRMI = NRMI_BLOCK(I)
          ENDIF
        ENDDO
      ENDDO

      DJADNR = NRMI
      END
