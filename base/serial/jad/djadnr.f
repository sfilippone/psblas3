C
C             Parallel Sparse BLAS  v2.0
C   (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
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
C     ... Compute infinity norma for sparse matrix in CSR Format ...
      DOUBLE PRECISION FUNCTION DJADNR(TRANS,M,N,NG,A,KA,JA,IA,
     +  INFOA,IERROR)
      use psb_const_mod
      use psb_spmat_type
      IMPLICIT NONE
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
