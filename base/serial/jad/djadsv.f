C
C             Parallel Sparse BLAS  version 2.2
C   (C) Copyright 2006/2007/2008
C                      Salvatore Filippone    University of Rome Tor Vergata
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
C     This routine compute only the casew with Unit diagonal ...
      SUBROUTINE DJADSV(UNITD,NROW,NG,A,KA,IA,JA,X,Y,IERROR)
      use psb_const_mod
      implicit none 
      real(psb_dpk_) A(*),X(*),Y(*)
      integer  IA(3,*),JA(*),KA(*),nrow,ng,i,j,k,ipg,ip2,ipx
      CHARACTER UNITD
      INTEGER   IERROR

      IERROR=0

      IF (UNITD.EQ.'U') THEN
        
        IF (NG.EQ.0) THEN
          DO I = 1, NROW
            Y(I) = X(I)
          ENDDO
        ENDIF

        DO IPG=1,NG
          DO I = IA(1,IPG),IA(1,IPG+1)-1
            Y(I) = X(I)
          END DO
*
*        LOOP ON COLUMNS
*        ---------------
*
          IP2 = IA(2,IPG)
          DO K = IP2, IA(3,IPG)-1
            IPX = IA(1,IPG)
            DO I = JA(K), JA(K+1)-1
              Y(IPX) = Y(IPX) - A(I)*Y(KA(I))
              IPX = IPX+1
            ENDDO
          ENDDO
*
*
*       LOOP ON ROWS
*       ---------------
*
          IPX = IA(1,IPG)
          DO K = IA(3,IPG), IA(2,IPG+1)-1
            DO I = JA(K), JA(K+1)-1
              Y(IPX) = Y(IPX) - A(I)*Y(KA(I))
            ENDDO
            IPX = IPX + 1
          ENDDO

****************************************
        END DO                  !END LOOP ON IPG=1,NG
****************************************
      ELSE
        WRITE(0,*) 'ERROR in DJADSV'
      ENDIF
      RETURN
      END

