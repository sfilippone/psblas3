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
C     SUBROUTINE DCREL(TRANS,M,N,DESCRA,A,IA1,IA2,IP1,DESCRN,
C                      AN,IAN1,IAN2,IP2,LAN,LIAN1,LIAN2,
C                      IAUX,LIAUX,IERRV)
C
C     Purpose: CSR to ELL format conversion
C     =======
C
C     Parameter:
C     =========
C
C     ...
C     IAN2  -  Vector: first element is max number of columns in matrices
C              ARN,IAN1, elements to M+1 are column index of diagonal
C              in ARN,IAN1 (in future releases)
C     ...
C
C
      SUBROUTINE DCREL(TRANS,M,N,DESCRA,A,IA1,IA2,IP1,DESCRN,
     *                 AN,IAN1,IAN2,IP2,LAN,LIAN1,LIAN2,
     *                 IAUX,LIAUX,IERRV)
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      INTEGER          LAN, LIAUX, LIAN1, LIAN2, M, N
      CHARACTER        TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION A(*), AN(*)
      INTEGER          IA1(*), IA2(*), IAN1(*), IAN2(*), IP1(*), IP2(*),
     *                 IAUX(LIAUX), IERRV(*)
      CHARACTER        DESCRA*11, DESCRN*11
C     .. Local Scalars ..
      INTEGER          I, J, LWORKR
C     .. External Subroutines ..
      EXTERNAL XSPERR
C     .. Executable Statements ..
C

C
C       Check for argument errors
C
      IF (TRANS.NE.'T' .AND. TRANS.NE.'N') THEN
         CALL XSPERR('TRANS   ',ICHAR(TRANS),1,'DCREL',IERRV)
      ENDIF
      IF (M.LE.0) THEN
         CALL XSPERR('MATDIM  ',M,2,'DCREL',IERRV)
      ENDIF
      IF (N.LE.0) THEN
         CALL XSPERR('MATDIM  ',N,3,'DCREL',IERRV)
      ENDIF
      IF(LIAN2.LT.1) THEN
         LIAN2 = 1
         CALL XSPERR('MATST   ',LIAN2,16,'DCREL',IERRV)
      ENDIF
      IF (TRANS.EQ.'N') THEN
         LWORKR = 0
      ELSE IF (TRANS.EQ.'T') THEN
         LWORKR =  N
      ENDIF
      IF (LIAUX.LT.LWORKR) THEN
         CALL XSPERR('LWORK   ',LIAUX,18,'DCREL',IERRV)
         LIAUX = LWORKR
      ENDIF
      IF (IERRV(1).NE.0)    RETURN


      descrn(1:3) = descra(1:3)
      IP1(1)=0
      IP2(1)=0

      IF(TRANS.EQ.'N') THEN
C
C       Input matrix need not be permuted
C
         IAN2(1)=IA2(2)-IA2(1)
         DO I = 2, M
            IAN2(1) = MAX0(IAN2(1),IA2(I+1)-IA2(I))
         ENDDO

         IF (LAN.LT.(M*IAN2(1))) THEN
            CALL XSPERR('MATST   ',LAN,14,'DCREL',IERRV)
            LAN = M * IAN2(1)
         ENDIF
         IF (LIAN1.LT.(M*IAN2(1))) THEN
            CALL XSPERR('MATST   ',LIAN1,15,'DCREL',IERRV)
            LIAN1 = M * IAN2(1)
         ENDIF
         IF (IERRV(1).NE.0)    RETURN

         DO I=1,M
            DO J=IA2(I),IA2(I+1)-1
               AN(M*(J-IA2(I))+I)=A(J)
               IAN1(M*(J-IA2(I))+I)=IA1(J)
            ENDDO
            DO J=IA2(I+1)-IA2(I)+1,IAN2(1)
               AN(M*(J-1)+I)=0.D0
               IAN1(M*(J-1)+I)=IAN2((J-2)*M+I)
            ENDDO
         ENDDO

      ELSE
C
C       Input matrix has to be permuted
C

         DO J=1,N
            IAUX(I)=0
         ENDDO
         DO I=1,M
            DO J=IA2(I),IA2(I+1)-1
               IAUX(IA1(J))=IAUX(IA1(J))+1
            ENDDO
         ENDDO
         IAN2(1)=IAUX(1)
         DO I = 2, M
            IAN2(1) = MAX0(IAN2(1),IAUX(I))
         ENDDO

         IF (LAN.LT.(N*IAN2(1))) THEN
            CALL XSPERR('MATST   ',LAN,14,'DCREL',IERRV)
            LAN = N * IAN2(1)
         ENDIF
         IF (LIAN1.LT.(N*IAN2(1))) THEN
            CALL XSPERR('MATST   ',LIAN1,15,'DCREL',IERRV)
            LIAN1 = N * IAN2(1)
         ENDIF
         IF (IERRV(1).NE.0)    RETURN

         DO J=1,N
            IAUX(I)=0
         ENDDO
         DO I=1,M
            DO J=IA2(I),IA2(I+1)-1
               IAUX(IA1(J))=IAUX(IA1(J))+1
               AN  (N*(IAUX(IA1(J)))+IA1(J))=A(J)
               IAN1(N*(IAUX(IA1(J)))+IA1(J))=I
            ENDDO
         ENDDO
         DO I=1,N
            DO J=IAUX(I)+1,IAN2(I)
               AN  (N*(J-1)+I)=0.D0
               IAN1(N*(J-1)+I)=IAN1(N*IAUX(I)+I)
            ENDDO
         ENDDO

      ENDIF


      RETURN
      END

