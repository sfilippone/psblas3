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
C     SUBROUTINE DCRDI(TRANS,M,N,DESCRA,A,IA1,IA2,IP1,DESCRN,
C                      AN,IAN1,IAN2,IP2,LAN,LIAN1,LIAN2,
C                      IAUX,LIAUX,IERRV)
C
C     Purpose: CSR to DIA format conversion
C     =======
C
C     Parameter:
C     =========
C
C     TRANS    - CHARACTER*1
C             On entry TRANS specifies whether A needs to be transposed
C             Unchanged on exit.
C
C     M        - INTEGER
C             On entry M specifies the number of rows of the matrix A.
C             M must be greater than zero.
C             Unchanged on exit.
C
C     N        - INTEGER
C             On entry M specifies the number of columns of the matrix A.
C             N must be equal to M (for the time being).
C             Unchanged on exit.
C
C     K        - INTEGER
C             On entry K specifies the number of columns of the matrix A.
C             K must be greater than or equal to zero.
C             Not used, because matrix supposed to be square.
C             Unchanged on exit.
C
C     DESCRA   - CHARACTER*5 array of DIMENSION (9)
C             On entry DESCRA defines the format of the sparse matrix.
C             Unchanged on exit.
C
C     A        - DOUBLE PRECISION array of DIMENSION (*)
C             On entry A specifies the values of the input sparse
C             matrix in CSR storage.
C             Unchanged on exit.
C
C     IA1      - INTEGER array of dimension (*)
C             On entry IA1 holds integer information on columns of input
C             sparse matrix A, i.e. which column corresponding element in
C             A belongs to.
C             Unchanged on exit.
C
C     IA2      - INTEGER array of dimension (*)
C             On entry IA2 holds rows pointers
C             Unchanged on exit.
C
C     DESCRN   - CHARACTER*5 array of DIMENSION (9)
C             On entry DESCRN defines the new format of the sparse matrix.
C             Unchanged on exit.
C
C     AN       - DOUBLE PRECISION array of DIMENSION (*)
C             On exit AN specifies the values of the input sparse
C             matrix in DIA storage (by diagonals).
C
C     IAN1     - INTEGER array of dimension (*)
C             On exit IAN1 holds integer information on columns of output
C             sparse matrix A, i.e. which diagonal is stored in each column.
C
C     IAN2     - INTEGER array of dimension (*)
C             On exit IAN2 holds in the first element the number of diagonals
C             of the matrix, i.e. the number of columns of output matrix AN.
C
C     IAUX     - INTEGER array of DIMENSION(LIAUX)
C             Work area.
C
C     LIAUX    - INTEGER
C             On entry LIAUX specifies the dimension of IAUX.
C             LIAUX must be greater than zero.
C             Unchanged on exit.
C
C     IERRV    - INTEGER array of dimension .....
C             On exit specifies if an error occur as follow:
C               IERRV(1) = 0    no error
C               IERRV(1) > 0    error
C
C
      SUBROUTINE DCRDI(TRANS,M,N,DESCRA,A,IA1,IA2,IP1,DESCRN,
     *                 AN,IAN1,IAN2,IP2,LAN,LIAN1,LIAN2,
     *                 IAUX,LIAUX,IERRV)
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      INTEGER          M, N, LAN, LIAN1, LIAN2, LIAUX
      CHARACTER        TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION A(*), AN(*)
      INTEGER          IA1(*), IA2(*), IAN1(*), IAN2(*), IP1(*), IP2(*),
     *                 IAUX(LIAUX), IERRV(*)
      CHARACTER        DESCRA*11, DESCRN*11
C     .. Local Scalars ..
      INTEGER          I, J, K, MAXND
C     .. External Subroutines ..
C     EXTERNAL XSPERR
C     .. Executable Statements ..
C

C
C       Check for argument errors
C
      IF (TRANS.NE.'T' .AND. TRANS.NE.'N') THEN
         CALL XSPERR('TRANS   ',ICHAR(TRANS),1,'DCRDI',IERRV)
      ENDIF
      IF (M.LE.0) THEN
         CALL XSPERR('MATDIM  ',M,2,'DCRDI',IERRV)
      ENDIF
      IF (N.LE.0) THEN
         CALL XSPERR('MATDIM  ',N,3,'DCRDI',IERRV)
      ENDIF
      IF(LIAN2.LT.1) THEN
         CALL XSPERR('MATST   ',LIAN2,16,'DCRDI',IERRV)
         LIAN2 = 1
      ENDIF
      MAXND = M + N - 1
      IF (LIAUX.LT.MAXND) THEN
         CALL XSPERR('LWORK   ',LIAUX,18,'DCRDI',IERRV)
         LIAUX = MAXND
      ENDIF
      IF (IERRV(1).NE.0)    RETURN


      DO J=1,MAXND
         IAUX(J)=0
      ENDDO
      DO I=1,M                          ! single out diagonals
         DO J=IA2(I),IA2(I+1)-1
            IAUX(M-I+IA1(J))=1
         ENDDO
      ENDDO

      IAN2(1)=0                         ! Computing number of diagonals
      DO J=1,MAXND
         IAN2(1)=IAN2(1)+IAUX(J)
      ENDDO

      IF (LAN.LT.(M*IAN2(1))) THEN
         CALL XSPERR('MATST   ',LAN,14,'DCRDI',IERRV)
         LAN = M * IAN2(1)
      ENDIF
      IF (LIAN1.LT.IAN2(1)) THEN
         CALL XSPERR('MATST   ',LIAN1,15,'DCRDI',IERRV)
         LIAN1 = IAN2(1)
      ENDIF
      IF (IERRV(1).NE.0)    RETURN

      DO I = 1, 3
         DESCRN(I:I)=DESCRA(I:I)
      ENDDO
      IP1(1)=0
      IP2(1)=0

      DO I=1,M                           ! zeroing AN
         DO J=1,IAN2(1)
            AN(M*(J-1)+I)=0.D0
         ENDDO
      ENDDO

      IF(TRANS.EQ.'N') THEN
C
C       Input matrix need not be permuted
C
         K=0
         DO J=M,MAXND                      ! main & upper diagonals
            IF(IAUX(J).EQ.1) THEN
               K=K+1
               IAN1(K)=J-M
            ENDIF
         ENDDO
         DO J=M-1,1,-1                     ! lower diagonals
            IF(IAUX(J).EQ.1) THEN
               K=K+1
               IAN1(K)=J-M
            ENDIF
         ENDDO

         DO I=1,M                           ! build AN (nonzeros only)
            DO J=IA2(I),IA2(I+1)-1
               DO K=1,IAN2(1)
                  IF((IA1(J)-I).EQ.IAN1(K))THEN
                     AN(M*(K-1)+I)=A(J)
                     GOTO 10
                  ENDIF
               ENDDO
   10       ENDDO
         ENDDO

      ELSE
C
C       Input matrix has to be permuted (square matrix only)
C
         K=0
         DO J=M,1,-1                       ! main & upper diagonals
            IF(IAUX(J).EQ.1) THEN
               K=K+1
               IAN1(K)=M-J
            ENDIF
         ENDDO
         DO J=M+1,MAXND                      ! lower diagonals
            IF(IAUX(J).EQ.1) THEN
               K=K+1
               IAN1(K)=M-J
            ENDIF
         ENDDO

         DO I=1,M                           ! build AN (nonzeros only)
            DO J=IA2(I),IA2(I+1)-1
               DO K=1,IAN2(1)
                  IF((I-IA1(J)).EQ.IAN1(K))THEN
                     AN(M*(K-1)+IA1(J))=A(J)
                     GOTO 20
                  ENDIF
               ENDDO
   20       ENDDO
         ENDDO

      ENDIF

      RETURN
      END

