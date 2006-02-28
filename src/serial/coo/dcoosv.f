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
C
C  Assumption: the triangular matrix has the diagonal element in the 
C  "right" place, i.e. the last in its row for Lower and the first 
C  for Upper. 
C
      SUBROUTINE DCOOSV (UPLO,TRANS,DIAG,N,AS,IA,JA,INFOA,B,X,IERROR)
      DOUBLE PRECISION   ZERO
      PARAMETER         (ZERO=0.0D0)
      LOGICAL           DEBUG
      PARAMETER         (DEBUG=.FALSE.)
      INTEGER           N,IERROR
      CHARACTER         DIAG, TRANS, UPLO
      DOUBLE PRECISION  AS(*), B(*), X(*)
      INTEGER           IA(*), JA(*),INFOA(*)
      DOUBLE PRECISION  ACC
      INTEGER           I, J, K, NNZ, II
      LOGICAL           LOW, TRA, UNI
      if (debug) write(*,*) 'DCOOSV ',n
      if (debug) write(*,*) 'DCOOSV ',n,nnz,diag,trans,uplo
      UNI = (DIAG.EQ.'U')
      TRA = (TRANS.EQ.'T')
      LOW = (UPLO.EQ.'L')
      NNZ = INFOA(1)
      if (debug) write(*,*) 'DCOOSV ',n,nnz,uni,tra,low,ia(1),ja(1)
      IERROR = 0
      if (debug) write(*,*) 'DCOOSV ierror ',ierror
      IF ( .NOT. TRA) THEN
        if (debug) write(*,*) 'DCOOSV  NOT TRA'
         IF (LOW) THEN
           if (debug) write(*,*) 'DCOOSV  LOW'
            IF ( .NOT. UNI) THEN              
              if (debug) write(*,*) 'DCOOSV  NOT UNI'
              I    = 1
              J    = I
              DO WHILE (I.LE.NNZ)              
                DO WHILE ((J.LE.NNZ).AND.(IA(J).EQ.IA(I)))
                  J = J+1
                ENDDO              
                ACC = ZERO
                IR = IA(I)
                DO K = I, J-2
                  JC = JA(K)
                  ACC = ACC + AS(K)*X(JC)
                ENDDO        
                X(IR) = (B(IR)-ACC)/AS(J-1)
                I = J 
              ENDDO

            ELSE IF (UNI) THEN
C
C    Bug warning: if UNI, some rows might be empty
C              
              I   = 1
              if (debug) write(*,*) 'DCOOSV UNILOW ',
     +           i,n,nnz,uni,tra,low
              DO II = 1, N 
                if (debug) write(*,*) 'Loop1 COOSV',i,j,ii,x(ii),b(ii)
                DO WHILE ((I.LE.NNZ).AND.(IA(I).LT.II))
                  I = I + 1 
c$$$                  if (debug) write(*,*) 'Loop2  COOSV',i,ia(i),ii
                ENDDO
                ACC = ZERO
                IF ((I.LE.NNZ).AND.(IA(I).EQ.II)) THEN
                  J  = I + 1
                  DO WHILE ((J.LE.NNZ).AND.(IA(J).EQ.IA(I)))
                    J = J+1
                  ENDDO              
                  DO K = I, J-1
                    JC = JA(K)
                    ACC = ACC + AS(K)*X(JC)
                  ENDDO        
                ELSE 
                  J = I 
                ENDIF
                X(II) = (B(II)-ACC)
                if (debug) write(*,*) 'Loop COOSV',i,j,ii,x(ii),b(ii)
                I = J                 
              ENDDO

            END IF

         ELSE IF ( .NOT. LOW) THEN
           if (debug) write(*,*) 'DCOOSV  NOT LOW'
            IF ( .NOT. UNI) THEN
              if (debug) write(*,*) 'DCOOSV  NOT UNI'
              I    = NNZ
              J    = NNZ
              DO WHILE (I.GT.0)              
                DO WHILE ((J.GT.0).AND.(IA(J).EQ.IA(I)))             
                  J = J-1
                ENDDO              
                ACC = ZERO
                IR = IA(I)
                DO K = I, J+2,-1
                  JC = JA(K)
                  ACC = ACC + AS(K)*X(JC)
                ENDDO        
                X(IR) = (B(IR)-ACC)/AS(J+1)
                I = J 
              ENDDO

            ELSE IF (UNI) THEN
              if (debug) write(*,*) 'DCOOSV  UNI'
              I   = NNZ
              DO II = N,1,-1
                DO WHILE ((I.GT.0).AND.(IA(I).GT.II))
                  I = I -1 
                ENDDO
                ACC = ZERO
                IF ((I.GT.0).AND.(IA(I).EQ.II)) THEN
                  J  = I - 1
                  DO WHILE ((J.GT.0).AND.(IA(J).EQ.IA(I)))
                    J = J-1
                  ENDDO              
                  DO K = I, J+1, -1
                    JC = JA(K)
                    ACC = ACC + AS(K)*X(JC)
                  ENDDO 
                ELSE 
                  J = I 
                ENDIF       
                X(II) = (B(II)-ACC)
                if (debug) write(*,*) 'Loop COOSV',i,j,ii,x(ii),b(ii)
                I = J 
              ENDDO

            END IF

         END IF

      ELSE IF (TRA) THEN
        IERROR = 3010
        return
CCCCCCCCCCCCCCCC
C
C  TBF
C
CCCCCCCCCCCCCCCC
         DO 180 I = 1, N
            X(I) = B(I)
  180    CONTINUE
         IF (LOW) THEN
            IF ( .NOT. UNI) THEN
               DO 220 I = N, 1, -1
                  X(I) = X(I)/AS(IA(I+1)-1)
                  ACC = X(I)
                  DO 200 J = IA(I), IA(I+1) - 2
                     K = JA(J)
                     X(K) = X(K) - AS(J)*ACC
  200             CONTINUE
  220          CONTINUE
            ELSE IF (UNI) THEN
               DO 260 I = N, 1, -1
                  ACC = X(I)
                  DO 240 J = IA(I), IA(I+1) - 1
                     K = JA(J)
                     X(K) = X(K) - AS(J)*ACC
  240             CONTINUE
  260          CONTINUE
            END IF
         ELSE IF ( .NOT. LOW) THEN
            IF ( .NOT. UNI) THEN
               DO 300 I = 1, N
                  X(I) = X(I)/AS(IA(I))
                  ACC = X(I)
                  DO 280 J = IA(I) + 1, IA(I+1) - 1
                     K = JA(J)
                     X(K) = X(K) - AS(J)*ACC
  280             CONTINUE
  300          CONTINUE
            ELSE IF (UNI) THEN
               DO 340 I = 1, N
                  ACC = X(I)
                  DO 320 J = IA(I), IA(I+1) - 1
                     K = JA(J)
                     X(K) = X(K) - AS(J)*ACC
  320             CONTINUE
  340          CONTINUE
            END IF
         END IF
      END IF
      RETURN
      END



