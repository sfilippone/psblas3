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
      SUBROUTINE ZSRSV (UPLO,TRANS,DIAG,N,AS,JA,IA,B,X)
      COMPLEX*16 ZERO
      PARAMETER  (ZERO = (0.0D0, 0.0D0))
      INTEGER    N
      CHARACTER  DIAG, TRANS, UPLO
      COMPLEX*16 AS(*), B(*), X(*)
      INTEGER    IA(*), JA(*)
      COMPLEX*16 ACC
      INTEGER    I, J, K
      LOGICAL    LOW, TRA, COTRA, UNI
      UNI = (DIAG.EQ.'U')
      TRA = (TRANS.EQ.'T')
      COTRA = (TRANS.EQ.'C')
      LOW = (UPLO.EQ.'L')
      IF ((.NOT.TRA).AND.(.NOT.COTRA)) THEN
        IF (LOW) THEN
          IF (.NOT.UNI) THEN
            DO 40 I = 1, N
              ACC = ZERO
              DO 20 J = IA(I), IA(I+1) - 2
                ACC = ACC + AS(J)*X(JA(J))
 20           CONTINUE
              X(I) = (B(I)-ACC)/AS(IA(I+1)-1)
 40         CONTINUE
          ELSE IF (UNI) THEN
            DO 80 I = 1, N
              ACC = ZERO
              DO 60 J = IA(I), IA(I+1) - 1
                ACC = ACC + AS(J)*X(JA(J))
 60           CONTINUE
              X(I) = B(I) - ACC
 80         CONTINUE
          END IF
        ELSE IF (.NOT.LOW) THEN
          IF (.NOT.UNI) THEN
            DO 120 I = N, 1, -1
              ACC = ZERO
              DO 100 J = IA(I) + 1, IA(I+1) - 1
                ACC = ACC + AS(J)*X(JA(J))
 100          CONTINUE
              X(I) = (B(I)-ACC)/AS(IA(I))
 120        CONTINUE
          ELSE IF (UNI) THEN
            DO 160 I = N, 1, -1
              ACC = ZERO
              DO 140 J = IA(I), IA(I+1) - 1
                ACC = ACC + AS(J)*X(JA(J))
 140          CONTINUE
              X(I) = B(I) - ACC
 160        CONTINUE
          END IF
        END IF
      ELSE IF (TRA) THEN
        DO 180 I = 1, N
          X(I) = B(I)
 180    CONTINUE
        IF (LOW) THEN
          IF (.NOT.UNI) THEN
            DO 220 I = N, 1, -1
              X(I) = X(I)/AS(IA(I+1)-1)
              ACC = X(I)
              DO 200 J = IA(I), IA(I+1) - 2
                K = JA(J)
                X(K) = X(K) - AS(J)*ACC
 200          CONTINUE
 220        CONTINUE
          ELSE IF (UNI) THEN
            DO 260 I = N, 1, -1
              ACC = X(I)
              DO 240 J = IA(I), IA(I+1) - 1
                K = JA(J)
                X(K) = X(K) - AS(J)*ACC
 240          CONTINUE
 260        CONTINUE
          END IF
        ELSE IF (.NOT.LOW) THEN
          IF (.NOT.UNI) THEN
            DO 300 I = 1, N
              X(I) = X(I)/AS(IA(I))
              ACC = X(I)
              DO 280 J = IA(I) + 1, IA(I+1) - 1
                K = JA(J)
                X(K) = X(K) - AS(J)*ACC
 280          CONTINUE
 300        CONTINUE
          ELSE IF (UNI) THEN
            DO 340 I = 1, N
              ACC = X(I)
              DO 320 J = IA(I), IA(I+1) - 1
                K = JA(J)
                X(K) = X(K) - AS(J)*ACC
 320          CONTINUE
 340        CONTINUE
          END IF
        END IF
      ELSE IF (COTRA) THEN
        DO 580 I = 1, N
          X(I) = B(I)
 580    CONTINUE
        IF (LOW) THEN
          IF (.NOT.UNI) THEN
            DO 620 I = N, 1, -1
              X(I) = X(I)/CONJG(AS(IA(I+1)-1))
              ACC = X(I)
              DO 600 J = IA(I), IA(I+1) - 2
                K = JA(J)
                X(K) = X(K) - CONJG(AS(J))*ACC
 600          CONTINUE
 620        CONTINUE
          ELSE IF (UNI) THEN
            DO 660 I = N, 1, -1
              ACC = X(I)
              DO 640 J = IA(I), IA(I+1) - 1
                K = JA(J)
                X(K) = X(K) - CONJG(AS(J))*ACC
 640          CONTINUE
 660        CONTINUE
          END IF
        ELSE IF (.NOT.LOW) THEN
          IF (.NOT.UNI) THEN
            DO 700 I = 1, N
              X(I) = X(I)/CONJG(AS(IA(I)))
              ACC = X(I)
              DO 680 J = IA(I) + 1, IA(I+1) - 1
                K = JA(J)
                X(K) = X(K) - CONJG(AS(J))*ACC
 680          CONTINUE
 700        CONTINUE
          ELSE IF (UNI) THEN
            DO 740 I = 1, N
              ACC = X(I)
              DO 720 J = IA(I), IA(I+1) - 1
                K = JA(J)
                X(K) = X(K) - CONJG(AS(J))*ACC
 720          CONTINUE
 740        CONTINUE
          END IF
        END IF
      END IF
      RETURN
      END



