      SUBROUTINE DCSRSV(UPLO,TRANS,DIAG,N,AS,JA,IA,B,X,IERROR)
      DOUBLE PRECISION   ZERO
      PARAMETER         (ZERO=0.0D0)
      INTEGER           N, IERROR
      CHARACTER         DIAG, TRANS, UPLO
      DOUBLE PRECISION  AS(*), B(*), X(*)
      INTEGER           IA(*), JA(*)
      DOUBLE PRECISION  ACC
      INTEGER           I, J, K
      LOGICAL           LOW, TRA, UNI
      UNI = (DIAG.EQ.'U')
      TRA = (TRANS.EQ.'T')
      LOW = (UPLO.EQ.'L')
      IF ( .NOT. TRA) THEN
         IF (LOW) THEN
            IF ( .NOT. UNI) THEN
               DO 40 I = 1, N
                  ACC = ZERO
                  DO 20 J = IA(I), IA(I+1) - 2
                     ACC = ACC + AS(J)*X(JA(J))
   20             CONTINUE
                  X(I) = (B(I)-ACC)/AS(IA(I+1)-1)
   40          CONTINUE
            ELSE IF (UNI) THEN
               DO 80 I = 1, N
                  ACC = ZERO
                  DO 60 J = IA(I), IA(I+1) - 1
                     ACC = ACC + AS(J)*X(JA(J))
   60             CONTINUE
                  X(I) = B(I) - ACC
   80          CONTINUE
            END IF
         ELSE IF ( .NOT. LOW) THEN
            IF ( .NOT. UNI) THEN
               DO 120 I = N, 1, -1
                  ACC = ZERO
                  DO 100 J = IA(I) + 1, IA(I+1) - 1
                     ACC = ACC + AS(J)*X(JA(J))
  100             CONTINUE
                  X(I) = (B(I)-ACC)/AS(IA(I))
  120          CONTINUE
            ELSE IF (UNI) THEN
               DO 160 I = N, 1, -1
                  ACC = ZERO
                  DO 140 J = IA(I), IA(I+1) - 1
                     ACC = ACC + AS(J)*X(JA(J))
  140             CONTINUE
                  X(I) = B(I) - ACC
  160          CONTINUE
            END IF
         END IF
      ELSE IF (TRA) THEN
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



