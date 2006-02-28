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
      SUBROUTINE DCRUPDATE(M, N, DESCRA, A, IA1,
     +   IA2, INFOA, IA, JA, DESCRH, H, IH1, IH2,
     +   INFOH, IH, JH, FLAG, GLOB_TO_LOC,
     +   IWORK, LIWORK, IERROR)
C
C     .. Matrix A to be updated is required to be stored with
C     .. column indices belonging to the same row ordered.
C     .. Block H to be inserted don't need to be stored in such way.
C
C     Flag = 0: put elements to 0.0D0;
C     Flag = 1: replace elements with new value;
C     Flag = 2: sum block value to elements;
C
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           IA, JA, IH, JH, M, N,
     +                  IERROR, FLAG, LIWORK 
C     .. Array Arguments ..
      INTEGER           IA1(*),IA2(*),IH1(*),IH2(*),
     +                  INFOA(*),INFOH(*),IWORK(*),
     +                  GLOB_TO_LOC(*)
      CHARACTER         DESCRA*11,DESCRH*11
      DOUBLE PRECISION  A(*),H(*)
C     .. Local scalars ..
      INTEGER           I, J, K, XBLCK, XMATR,
     +                  AUX, AUX1, AUX2, AUX3, 
     +                  LOC_COLUMN, LOC_POINTER, SHIFT1,
     +                  SHIFT2
C     .. Local arrays ..
      IERROR = 0

      DO I = 1, M
         XBLCK = IH + I - 1
         XMATR = IA + I - 1
         SHIFT1 = IA2(XMATR + 1) - IA2(XMATR)
         SHIFT2 = 2 * SHIFT1
C        If columns are already sorted, return point is 100
C         CALL MRGSRT(IA2(XMATR + 1) - IA2(XMATR),
C     +               IA1(IA2(XMATR)),
C     +               IWORK(SHIFT2 + 1),
C     +               *100)
         GOTO 100 
         K = IWORK(SHIFT2 + 1)
C     If columns have been sorted by mrgsrt
         DO J = 1, IA2(XMATR + 1) - IA2(XMATR)
            IWORK(J) = IA1(IA2(XMATR) - 1 + K)
            IWORK(SHIFT1 + J) = IA2(XMATR) - 1 + K
            K = IWORK(SHIFT2 + 1 + K)
         ENDDO
         GOTO 101
C     Else
 100     CONTINUE
         DO J = IA2(XMATR), IA2(XMATR + 1) - 1
            AUX = J - IA2(XMATR) + 1
            IWORK(AUX) = IA1(J)
            IWORK(SHIFT1 + AUX) = J
         ENDDO
C     End If 
C        Now IWORK(1: .. ) contains ordered column indices
C        and IWORK(SHIFT1 + 1: .. ) contains position of those
C        indices in the stored matrix data structures.
 101     CONTINUE
         DO J = IH2(XBLCK), IH2(XBLCK + 1) - 1
            IF ((JH .LE. IH1(J)) .AND.
     +          (IH1(J) .LE. (JH + N - 1))) THEN
               LOC_COLUMN = GLOB_TO_LOC(JA - JH + IH1(J))             
C              Binary search 
               AUX1 = 1
               AUX2 = IA2(XMATR + 1) - IA2(XMATR)
               DO K = 1, IA2(XMATR + 1) - IA2(XMATR)
                  AUX = (AUX1 + AUX2) / 2
                  IF (LOC_COLUMN .GT. IWORK(AUX)) THEN
                     AUX1 = AUX + 1
                  ELSE
                     AUX2 = AUX - 1
                  ENDIF
                  IF ((LOC_COLUMN .EQ. IWORK(AUX)) .OR.
     +                (AUX1 .GT. AUX2))
     +               EXIT
               ENDDO
               IF (LOC_COLUMN .EQ. IWORK(AUX)) THEN
                  LOC_POINTER = IWORK(SHIFT1 + AUX)
               ELSE
                  IERROR = 1
                  GOTO 9999
               ENDIF
               IF (FLAG .EQ. 0) THEN                          
                  A(LOC_POINTER) = 0.0D0
               ELSE IF (FLAG .EQ. 1) THEN
                  A(LOC_POINTER) = H(J)
               ELSE IF (FLAG .EQ. 2) THEN
                  A(LOC_POINTER) = A(LOC_POINTER) + H(J)
               ELSE
                  IERROR = 1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
 9999 RETURN
      END




