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
      SUBROUTINE DCOOZERO(M,N,DESCRA,A,IA1,IA2,
     +   INFOA,IA,JA,MZ,NZ,IERROR)
C
C     This subroutione performs the operation:
C
C     A(IA : IA + MZ - 1, JA : JA + NZ - 1) = 0
C
C     This isn't accomplished by removing elements 
C     from sparse matrix representation, but assigning them
C     the zero value.
C     Columns are supposed to be ordered 
C     into the same row. This subroutine will
C     not work properly otherwise.
C
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           M,N,IA,JA,MZ,NZ,IERROR
C     .. Array Arguments ..
      INTEGER           IA1(*),IA2(*),INFOA(*)
      CHARACTER         DESCRA*11
      DOUBLE PRECISION  A(*)
C     .. Local scalars ..
      INTEGER I, J, JBEGIN, JEND, AUX, NNZ
      DOUBLE PRECISION 

      IERROR=0
      IF (((JA + NZ - 1) .GT. N) .OR.
     +    ((IA + MZ - 1) .GT. M) .OR.
     +    (IA .LT. 1) .OR. (JA .LT. 1)) THEN
         IERROR = 1
         GOTO 9999
      ENDIF   
      NNZ = INFOA(1)
      I   = 1
      DO WHILE ((IA1(I).LT.IA).AND.(I.LE.NNZ))
        I = I + 1
      ENDDO
      DO WHILE ((IA1(I).LE.(IA+MZ-1)).AND.(I.LE.NNZ))
        IF ((JA.LE.IA2(I)).AND.(IA2(I).LE.(JA+NZ-1))) THEN
          A(I) = 0.0D0
        ENDIF
        I = I + 1 
      ENDDO

      RETURN
      END
