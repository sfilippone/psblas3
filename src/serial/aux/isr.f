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
      SUBROUTINE ISR(N,X)
C
C  Quicksort.
C  Adapted from a number of sources, including Don Knuth's TAOCP.
C
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER X(N)
C     ..
C     .. Local Scalars ..
      INTEGER I, J, XX, ILX, IUX, ISTP, PIV, LPIV
      INTEGER IT1, N1, N2
      INTEGER MAXSTACK,NPARMS,ITHRS
      PARAMETER (MAXSTACK=64,NPARMS=3,ITHRS=16)
      INTEGER ISTACK(NPARMS,MAXSTACK)
C     ..

C
C     Small inputs will only get through insertion sort. 
C
      IF (N.GT.ITHRS) THEN          
C
C     Init stack pointer
C
        ISTP = 1
        ISTACK(1,ISTP) = 1
        ISTACK(2,ISTP) = N
        
        DO WHILE (ISTP.GT.0)             
          ILX  = ISTACK(1,ISTP)
          IUX  = ISTACK(2,ISTP)
          ISTP = ISTP - 1
c$$$            write(0,*) 'Debug 1: ',ilx,iux
C
C       Choose a pivot with median-of-three heuristics, leave it 
C       in the LPIV location
C            
          I = ILX
          J = IUX 
          LPIV = (I+J)/2
          PIV  = X(LPIV)
          IF (PIV.LT.X(I)) THEN
            IT1 = X(I)
            X(I) = X(LPIV)
            X(LPIV) = IT1
            PIV = X(LPIV)
          ENDIF
          IF (PIV.GT.X(J)) THEN
            IT1 = X(J)
            X(J) = X(LPIV)
            X(LPIV) = IT1
            PIV = X(LPIV)
          ENDIF
          IF (PIV.LT.X(I)) THEN
            IT1 = X(I)
            X(I) = X(LPIV)
            X(LPIV) = IT1
            PIV = X(LPIV)
          ENDIF
C
C     Now PIV is correct;  place it into first location

          IT1 = X(I)
          X(I) = X(LPIV)
          X(LPIV) = IT1
          
          I = ILX - 1 
          J = IUX + 1 
          
 130      CONTINUE
          I = I + 1
          XK = X(I)
          IF (XK.LT.PIV) GOTO 130
C
C     Ensure finite termination for next loop
C
          IT1  = XK
          X(I) = PIV
 140      CONTINUE
          J = J - 1
          XK = X(J)
          IF (XK.GT.PIV) GOTO 140
          X(I) = IT1  
 150      CONTINUE
          
          IF (J.GT.I) THEN
            IT1  = X(I)
            X(I) = X(J)
            X(J) = IT1 
            GO TO 130
          END IF
          
          if (i.eq.ilx) then 
            if (x(i).ne.piv) then
              write(0,*) 'Should never ever get here????!!!!'
              stop
            endif
            i = i + 1 
          endif
          
          N1 = (I-1)-ILX+1
          N2 = IUX-(I)+1
          IF (N1.GT.N2) THEN
            if (n1.gt.ithrs) then 
              ISTP = ISTP + 1
              ISTACK(1,ISTP) = ILX
              ISTACK(2,ISTP) = I-1
            endif
            if (n2.gt.ithrs) then
              ISTP = ISTP + 1
              ISTACK(1,ISTP) = I
              ISTACK(2,ISTP) = IUX
            endif
          ELSE
            if (n2.gt.ithrs) then
              ISTP = ISTP + 1
              ISTACK(1,ISTP) = I
              ISTACK(2,ISTP) = IUX
            endif
            if (n1.gt.ithrs) then 
              ISTP = ISTP + 1
              ISTACK(1,ISTP) = ILX
              ISTACK(2,ISTP) = I-1
            endif
          ENDIF               
        ENDDO
      ENDIF
      
      DO J=N-1,1,-1
        IF (X(J+1).LT.X(J)) THEN
          XX = X(J)
          I=J+1
 100      CONTINUE
          X(I-1) = X(I)
          I = I+1
          IF ((I.LE.N)) then 
            if (X(I).LT.XX) GOTO 100
          endif 
          X(I-1) = XX
        ENDIF
      ENDDO
      
      RETURN

      END
