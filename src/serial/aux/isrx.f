      SUBROUTINE ISRX(N,X,INDX)
C
C  Quicksort with indices into original positions.
C  Adapted from a number of sources, including Don Knuth's TAOCP.
C
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER INDX(N),X(N)
C     ..
C     .. Local Scalars ..
      INTEGER I, J, II, XX, ILX, IUX, ISTP, PIV, LPIV
      INTEGER IT1, IT2, N1, N2
      INTEGER MAXSTACK,NPARMS,ITHRS
      PARAMETER (MAXSTACK=64,NPARMS=3,ITHRS=16)
      INTEGER ISTACK(NPARMS,MAXSTACK)
C     ..

      DO I=1, N
        INDX(I) = I
      ENDDO      
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
            IT2 = INDX(I)
            X(I) = X(LPIV)
            INDX(I) = INDX(LPIV)
            X(LPIV) = IT1
            INDX(LPIV) = IT2
            PIV = X(LPIV)
          ENDIF
          IF (PIV.GT.X(J)) THEN
            IT1 = X(J)
            IT2 = INDX(J)
            X(J) = X(LPIV)
            INDX(J) = INDX(LPIV)
            X(LPIV) = IT1
            INDX(LPIV) = IT2
            PIV = X(LPIV)
          ENDIF
          IF (PIV.LT.X(I)) THEN
            IT1 = X(I)
            IT2 = INDX(I)
            X(I) = X(LPIV)
            INDX(I) = INDX(LPIV)
            X(LPIV) = IT1
            INDX(LPIV) = IT2
            PIV = X(LPIV)
          ENDIF
C
C     Now PIV is correct;  place it into first location
C            
          IT1 = X(I)
          IT2 = INDX(I)
          X(I) = X(LPIV)
          INDX(I) = INDX(LPIV)
          X(LPIV) = IT1
          INDX(LPIV) = IT2
          
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
            IT2   = INDX(I)
            X(I) = X(J)
            INDX(I) = INDX(J)
            X(J) = IT1 
            INDX(J) = IT2  
            GO TO 130
          END IF
          
          if (i.eq.ilx) then 
            if (x(i).ne.piv) then
              write(0,*)
     +          'ISRX:: Should never ever get here????!!!!'
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
          II = INDX(J)
          I=J+1
 100      CONTINUE
          X(I-1) = X(I)
          INDX(I-1) = INDX(I)
          I = I+1
          IF ((I.LE.N)) then 
            if (X(I).LT.XX) GOTO 100
          endif
          X(I-1) = XX
          INDX(I-1) = II
        ENDIF
      ENDDO
      
      RETURN

      END
