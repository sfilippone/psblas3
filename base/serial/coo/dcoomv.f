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
***********************************************************************
*  DCOOMV. Prolog to be updated.                                      *
*                                                                     *
*    FUNCTION: Driver for routines performing one of the sparse       *
*              matrix vector operations                               *
*                                                                     *
*                   y = alpha*op(A)*x + beta*y                        *
*                                                                     *
*              where op(A) is one of:                                 *
*                                                                     *
*                  op(A) = A or op(A) = A'  or                        *
*                  op(A) = lower or upper part of A                   *
*                                                                     *
*              alpha and beta are scalars.                            *
*              The data structure of the matrix is related            *
*              to the scalar computer.                                *
*              This is an internal routine called by:                 *
*              DSMMV                                                  *
*                                                                     *
*    ENTRY-POINT = DSRMV                                              *
*   INPUT =                                                           *
*                                                                     *
*                                                                     *
*      SYMBOLIC NAME: TRANS                                           *
*      POSITION:      PARAMETER NO 1.                                 *
*      ATTRIBUTES:    CHARACTER*1                                     *
*      VALUES:        'N' 'T' 'L' 'U'                                 *
*      DESCRIPTION:   Specifies the form of op(A) to be used in the   *
*                     matrix vector multiplications as follows:       *
*                                                                     *
*                     TRANS = 'N'       ,  op( A ) = A.               *
*                                                                     *
*                     TRANS = 'T'       ,  op( A ) = A'.              *
*                                                                     *
*                     TRANS = 'L' or 'U',  op( A ) = lower or         *
*                             upper part of A                         *
*                                                                     *
*      SYMBOLIC NAME: DIAG                                            *
*      POSITION:      PARAMETER NO 2.                                 *
*      ATTRIBUTES:    CHARACTER*1                                     *
*      VALUES:        'N' 'U'                                         *
*      DESCRIPTION:                                                   *
*                     Specifies whether or not the matrix A has       *
*                     unit diagonal as follows:                       *
*                                                                     *
*                     DIAG = 'N'  A is not assumed                    *
*                            to have unit diagonal                    *
*                                                                     *
*                     DIAG = 'U'  A is assumed                        *
*                            to have unit diagonal.                   *
*                                                                     *
*                     WARNING: it is the caller's responsibility      *
*                     to ensure that if the matrix has unit           *
*                     diagonal, there are no elements of the          *
*                     diagonal are stored in the arrays AS and JA.    *
*                                                                     *
*       SYMBOLIC NAME: M                                              *
*       POSITION:      PARAMETER NO 3.                                *
*       ATTRIBUTES:    INTEGER*4.                                     *
*       VALUES:        M >= 0                                         *
*       DESCRIPTION:   Number of rows of the matrix op(A).            *
*                                                                     *
*       SYMBOLIC NAME: N                                              *
*       POSITION:      PARAMETER NO 4.                                *
*       ATTRIBUTES:    INTEGER*4.                                     *
*       VALUES:        N >= 0                                         *
*       DESCRIPTION:   Number of columns of the matrix op(A)          *
*                                                                     *
*       SYMBOLIC NAME: ALPHA                                          *
*       POSITION:      PARAMETER NO 5.                                *
*       ATTRIBUTES:    REAL*8.                                        *
*       VALUES:        any.                                           *
*       DESCRIPTION:   Specifies the scalar alpha.                    *
*                                                                     *
*                                                                     *
*       SYMBOLIC NAME: AS                                             *
*       POSITION:      PARAMETER NO 6.                                *
*       ATTRIBUTES:    REAL*8: ARRAY(IA(M+1)-1)                       *
*       VALUES:        ANY                                            *
*       DESCRIPTION:   Array containing the non zero coefficients of  *
*                      the sparse matrix op(A).                       *
*                                                                     *
*       SYMBOLIC NAME: JA                                             *
*       POSITION:      PARAMETER NO 7.                                *
*       ATTRIBUTES:    INTEGER*4: ARRAY(IA(M+1)-1)                    *
*       VALUES:        0 < JA(I) <= M                                 *
*       DESCRIPTION:   Array containing the column number of the      *
*                      nonzero coefficients stored in array AS.       *
*                                                                     *
*       SYMBOLIC NAME: IA                                             *
*       POSITION:      PARAMETER NO 8.                                *
*       ATTRIBUTES:    INTEGER*4: ARRAY(*)                            *
*       VALUES:        IA(I) > 0                                      *
*       DESCRIPTION:   Contains the pointers for the beginning of     *
*                      each rows.                                     *
*                                                                     *
*                                                                     *
*       SYMBOLIC NAME: X                                              *
*       POSITION:      PARAMETER NO 9.                                *
*       ATTRIBUTES:    REAL*8 ARRAY(N) (or ARRAY(M) when op(A) = A')  *
*       VALUES:        any.                                           *
*       DESCRIPTION:   Contains the values of the vector to be        *
*                      multiplied by the matrix A.                    *
*                                                                     *
*       SYMBOLIC NAME: BETA                                           *
*       POSITION:      PARAMETER NO 10.                               *
*       ATTRIBUTES:    REAL*8.                                        *
*       VALUES:        any.                                           *
*       DESCRIPTION:   Specifies the scalar beta.                     *
*                                                                     *
*       SYMBOLIC NAME: Y                                              *
*       POSITION:      PARAMETER NO 11.                               *
*       ATTRIBUTES:    REAL*8 ARRAY(M) (or ARRAY(N) when op(A) = A')  *
*       VALUES:        any.                                           *
*       DESCRIPTION:   Contains the values of the vector to be        *
*                      updated by the matrix-vector multiplication.   *
*                                                                     *
*       SYMBOLIC NAME: WORK                                           *
*       POSITION:      PARAMETER NO 12.                               *
*       ATTRIBUTES:    REAL*8 ARRAY(M) (or ARRAY(N) when op(A) = A')  *
*       VALUES:        any.                                           *
*       DESCRIPTION:   Work area available to the program. It is used *
*                      only when TRANS = 'T'.                         *
*                                                                     *
*  OUTPUT =                                                           *
*                                                                     *
*                                                                     *
*       SYMBOLIC NAME: Y                                              *
*       POSITION:      PARAMETER NO 11.                               *
*       ATTRIBUTES:    REAL*8 ARRAY(M) (or ARRAY(N) when op(A) = A')  *
*       VALUES:        any.                                           *
*       DESCRIPTION:   Contains the values of the vector              *
*                      updated by the matrix-vector multiplication.   *
*                                                                     *
*                                                                     *
***********************************************************************
      SUBROUTINE DCOOMV (TRANS,DIAG,M,N,ALPHA,AS,IA,JA,INFOA,X,
     +  BETA,Y,WORK,IERROR)
      use psb_const_mod
C     .. Parameters ..
      real(psb_dpk_)  ONE, ZERO
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      real(psb_dpk_)  ALPHA, BETA
      INTEGER           M, N, IERROR
      CHARACTER         DIAG, TRANS
C     .. Array Arguments ..
      real(psb_dpk_)  AS(*), WORK(*), X(*), Y(*)
      INTEGER           IA(*), JA(*),infoa(*)
C     .. Local Scalars ..
      real(psb_dpk_)  ACC,  TX
      INTEGER           I, J, K, NNZ, IR, JC
      LOGICAL           SYM, TRA, UNI
C     .. Executable Statements ..
C
      IERROR=0
      UNI = (DIAG.EQ.'U')
      TRA = (TRANS.EQ.'T')

C     Symmetric matrix upper or lower 
      SYM = ((TRANS.EQ.'L').OR.(TRANS.EQ.'U'))
C

      IF (ALPHA.EQ.ZERO) THEN
        IF (BETA.EQ.ZERO) THEN
          DO I = 1, M
            Y(I) = ZERO
          ENDDO
        ELSE
          DO 20 I = 1, M
            Y(I) = BETA*Y(I)
 20       CONTINUE
        ENDIF
        RETURN
      END IF

      NNZ = INFOA(1)
C
      IF (SYM) THEN
        IF (UNI) THEN
C
C              ......Symmetric with unitary diagonal.......
C              ....OK!!
C              To be optimized
          
          IF (BETA.NE.ZERO) THEN
            DO  I = 1, M
C
C                 Product for diagonal elements
C
              Y(I) = BETA*Y(I) + ALPHA*X(I)
            ENDDO
          ELSE
            DO I = 1, M
              Y(I) = ALPHA*X(I)
            ENDDO
          ENDIF

C              Product for other elements


          I    = 1
          J    = I
          DO WHILE (I.LE.NNZ)
            
            DO WHILE ((IA(J).EQ.IA(I)).AND.
     +        (J.LE.NNZ))
              J = J+1
            ENDDO
            
            ACC = ZERO
            IR = IA(I)
            TX = X(IR)
            DO K = I, J-1
              JC = JA(K)
              ACC = ACC + AS(K)*X(JC)
              Y(JC) = Y(JC) + ALPHA * AS(K)*TX
            ENDDO        
            Y(IR) = Y(IR) + ALPHA * ACC
            I = J 
          ENDDO
C
        ELSE IF ( .NOT. UNI) THEN
C
C            Check if matrix is lower or upper
C
          IF (TRANS.EQ.'L') THEN
C
C               LOWER CASE: diagonal element is the last element of row
C               ....OK!

            IF (BETA.NE.ZERO) THEN
              DO I = 1, M
                Y(I) = BETA*Y(I)
              ENDDO
            ELSE
              DO I = 1, M
                Y(I) = ZERO
              ENDDO
            ENDIF

            I    = 1
            J    = I
            DO WHILE (I.LE.NNZ)              
              DO WHILE ((IA(J).EQ.IA(I)).AND.
     +          (J.LE.NNZ))
                J = J+1
              ENDDO              
              ACC = ZERO
              IR = IA(I)
              TX = X(IR)
              DO K = I, J-1
                JC = JA(K)
                ACC = ACC + AS(K)*X(JC)
                IF (IR.NE.JC) THEN 
                  Y(JC) = Y(JC) + ALPHA * AS(K)*TX
                ENDIF
              ENDDO        
              Y(IR) = Y(IR) + ALPHA * ACC
              I = J 
            ENDDO


          ELSE                  ! ....Trans<>L
C
C              UPPER CASE
C              ....OK!! (Actually it's just the same as above!)
C

            IF (BETA.NE.ZERO) THEN
              DO I = 1, M
                Y(I) = BETA*Y(I)
              ENDDO
            ELSE
              DO I = 1, M
                Y(I) = ZERO
              ENDDO
            ENDIF

            I    = 1
            J    = I
            DO WHILE (I.LE.NNZ)              
              DO WHILE ((IA(J).EQ.IA(I)).AND.
     +          (J.LE.NNZ))
                J = J+1
              ENDDO              
              ACC = ZERO
              IR = IA(I)
              TX = X(IR)
              DO K = I, J-1
                JC = JA(K)
                ACC = ACC + AS(K)*X(JC)
                IF (IR.NE.JC) THEN 
                  Y(JC) = Y(JC) + ALPHA * AS(K)*TX
                ENDIF
              ENDDO        
              Y(IR) = Y(IR) + ALPHA * ACC
              I = J 
            ENDDO

          END IF                ! ......TRANS=='L'

        END IF                  ! ......Not UNI
C
      ELSE IF ( .NOT. TRA) THEN !......NOT SYM

        IF ( .NOT. UNI) THEN      
C
C          .......General Not Unit, No Traspose
C
          IF (BETA.NE.ZERO) THEN
            DO I = 1, M
              Y(I) = BETA*Y(I)
            ENDDO
          ELSE
            DO I = 1, M
              Y(I) = ZERO
            ENDDO
          ENDIF

          I    = 1
          J    = I
          IF (nnz > 0) then 
            IR   = IA(1) 
            ACC  = zero
            DO 
              if (i>nnz) then 
                Y(IR) = Y(IR) + ALPHA * ACC
                exit
              endif
              IF (IA(I) /= IR) THEN 
                Y(IR) = Y(IR) + ALPHA * ACC
                IR    = IA(I) 
                ACC   = ZERO
              ENDIF
              ACC     = ACC + AS(I) * X(JA(I))
              I       = I + 1               
          ENDDO
        endif
C
        ELSE IF (UNI) THEN
C

          IF (BETA.NE.ZERO) THEN
            DO I = 1, M
              Y(I) = BETA*Y(I)+ALPHA*X(I)
            ENDDO
          ELSE
            DO I = 1, M
              Y(I) = ALPHA*X(I)
            ENDDO
          ENDIF

          I    = 1
          J    = I
          DO WHILE (I.LE.NNZ)              
            DO WHILE ((IA(J).EQ.IA(I)).AND.
     +        (J.LE.NNZ))
              J = J+1
            ENDDO              
            ACC = ZERO
            IR = IA(I)
            DO K = I, J-1
              JC = JA(K)
              ACC = ACC + AS(K)*X(JC)
            ENDDO        
            Y(IR) = Y(IR) + ALPHA * ACC
            I = J 
          ENDDO
          
        END IF                  !....End Testing on UNI
C
      ELSE IF (TRA) THEN        !....Else on SYM (swapped M and N)
C
        IF ( .NOT. UNI) THEN
C
          IF (BETA.NE.ZERO) THEN
            DO I = 1, M
              Y(I) = BETA*Y(I)
            ENDDO
          ELSE
            DO I = 1, M
              Y(I) = ZERO
            ENDDO
          ENDIF
C
        ELSE IF (UNI) THEN
C
          
          IF (BETA.NE.ZERO) THEN
            DO I = 1, M
              Y(I) = BETA*Y(I)+ALPHA*X(I)
            ENDDO
          ELSE
            DO I = 1, M
              Y(I) = ALPHA*X(I)
            ENDDO
          ENDIF          
C
        END IF                  !....UNI
C
        IF (ALPHA.EQ.ONE) THEN
C
          I    = 1
          DO I=1,NNZ
            IR = JA(I)
            JC = IA(I)
            Y(IR) = Y(IR) +  AS(I)*X(JC)
          ENDDO
C
        ELSE IF (ALPHA.EQ.-ONE) THEN
C

          DO I=1,NNZ
            IR = JA(I)
            JC = IA(I)
            Y(IR) = Y(IR) - AS(I)*X(JC)
          ENDDO
C
        ELSE                    !.....Else on TRA

          DO I=1,M
            WORK(I) =  ALPHA*X(I)
          ENDDO

          DO I=1,NNZ
            IR = JA(I)
            JC = IA(I)
            Y(IR) = Y(IR) + AS(I)*WORK(JC)
          ENDDO

        END IF                  !.....End testing on ALPHA

      END IF                    !.....End testing on SYM
C
      RETURN
C
C     END OF DSRMV
C
      END
     
