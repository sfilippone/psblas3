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
***********************************************************************
*    DSRMV modified for SPARKER
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
      SUBROUTINE DCSRMV(TRANS,DIAG,M,N,ALPHA,AS,JA,IA,X,BETA,Y,
     +  WORK,LWORK,IERROR)
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, BETA
      INTEGER           M, N,LWORK,IERROR
      CHARACTER         DIAG, TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION  AS(*), WORK(*), X(*), Y(*)
      INTEGER           IA(*), JA(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC
      INTEGER           I, J, K, NCOLA, NROWA
      LOGICAL           SYM, TRA, UNI
C     .. Executable Statements ..
C
      IERROR = 0
      UNI = (DIAG.EQ.'U')
      TRA = (TRANS.EQ.'T')

C     Symmetric matrix upper or lower 
      SYM = ((TRANS.EQ.'L').OR.(TRANS.EQ.'U'))
C
      IF ( .NOT. TRA) THEN
        NROWA = M
        NCOLA = N
      ELSE IF (TRA) THEN
        NROWA = N
        NCOLA = M
      END IF                    !(....TRA)

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

C
      IF (SYM) THEN
        IF (UNI) THEN
C
C              ......Symmetric with unitary diagonal.......
C              ....OK!!
C              To be optimized
          
          IF (BETA.NE.ZERO) THEN
            DO 40 I = 1, M
C
C                 Product for diagonal elements
C
              Y(I) = BETA*Y(I) + ALPHA*X(I)
 40         CONTINUE
          ELSE
            DO I = 1, M
              Y(I) = ALPHA*X(I)
            ENDDO
          ENDIF

C              Product for other elements
          DO 80 I = 1, M
            ACC = ZERO
            DO 60 J = IA(I), IA(I+1) - 1
              K = JA(J)
              Y(K) = Y(K) + ALPHA*AS(J)*X(I)
              ACC = ACC + AS(J)*X(K)
 60         CONTINUE
            Y(I) = Y(I) + ALPHA*ACC
 80       CONTINUE
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
              DO 100 I = 1, M
                Y(I) = BETA*Y(I)
 100          CONTINUE
            ELSE
              DO I = 1, M
                Y(I) = ZERO
              ENDDO
            ENDIF

            DO 140 I = 1, M
              ACC = ZERO
              DO 120 J = IA(I), IA(I+1) - 1 ! it was -2
                K = JA(J)
C
C                 To be optimized
C
                IF (K.NE.I) THEN !for symmetric element 
                  Y(K) = Y(K) + ALPHA*AS(J)*X(I)
                ENDIF
                ACC = ACC + AS(J)*X(K)
 120          CONTINUE

              Y(I) = Y(I) + ALPHA*ACC 
 140        CONTINUE
          ELSE                  ! ....Trans<>L
C
C              UPPER CASE
C              ....OK!!
C
            IF (BETA.NE.ZERO) THEN
              DO 160 I = 1, M
                Y(I) = BETA*Y(I)
 160          CONTINUE
            ELSE
              DO I = 1, M
                Y(I) = ZERO
              ENDDO
            ENDIF

            DO 200 I = 1, M
              ACC = ZERO
              DO 180 J = IA(I) , IA(I+1) - 1 ! removed +1
                K = JA(J)
C
C                 To be optimized
C
                IF(K.NE.I) THEN
                  Y(K) = Y(K) + ALPHA*AS(J)*X(I)
                ENDIF
                ACC = ACC + AS(J)*X(K)
 180          CONTINUE
              Y(I) = Y(I) + ALPHA*ACC
 200        CONTINUE
          END IF                ! ......TRANS=='L'
        END IF                  ! ......Not UNI
C
      ELSE IF ( .NOT. TRA) THEN !......NOT SYM

        IF ( .NOT. UNI) THEN      
C
C          .......General Not Unit, No Traspose
C
          if (beta  == zero) then
            if (alpha==one) then 
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j))
                enddo
                y(i) = acc
              enddo
            else if (alpha==-one) then 
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc - as(j)*x(ja(j))
                enddo
                y(i) = acc
              enddo
            else
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j))
                enddo
                y(i) = alpha*acc
              enddo
              
            endif 

          else if (beta==one) then 

            if (alpha==one) then 
              do i = 1, m
                acc = y(i)
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j))
                enddo
                y(i) = acc
              enddo
            else if (alpha==-one) then 
              do i = 1, m
                acc = y(i)
                do j = ia(i), ia(i+1) - 1
                  acc = acc - as(j)*x(ja(j))
                enddo
                y(i) = acc
              enddo
            else
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j))
                enddo
                y(i) = alpha*acc + y(i)
              enddo
            endif 

          else if (beta==-one) then 

            if (alpha==one) then 
              do i = 1, m
                acc = -y(i)
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j))
                enddo
                y(i) = acc
              enddo
            else if (alpha==-one) then 
              do i = 1, m
                acc = -y(i)
                do j = ia(i), ia(i+1) - 1
                  acc = acc - as(j)*x(ja(j))
                enddo
                y(i) = acc
              enddo
            else
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j))
                enddo
                y(i) = alpha*acc - y(i)
              enddo
            endif 
          else  
            if (alpha==one) then 
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j))
                enddo
                y(i) = acc + beta*y(i)
              enddo
            else if (alpha==-one) then 
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc - as(j)*x(ja(j))
                enddo
                y(i) = acc + beta*y(i)
              enddo
            else
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j))
                enddo
                y(i) = alpha*acc + beta*y(i)
              enddo
            endif 
          end if 
C
        ELSE IF (UNI) THEN
C
          IF (BETA.NE.ZERO) THEN
            DO 280 I = 1, M
              ACC = ZERO
              DO 260 J = IA(I), IA(I+1) - 1
                ACC = ACC + AS(J)*X(JA(J))
 260          CONTINUE
              Y(I) = ALPHA*(ACC+X(I)) + BETA*Y(I)
 280        CONTINUE
          ELSE                  !(BETA.EQ.ZERO)
            DO I = 1, M
              ACC = ZERO
              DO J = IA(I), IA(I+1) - 1
                ACC = ACC + AS(J)*X(JA(J))
              ENDDO
              Y(I) = ALPHA*(ACC+X(I))
            ENDDO
          ENDIF
        END IF                  !....End Testing on UNI
C
      ELSE IF (TRA) THEN        !....Else on SYM (swapped M and N)
C
        IF ( .NOT. UNI) THEN
C
          IF (BETA.NE.ZERO) THEN
            DO 300 I = 1, M
              Y(I) = BETA*Y(I)
 300        CONTINUE
          ELSE                  !(BETA.EQ.ZERO)
            DO I = 1, M
              Y(I) = ZERO
            ENDDO
          ENDIF
C
        ELSE IF (UNI) THEN
C

          IF (BETA.NE.ZERO) THEN
            DO 320 I = 1, M
              Y(I) = BETA*Y(I) + ALPHA*X(I)
 320        CONTINUE
          ELSE                  !(BETA.EQ.ZERO)
            DO I = 1, M
              Y(I) = ALPHA*X(I)
            ENDDO
          ENDIF

C
        END IF                  !....UNI
C
        IF (ALPHA.EQ.ONE) THEN
C
          DO 360 I = 1, N
            DO 340 J = IA(I), IA(I+1) - 1
              K = JA(J)
              Y(K) = Y(K) + AS(J)*X(I)
 340        CONTINUE
 360      CONTINUE
C
        ELSE IF (ALPHA.EQ.-ONE) THEN
C
          DO 400 I = 1, n
            DO 380 J = IA(I), IA(I+1) - 1
              K = JA(J)
              Y(K) = Y(K) - AS(J)*X(I)
 380        CONTINUE
 400      CONTINUE
C
        ELSE                    !.....Else on TRA
C
C           This work array is used for efficiency
C
          IF (LWORK.LT.N) THEN
            IERROR = 60
            WORK(1) = DBLE(N)
            RETURN
          ENDIF
          DO 420 I = 1, N
            WORK(I) = ALPHA*X(I)
 420      CONTINUE
C
          DO 460 I = 1, n
            DO 440 J = IA(I), IA(I+1) - 1
              K = JA(J)
              Y(K) = Y(K) + AS(J)*WORK(I)
 440        CONTINUE
 460      CONTINUE
C
        END IF                  !.....End testing on ALPHA

      END IF                    !.....End testing on SYM
C
      RETURN
C
C     END OF DSRMV
C
      END
     
