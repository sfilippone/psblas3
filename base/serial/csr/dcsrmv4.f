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
      SUBROUTINE DCSRMV4(TRANS,DIAG,M,N,ALPHA,AS,JA,IA,X,LDX,
     +  BETA,Y,LDY, WORK,LWORK,IERROR)
      integer  nb
      parameter (nb=4)
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, BETA
      INTEGER           M, N,LWORK,IERROR,ldx,ldy
      CHARACTER         DIAG, TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION  AS(*), WORK(*), X(LDX,NB), Y(LDY,NB)
      INTEGER           IA(*), JA(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC(nb)
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
      if ( .not. tra) then
        nrowa = m
        ncola = n
      else if (tra) then
        nrowa = n
        ncola = m
      end if                    !(....tra)

      if (alpha.eq.zero) then
        if (beta.eq.zero) then
          do i = 1, m
            y(i,1:nb) = zero
          enddo
        else
          do 20 i = 1, m
            y(i,1:nb) = beta*y(i,1:nb)
 20       continue
        endif
        return
      end if

c
      if (sym) then
        if (uni) then
c
c              ......Symmetric with unitary diagonal.......
C              ....OK!!
C              To be optimized
          
          if (beta.ne.zero) then
            do i = 1, m
C
C                 Product for diagonal elements
c
              y(i,1:nb) = beta*y(i,1:nb) + alpha*x(i,1:nb)
            enddo
          else
            do i = 1, m
              y(i,1:nb) = alpha*x(i,1:nb)
            enddo
          endif

C              Product for other elements
          do 80 i = 1, m
            acc = zero
            do 60 j = ia(i), ia(i+1) - 1
              k = ja(j)
              y(k,1:nb) = y(k,1:nb) + alpha*as(j)*x(i,1:nb)
              acc(1:nb) = acc(1:nb) + as(j)*x(k,1:nb)
 60         continue
            y(i,1:nb) = y(i,1:nb) + alpha*acc(1:nb)
 80       continue
C
        else if ( .not. uni) then
C
C            Check if matrix is lower or upper
C
          if (trans.eq.'L') then
C
C               LOWER CASE: diagonal element is the last element of row
C               ....OK!

            if (beta.ne.zero) then
              do 100 i = 1, m
                y(i,1:nb) = beta*y(i,1:nb)
 100          continue
            else
              do i = 1, m
                y(i,1:nb) = zero
              enddo
            endif

            do 140 i = 1, m
              acc = zero
              do 120 j = ia(i), ia(i+1) - 1 ! it was -2
                K = ja(j)
C
C                 To be optimized
C
                if (k.ne.i) then !for symmetric element 
                  y(k,1:nb) = y(k,1:nb) + alpha*as(j)*x(i,1:nb)
                endif
                acc(1:nb) = acc(1:nb) + as(j)*x(k,1:nb)
 120          continue

              y(i,1:nb) = y(i,1:nb) + alpha*acc(1:nb)
 140        continue
          else                  ! ....Trans<>L
C
C              UPPER CASE
C              ....OK!!
C
            if (beta.ne.zero) then
              do 160 i = 1, m
                y(i,1:nb) = beta*y(i,1:nb)
 160          continue
            else
              do i = 1, m
                y(i,1:nb) = zero
              enddo
            endif

            do 200 i = 1, m
              acc = zero
              do 180 j = ia(i) , ia(i+1) - 1 ! removed +1
                k = ja(j)
C
C                 To be optimized
C
                if (k.ne.i) then
                  y(k,1:nb) = y(k,1:nb) + alpha*as(j)*x(i,1:nb)
                endif
                acc(1:nb) = acc(1:nb) + as(j)*x(k,1:nb)
 180          continue
              y(i,1:nb) = y(i,1:nb) + alpha*acc(1:nb)
 200        continue
          end if                ! ......TRANS=='L'
        end if                  ! ......Not UNI
c
      else if ( .not. tra) then !......NOT SYM

        if ( .not. uni) then      
C
C          .......General Not Unit, No Traspose
C
          if (beta  == zero) then
            if (alpha==one) then 
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j),1:nb)
                enddo
                y(i,1:nb) = acc
              enddo
            else if (alpha==-one) then 
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc - as(j)*x(ja(j),1:nb)
                enddo
                y(i,1:nb) = acc
              enddo
            else
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j),1:nb)
                enddo
                y(i,1:nb) = alpha*acc
              enddo
              
            endif 

          else if (beta==one) then 

            if (alpha==one) then 
              do i = 1, m
                acc = y(i,1:nb)
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j),1:nb)
                enddo
                y(i,1:nb) = acc
              enddo
            else if (alpha==-one) then 
              do i = 1, m
                acc = y(i,1:nb)
                do j = ia(i), ia(i+1) - 1
                  acc = acc - as(j)*x(ja(j),1:nb)
                enddo
                y(i,1:nb) = acc
              enddo
            else
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j),1:nb)
                enddo
                y(i,1:nb) = alpha*acc + y(i,1:nb)
              enddo
            endif 

          else if (beta==-one) then 

            if (alpha==one) then 
              do i = 1, m
                acc = -y(i,1:nb)
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j),1:nb)
                enddo
                y(i,1:nb) = acc
              enddo
            else if (alpha==-one) then 
              do i = 1, m
                acc = -y(i,1:nb)
                do j = ia(i), ia(i+1) - 1
                  acc = acc - as(j)*x(ja(j),1:nb)
                enddo
                y(i,1:nb) = acc
              enddo
            else
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j),1:nb)
                enddo
                y(i,1:nb) = alpha*acc - y(i,1:nb)
              enddo
            endif 
          else  
            if (alpha==one) then 
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j),1:nb)
                enddo
                y(i,1:nb) = acc + beta*y(i,1:nb)
              enddo
            else if (alpha==-one) then 
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc - as(j)*x(ja(j),1:nb)
                enddo
                y(i,1:nb) = acc + beta*y(i,1:nb)
              enddo
            else
              do i = 1, m
                acc = zero
                do j = ia(i), ia(i+1) - 1
                  acc = acc + as(j)*x(ja(j),1:nb)
                enddo
                y(i,1:nb) = alpha*acc + beta*y(i,1:nb)
              enddo
            endif 
          end if 

c
        else if (uni) then
c
          if (beta.ne.zero) then
            do 280 i = 1, m
              acc(1:nb) = zero
              do 260 j = ia(i), ia(i+1) - 1
                acc(1:nb) = acc(1:nb) + as(j)*x(ja(j),1:nb)
 260          continue
              y(i,1:nb) = alpha*(acc(1:nb)+x(i,1:nb)) + beta*y(i,1:nb)
 280        continue
          else                  !(beta.eq.zero)
            do i = 1, m
              acc(1:nb) = zero
              do j = ia(i), ia(i+1) - 1
                acc(1:nb) = acc(1:nb) + as(j)*x(ja(j),1:nb)
              enddo
              y(i,1:nb) = alpha*(acc(1:nb)+x(i,1:nb))
            enddo
          endif
        end if                  !....End Testing on UNI
C
      else if (tra) then        !....Else on SYM (swapped M and N)
C
        if ( .not. uni) then
c
          if (beta.ne.zero) then
            do 300 i = 1, m
              y(i,1:nb) = beta*y(i,1:nb)
 300        continue
          else                  !(BETA.EQ.ZERO)
            do i = 1, m
              y(i,1:nb) = zero
            enddo
          endif
c
        else if (uni) then
c

          if (beta.ne.zero) then
            do 320 i = 1, m
              y(i,1:nb) = beta*y(i,1:nb) + alpha*x(i,1:nb)
 320        continue
          else                  !(BETA.EQ.ZERO)
            do i = 1, m
              y(i,1:nb) = alpha*x(i,1:nb)
            enddo
          endif

c
        end if                  !....UNI
C
        if (alpha.eq.one) then
c
          do 360 i = 1, n
            do 340 j = ia(i), ia(i+1) - 1
              k = ja(j)
              y(k,1:nb) = y(k,1:nb) + as(j)*x(i,1:nb)
 340        continue
 360      continue
c
        else if (alpha.eq.-one) then
c
          do 400 i = 1, n
            do 380 j = ia(i), ia(i+1) - 1
              k = ja(j)
              y(k,1:nb) = y(k,1:nb) - as(j)*x(i,1:nb)
 380        continue
 400      continue
c
        else                    !.....Else on TRA
C
C           This work array is used for efficiency
C
          if (lwork.lt.n) then
            ierror = 60
            work(1) = dble(n)
            return
          endif
c$$$          do 420 i = 1, n
c$$$            work(i) = alpha*x(i,1:4)
c$$$ 420      continue
c$$$C
c$$$          DO 460 I = 1, n
c$$$            DO 440 J = IA(I), IA(I+1) - 1
c$$$              K = JA(J)
c$$$              Y(K) = Y(K) + AS(J)*WORK(I)
c$$$ 440        CONTINUE
c$$$ 460      CONTINUE
c
        end if                  !.....end testing on alpha

      end if                    !.....end testing on sym
c
      return
c
c     end of dsrmv
c
      end
     
