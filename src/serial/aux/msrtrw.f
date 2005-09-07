C  msrtrw.f
C  Author: Carlo Vittoli 
C  Date:   May 19, 1994
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Subroutine msrtrw                                                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Purpose: Sort rows by column indices (merge sorting)                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Parameters:                                                         C
C     Input:                                                           C
C                                                                      C
C     Output:                                                          C
C                                                                      C
C     Others:                                                          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Algorithm:                                                          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  References:                                                         C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE  msrtrw(m,a,ia1,ia2,work,lwork,awork,lawork,ierrv)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           m, lwork, lawork
C     .. Array Arguments ..
      DOUBLE PRECISION  a(*), awork(m)
      INTEGER           ia1(*), ia2(*), work(2*m), ierrv(*)
C     .. Local Scalars ..
      INTEGER           i1, i, j2, lrow, jend, jstart, irow, iret
C     .. External Subroutines ..
      EXTERNAL          xsperr
C     .. Executable Statements ..
      if(lawork.lt.m) then
         call xsperr('LWORK   ',lawork,8,'MSRTRW',IERRV)
         goto 9999
      endif
      if(lwork.lt.2*m) then
         call xsperr('LWORK   ',lwork,6,'MSRTRW',IERRV)
         goto 9999
      endif
C     Start
      DO 160 IROW = 1, M
         JSTART = ia2(IROW)
         JEND = ia2(IROW+1) - 1
         LROW = JEND - JSTART + 1
         call mrgsrt(lrow,ia1(jstart),work,iret)
         if (iret.eq.0) then
           I1 = WORK(1)
           DO 20 I = 1, LROW
             work(m+I) = I1 + JSTART - 1
             I1 = WORK(I1+1)
c            if (i1+1.gt.m) stop 'AAAAAAAArgh!!!!!!!!!'
 20        CONTINUE
           DO 40 I = 1, LROW
             WORK(I) = ia1(work(m+I))
             AWORK(I) = A(work(m+I))
 40        CONTINUE
           DO 60 I = 1, LROW
             J2 = I + JSTART - 1
             A(J2) = AWORK(I)
             ia1(J2) = WORK(I)
 60        CONTINUE
         endif
  160 continue
 9999 continue
      RETURN
      END
