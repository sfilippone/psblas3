      subroutine reordvn(nnz,ar,ia1,ia2,idx)
      integer nnz
      integer ia1(*),ia2(*),idx(0:*)
      double precision ar(*)
      integer lp, kk, swapia1, swapia2, lswap
      double precision swapar

      LP = IDX(0)
      KK = 1
 500  CONTINUE
      IF ((LP.EQ.0).OR.(KK.GT.NNZ)) GOTO 800
 600  CONTINUE
      IF (LP.GE.KK) GOTO 700
      LP = IDX(LP)
      GOTO 600
 700  CONTINUE
C        ... Swap of vectors IA2, IA1, AR ...
      SWAPIA2 = IA2(KK)
      SWAPIA1 = IA1(KK)
      SWAPAR  = AR(KK)
      IA2(KK) = IA2(LP)
      IA1(KK) = IA1(LP)
      AR(KK)  = AR(LP)
      IA2(LP) = SWAPIA2
      IA1(LP) = SWAPIA1
      AR(LP)  = SWAPAR
      LSWAP   = IDX(LP)
      IDX(LP) = IDX(KK)
      IDX(KK) = LP
      LP      = LSWAP
      KK      = KK+1
      GOTO 500
 800  CONTINUE
      return
      end
      subroutine ireordv1(nnz,ia1,idx)
      integer nnz
      integer ia1(*),idx(0:*)
      integer lp, kk, swapia1, lswap

      LP = IDX(0)
      KK = 1
 500  CONTINUE
      IF ((LP.EQ.0).OR.(KK.GT.NNZ)) GOTO 800
 600  CONTINUE
      IF (LP.GE.KK) GOTO 700
      LP = IDX(LP)
      GOTO 600
 700  CONTINUE
C        ... Swap of vectors IA2, IA1, AR ...
      SWAPIA1 = IA1(KK)
      IA1(KK) = IA1(LP)
      IA1(LP) = SWAPIA1
      LSWAP   = IDX(LP)
      IDX(LP) = IDX(KK)
      IDX(KK) = LP
      LP      = LSWAP
      KK      = KK+1
      GOTO 500
 800  CONTINUE
      return
      end
      subroutine reordvn3(nnz,ar,ia1,ia2,ia3,idx)
      integer nnz
      integer ia1(*),ia2(*),ia3(*),idx(0:*)
      double precision ar(*)
      integer lp, kk, swapia1, swapia2, swapia3,lswap
      double precision swapar

      LP = IDX(0)
      KK = 1
 500  CONTINUE
      IF ((LP.EQ.0).OR.(KK.GT.NNZ)) GOTO 800
 600  CONTINUE
      IF (LP.GE.KK) GOTO 700
      LP = IDX(LP)
      GOTO 600
 700  CONTINUE
C        ... Swap of vectors IA2, IA1, AR ...
      SWAPIA3 = IA3(KK)
      SWAPIA2 = IA2(KK)
      SWAPIA1 = IA1(KK)
      SWAPAR = AR(KK)
      IA3(KK) = IA3(LP)
      IA2(KK) = IA2(LP)
      IA1(KK) = IA1(LP)
      AR(KK) = AR(LP)
      IA3(LP) = SWAPIA3
      IA2(LP) = SWAPIA2
      IA1(LP) = SWAPIA1
      AR(LP) = SWAPAR
      LSWAP = IDX(LP)
      IDX(LP) = IDX(KK)
      IDX(KK) = LP
      LP    = LSWAP
      KK = KK+1
      GOTO 500
 800  CONTINUE
      return
      end

      subroutine zreordvn(nnz,ar,ia1,ia2,idx)
      integer nnz
      integer ia1(*),ia2(*),idx(0:*)
      complex*16 ar(*)
      integer lp, kk, swapia1, swapia2, lswap
      complex*16 swapar

      LP = IDX(0)
      KK = 1
 500  CONTINUE
      IF ((LP.EQ.0).OR.(KK.GT.NNZ)) GOTO 800
 600  CONTINUE
      IF (LP.GE.KK) GOTO 700
      LP = IDX(LP)
      GOTO 600
 700  CONTINUE
C        ... Swap of vectors IA2, IA1, AR ...
      SWAPIA2 = IA2(KK)
      SWAPIA1 = IA1(KK)
      SWAPAR  = AR(KK)
      IA2(KK) = IA2(LP)
      IA1(KK) = IA1(LP)
      AR(KK)  = AR(LP)
      IA2(LP) = SWAPIA2
      IA1(LP) = SWAPIA1
      AR(LP)  = SWAPAR
      LSWAP   = IDX(LP)
      IDX(LP) = IDX(KK)
      IDX(KK) = LP
      LP      = LSWAP
      KK      = KK+1
      GOTO 500
 800  CONTINUE
      return
      end
      subroutine zreordvn3(nnz,ar,ia1,ia2,ia3,idx)
      integer nnz
      integer ia1(*),ia2(*),ia3(*),idx(0:*)
      complex*16 ar(*)
      integer lp, kk, swapia1, swapia2, swapia3,lswap
      complex*16 swapar

      LP = IDX(0)
      KK = 1
 500  CONTINUE
      IF ((LP.EQ.0).OR.(KK.GT.NNZ)) GOTO 800
 600  CONTINUE
      IF (LP.GE.KK) GOTO 700
      LP = IDX(LP)
      GOTO 600
 700  CONTINUE
C        ... Swap of vectors IA2, IA1, AR ...
      SWAPIA3 = IA3(KK)
      SWAPIA2 = IA2(KK)
      SWAPIA1 = IA1(KK)
      SWAPAR = AR(KK)
      IA3(KK) = IA3(LP)
      IA2(KK) = IA2(LP)
      IA1(KK) = IA1(LP)
      AR(KK) = AR(LP)
      IA3(LP) = SWAPIA3
      IA2(LP) = SWAPIA2
      IA1(LP) = SWAPIA1
      AR(LP) = SWAPAR
      LSWAP = IDX(LP)
      IDX(LP) = IDX(KK)
      IDX(KK) = LP
      LP    = LSWAP
      KK = KK+1
      GOTO 500
 800  CONTINUE
      return
      end

