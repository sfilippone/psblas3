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
      subroutine reordvn(nnz,ar,ia1,ia2,idx)
      use psb_const_mod
      integer nnz
      integer ia1(*),ia2(*),idx(0:*)
      real(psb_dpk_) ar(*)
      integer lp, kk, swapia1, swapia2, lswap
      real(psb_dpk_) swapar

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
      subroutine ireordv2(nnz,ia1,ia2,idx)
      use psb_const_mod
      integer nnz
      integer ia1(*),ia2(*),idx(0:*)
      integer lp, kk, swapia1, swapia2, lswap

      LP = IDX(0)
      KK = 1
 500  CONTINUE
      IF ((LP.EQ.0).OR.(KK.GT.NNZ)) GOTO 800
 600  CONTINUE
      IF (LP.GE.KK) GOTO 700
      LP = IDX(LP)
      GOTO 600
 700  CONTINUE
C        ... Swap of vectors IA2, IA1 ..
      SWAPIA2 = IA2(KK)
      SWAPIA1 = IA1(KK)
      IA2(KK) = IA2(LP)
      IA1(KK) = IA1(LP)
      IA2(LP) = SWAPIA2
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
      subroutine ireordv1(nnz,ia1,idx)
      use psb_const_mod
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
      use psb_const_mod
      integer nnz
      integer ia1(*),ia2(*),ia3(*),idx(0:*)
      real(psb_dpk_) ar(*)
      integer lp, kk, swapia1, swapia2, swapia3,lswap
      real(psb_dpk_) swapar

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
      use psb_const_mod
      integer nnz
      integer ia1(*),ia2(*),idx(0:*)
      complex(psb_dpk_) ar(*)
      integer lp, kk, swapia1, swapia2, lswap
      complex(psb_dpk_) swapar

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
      use psb_const_mod
      integer nnz
      integer ia1(*),ia2(*),ia3(*),idx(0:*)
      complex(psb_dpk_) ar(*)
      integer lp, kk, swapia1, swapia2, swapia3,lswap
      complex(psb_dpk_) swapar

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

