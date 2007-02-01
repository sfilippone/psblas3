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
C
C   This file imported from ScaLAPACK. 
C
C
      SUBROUTINE PDTREECOMB( ICTXT, SCOPE, N, MINE, RDEST0, CDEST0,
     $                       SUBPTR )
*
*  -- ScaLAPACK tools routine (version 1.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     February 28, 1995
*
*     .. Scalar Arguments ..
      CHARACTER          SCOPE
      INTEGER            CDEST0, ICTXT, N, RDEST0
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   MINE( * )
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL           SUBPTR
*     ..
*
*  Purpose
*  =======
*
*  PDTREECOMB does a 1-tree parallel combine operation on scalars,
*  using the subroutine indicated by SUBPTR to perform the required
*  computation.
*
*  Arguments
*  =========
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle, indicating the global context of
*          the operation. The context itself is global.
*
*  SCOPE   (global input) CHARACTER
*          The scope of the operation:  'Rowwise', 'Columnwise', or
*          'All'.
*
*  N       (global input) INTEGER
*          The number of elements in MINE.  N = 1 for the norm-2
*          computation and 2 for the sum of square.
*
*  MINE    (local input/global output) DOUBLE PRECISION array of
*          dimension at least equal to N. The local data to use in the
*          combine.
*
*  RDEST0  (global input) INTEGER
*          The process row to receive the answer. If RDEST0 = -1,
*          every process in the scope gets the answer.
*
*  CDEST0  (global input) INTEGER
*          The process column to receive the answer. If CDEST0 = -1,
*          every process in the scope gets the answer.
*
*  SUBPTR  (local input) Pointer to the subroutine to call to perform
*          the required combine.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            BCAST, RSCOPE, CSCOPE
      INTEGER            CMSSG, DEST, DIST, HISDIST, I, IAM, MYCOL,
     $                   MYROW, MYDIST, MYDIST2, NP, NPCOL, NPROW,
     $                   RMSSG, TCDEST, TRDEST
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   HIS( 2 )
*     ..
*     .. External Subroutines ..
#if !defined(SERIAL_MPI)
      EXTERNAL           BLACS_GRIDINFO, DGEBR2D, DGEBS2D,
     $                   DGERV2D, DGESD2D
#endif    
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     See if everyone wants the answer (need to broadcast the answer)
*
      BCAST = ( ( RDEST0.EQ.-1 ).OR.( CDEST0.EQ.-1 ) )
      IF( BCAST ) THEN
         TRDEST = 0
         TCDEST = 0
      ELSE
         TRDEST = RDEST0
         TCDEST = CDEST0
      END IF
#if  !defined(SERIAL_MPI)

*
*     Get grid parameters.
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Figure scope-dependant variables, or report illegal scope
*
      RSCOPE = LSAME( SCOPE, 'R' )
      CSCOPE = LSAME( SCOPE, 'C' )
*
      IF( RSCOPE ) THEN
         IF( BCAST ) THEN
            TRDEST = MYROW
         ELSE IF( MYROW.NE.TRDEST ) THEN
            RETURN
         END IF
         NP = NPCOL
         MYDIST = MOD( NPCOL + MYCOL - TCDEST, NPCOL )
      ELSE IF( CSCOPE ) THEN
         IF( BCAST ) THEN
            TCDEST = MYCOL
         ELSE IF( MYCOL.NE.TCDEST ) THEN
            RETURN
         END IF
         NP = NPROW
         MYDIST = MOD( NPROW + MYROW - TRDEST, NPROW )
      ELSE IF( LSAME( SCOPE, 'A' ) ) THEN
         NP = NPROW * NPCOL
         IAM = MYROW*NPCOL + MYCOL
         DEST = TRDEST*NPCOL + TCDEST
         MYDIST = MOD( NP + IAM - DEST, NP )
      ELSE
         RETURN
      END IF
*
      IF( NP.LT.2 )
     $   RETURN
*
      MYDIST2 = MYDIST
      RMSSG = MYROW
      CMSSG = MYCOL
      I = 1
*
   10 CONTINUE
*
         IF( MOD( MYDIST, 2 ).NE.0 ) THEN
*
*           If I am process that sends information
*
            DIST = I * ( MYDIST - MOD( MYDIST, 2 ) )
*
*           Figure coordinates of dest of message
*
            IF( RSCOPE ) THEN
               CMSSG = MOD( TCDEST + DIST, NP )
            ELSE IF( CSCOPE ) THEN
               RMSSG = MOD( TRDEST + DIST, NP )
            ELSE
               CMSSG = MOD( DEST + DIST, NP )
               RMSSG = CMSSG / NPCOL
               CMSSG = MOD( CMSSG, NPCOL )
            END IF
*
            CALL DGESD2D( ICTXT, N, 1, MINE, N, RMSSG, CMSSG )
*
            GO TO 20
*
         ELSE
*
*           If I am a process receiving information, figure coordinates
*           of source of message
*
            DIST = MYDIST2 + I
            IF( RSCOPE ) THEN
               CMSSG = MOD( TCDEST + DIST, NP )
               HISDIST = MOD( NP + CMSSG - TCDEST, NP )
            ELSE IF( CSCOPE ) THEN
               RMSSG = MOD( TRDEST + DIST, NP )
               HISDIST = MOD( NP + RMSSG - TRDEST, NP )
            ELSE
               CMSSG = MOD( DEST + DIST, NP )
               RMSSG = CMSSG / NPCOL
               CMSSG = MOD( CMSSG, NPCOL )
               HISDIST = MOD( NP + RMSSG*NPCOL+CMSSG - DEST, NP )
            END IF
*
            IF( MYDIST2.LT.HISDIST ) THEN
*
*              If I have anyone sending to me
*
               CALL DGERV2D( ICTXT, N, 1, HIS, N, RMSSG, CMSSG )
               CALL SUBPTR( MINE, HIS )
*
            END IF
            MYDIST = MYDIST / 2
*
         END IF
         I = I * 2
*
      IF( I.LT.NP )
     $   GO TO 10
*
   20 CONTINUE
*
      IF( BCAST ) THEN
         IF( MYDIST2.EQ.0 ) THEN
            CALL DGEBS2D( ICTXT, SCOPE, ' ', N, 1, MINE, N )
         ELSE
            CALL DGEBR2D( ICTXT, SCOPE, ' ', N, 1, MINE, N,
     $                    TRDEST, TCDEST )
         END IF
      END IF
#endif    
*
      RETURN
*
*     End of PDTREECOMB
*
      END
*
      SUBROUTINE DCOMBAMAX( V1, V2 )
*
*  -- ScaLAPACK tools routine (version 1.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     February 28, 1995
*
*     .. Array Arguments ..
      DOUBLE PRECISION   V1( 2 ), V2( 2 )
*     ..
*
*  Purpose
*  =======
*
*  DCOMBAMAX finds the element having max. absolute value as well
*  as its corresponding globl index.
*
*  Arguments
*  =========
*
*  V1        (local input/local output) DOUBLE PRECISION array of
*            dimension 2.  The first maximum absolute value element and
*            its global index. V1(1) = AMAX, V1(2) = INDX.
*
*  V2        (local input) DOUBLE PRECISION array of dimension 2.
*            The second maximum absolute value element and its global
*            index. V2(1) = AMAX, V2(2) = INDX.
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( ABS( V1( 1 ) ).LT.ABS( V2( 1 ) ) ) THEN
         V1( 1 ) = V2( 1 )
         V1( 2 ) = V2( 2 )
      END IF
*
      RETURN
*
*     End of DCOMBAMAX
*
      END
*
      SUBROUTINE DCOMBSSQ( V1, V2 )
*
*  -- ScaLAPACK tools routine (version 1.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     February 28, 1995
*
*     .. Array Arguments ..
      DOUBLE PRECISION   V1( 2 ), V2( 2 )
*     ..
*
*  Purpose
*  =======
*
*  DCOMBSSQ does a scaled sum of squares on two scalars.
*
*  Arguments
*  =========
*
*  V1        (local input/local output) DOUBLE PRECISION array of
*            dimension 2.  The first scaled sum. V1(1) = SCALE,
*            V1(2) = SUMSQ.
*
*  V2        (local input) DOUBLE PRECISION array of dimension 2.
*            The second scaled sum. V2(1) = SCALE, V2(2) = SUMSQ.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
      IF( V1( 1 ).GE.V2( 1 ) ) THEN
         IF( V1( 1 ).NE.ZERO )
     $      V1( 2 ) = V1( 2 ) + ( V2( 1 ) / V1( 1 ) )**2 * V2( 2 )
      ELSE
         V1( 2 ) = V2( 2 ) + ( V1( 1 ) / V2( 1 ) )**2 * V1( 2 )
         V1( 1 ) = V2( 1 )
      END IF
*
      RETURN
*
*     End of DCOMBSSQ
*
      END
*
      SUBROUTINE DCOMBNRM2( X, Y )
*
*  -- ScaLAPACK tools routine (version 1.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     February 28, 1995
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
*     ..
*
*  Purpose
*  =======
*
*  DCOMBNRM2 combines local norm 2 results, taking care not to cause
*  unnecessary overflow.
*
*  Arguments
*  =========
*
*  X       (local input) DOUBLE PRECISION
*  Y       (local input) DOUBLE PRECISION
*          X and Y specify the values x and y. X and Y are supposed to
*          be >= 0.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   W, Z
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      W = MAX( X, Y )
      Z = MIN( X, Y )
*
      IF( Z.EQ.ZERO ) THEN
         X = W
      ELSE
         X = W*SQRT( ONE+( Z / W )**2 )
      END IF
*
      RETURN
*
*     End of DCOMBNRM2
*
      END
