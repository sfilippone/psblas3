***********************************************************************
*                                                                     *
*    FUNCTION = This subroutine returns an array of pointers, L,      *
*               to be used to sort the integer input vector K;        *
*               the routine implements a list merge-sort              *
*                                                                     *
***********************************************************************
*                                                                     *
*                CALL MRGSRT(N,K,L,IRET)                              *
*                                                                     *
*    INPUT =                                                          *
*                                                                     *
*       SYMBOLIC NAME: N                                              *
*       POSITION:      First parameter.                               *
*       ATTRIBUTES:    INTEGER                                        *
*       VALUES:        >= 0                                           *
*       DESCRIPTION:   Dimension of the array to be sorted            *
*                                                                     *
*       SYMBOLIC NAME: K                                              *
*       POSITION:      Second parameter                               *
*       ATTRIBUTES:    INTEGER  ARRAY(N)                              *
*       VALUES:        Any                                            *
*       DESCRIPTION:   Input array containing the keys, i.e., values  *
*                      to be sorted                                   *
*                                                                     *
*                                                                     *
*                                                                     *
*    OUTPUT =                                                         *
*                                                                     *
*       SYMBOLIC NAME: L                                              *
*       POSITION:      Third parameter                                *
*       ATTRIBUTES:    INTEGER   ARRAY(N+2)                           *
*       VALUES:        >= 0                                           *
*       DESCRIPTION:   On exit, this array contains pointers to       *
*                      the keys array.                                *
*                                                                     *
***********************************************************************
***********************************************************************
*                                                                     *
***********************************************************************
***********************************************************************
*                ALGORITHM DESCRIPTION                                *
*                                                                     *
* REFERENCES  = (1) D. E. Knuth                                       *
*                   The Art of Computer Programming,                  *
*                     vol.3: Sorting and Searching                    *
*                   Addison-Wesley, 1973                              *
*                                                                     *
* FUNCTION    = This subroutine is based on the well-known merge-sort *
*               algorithm; according to (1) we are sorting 'records'  *
*               R(I) with respect to keys K(I), and to this purpose   *
*               we use 'links' L(I); at the end of the subroutine,    *
*               L(0) is the index of the first record in the sorted   *
*               sequence, then for every record R(I), we have into    *
*               L(I) the index of the next one in the sequence. A     *
*               value L(I)=0 signals the end of the sequence.         *
*               The sorting is stable, i.e., if K(I)=K(J) and I<J,    *
*               then in the sorted sequence R(I) precedes R(J); many  *
*               sorting algorithms, e.g. quicksort, are not stable.   *
*               The list merge-sort is one of the fastest stable      *
*               sortings available; it is guaranteed to run in        *
*               O(N log N) time on both the average and worst cases.  *
*                                                                     *
*                                                                     *
***********************************************************************
***********************************************************************
*                ALGORITHM EXAMPLE(S)                                 *
*                                                                     *
* EXAMPLE: Construct a sorted array of records RS from a vector R     *
*          according to the keys stored in K                          *
*                                                                     *
*          CALL MRGSRT(N,K,L,*100)                                    *
*          I = L(0)                                                   *
*          DO 100 J = 1, N                                            *
*             RS(J) = R(I)                                            *
*             I = L(I)                                                *
* 100      CONTINUE     ! RETURN POINT IF ARRAY ALREADY SORTED        *
*                                                                     *
*                                                                     *
* EXAMPLE: Sort in place array R                                      *
*                                                                     *
*          CALL MRGSRT(N,K,L,*400)                                    *
*          LP = L(0)                                                  *
*          KK = 1                                                     *
* 100      CONTINUE                                                   *
*          IF ((LP.EQ.0).OR.(KK.GT.N)) GOTO 400                       *
* 200        CONTINUE                                                 *
*            IF (LP.GE.KK) GOTO 300                                   *
*            LP = L(LP)                                               *
*            GOTO 200                                                 *
* 300        CONTINUE                                                 *
*            SWAP = R(KK)                                             *
*            R(KK) = R(LP)                                            *
*            R(LP) = SWAP                                             *
*            LSWAP = L(LP)                                            *
*            L(LP) = L(KK)                                            *
*            L(KK) = LP                                               *
*            LP    = LSWAP                                            *
*            KK = KK+1                                                *
*          GOTO 100                                                   *
* 400      CONTINUE                                                   *
*                                                                     *
*                                                                     *
***********************************************************************
      SUBROUTINE MRGSRT(N,K,L,IRET)

C     .. Scalar Arguments ..
      INTEGER N, IRET
C     ..
C     .. Array Arguments ..
      INTEGER K(N),L(0:N+1)
C     ..
C     .. Local Scalars ..
      INTEGER P,Q,S,T
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS,ISIGN
C     ..
      IRET = 0
C  First step: we are preparing ordered sublists, exploiting
C  what order was already in the input data; negative links
C  mark the end of the sublists
      L(0) = 1
      T = N + 1
      DO 100 P = 1,N - 1
        IF (K(P).LE.K(P+1)) THEN
          L(P) = P + 1
        ELSE
          L(T) = - (P+1)
          T = P
        END IF
 100  CONTINUE
      L(T) = 0
      L(N) = 0
C See if the input was already sorted
      IF (L(N+1).EQ.0) THEN
        IRET = 1
        RETURN 
      ELSE
        L(N+1) = IABS(L(N+1))
      END IF
 200  CONTINUE
C Otherwise, begin a pass through the list.
C Throughout all the subroutine we have:
C  P, Q: pointing to the sublists being merged
C  S: pointing to the most recently processed record
C  T: pointing to the end of previously completed sublist
      S = 0
      T = N + 1
      P = L(S)
      Q = L(T)
      IF (Q.EQ.0) RETURN
 300  CONTINUE
      IF (K(P).GT.K(Q)) GO TO 600
 400  CONTINUE
      L(S) = ISIGN(P,L(S))
      S = P
      P = L(P)
      IF (P.GT.0)  GO TO 3100
C  Otherwise, one sublist ended, and we append to it the rest
C  of the other one.
 500  CONTINUE
      L(S) = Q
      S = T
 550  CONTINUE
      T = Q
      Q = L(Q)
      IF (Q.GT.0) GO TO 550
      GO TO 800
 600  CONTINUE
      L(S) = ISIGN(Q,L(S))
      S = Q
      Q = L(Q)
      IF (Q.GT.0)  GO TO 3200
 700  CONTINUE
      L(S) = P
      S = T
 750  CONTINUE
      T = P
      P = L(P)
      IF (P.GT.0) GO TO 750
 800  CONTINUE
      P = -P
      Q = -Q
      IF (Q.EQ.0) THEN
        L(S) = ISIGN(P,L(S))
        L(T) = 0
        GO TO 200
      ELSE
        GO TO 300
      END IF
 3100 CONTINUE
      IF (K(P).GT.K(Q)) GO TO 600
      S = P
      P = L(P)
      IF (P.GT.0) GO TO 3100
      GO TO 500
 3200 CONTINUE
      IF (K(P).LE.K(Q)) GO TO 400
      S = Q
      Q = L(Q)
      IF (Q.GT.0) GO TO 3200
      GO TO 700
      END
