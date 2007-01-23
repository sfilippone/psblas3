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
      SUBROUTINE DCRJD(TRANS,M,N,UNITD,D,DESCRA,AR,IA1,IA2,INFO,
     *  IP1,DESCRN,ARN,IAN1,IAN2,INFON,IP2,LARN,LIAN1,
     *  LIAN2,AUX,LAUX,SIZE_REQ,IERROR)
C
C     Purpose
C     =======
C
C     DCRJD converts a CSR matrix into a Jagged Diagonal.
C
C  
C     Notes
C     =====
C   
C     Parameters
C     ==========
C   
C     TRANS   Whether the transpose should be converted.
C   
C     M,N     Size of input matrix A                
C   
C     UNITD   Scaling by diagonal D: 'U'nit, 'L'eft, 'R'ight 
C     D(*)    
C   
C     DESCRA  Input matrix A.  
C     AR,IA1, 
C     IA2,INFO
C   
C     DESCRN  Output matrix in JAD format
C     ARN,IAN1
C     IAN2,INFON, IP1, IP2
C   
      use psb_const_mod
      use psb_spmat_type
      IMPLICIT NONE

C
C     .. Scalar Arguments ..
      INTEGER            LARN, LAUX, LAUX2, LIAN1, LIAN2, M, N,
     *     SIZE_REQ, IERROR
      CHARACTER          TRANS,UNITD
C     .. Array Arguments ..
      DOUBLE PRECISION   AR(*), ARN(*), D(*), AUX(LAUX)
      INTEGER            IA1(*), IA2(*), INFO(*), IAN1(*), IAN2(*),
     *  INFON(*), IP1(*), IP2(*)
      CHARACTER          DESCRA*11, DESCRN*11
C     .. Local Scalars ..
      INTEGER            IOFF, ISTROW, NJA, NZ, PIA,
     +  PJA, PNG, K, MAX_NG, NG, LJA, ERR_ACT
      LOGICAL            SCALE
      logical  debug
      parameter (debug=.false.)
      CHARACTER          UPLO
      INTEGER MAX_NNZERO
c     .. Local Arrays ..
      CHARACTER*20       NAME
      INTEGER            INT_VAL(5), IERRV(5)

C     .. External Subroutines ..
      EXTERNAL           DVTFG
      EXTERNAL           MAX_NNZERO
C     .. Executable Statements ..
C
      NAME = 'DCRJD\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

      IF (LAUX.LT.4) THEN
         IERROR = 60
         INT_VAL(1) = 22
         INT_VAL(2) = 4
         INT_VAL(3) = LAUX
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      ENDIF

      IF (TRANS.EQ.'N') THEN
C
        NJA    = 3*M
        SCALE  = (UNITD.EQ.'L') ! meaningless
        IOFF   = 5
C
C        SET THE VALUES OF POINTERS TO VECTOR IAN2 AND AUX
C
        PNG    = IOFF
        PIA    = PNG + 1
        PJA    = PIA + 3*(M+2)

        IF (DESCRA(1:1).EQ.'G') THEN

C
C        CHECK ON DIMENSION OF IAN2 AND AUX
C
          MAX_NG = M/PSB_MINJDROWS_+1

          IF ((PIA+3*(MAX_NG+1).GT.LIAN2).OR.(M+1 .GT. LAUX)) THEN
C              ... If I haven't sufficent memory to compute NG in IAN2 ...
            IF (M+1+3*(MAX_NG+1)/PSB_DBLEINT_+1.GT.LAUX) THEN
C              ... If I haven't sufficent memory to compute NG in AUX ...
               IERROR = 60
               INT_VAL(1) = 22
               INT_VAL(2) = M+1+3*(MAX_NG+1)/PSB_DBLEINT_+1
               INT_VAL(3) = LAUX
               CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
               GOTO 9999
            ELSE
C                 ... I have sufficent memory to compute NG in AUX ...
              CALL DGBLOCK(M,IA2,IP1,AUX(M+2),NG, AUX, LAUX*2)
              CALL CHECK_DIM(M,N,AUX(M+2),NG,IA2,
     +          NZ,LARN,LIAN1,LIAN2,IERRV)
              IF (IERRV(1).NE.0) THEN 
                SIZE_REQ = MAX(IERRV(2),IERRV(3),IERRV(4))

                GOTO 9998
              ENDIF
            ENDIF
          END IF
          
          NZ     = IA2(M+1) - 1
C
C           ... Initialize Permutation Matrix ...
C
          DO 10 K = 1, M
            IP1(K) = K
 10       CONTINUE

          IP2(1) = 0

          CALL DGBLOCK(M,IA2,IP1,IAN2(PIA),IAN2(PNG), AUX, LAUX*2)
          
          PJA = PIA + 3*(IAN2(PNG)+1)
C
C           CHECK FOR ARRAY DIMENSIONS
C
          CALL CHECK_DIM(M,N,IAN2(PIA),IAN2(PNG),IA2,
     +      NZ,LARN,LIAN1,LIAN2,IERRV)
          IF (IERRV(1) .NE.0) THEN 
            SIZE_REQ = MAX(IERRV(2),IERRV(3),IERRV(4))
            GOTO 9998
          ENDIF

          LJA = LIAN2-PJA
          CALL DGINDEX(M,IAN2(PNG),AR,IA1,IA2,ARN,IAN1,IAN2(PIA), 
     +         IAN2(PJA), INFON, LARN,LIAN1,
     +         LJA,IP1, AUX, LAUX*2, SIZE_REQ,IERROR)

          IF (IERROR.NE.0) THEN
             CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
             GOTO 9999
          ENDIF

          DESCRN(1:1) = 'G'
          DESCRN(2:2) = 'U'
          DESCRN(3:3) = 'N'

        ELSE IF (DESCRA(1:1).EQ.'S' .AND. DESCRA(2:2).EQ.'U') THEN
C
          ISTROW = 1
          NZ     = 2*(IA2(M+1)-1) - M
C
C           CHECK ON DIMENSION OF IAN1 AND ARN
C
          IF (NZ .GT. LIAN1) THEN
             IERROR = 60
             INT_VAL(1) = 19
             INT_VAL(2) = NZ
             INT_VAL(3) = LAUX
             LIAN1  = NZ
             CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
             GOTO 9999
          END IF
          IF (NZ .GT. LARN) THEN
             IERROR = 60
             INT_VAL(1) = 18
             INT_VAL(2) = NZ
             INT_VAL(3) = LAUX
             LIAN1  = NZ
             CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
             GOTO 9999
          END IF

          DO 20 K = 1, M
            IP2(K) = K
 20       CONTINUE

c$$$            CALL DVSSG(M,IA1,IA2,IP2,IAN2(PNG),IP1,IP2,AUX(IWLEN),
c$$$     *                 AUX(IWORK1))
c$$$            CALL DVSMR(M,AR,IA1,IA2,IAN2(PNG),AUX(IWLEN),IP1,IP2,
c$$$     *                 IAN2(PIA),IAN2(PJA),IAN1,ARN,AUX(IWORK1),
c$$$     *                 AUX(IWORK2),NJA,IER,SCALE)
C
        ELSE IF (DESCRA(1:1).EQ.'T') THEN
C
C  Only unit diagonal so far for triangular matrices. 
C


          IF (DESCRA(3:3).NE.'U') THEN 
            IERROR=3022
            CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
            GOTO 9999
          ENDIF          
          
          UPLO = DESCRA(2:2)
          NZ     = IA2(M+1) - 1
C
C           ...Compute levels...
C           Each level correspond to a block
C           IAN1 is used as a work area

          CALL DVTFG(UPLO,M,IA1,IA2,IAN2(PNG),IP2,IP1,IAN1,
     +      AUX,AUX(M+1),AUX(2*(M+1)))

C           Generate IA(1,*)
          DO K = 1, IAN2(PNG)+1
            IAN2(PIA+3*(K-1)) = IAN1(K)
          ENDDO

          CALL GEN_BLOCK(M,IAN2(PNG),IAN2(PIA),AUX)

          PJA = PIA + 3*(IAN2(PNG)+1)

C
C           CHECK FOR ARRAY DIMENSIONS
C

          CALL CHECK_DIM(M,N,IAN2(PIA),IAN2(PNG),IA2,
     +      NZ,LARN,LIAN1,LIAN2,IERRV)
          
          IF (IERRV(1).NE.0) THEN 
            size_req = max(ierrv(2),ierrv(3),ierrv(4))
c$$$                write(0,*) "error 2",ierrv(1)
            GOTO 9998
          endif
          LJA = LIAN2-PJA

          CALL DGIND_TRI(M,IAN2(PNG),AR,IA1,IA2,ARN,IAN1,IAN2(PIA), 
     +      IAN2(PJA),LARN,LIAN1,LJA,IP1,AUX, LAUX*2, IERROR)

          IF (IERROR.NE.0) THEN
             IERROR=4011
           CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
           GOTO 9999
          ENDIF
          
          DESCRN(1:1) = 'T'
          DESCRN(2:2) = DESCRA(2:2)
          DESCRN(3:3) = DESCRA(3:3)

        END IF
C
C        SET THE OUTPUT PARAMETER
C
        IAN2(1) = PNG
        IAN2(2) = PIA
        IAN2(3) = PJA
        LARN    = NZ
        LIAN1   = NZ
        LIAN2   = 3*M + 10
        LAUX2   = 4*M + 2
C
      ELSE IF (TRANS.NE.'N') THEN
C
C           TO BE DONE
C
         IERROR = 3021
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999

      END IF
      
 9998 CONTINUE 
      CALL FCPSB_ERRACTIONRESTORE(ERR_ACT)
      RETURN

 9999 CONTINUE
      CALL FCPSB_ERRACTIONRESTORE(ERR_ACT)

      IF ( ERR_ACT .NE. 0 ) THEN 
         CALL FCPSB_SERROR()
         RETURN
      ENDIF

      RETURN
      END
