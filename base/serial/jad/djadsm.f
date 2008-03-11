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
      SUBROUTINE DJADSM(TRANST,M,N,VDIAG,TDIAG,PERMQ,ALPHA,DESCRA,
     +   AR,JA,IA,PERMP,B,LDB,BETA,C,LDC,WORK)
      use psb_const_mod
      use psb_error_mod
      use psb_string_mod
      implicit none
C
C
C     .. Scalar Arguments ..
      INTEGER           LDB, LDC, M, N
      CHARACTER         TDIAG, TRANST
      real(psb_dpk_)  ALPHA, BETA
C     .. Array Arguments ..
      real(psb_dpk_)  AR(*), B(LDB,*), C(LDC,*), VDIAG(*), WORK(*)
      INTEGER           IA(*), JA(*), PERMP(*), PERMQ(*)
      CHARACTER         DESCRA*11
C     .. Local Scalars ..
      INTEGER           PIA, PJA, PNG
      INTEGER           I, K, ERR_ACT
      CHARACTER         UPLO,UNITD
      integer        :: debug_level, debug_unit
      CHARACTER*20      NAME
      INTEGER           INT_VAL(5),ierror
C     .. Executable Statements ..
C
      NAME = 'DJADSM'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)
      debug_unit  = psb_get_debug_unit()
      debug_level = psb_get_debug_level()

      IF((ALPHA.NE.1.D0) .OR. (BETA.NE.0.D0))then
         IERROR=5
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      ENDIF
      UPLO = '?'
      IF (psb_toupper(DESCRA(1:1)).EQ.'T' .AND.
     +  psb_toupper(DESCRA(2:2)).EQ.'U') UPLO = 'U'
      IF (psb_toupper(DESCRA(1:1)).EQ.'T' .AND.
     +  psb_toupper(DESCRA(2:2)).EQ.'L') UPLO = 'L'
C
      IF (UPLO.EQ.'?') THEN
         IERROR=5
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      END IF

      IF (psb_toupper(DESCRA(3:3)).NE.'U') THEN
         IERROR=5
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      END IF
      UNITD=DESCRA(3:3)
C
C        B = INV(A)*B  OR B=INV(A')*B
C
      if (debug_level >= psb_debug_serial_comp_)
     +  write(debug_unit,*) trim(name),': entry',m,n

      IF (psb_toupper(TDIAG).EQ.'R') THEN
        if (debug_level >= psb_debug_serial_comp_)
     +    write(debug_unit,*) trim(name),': Right Scale'
        DO  I = 1, N
          DO  K = 1, M
            B(K,I) = B(K,I)*VDIAG(K)
          ENDDO
        ENDDO
      END IF
      
      PNG = IA(1)
      PIA = IA(2)
      PJA = IA(3)

      DO I = 1, N
         CALL DJADSV(UNITD,M,IA(PNG),
     +      AR,JA,IA(PIA),IA(PJA),B(1,I),C(1,I),IERROR)
      ENDDO
      IF(IERROR.NE.0) THEN
         INT_VAL(1)=IERROR
         CALL FCPSB_ERRPUSH(4012,NAME,INT_VAL)
         GOTO 9999
      END IF


      IF (psb_toupper(TDIAG).EQ.'L') THEN
        if (debug_level >= psb_debug_serial_comp_)
     +    write(debug_unit,*) trim(name),': Left Scale'
         DO I = 1, N
            DO K = 1, M
               C(K,I) = C(K,I)*VDIAG(K)
            ENDDO
         ENDDO
      END IF

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
