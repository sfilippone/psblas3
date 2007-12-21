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
      SUBROUTINE DJDCO(TRANS,M,N,DESCRA,AR,IA1,IA2,IPERM,INFO,
     *     IP1,DESCRN,ARN,IA1N,IA2N,INFON,IP2,LARN,LIA1N,
     *     LIA2N,AUX,LAUX,IERROR)      
      use psb_const_mod
      use psb_error_mod
      IMPLICIT NONE
C     
C     .. Scalar Arguments ..
      INTEGER            LARN, LAUX, LIA1N, LIA2N, M, N, IERROR
      CHARACTER          TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION   AR(*), ARN(*)
      INTEGER            AUX(0:LAUX-1),IPERM(*)
      INTEGER            IA1(*), IA2(*), INFO(*), IA1N(*), 
     *     IA2N(*), INFON(*), IP1(*), IP2(*)
      CHARACTER          DESCRA*11, DESCRN*11
C     .. Local Scalars .. 
      INTEGER            PIA, PJA, PNG, ERR_ACT
      integer         :: debug_level, debug_unit
c     .. Local Arrays ..
      CHARACTER*20       NAME
      INTEGER            INT_VAL(5)

      NAME = 'DJDCO'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)
      debug_unit  = psb_get_debug_unit()
      debug_level = psb_get_debug_level()
            
      PNG = IA2(1)
      PIA = IA2(2)
      PJA = IA2(3)

      if (debug_level >= psb_debug_serial_)
     +  write(debug_unit,*)  trim(name),': On entry NNZ LAUX ',
     +     info(1),laux,larn,lia1n,lia2n
      
      CALL DJDCOX(TRANS,M,N,DESCRA,AR,IA2(PIA),IA2(PJA),
     *  IA1,IA2(PNG),IPERM, INFO, IP1,DESCRN,ARN,IA1N,IA2N,INFON,
     *  IP2,LARN,LIA1N, LIA2N,AUX,LAUX,IERROR)
        IF(IERROR.NE.0) THEN
           IERROR=4011
           CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
           GOTO 9999
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
      
