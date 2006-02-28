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
      SUBROUTINE DCSUPD(M, N, FIDA, DESCRA, A, IA1,
     +  IA2, INFOA, IA, JA, FIDH, DESCRH, H, IH1, IH2,
     +  INFOH, IH, JH, FLAG, GLOB_TO_LOC,
     +  IWORK, LIWORK, IERROR)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER           IA, JA, IH, JH, M, N,
     +  IERROR, FLAG, LIWORK 
C     .. Array Arguments ..
      INTEGER           IA1(*),IA2(*),IH1(*),IH2(*),
     +  INFOA(*),INFOH(*),IWORK(*),
     +  GLOB_TO_LOC(*)
      CHARACTER         DESCRA*11,DESCRH*11, FIDA*5, FIDH*5
      DOUBLE PRECISION  A(*),H(*)    
C     .. Local Scalars..
      INTEGER          ERR_ACT
C     .. Local Array..
      integer           int_val(5)
      double precision  real_val(5)
      character*20      name, strings(2)

C     .. Executable Statements ..                                        
C                                                                        
C                                                                        
C     Check parameters                                                   
C                                                                        
      IERROR = 0
      NAME = 'DCSUPD\0'
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

      IF (M.LT.0) THEN
        IERROR = 10
        INT_VAL(1) = 1
        INT_VAL(2) = M
      ELSE IF (N.LT.0) THEN
        IERROR = 10
        INT_VAL(1) = 2
        INT_VAL(2) = N
      ENDIF
C                                                                        
C     Error handling                                                     
C                                                                        
      IF(IERROR.NE.0) THEN
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      ENDIF
C                                                                        
C     Check for M, N                                                     
C                                                                        
      IF(M.LE.0 .OR. N.LE.0) THEN                                        
        GOTO 9999                                                       
      ENDIF      

C                                                                        
C     Switching on FIDA                                               
C                                                                        
      IF (FIDA(1:3).EQ.'CSR') THEN                              
        IF (FIDH(1:3).EQ.'CSR') THEN
C
C           Submatrix H in CSR format into A matrix in CSR format
C
          CALL DCRCRUPD(M, N, DESCRA, A, IA1,
     +      IA2, INFOA, IA, JA, DESCRH, H, IH1, IH2,
     +      INFOH, IH, JH, FLAG, GLOB_TO_LOC,
     +      IWORK, LIWORK, IERROR)         
        ELSE IF (FIDH(1:3).EQ.'COO') THEN
C
C           Submatrix H in COO format into A matrix in CSR format
C
          CALL DCOCRUPD(M, N, DESCRA, A, IA1,
     +      IA2, INFOA, IA, JA, DESCRH, H, IH1, IH2,
     +      INFOH, IH, JH, FLAG, GLOB_TO_LOC,
     +      IWORK, LIWORK, IERROR)         
        ELSE
          write(*,*) 'FIDA', FIDA(1:3), 'FIDH', FIDH(1:3)
          IERROR = 3010
          STRINGS(1) = FIDH(1:3)
          NAME = 'DCSUPD\0'
        ENDIF
      ELSE IF (FIDA(1:3).EQ.'COO') THEN
        IF (FIDH(1:3).EQ.'COO')THEN
          CALL DCOCRUPD(M, N, DESCRA, A, IA1,
     +      IA2, INFOA, IA, JA, DESCRH, H, IH1, IH2,
     +      INFOH, IH, JH, FLAG, GLOB_TO_LOC,
     +      IWORK, LIWORK, IERROR)         
        ENDIF
      ELSE IF (FIDA(1:3).EQ.'JAD') THEN
        IF (FIDH(1:3).EQ.'COO')THEN
          CALL DCOJDUPD(M, N, DESCRA, A, IA1,
     +      IA2, INFOA, IA, JA, DESCRH, H, IH1, IH2,
     +      INFOH, IH, JH, FLAG, GLOB_TO_LOC,
     +      IWORK, LIWORK, IERROR)         
        ENDIF
C     IERROR = 3010
C     STRINGS(1) = FIDH(1:3)
C     NAME = 'DCSUPD\0'
      ENDIF

      IF(IERROR.NE.0) THEN
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



