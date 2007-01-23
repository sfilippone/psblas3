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
C     Covert matrix from COO format to CSR Format
C     Note: this never sets IP1 and P2!
C
      SUBROUTINE ZCOCR(TRANS,M,N,UNITD,D,DESCRA,AR,JA,IA,INFO,
     *  P1,DESCRN,ARN,IAN1,IAN2,INFON,P2,LARN,LIAN1,
     *  LIAN2,AUX,LAUX,IERROR)

      use psb_const_mod
      use psb_spmat_type
      IMPLICIT NONE

C
C     .. Scalar Arguments ..
      INTEGER            LARN, LAUX, LAUX2, LIAN1, LIAN2, M, 
     +  N, IUPDUP, IERROR
      CHARACTER          TRANS,UNITD
C     .. Array Arguments ..
      complex(kind(1.d0))   AR(*), ARN(*), D(*)
      INTEGER            AUX(0:LAUX-1)
      INTEGER            JA(*), IA(*), INFO(*), IAN1(*), IAN2(*),
     *  INFON(*), P1(*), P2(*)
      CHARACTER          DESCRA*11, DESCRN*11
C     .. Local Scalars ..
      INTEGER            NNZ, K, ROW, I, J, NZL, IRET
      integer            ipx, ip1, ip2, CHECK_FLAG, err_act
      INTEGER            ELEM, ELEM_CSR
      LOGICAL            SCALE
      INTEGER MAX_NNZERO
      logical     debug
      parameter   (debug=.false.)
c     .. Local Arrays ..
      CHARACTER*20       NAME
      INTEGER            INT_VAL(5)

C
C     ...Common variables...

C     .. External Subroutines ..
      EXTERNAL           MAX_NNZERO
C     .. Executable Statements ..
C

      NAME = 'ZCOCR\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

      call psb_getifield(check_flag,psb_dupl_,infon,psb_ifasize_,ierror)

      IF ((TRANS.EQ.'N').or.(TRANS.EQ.'n')) THEN

        SCALE  = (UNITD.EQ.'L') ! meaningless
        P1(1) = 0
        P2(1) = 0

        NNZ = INFO(1)
        if (debug) then 
          write(0,*) 'On entry to ZCOCR: NNZ LAUX ',
     +      nnz,laux,larn,lian1,lian2
        endif
        IF (LAUX.LT.NNZ+2) THEN
          IERROR = 60
          INT_VAL(1) = 22
          INT_VAL(2) = NNZ+2
          INT_VAL(3) = LAUX
        ELSE IF (LARN.LT.NNZ) THEN
          IERROR = 60
          INT_VAL(1) = 18
          INT_VAL(2) = NNZ+2
          INT_VAL(3) = LAUX
        ELSE IF (LIAN1.LT.NNZ) THEN
          IERROR = 60
          INT_VAL(1) = 19
          INT_VAL(2) = NNZ+2
          INT_VAL(3) = LAUX
        ELSE IF (LIAN2.LT.M+1) THEN
          IERROR = 60
          INT_VAL(1) = 20
          INT_VAL(2) = NNZ+2
          INT_VAL(3) = LAUX
        ENDIF
        
C
C     Error handling
C
        IF(IERROR.NE.0) THEN
          CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
          GOTO 9999
        END IF

        IF (DESCRA(1:1).EQ.'G') THEN
C
C        Sort COO data structure
C
          if (debug) write(0,*)'First sort',nnz
c$$$          if (debug) then
c$$$            do k=1,nnz
c$$$              write(*,*) k,ia(k),ja(k),ar(k)
c$$$            enddo
c$$$          endif
c$$$          write(0,*) 'ZCOCR Sizes ',lian2,((m+1)+nnz+psb_ireg_flgs_+1),
c$$$     +      m+1,nnz,psb_ireg_flgs_,
c$$$     +      laux,2*(2+nnz)
          if ((lian2.ge.((m+1)+nnz+psb_ireg_flgs_+1))
     +      .and.(laux.ge.2*(2+nnz))) then 
C
C       Prepare for smart regeneration
c             
            
            ipx = nnz+3            
            do i=1, nnz
              aux(ipx+i-1) = i
            enddo
            ip1                  = m+2
            infon(psb_upd_pnt_)  = ip1
            ip2                  = ip1+psb_ireg_flgs_
            ian2(ip1+psb_ip2_)   = ip2
            ian2(ip1+psb_iflag_) = check_flag
            ian2(ip1+psb_nnzt_)  = nnz
            ian2(ip1+psb_nnz_)   = 0
            ian2(ip1+psb_ichk_)  = nnz+check_flag

c$$$            write(0,*)'ZCOCR Check: ',ip2,ian2(ip1+psb_iflag_),
c$$$     +        ian2(ip1+psb_nnzt_), ian2(ip1+psb_nnz_),
c$$$     +        ian2(ip1+psb_ichk_)

c$$$     +        ip1,ip2,nnz,ian2(ip1+nnzt_) 

            if (debug) write(0,*) 'Build check :',ian2(ip1+psb_nnzt_) 
C       .... Order with key IA ...
            CALL MRGSRT(NNZ,IA,AUX,IRET)
            IF (IRET.EQ.0) CALL ZREORDVN3(NNZ,AR,IA,JA,AUX(IPX),AUX)
            if (debug) then 
              do i=1, nnz-1
                if (ia(i).gt.ia(i+1)) then 
                  write(0,*) 'Sorting error:',i,ia(i),ia(i+1)
                endif
              enddo
              write(0,*) 'nnz :',m,nnz,ia(nnz),ja(nnz)
            endif

C       .... Order with key IA2N ...
            
            I    = 1
            J    = I
c$$$            DO WHILE (I.LE.NNZ)
c$$$              DO WHILE ((IA(J).EQ.IA(I)).AND.
c$$$     +          (J.LE.NNZ))
            DO 
              if (I>NNZ) exit
              DO
                if (j>nnz) exit
                if (ia(j) /= ia(i)) exit
                J = J+1
              ENDDO
              NZL = J - I
              CALL MRGSRT(NZL,JA(I),AUX,IRET)
              IF (IRET.EQ.0) CALL ZREORDVN3(NZL,AR(I),IA(I),JA(I),
     +          AUX(IPX+I-1),AUX)
              I = J
            ENDDO

C        ... Construct CSR Representation...
            ELEM = 1
            ELEM_CSR = 1
C        ... Insert first element ...
            DO ROW = 1, IA(1)
              IAN2(ROW) = 1
            ENDDO
            if (debug) write(0,*)'Rebuild CSR',ia(1),elem_csr
            IAN1(ELEM_CSR) = JA(ELEM)
            ARN(ELEM_CSR)  = AR(ELEM)
            ian2(ip2+aux(ipx+elem-1)-1) = elem_csr
            ELEM           = ELEM+1
            ELEM_CSR       = ELEM_CSR+1
C        ... Insert remaining element ...
            DO ROW = IA(1), M
c$$$              if (debug) write(*,*)'CSR Loop:',row,m,elem_csr
c$$$              DO WHILE ((IA(ELEM).EQ.ROW).AND.(ELEM.LE.NNZ))
              DO 
                if (elem > nnz) exit
                if (ia(elem) /= row) exit
                IF (IA(ELEM).NE.IA(ELEM-1)) THEN
C                 ... Insert first element of a row ...
                  IAN1(ELEM_CSR) = JA(ELEM)
                  ARN(ELEM_CSR)  = AR(ELEM)
                  ian2(ip2+aux(ipx+elem-1)-1) = elem_csr
                  ELEM_CSR       = ELEM_CSR+1
                ELSE IF (JA(ELEM).NE.JA(ELEM-1)) THEN
C                 ... Insert other element of row ...
                  IAN1(ELEM_CSR) = JA(ELEM)
                  ARN(ELEM_CSR)  = AR(ELEM)
                  ian2(ip2+aux(ipx+elem-1)-1) = elem_csr
                  ELEM_CSR = ELEM_CSR+1
                ELSE
                  IF (CHECK_FLAG.EQ.psb_dupl_err_) THEN
C                    ... Error, there are duplicated elements ...
                    IERROR = 130
                    CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                    GOTO 9999
                  ELSE IF (CHECK_FLAG.EQ.psb_dupl_ovwrt_) THEN
C                    ... Insert only the last duplicated element ...
                    ARN(ELEM_CSR-1) = AR(ELEM)
                    ian2(ip2+aux(ipx+elem-1)-1) = elem_csr-1
                  ELSE IF (CHECK_FLAG.EQ.psb_dupl_add_) THEN 
C                    ... Sum the duplicated element ...
                    ARN(ELEM_CSR-1) = ARN(ELEM_CSR-1) + AR(ELEM)
                    ian2(ip2+aux(ipx+elem-1)-1) = elem_csr-1
                  END IF
                ENDIF
                ELEM = ELEM + 1
              ENDDO
              IAN2(ROW+1) = ELEM_CSR
            ENDDO
          ELSE
c$$$            write(0,*) 'Going for ELSE !!!?'
C       .... Order with key IA ...
            CALL MRGSRT(NNZ,IA,AUX,IRET)
            IF (IRET.EQ.0) CALL ZREORDVN(NNZ,AR,IA,JA,AUX)
C       .... Order with key IA2N ...
            I    = 1
            J    = I
c$$$            DO WHILE (I.LE.NNZ)
c$$$              DO WHILE ((IA(J).EQ.IA(I)).AND.
c$$$     +          (J.LE.NNZ))
            DO 
              if (I>NNZ) exit
              DO
                if (j>nnz) exit
                if (ia(j) /= ia(i)) exit
                J = J+1
              ENDDO
              NZL = J - I
              CALL MRGSRT(NZL,JA(I),AUX,IRET)
              IF (IRET.EQ.0) CALL ZREORDVN(NZL,AR(I),IA(I),JA(I),AUX)
              I = J
            ENDDO



C        ... Construct CSR Representation...
            ELEM = 1
            ELEM_CSR = 1
C        ... Insert first element ...
            DO ROW = 1, IA(1)
              IAN2(ROW) = 1
            ENDDO
            if (debug) write(0,*)'Rebuild CSR',ia(1),elem_csr
            IAN1(ELEM_CSR) = JA(ELEM)
            ARN(ELEM_CSR) = AR(ELEM)
            ELEM = ELEM+1
            ELEM_CSR = ELEM_CSR+1
C        ... Insert remaining element ...
            DO ROW = IA(1), M
c$$$              if (debug) write(*,*)'CSR Loop:',row,m,elem_csr
c$$$              DO WHILE ((IA(ELEM).EQ.ROW).AND.(ELEM.LE.NNZ))
              DO 
                if (elem > nnz) exit
                if (ia(elem) /= row) exit
                IF (IA(ELEM).NE.IA(ELEM-1)) THEN
C                 ... Insert first element of a row ...
                  IAN1(ELEM_CSR) = JA(ELEM)
                  ARN(ELEM_CSR) = AR(ELEM)
                  ELEM_CSR = ELEM_CSR+1
                ELSE IF (JA(ELEM).NE.JA(ELEM-1)) THEN
C                 ... Insert other element of row ...
                  IAN1(ELEM_CSR) = JA(ELEM)
                  ARN(ELEM_CSR) = AR(ELEM)
                  ELEM_CSR = ELEM_CSR+1
                ELSE
                  IF (CHECK_FLAG.EQ.psb_dupl_err_) THEN
C     ... Error, there are duplicated elements ...
                    IERROR = 130
                    CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                    GOTO 9999
                  ELSE IF (CHECK_FLAG.EQ.psb_dupl_ovwrt_) THEN
C                    ... Insert only the last duplicated element ...
                    ARN(ELEM_CSR-1) = AR(ELEM)
                    if (debug) write(0,*) 'Duplicated overwrite',
     +                elem_csr-1,elem
                  ELSE IF (CHECK_FLAG.EQ.psb_dupl_add_) THEN 
C                    ... Sum the duplicated element ...
                    ARN(ELEM_CSR-1) = ARN(ELEM_CSR-1) + AR(ELEM)
                    if (debug) write(0,*) 'Duplicated add',
     +                elem_csr-1,elem
                  END IF
                ENDIF
                ELEM = ELEM + 1
              ENDDO
              IAN2(ROW+1) = ELEM_CSR
            ENDDO
          ENDIF

          if (debug) write(0,*)'Done Rebuild CSR',
     +      ian2(m+1),ia(elem)
          if (debug) then 
            do i=ian2(m+1), nnz
              write(0,*) 'Overflow check :',ia(i),ja(i),ar(i)
            enddo
          endif

        ELSE IF (DESCRA(1:1).EQ.'S' .AND. DESCRA(2:2).EQ.'U') THEN

          DO 20 K = 1, M
            P2(K) = K
 20       CONTINUE

        ELSE IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'U') THEN

C       .... Order with key IA ...
          CALL MRGSRT(NNZ,IA,AUX,IRET)
          IF (IRET.EQ.0) CALL ZREORDVN(NNZ,AR,IA,JA,AUX)
C       .... Order with key IA2N ...
          I    = 1
          J    = I
c$$$            DO WHILE (I.LE.NNZ)
c$$$              DO WHILE ((IA(J).EQ.IA(I)).AND.
c$$$     +          (J.LE.NNZ))
          DO 
            if (I>NNZ) exit
            DO
              if (j>nnz) exit
              if (ia(j) /= ia(i)) exit
              J = J+1
            ENDDO
            NZL = J - I
            CALL MRGSRT(NZL,JA(I),AUX,IRET)
            IF (IRET.EQ.0) CALL ZREORDVN(NZL,AR(I),IA(I),JA(I),AUX)
            I = J
          ENDDO


C        ... Construct CSR Representation...
          ELEM = 1
          ELEM_CSR = 1
C        ... Insert first element ...
          DO ROW = 1, IA(1)
            IAN2(ROW) = 1
          ENDDO
          if (debug) write(0,*)'Rebuild CSR',ia(1),elem_csr
          IF(JA(ELEM).GT.IA(ELEM)) THEN
            IAN1(ELEM_CSR) = JA(ELEM)
            ARN(ELEM_CSR) = AR(ELEM)
            ELEM_CSR = ELEM_CSR+1
          ENDIF

          ELEM = ELEM+1

C        ... Insert remaining element ...
          DO ROW = IA(1), M
c$$$  if (debug) write(*,*)'CSR Loop:',row,m,elem_csr
c$$$              DO WHILE ((IA(ELEM).EQ.ROW).AND.(ELEM.LE.NNZ))
            DO 
              if (elem > nnz) exit
              if (ia(elem) /= row) exit
              IF (IA(ELEM).NE.IA(ELEM-1)) THEN
C     ... Insert first element of a row ...
                IF(JA(ELEM).GT.IA(ELEM)) THEN                   
                  IAN1(ELEM_CSR) = JA(ELEM)
                  ARN(ELEM_CSR) = AR(ELEM)
                  ELEM_CSR = ELEM_CSR+1
                ENDIF
              ELSE IF (JA(ELEM).NE.JA(ELEM-1)) THEN
C     ... Insert other element of row ...
                IF(JA(ELEM).GT.IA(ELEM)) THEN                   
                  IAN1(ELEM_CSR) = JA(ELEM)
                  ARN(ELEM_CSR) = AR(ELEM)
                  ELEM_CSR = ELEM_CSR+1
                ENDIF
              ELSE
                IF (CHECK_FLAG.EQ.psb_dupl_err_) THEN
C     ... Error, there are duplicated elements ...
                  IERROR = 130
                  CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                  GOTO 9999
                ELSE IF (CHECK_FLAG.EQ.psb_dupl_ovwrt_) THEN
C     ... Insert only the last duplicated element ...
                  IF(JA(ELEM).GT.IA(ELEM)) THEN                   
                    ARN(ELEM_CSR-1) = AR(ELEM)
                  ENDIF
                  if (debug) write(0,*) 'Duplicated overwrite',
     +              elem_csr-1,elem
                ELSE IF (CHECK_FLAG.EQ.psb_dupl_add_) THEN 
C     ... Sum the duplicated element ...
                  IF(JA(ELEM).GT.IA(ELEM)) THEN                   
                    ARN(ELEM_CSR-1) = ARN(ELEM_CSR-1) + AR(ELEM)
                  ENDIF
                  if (debug) write(0,*) 'Duplicated add',
     +              elem_csr-1,elem
                END IF
              ENDIF
              ELEM = ELEM + 1
            ENDDO
            IAN2(ROW+1) = ELEM_CSR
          ENDDO

          
          if (debug) write(0,*)'Done Rebuild CSR',
     +      ian2(m+1),ia(elem)
          if (debug) then 
            do i=ian2(m+1), nnz
              write(0,*) 'Overflow check :',ia(i),ja(i),ar(i)
            enddo
          endif
          
          
          
        ELSE IF (DESCRA(1:1).EQ.'T' .AND. DESCRA(2:2).EQ.'L') THEN

C       .... Order with key IA ...
          CALL MRGSRT(NNZ,IA,AUX,IRET)
          IF (IRET.EQ.0) CALL ZREORDVN(NNZ,AR,IA,JA,AUX)
C       .... Order with key IA2N ...
          I    = 1
          J    = I
c$$$            DO WHILE (I.LE.NNZ)
c$$$              DO WHILE ((IA(J).EQ.IA(I)).AND.
c$$$     +          (J.LE.NNZ))
          DO 
            if (I>NNZ) exit
            DO
              if (j>nnz) exit
              if (ia(j) /= ia(i)) exit
              J = J+1
            ENDDO
            NZL = J - I
            CALL MRGSRT(NZL,JA(I),AUX,IRET)
            IF (IRET.EQ.0) CALL ZREORDVN(NZL,AR(I),IA(I),JA(I),AUX)
            I = J
          ENDDO

C        ... Construct CSR Representation...
          ELEM = 1
          ELEM_CSR = 1
C        ... Insert first element ...
          DO ROW = 1, IA(1)
            IAN2(ROW) = 1
          ENDDO
          if (debug) write(0,*)'Rebuild CSR',ia(1),elem_csr
          IF(JA(ELEM).LT.IA(ELEM)) THEN                   
            IAN1(ELEM_CSR) = JA(ELEM)
            ARN(ELEM_CSR) = AR(ELEM)
            ELEM_CSR = ELEM_CSR+1
          ENDIF
          ELEM = ELEM+1
          
C     ... Insert remaining element ...
          DO ROW = IA(1), M
c$$$              if (debug) write(*,*)'CSR Loop:',row,m,elem_csr
c$$$              DO WHILE ((IA(ELEM).EQ.ROW).AND.(ELEM.LE.NNZ))
            DO 
              if (elem > nnz) exit
              if (ia(elem) /= row) exit
              IF (IA(ELEM).NE.IA(ELEM-1)) THEN
C     ... Insert first element of a row ...
                IF(JA(ELEM).LT.IA(ELEM)) THEN                   
                  IAN1(ELEM_CSR) = JA(ELEM)
                  ARN(ELEM_CSR) = AR(ELEM)
                  ELEM_CSR = ELEM_CSR+1
                ENDIF
              ELSE IF (JA(ELEM).NE.JA(ELEM-1)) THEN
C     ... Insert other element of row ...
                IF(JA(ELEM).LT.IA(ELEM)) THEN                               
                  IAN1(ELEM_CSR) = JA(ELEM)
                  ARN(ELEM_CSR) = AR(ELEM)
                  ELEM_CSR = ELEM_CSR+1
                ENDIF
              ELSE
                IF (CHECK_FLAG.EQ.psb_dupl_err_) THEN
C                    ... Error, there are duplicated elements ...
                  IERROR = 130
                  CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
                  GOTO 9999
                ELSE IF (CHECK_FLAG.EQ.psb_dupl_ovwrt_) THEN
C                    ... Insert only the last duplicated element ...
                  IF(JA(ELEM).LT.IA(ELEM)) THEN                   
                    ARN(ELEM_CSR-1) = AR(ELEM)
                  ENDIF
                  if (debug) write(0,*) 'Duplicated overwrite',
     +              elem_csr-1,elem
                ELSE IF (CHECK_FLAG.EQ.psb_dupl_add_) THEN 
C                    ... Sum the duplicated element ...
                  IF(JA(ELEM).LT.IA(ELEM)) THEN                   
                    ARN(ELEM_CSR-1) = ARN(ELEM_CSR-1) + AR(ELEM)
                  ENDIF
                  if (debug) write(0,*) 'Duplicated add',
     +              elem_csr-1,elem
                END IF
              ENDIF
              ELEM = ELEM + 1
            ENDDO
            IAN2(ROW+1) = ELEM_CSR
          ENDDO


          if (debug) write(0,*)'Done Rebuild CSR',
     +      ian2(m+1),ia(elem)
          if (debug) then 
            do i=ian2(m+1), nnz
              write(0,*) 'Overflow check :',ia(i),ja(i),ar(i)
            enddo
          endif


        END IF
C
      ELSE IF (TRANS.NE.'N') THEN
C
C           TO DO
C
        IERROR = 3021
        CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
        GOTO 9999

      END IF
      INFON(1)=ELEM_CSR-1

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
