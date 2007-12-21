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
      use psb_error_mod
      use psb_spmat_type
      use psb_string_mod
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
      integer            nnz, k, row, i, j, nzl, iret
      integer            ipx, ip1, ip2, check_flag, err_act
      integer            elem, elem_csr,regen_flag
      logical            scale
      integer            max_nnzero
      integer, allocatable :: itmp(:)
c     .. local arrays ..
      character*20       name
      integer            int_val(5)
      integer         :: debug_level, debug_unit

C
C     ...Common variables...

C     .. External Subroutines ..
      EXTERNAL           MAX_NNZERO
C     .. Executable Statements ..
C
      NAME = 'ZCOCR'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)
      debug_unit  = psb_get_debug_unit()
      debug_level = psb_get_debug_level()

      call psb_getifield(check_flag,psb_dupl_,infon,psb_ifasize_,ierror)
      call psb_getifield(regen_flag,psb_upd_,infon,psb_ifasize_,ierror)

      IF (toupper(TRANS).EQ.'N') THEN

        SCALE  = (toupper(UNITD).EQ.'L') ! meaningless
        P1(1) = 0
        P2(1) = 0
        nnz = info(1)
        if (debug_level >= psb_debug_serial_) then 
          write(debug_unit,*) trim(name),': On entry  NNZ LAUX ',
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
          INT_VAL(2) = NNZ
          INT_VAL(3) = LARN
        ELSE IF (LIAN1.LT.NNZ) THEN
          IERROR = 60
          INT_VAL(1) = 19
          INT_VAL(2) = NNZ
          INT_VAL(3) = LIAN1
        ELSE IF (LIAN2.LT.M+1) THEN
          IERROR = 60
          INT_VAL(1) = 20
          INT_VAL(2) = M+1
          INT_VAL(3) = LIAN2
        ENDIF
        
C
C     Error handling
C
        IF(IERROR.NE.0) THEN
          CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
          GOTO 9999
        END IF

        allocate(itmp(nnz),stat=iret)
        if (iret /= 0) then 
          call fcpsb_errpush(4010,name,int_val)
          goto 9999
        end if
        
        do k=1, nnz
          arn(k)  = ar(k)
          ian1(k) = ja(k)
          itmp(k) = ia(k)
        enddo
        ! Mark as unavailable by default.
        infon(psb_upd_pnt_)  = 0

        
        IF (toupper(descra(1:1)).EQ.'G') THEN
C
C        Sort COO data structure
C
          if (debug_level >= psb_debug_serial_)
     +      write(debug_unit,*) trim(name),': First sort',nnz

          if ((regen_flag==psb_upd_perm_).and.
     +      (lian2.ge.((m+1)+nnz+psb_ireg_flgs_+1))
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

            if (debug_level >= psb_debug_serial_)
     +        write(debug_unit,*)  trim(name),
     +        ': Build check :',ian2(ip1+psb_nnzt_) 

C       .... Order with key IA ...
            call msort_up(nnz,itmp,aux,iret)
            if (iret.eq.0)
     +         call zreordvn3(nnz,arn,itmp,ian1,aux(ipx),aux)
            if (debug_level >= psb_debug_serial_) then 
              do i=1, nnz-1
                if (itmp(i).gt.itmp(i+1)) then 
                  write(debug_unit,*)  trim(name),
     +              'Sorting error:',i,itmp(i),itmp(i+1)
                endif
              enddo
              write(debug_unit,*)  trim(name),
     +          'nnz :',m,nnz,itmp(nnz),ian1(nnz)
            endif

C       .... Order with key JA ...
            
            I    = 1
            J    = I
            do 
              if (i>nnz) exit
              do
                if (j>nnz) exit
                if (itmp(j) /= itmp(i)) exit
                j = j+1
              enddo
              nzl = j - i
              call msort_up(nzl,ian1(i),aux,iret)
              if (iret.eq.0) call zreordvn3(nzl,arn(i),itmp(i),ian1(i),
     +          aux(ipx+i-1),aux)
              i = j
            enddo

c        ... Construct CSR Representation...
            elem = 1
            elem_csr = 1
c        ... Insert first element ...
            do row = 1, itmp(1)
              ian2(row) = 1
            enddo
            if (debug_level >= psb_debug_serial_)
     +        write(debug_unit,*)  trim(name),
     +        ': Rebuild CSR',ia(1),elem_csr

            ian1(elem_csr) = ian1(elem)
            arn(elem_csr)  = arn(elem)
            ian2(ip2+aux(ipx+elem-1)-1) = elem_csr
            elem           = elem+1
            elem_csr       = elem_csr+1
c        ... insert remaining element ...
            do row = itmp(1), m
              do 
                if (elem > nnz) exit
                if (itmp(elem) /= row) exit
                if (itmp(elem).ne.itmp(elem-1)) then
c                 ... insert first element of a row ...
                  ian1(elem_csr) = ian1(elem)
                  arn(elem_csr)  = arn(elem)
                  ian2(ip2+aux(ipx+elem-1)-1) = elem_csr
                  elem_csr       = elem_csr+1
                else if (ian1(elem).ne.ian1(elem-1)) then
c                 ... insert other element of row ...
                  ian1(elem_csr) = ian1(elem)
                  arn(elem_csr)  = arn(elem)
                  ian2(ip2+aux(ipx+elem-1)-1) = elem_csr
                  elem_csr = elem_csr+1
                else
                  if (check_flag.eq.psb_dupl_err_) then
c                    ... error, there are duplicated elements ...
                    ierror = 130
                    call fcpsb_errpush(ierror,name,int_val)
                    goto 9999
                  else if (check_flag.eq.psb_dupl_ovwrt_) then
c                    ... insert only the last duplicated element ...
                    arn(elem_csr-1) = arn(elem)
                    ian2(ip2+aux(ipx+elem-1)-1) = elem_csr-1
                  else if (check_flag.eq.psb_dupl_add_) then 
c                    ... sum the duplicated element ...
                    arn(elem_csr-1) = arn(elem_csr-1) + arn(elem)
                    ian2(ip2+aux(ipx+elem-1)-1) = elem_csr-1
                  end if
                endif
                elem = elem + 1
              enddo
              ian2(row+1) = elem_csr
            enddo

            
          else
C       .... Order with key IA ...

            call msort_up(nnz,itmp,aux,iret)
            if (iret.eq.0) call zreordvn(nnz,arn,itmp,ian1,aux)
C       .... Order with key JA ...
            i    = 1
            j    = i
            do 
              if (i>nnz) exit
              do
                if (j>nnz) exit
                if (itmp(j) /= itmp(i)) exit
                j = j+1
              enddo
              nzl = j - i
              call msort_up(nzl,ian1(i),aux,iret)
              if (iret.eq.0)
     +          call zreordvn(nzl,arn(i),itmp(i),ian1(i),aux)
              i = j
            enddo



C        ... Construct CSR Representation...
            elem = 1
            elem_csr = 1
C        ... Insert first element ...
            do row = 1, itmp(1)
              ian2(row) = 1
            enddo
            if (debug_level >= psb_debug_serial_)
     +        write(debug_unit,*)  trim(name),
     +        ': Rebuild CSR',ia(1),elem_csr

            ian1(elem_csr) = ian1(elem)
            arn(elem_csr)  = arn(elem)
            elem = elem+1
            elem_csr = elem_csr+1
C        ... Insert remaining element ...
            do row = itmp(1), m
              do 
                if (elem > nnz) exit
                if (itmp(elem) /= row) exit
                if (itmp(elem).ne.itmp(elem-1)) then
c                 ... insert first element of a row ...
                  ian1(elem_csr) = ian1(elem)
                  arn(elem_csr)  = arn(elem)
                  elem_csr = elem_csr+1
                else if (ian1(elem).ne.ian1(elem-1)) then
C                 ... Insert other element of row ...
                  ian1(elem_csr) = ian1(elem)
                  arn(elem_csr)  = arn(elem)
                  elem_csr = elem_csr+1
                else
                  if (check_flag.eq.psb_dupl_err_) then
c     ... error, there are duplicated elements ...
                    ierror = 130
                    call fcpsb_errpush(ierror,name,int_val)
                    goto 9999
                  else if (check_flag.eq.psb_dupl_ovwrt_) then
c                    ... insert only the last duplicated element ...
                    arn(elem_csr-1) = arn(elem)
                  else if (check_flag.eq.psb_dupl_add_) then 
c                    ... sum the duplicated element ...
                    arn(elem_csr-1) = arn(elem_csr-1) + arn(elem)
                  end if
                endif
                elem = elem + 1
              enddo
              ian2(row+1) = elem_csr
            enddo
          endif

          if (debug_level >= psb_debug_serial_)
     +      write(debug_unit,*)  trim(name),': Done Rebuild CSR',
     +      ian2(m+1),ia(elem)

        ELSE IF (toupper(DESCRA(1:1)).EQ.'S' .AND.
     +      toupper(DESCRA(2:2)).EQ.'U') THEN

          do 20 k = 1, m
            p2(k) = k
 20       continue

        else if (toupper(DESCRA(1:1)).EQ.'T' .AND.
     +      toupper(DESCRA(2:2)).EQ.'U') THEN

            call msort_up(nnz,itmp,aux,iret)
            if (iret.eq.0) call zreordvn(nnz,arn,itmp,ian1,aux)
C       .... Order with key JA ...
            i    = 1
            j    = i
            do 
              if (i>nnz) exit
              do
                if (j>nnz) exit
                if (itmp(j) /= itmp(i)) exit
                j = j+1
              enddo
              nzl = j - i
              call msort_up(nzl,ian1(i),aux,iret)
              if (iret.eq.0)
     +          call zreordvn(nzl,arn(i),itmp(i),ian1(i),aux)
              i = j
            enddo



C        ... Construct CSR Representation...
            elem = 1
            elem_csr = 1
C        ... Insert first element ...
            do row = 1, itmp(1)
              ian2(row) = 1
            enddo
            if (debug_level >= psb_debug_serial_)
     +        write(debug_unit,*)  trim(name),
     +        ': Rebuild CSR',ia(1),elem_csr

            ian1(elem_csr) = ian1(elem)
            arn(elem_csr)  = arn(elem)
            elem = elem+1
            elem_csr = elem_csr+1
C        ... Insert remaining element ...
            do row = itmp(1), m
              do 
                if (elem > nnz) exit
                if (itmp(elem) /= row) exit
                if (itmp(elem).ne.itmp(elem-1)) then
c                 ... insert first element of a row ...
                  ian1(elem_csr) = ian1(elem)
                  arn(elem_csr)  = arn(elem)
                  elem_csr = elem_csr+1
                else if (ian1(elem).ne.ian1(elem-1)) then
C                 ... Insert other element of row ...
                  ian1(elem_csr) = ian1(elem)
                  arn(elem_csr)  = arn(elem)
                  elem_csr = elem_csr+1
                else
                  if (check_flag.eq.psb_dupl_err_) then
c     ... error, there are duplicated elements ...
                    ierror = 130
                    call fcpsb_errpush(ierror,name,int_val)
                    goto 9999
                  else if (check_flag.eq.psb_dupl_ovwrt_) then
c                    ... insert only the last duplicated element ...
                    arn(elem_csr-1) = arn(elem)
                  else if (check_flag.eq.psb_dupl_add_) then 
c                    ... sum the duplicated element ...
                    arn(elem_csr-1) = arn(elem_csr-1) + arn(elem)
                  end if
                endif
                elem = elem + 1
              enddo
              ian2(row+1) = elem_csr
            enddo
          
        else if (toupper(descra(1:1)).EQ.'T' .AND.
     +        toupper(DESCRA(2:2)).EQ.'L') THEN

            call msort_up(nnz,itmp,aux,iret)
            if (iret.eq.0) call zreordvn(nnz,arn,itmp,ian1,aux)
C       .... Order with key JA ...
            i    = 1
            j    = i
            do 
              if (i>nnz) exit
              do
                if (j>nnz) exit
                if (itmp(j) /= itmp(i)) exit
                j = j+1
              enddo
              nzl = j - i
              call msort_up(nzl,ian1(i),aux,iret)
              if (iret.eq.0)
     +          call zreordvn(nzl,arn(i),itmp(i),ian1(i),aux)
              i = j
            enddo

C        ... Construct CSR Representation...
            elem = 1
            elem_csr = 1
C        ... Insert first element ...
            do row = 1, itmp(1)
              ian2(row) = 1
            enddo
            if (debug_level >= psb_debug_serial_)
     +        write(debug_unit,*)  trim(name),
     +        ': Rebuild CSR',ia(1),elem_csr

            ian1(elem_csr) = ian1(elem)
            arn(elem_csr)  = arn(elem)
            elem = elem+1
            elem_csr = elem_csr+1
C        ... Insert remaining element ...
            do row = itmp(1), m
              do 
                if (elem > nnz) exit
                if (itmp(elem) /= row) exit
                if (itmp(elem).ne.itmp(elem-1)) then
c                 ... insert first element of a row ...
                  ian1(elem_csr) = ian1(elem)
                  arn(elem_csr)  = arn(elem)
                  elem_csr = elem_csr+1
                else if (ian1(elem).ne.ian1(elem-1)) then
C                 ... Insert other element of row ...
                  ian1(elem_csr) = ian1(elem)
                  arn(elem_csr)  = arn(elem)
                  elem_csr = elem_csr+1
                else
                  if (check_flag.eq.psb_dupl_err_) then
c     ... error, there are duplicated elements ...
                    ierror = 130
                    call fcpsb_errpush(ierror,name,int_val)
                    goto 9999
                  else if (check_flag.eq.psb_dupl_ovwrt_) then
c                    ... insert only the last duplicated element ...
                    arn(elem_csr-1) = arn(elem)
                  else if (check_flag.eq.psb_dupl_add_) then 
c                    ... sum the duplicated element ...
                    arn(elem_csr-1) = arn(elem_csr-1) + arn(elem)
                  end if
                endif
                elem = elem + 1
              enddo
              ian2(row+1) = elem_csr
            enddo
          
          if (debug_level >= psb_debug_serial_)
     +        write(debug_unit,*)  trim(name),': Done Rebuild CSR',
     +        ian2(m+1),ia(elem)


        end if
c
      else if (toupper(TRANS).NE.'N') then
c
c           to do
c
        ierror = 3021
        call fcpsb_errpush(ierror,name,int_val)
        goto 9999

      end if
      infon(1)=elem_csr-1

      call fcpsb_erractionrestore(err_act)
      return

 9999 continue
      call fcpsb_erractionrestore(err_act)

      if ( err_act .ne. 0 ) then 
        call fcpsb_serror()
        return
      endif

      return
      end
