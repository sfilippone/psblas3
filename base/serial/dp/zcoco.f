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
c     covert matrix from COO format to COO format
c
      subroutine zcoco(trans,m,n,unitd,d,descra,ar,ia1,ia2,info,
     *  p1,descrn,arn,ia1n,ia2n,infon,p2,larn,lia1n,
     *  lia2n,aux,laux,ierror)

      use psb_const_mod
      use psb_spmat_type
      implicit none

c     .. scalar arguments ..
      integer            larn, laux, lia1n, lia2n, 
     +  m, n, ierror
      character          trans,unitd
c     .. array arguments ..
      complex(kind(1.d0)) ::  ar(*), arn(*), d(*)
      integer            aux(0:laux-1)
      integer            ia1(*), ia2(*), info(*), ia1n(*), ia2n(*),
     *  infon(*), p1(*), p2(*)
      character          descra*11, descrn*11
c     .. local scalars ..
      integer            ipx, ip1, ip2, check_flag, err_act
      integer            nnz, k, i, j, nzl, iret
      integer            elem_in, elem_out
      logical            scale
      integer max_nnzero
      logical     debug
      parameter   (debug=.false.)
c     .. local arrays ..
      character*20       name
      integer            int_val(5)
c
c     ...common variables...
c     this flag describe the action to do
      
c     .. external subroutines ..
      external           max_nnzero
c     .. executable statements ..
c

      name = 'zcoco'
      ierror = 0
      call fcpsb_erractionsave(err_act)

      call psb_getifield(check_flag,psb_dupl_,infon,psb_ifasize_,ierror)
      if (trans.eq.'N') then
        scale  = (unitd.eq.'L') ! meaningless
        p1(1) = 0
        p2(1) = 0

        call psb_getifield(nnz,psb_nnz_,info,psb_ifasize_,ierror) 
        if (debug) then 
          write(*,*) 'on entry to dcoco: nnz laux ',
     +      nnz,laux,larn,lia1n,lia2n
        endif
        if (laux.lt.nnz+2) then
          ierror = 60
          int_val(1) = 22
          int_val(2) = nnz+2
          int_val(3) = laux
        else if (larn.lt.nnz) then
          ierror = 60
          int_val(1) = 18
          int_val(2) = nnz+2
          int_val(3) = laux
        else if (lia1n.lt.nnz) then
          ierror = 60
          int_val(1) = 19
          int_val(2) = nnz+2
          int_val(3) = laux
        else if (lia2n.lt.m+1) then
          ierror = 60
          int_val(1) = 20
          int_val(2) = nnz+2
          int_val(3) = laux
        endif
        
c
c     error handling
c
        if(ierror.ne.0) then
          call fcpsb_errpush(ierror,name,int_val)
          goto 9999
        end if

        if (descra(1:1).eq.'G') then
c
c     sort COO data structure
c     
          if (debug) write(*,*)'first sort',nnz
          do k=1, nnz
            arn(k)  = ar(k)
            ia1n(k) = ia1(k)
            ia2n(k) = ia2(k)
          enddo
          
          if (debug) write(*,*)'second sort'            

          if ((lia2n.ge.(2*nnz+psb_ireg_flgs_+1))
     +      .and.(laux.ge.2*(2+nnz))) then 
c     
c     prepare for smart regeneration
c     
            ipx = nnz+3            
            do i=1, nnz
              aux(ipx+i-1) = i
            enddo
            ip1                  = nnz+2
            infon(psb_upd_pnt_)  = ip1
            ip2                  = ip1+psb_ireg_flgs_
            ia2n(ip1+psb_ip2_)   = ip2
            ia2n(ip1+psb_iflag_) = check_flag
            ia2n(ip1+psb_nnzt_)  = nnz
            ia2n(ip1+psb_nnz_)   = 0
            ia2n(ip1+psb_ichk_)  = nnz+check_flag
            if (debug) write(0,*) 'build check :',ia2n(ip1+psb_nnzt_) 
            
c     .... order with key ia1n ...
            call mrgsrt(nnz,ia1n,aux,iret)
            if (iret.eq.0)
     +        call zreordvn3(nnz,arn,ia1n,ia2n,aux(ipx),aux)
c     .... order with key ia2n ...
            
            i    = 1
            j    = i
            do while (i.le.nnz)
              do while ((ia1n(j).eq.ia1n(i)).and.
     +          (j.le.nnz))
                j = j+1
              enddo
              nzl = j - i
              call mrgsrt(nzl,ia2n(i),aux,iret)
              if (iret.eq.0) call zreordvn3(nzl,arn(i),ia1n(i),ia2n(i),
     +          aux(ipx+i-1),aux)
              i = j
            enddo
            
            ia2n(ip2+aux(ipx+1-1)-1) = 1

c     ... construct final COO  representation...
            elem_out = 1
c     ... insert remaining element ...
            do elem_in  = 2, nnz
              if ((ia1n(elem_in).eq.ia1n(elem_out)).and.
     +          (ia2n(elem_in).eq.ia2n(elem_out))) then 
                if (check_flag.eq.psb_dupl_err_) then
c     ... error, there are duplicated elements ...
                  ierror = 130
                  call fcpsb_errpush(ierror,name,int_val)
                  goto 9999
                else if (check_flag.eq.psb_dupl_ovwrt_) then
c     ... insert only the first duplicated element ...
                  ia2n(ip2+aux(ipx+elem_in-1)-1) = elem_out
                else if (check_flag.eq.psb_dupl_add_) then
c     ... sum the duplicated element ...
                  arn(elem_out) = arn(elem_out) + arn(elem_in)
                  ia2n(ip2+aux(ipx+elem_in-1)-1) = elem_out
                end if
              else
                elem_out = elem_out + 1
                arn(elem_out)  = arn(elem_in)
                ia2n(ip2+aux(ipx+elem_in-1)-1) = elem_out
                ia1n(elem_out) = ia1n(elem_in)
                ia2n(elem_out) = ia2n(elem_in)
              endif
            enddo
            
          else
            
c     .... order with key ia1n ...
            call mrgsrt(nnz,ia1n,aux,iret)
            if (iret.eq.0) call zreordvn(nnz,arn,ia1n,ia2n,aux)
c     .... order with key ia2n ...
            
            i    = 1
            j    = i
            do while (i.le.nnz)
              do while ((ia1n(j).eq.ia1n(i)).and.
     +          (j.le.nnz))
                j = j+1
              enddo
              nzl = j - i
              call mrgsrt(nzl,ia2n(i),aux,iret)
              if (iret.eq.0) call zreordvn(nzl,arn(i),ia1n(i),ia2n(i),
     +          aux)
              i = j
            enddo
c     ... construct final COO  representation...
            elem_out = 1
c     ... insert remaining element ...
            do elem_in  = 2, nnz
              if ((ia1n(elem_in).eq.ia1n(elem_out)).and.
     +          (ia2n(elem_in).eq.ia2n(elem_out))) then 
                if (check_flag.eq.psb_dupl_err_) then
c     ... error, there are duplicated elements ...
                  ierror = 130
                  call fcpsb_errpush(ierror,name,int_val)
                  goto 9999
                else if (check_flag.eq.psb_dupl_ovwrt_) then
c     ... insert only the first duplicated element ...
                else if (check_flag.eq.psb_dupl_add_) then
c     ... sum the duplicated element ...
                  arn(elem_out) = arn(elem_out) + arn(elem_in)
                end if
              else
                elem_out = elem_out + 1
                arn(elem_out)  = arn(elem_in)
                ia1n(elem_out) = ia1n(elem_in)
                ia2n(elem_out) = ia2n(elem_in)
              endif
            enddo
          endif
          infon(psb_nnz_)  = elem_out
          infon(psb_srtd_) = psb_isrtdcoo_
          
          if (debug) write(*,*)'done rebuild COO',infon(1)
          
        else if (descra(1:1).eq.'S' .and. descra(2:2).eq.'U') then

          do 20 k = 1, m
            p2(k) = k
 20       continue

        else if (descra(1:1).eq.'T' .and. descra(2:2).eq.'U') then
          ierror = 3021
          call fcpsb_errpush(ierror,name,int_val)
          goto 9999

          
        else if (descra(1:1).eq.'T' .and. descra(2:2).eq.'L') then
          ierror = 3021
          call fcpsb_errpush(ierror,name,int_val)
          goto 9999
          
        end if
c
      else if (trans.ne.'N') then
c
c           to do
c
        ierror = 3021
        call fcpsb_errpush(ierror,name,int_val)
        goto 9999

      end if

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
