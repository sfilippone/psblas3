C
C             Parallel Sparse BLAS  version 3.0
C   (C) Copyright 2006, 2007, 2008, 2009, 2010
C                      Salvatore Filippone    University of Rome Tor Vergata
C                      Alfredo Buttari        CNRS-IRIT, Toulouse
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
      subroutine  caxpby(m, n, alpha, X, lldx, beta, Y, lldy, info)
      use psb_const_mod
      complex(psb_spk_), parameter :: one=(1.0,0.0)
      complex(psb_spk_), parameter :: zero=(0.0,0.0)
      integer n, m, lldx, lldy, info
      complex(psb_spk_) X(lldx,*), Y(lldy,*)
      complex(psb_spk_) alpha, beta
      integer i, j
      integer int_err(5)
      character  name*20
      name='caxpby'


C
C     Error handling
C
      info = psb_success_
      if (m.lt.0) then 
        info=psb_err_iarg_neg_
        int_err(1)=1
        int_err(2)=m
        call fcpsb_errpush(info,name,int_err)
        goto 9999
      else if (n.lt.0) then 
        info=psb_err_iarg_neg_
        int_err(1)=1
        int_err(2)=n
        call fcpsb_errpush(info,name,int_err)
        goto 9999
      else if (lldx.lt.max(1,m)) then 
        info=psb_err_iarg_not_gtia_ii_
        int_err(1)=5
        int_err(2)=1
        int_err(3)=lldx
        int_err(4)=m
        call fcpsb_errpush(info,name,int_err)
        goto 9999
      else if (lldy.lt.max(1,m)) then 
        info=psb_err_iarg_not_gtia_ii_
        int_err(1)=8
        int_err(2)=1
        int_err(3)=lldy
        int_err(4)=m
        call fcpsb_errpush(info,name,int_err)
        goto 9999
      endif

      if (alpha.eq.zero) then 
        if (beta.eq.zero) then 
          do j=1, n 
            do i=1,m 
              y(i,j) = zero
            enddo 
          enddo
        else if (beta.eq.one) then
c$$$
c$$$     Do nothing! 
c$$$            

        else if (beta.eq.-one) then 
          do j=1,n 
            do i=1,m 
              y(i,j) = - y(i,j)
            enddo
          enddo
        else  
          do j=1,n 
            do i=1,m 
              y(i,j) =  beta*y(i,j)
            enddo
          enddo            
        endif

      else if (alpha.eq.one) then

        if (beta.eq.zero) then 
          do j=1,n 
            do i=1,m 
              y(i,j) = x(i,j)
            enddo
          enddo
        else if (beta.eq.one) then
          do j=1,n 
            do i=1,m 
              y(i,j) = x(i,j) + y(i,j)
            enddo
          enddo

        else if (beta.eq.-one) then 
          do j=1,n 
            do i=1,m 
              y(i,j) = x(i,j) - y(i,j)
            enddo
          enddo
        else  
          do j=1,n 
            do i=1,m 
              y(i,j) = x(i,j) + beta*y(i,j)
            enddo
          enddo            
        endif

      else if (alpha.eq.-one) then 

        if (beta.eq.zero) then 
          do j=1,n 
            do i=1,m 
              y(i,j) = -x(i,j)
            enddo
          enddo
        else if (beta.eq.one) then
          do j=1,n 
            do i=1,m 
              y(i,j) = -x(i,j) + y(i,j)
            enddo
          enddo

        else if (beta.eq.-one) then 
          do j=1,n 
            do i=1,m 
              y(i,j) = -x(i,j) - y(i,j)
            enddo
          enddo
        else  
          do j=1,n 
            do i=1,m 
              y(i,j) = -x(i,j) + beta*y(i,j)
            enddo
          enddo            
        endif

      else  

        if (beta.eq.zero) then 
          do j=1,n 
            do i=1,m 
              y(i,j) = alpha*x(i,j)
            enddo
          enddo
        else if (beta.eq.one) then
          do j=1,n 
            do i=1,m 
              y(i,j) = alpha*x(i,j) + y(i,j)
            enddo
          enddo

        else if (beta.eq.-one) then 
          do j=1,n 
            do i=1,m 
              y(i,j) = alpha*x(i,j) - y(i,j)
            enddo
          enddo
        else  
          do j=1,n 
            do i=1,m 
              y(i,j) = alpha*x(i,j) + beta*y(i,j)
            enddo
          enddo            
        endif

      endif

      return

 9999 continue
      call fcpsb_serror()
      return

      end
