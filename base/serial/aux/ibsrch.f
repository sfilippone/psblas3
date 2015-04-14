C
C             Parallel Sparse BLAS  version 3.4
C   (C) Copyright 2006, 2010, 2015
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
      subroutine ibsrch(ipos,key,n,v)
      use psb_serial_mod
      implicit none
      integer(psb_ipk_) :: ipos, key, n
      integer(psb_ipk_) :: v(n)
      
      integer(psb_ipk_) :: lb, ub, m
      
      
      lb = 1 
      ub = n
      ipos = -1 
      
      do while (lb.le.ub) 
        m = (lb+ub)/2
        if (key.eq.v(m))  then
          ipos = m 
          lb   = ub + 1
        else if (key.lt.v(m))  then
          ub = m-1
        else 
          lb = m + 1
        end if
      enddo
      return
      end

