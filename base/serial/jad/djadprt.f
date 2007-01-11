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
c
c       What if a wrong DESCRA is passed?
c       WARNING: THIS CANNOT POSSIBLY WORK CORRECTLY BECAUSE 
c       IT DOES NOT ACCOUNT FOR ROW PERMUTATION.
*
*
      SUBROUTINE DJADPRT(NROW,NCOL,NG,A,KA,JA,IA,TITLE,IOUT)
C
C
C     .. Scalar Arguments ..
      INTEGER           IOUT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*)
      INTEGER           IA(3,*), JA(*), KA(*)
      CHARACTER         DESCRA*11, TITLE*(*)
C     .. Local Scalars ..
      INTEGER           I, K


C     .. External Subroutines ..
C
C

      nnzero = ja(ia(2,ng+1)-1+1)-1

      write(iout,fmt=998)
      
      write(iout,fmt=992)
      write(iout,fmt=996)
      write(iout,fmt=996) title
      write(iout,fmt=995) 'Number of rows:    ',nrow
      write(iout,fmt=995) 'Number of columns: ',ncol
      write(iout,fmt=995) 'Nonzero entries:   ',nnzero 
      write(iout,fmt=996)
      write(iout,fmt=992)
      write(iout,*) nrow,ncol,nnzero
 998  format('%%MatrixMarket matrix coordinate real general')
 997  format('%%MatrixMarket matrix coordinate real symmetric')
 992  format('%======================================== ')
 996  format('%  ',a)
 995  format('%  ',a,i9,a,i9,a,i9)
 994  format(i6,1x,i6,1x,e16.8)

      do ipg=1, ng 
        do k = ia(2,ipg), ia(3,ipg)-1                                  
          ipx = ia(1,ipg)                                                
          do i = ja(k), ja(k+1) - 1                                   
            write(iout,994) ipx,ka(i),a(i)
            ipx = ipx + 1                                               
          enddo 
        enddo
        
        ipx = ia(1,ipg)                                                   
        do k = ia(3,ipg), ia(2,ipg+1)-1                                
          do i = ja(k), ja(k+1) - 1       
            write(iout,994) ipx,ka(i),a(i)
          enddo
          ipx = ipx + 1                                                  
        enddo
      enddo

      return
      end
