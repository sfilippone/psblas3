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
      SUBROUTINE DCRINCO(M,N,DESCRA,A,IA1,IA2,
     +  INFOA,IA,JA,LATOT,LIA1TOT,LIA2TOT,
     +  DESCRH,H,IH1,IH2,INFOH,IH,JH,WORK,LWORK,IERRV)
      IMPLICIT NONE                                                      
      INCLUDE  'psb_const.fh'
C     .. Scalar Arguments ..                                             
      INTEGER          LWORK, M, N
      INTEGER          LATOT,LIA1TOT,LIA2TOT,IA,JA,IH,JH
C     .. Array Arguments ..                                              
      DOUBLE PRECISION A(*), H(*), WORK(LWORK)                     
      INTEGER          IA1(*), IA2(*), IH1(*), IH2(*),
     +  INFOA(*), INFOH(*), IERRV(*)
      CHARACTER        DESCRA*11, DESCRH*11
      
C     .. Local scalars ..
      INTEGER I, J, NZA, nzh
      logical debug
      parameter (debug=.false.)


      IERRV(1) = 0
      nza      = infoa(nnz_)
      nzh      = ih2(ih+m)-ih2(ih)
      
      if ((nza+nzh).le.min(latot,lia1tot,lia2tot)) then 
C
C     In this case we are (hopefully) safe
C        
        DO I = IH, IH+M-1
          if (debug) write(0,*) 'DCRINCO: loop ',i,ih2(i),ih2(i+1)
          DO J = IH2(I), IH2(I+1)-1
            IF ((IH1(J).GE.JH).AND.(IH1(J).LT.JH+N)) THEN
C              If current element belongs to submatrix to insert
              nza      = nza +1
              A(nza)   = H(J)
              IA1(nza) = I+IA-IH
              IA2(nza) = IH1(J)+JA-JH
            ELSE  
              if (debug) then 
                write(*,*) 'DCRINCO:  out of range',jh,ih1(j),jh+n
              endif
            ENDIF
          ENDDO
        ENDDO
        
      else
C
C    Slow but safe
C        
        if (debug) write(0,*) 'DCRINCO: ',m,ih,jh,infoa(nnz_)
C     Insert Element in COO Format
        DO I = IH, IH+M-1
          if (debug) write(0,*) 'DCRINCO: loop ',i,ih2(i),ih2(i+1)
          DO J = IH2(I), IH2(I+1)-1
            IF ((IH1(J).GE.JH).AND.(IH1(J).LT.JH+N)) THEN
C              If current element belongs to submatrix to insert
              nza      = nza +1
              IF ((nza.le.LATOT) .and.(nza.le.LIA1TOT)
     +          .and.(nza.le.LIA2TOT)) THEN
                A(nza)   = H(J)
                IA1(nza) = I+IA-IH
                IA2(nza) = IH1(J)+JA-JH
                if (debug) then 
                  write(*,*) 'DCRINCO: ',j,h(j),i+ia-ih,ih1(j)+ja-jh
                endif
              else
                IF (nza.GT.LATOT) THEN
                  IERRV(1) = 10
                  IERRV(2) = nza
                ELSE IF (nza.GT.LIA1TOT) THEN
                  IERRV(1) = 20
                  IERRV(2) = nza
                ELSE IF (nza.GT.LIA2TOT) THEN
                  IERRV(1) = 30
                  IERRV(2) = nza
                ENDIF
                RETURN
              endif
            ELSE  
              if (debug) then 
                write(*,*) 'DCRINCO:  out of range',jh,ih1(j),jh+n
              endif
            ENDIF
          ENDDO
        ENDDO

      endif

      infoa(nnz_) = nza
      return
      END
