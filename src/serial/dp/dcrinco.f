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
