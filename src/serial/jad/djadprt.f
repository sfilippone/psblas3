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
