c
c       What if a wrong DESCRA is passed?
c
c
*
*
      SUBROUTINE DCSRPRT(M,N,DESCRA,AR,JA,IA,TITLE,IOUT)
C
C
C     .. Scalar Arguments ..
      INTEGER           M, N, IOUT
C     .. Array Arguments ..
      DOUBLE PRECISION  AR(*)
      INTEGER           IA(*), JA(*)
      CHARACTER         DESCRA*11, TITLE*(*)
C     .. Local Scalars ..
      INTEGER           I, J, nnzero


C     .. External Subroutines ..
C
C
      if ((descra(1:1).eq.'g').or.(descra(1:1).eq.'G')) then 
        write(iout,fmt=998)
      else if ((descra(1:1).eq.'s').or.(descra(1:1).eq.'S')) then 
        write(iout,fmt=997)
      else
        write(iout,fmt=998)
      endif
      nnzero = ia(m+1) -1 
      write(iout,fmt=992)
      write(iout,fmt=996)
      write(iout,fmt=996) title
      write(iout,fmt=995) 'Number of rows:    ',m
      write(iout,fmt=995) 'Number of columns: ',n
      write(iout,fmt=995) 'Nonzero entries:   ',nnzero 
      write(iout,fmt=996)
      write(iout,fmt=992)
      write(iout,*) m,n,nnzero
 998  format('%%MatrixMarket matrix coordinate real general')
 997  format('%%MatrixMarket matrix coordinate real symmetric')
 992  format('%======================================== ')
 996  format('%  ',a)
 995  format('%  ',a,i9,a,i9,a,i9)

      do i=1, m
        do j=ia(i),ia(i+1)-1
          write(iout,fmt=994) i,ja(j),ar(j)
 994      format(i6,1x,i6,1x,e16.8)
        enddo
      enddo

      RETURN
      END
