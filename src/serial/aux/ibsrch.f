      subroutine ibsrch(ipos,key,n,v)
      integer ipos, key, n
      integer v(n)

      integer lb, ub, m
      

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

