      subroutine  daxpby(m, n, alpha, X, lldx, beta, Y, lldy, info)
      double precision one, zero
      parameter  (one=1.d0,zero=0.d0)
      integer n, m, lldx, lldy, info
      double precision X(lldx,*), Y(lldy,*)
      double precision alpha, beta
      integer i, j
      integer int_err(5)
      double precision real_err(5)
      character  name*20
      name='daxpby'


C
C     Error handling
C
      info = 0
      if (m.lt.0) then 
         info=10
         int_err(1)=1
         int_err(2)=m
         call fcpsb_errpush(info,name,int_err)
         goto 9999
      else if (n.lt.0) then 
         info=10
         int_err(1)=1
         int_err(2)=n
         call fcpsb_errpush(info,name,int_err)
         goto 9999
      else if (lldx.lt.max(1,m)) then 
         info=50
         int_err(1)=5
         int_err(2)=1
         int_err(3)=lldx
         int_err(4)=m
         call fcpsb_errpush(info,name,int_err)
         goto 9999
      else if (lldy.lt.max(1,m)) then 
         info=50
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
