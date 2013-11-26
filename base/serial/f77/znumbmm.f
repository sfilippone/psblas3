c == =====================================================================
c Sparse Matrix Multiplication Package
c
c Randolph E. Bank and Craig C. Douglas
c
c na.bank@na-net.ornl.gov and na.cdouglas@na-net.ornl.gov
c
c Compile this with the following command (or a similar one):
c
c     f77 -c -O smmp.f
c
c == =====================================================================
      subroutine znumbmm(n, m, l,
     *  ia, ja, diaga, a,
     *  ib, jb, diagb, b,
     *  ic, jc, diagc, c,
     *  temp)
c
      use psb_const_mod
      integer(psb_ipk_) :: ia(*), ja(*), diaga,
     *  ib(*), jb(*), diagb,
     *  ic(*), jc(*), diagc 
c
      complex(psb_dpk_) :: a(*), b(*), c(*), temp(*),ajj
c
c       numeric matrix multiply c=a*b
c
      maxlmn = max(l,m,n)
      do 10 i = 1,maxlmn
        temp(i) = 0.
 10   continue 
      minlm = min(l,m)
      minln = min(l,n)
      minmn = min(m,n)
c
c   c = a*b
c
      do 50 i = 1,n
        do 30 jj = ia(i),ia(i+1)
c    a = d + ...
          if (jj.eq.ia(i+1)) then
            if (diaga.eq.0 .or. i.gt.minmn) goto 30
            j = i
            ajj = a(i)
          else
            j=ja(jj)
            ajj = a(jj)
          endif
c    b = d + ...
          if (diagb.eq.1 .and. j.le.minlm) 
     *      temp(j) = temp(j) + ajj * b(j)
          if ((j<1).or.(j>m)) then 
            write(psb_err_unit,*)
     +        ' NUMBMM: Problem with A ',i,jj,j,m
          endif
          
          do 20 k = ib(j),ib(j+1)-1
            if((jb(k)<1).or. (jb(k) > maxlmn))  then 
              write(psb_err_unit,*)
     +          ' NUMBMM: jb  problem',j,k,jb(k),maxlmn
            else
              temp(jb(k)) = temp(jb(k)) + ajj * b(k)
            endif
 20       continue 
 30     continue
c    c = d + ...
        if (diagc.eq.1 .and. i.le.minln) then
          c(i) = temp(i)
          temp(i) = 0.
        endif
c$$$        if (mod(i,100) == 1)
c$$$     +    write(psb_err_unit,*)
c$$$     ' NUMBMM: Fixing row ',i,ic(i),ic(i+1)-1
        do 40 j = ic(i),ic(i+1)-1
          if((jc(j)<1).or. (jc(j) > maxlmn))  then 
            write(psb_err_unit,*)
     +        ' NUMBMM: output problem',i,j,jc(j),maxlmn
          else
            c(j) = temp(jc(j))
            temp(jc(j)) = 0.
          endif
 40     continue 
 50   continue
      return
      end
