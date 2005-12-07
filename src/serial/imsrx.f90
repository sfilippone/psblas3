! File:  imsrx.f90 
 ! Subroutine: 
 ! Parameters:subroutine imsrx(n,x,indx)
  integer :: n
  integer :: x(n)
  integer :: indx(n)
  
  integer, allocatable :: iaux(:)
  
  integer :: iswap, iret, info, lp, k
  integer :: lswap, ixswap

  if (n<0) then 
    write(0,*) 'Error: IMSRX: N<0'
    return
  endif
  
  if (n==0) return
  if (n==1) then 
    indx(1) = 1
    return
  endif

  allocate(iaux(0:n+1),stat=info)
  if (info/=0) then 
    write(0,*) 'IMSRX: memory allocation failed',info
    return
  endif

  do k=1,n
    indx(k) = k
  enddo

  call mrgsrt(n,x(1),iaux(1),iret)
  
  if (iret /= 1) then 
    lp = iaux(0)
    k  = 1
    do 
      if ((lp==0).or.(k>n)) exit
      do 
        if (lp >= k) exit
        lp = iaux(lp)
      end do
      iswap    = x(lp)
      x(lp)    = x(k)
      x(k)     = iswap
      ixswap   = indx(lp)
      indx(lp) = indx(k)
      indx(k)  = ixswap
      lswap    = iaux(lp)
      iaux(lp) = iaux(k)
      iaux(k)  = lswap
      lp = lswap 
      k  = k + 1
    enddo
  end if

  deallocate(iaux,stat=info)
  if (info/=0) then 
    write(0,*) 'IMSRX: memory deallocation failed',info
  endif
  return
end subroutine imsrx
