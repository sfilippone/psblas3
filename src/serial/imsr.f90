! File:  imsr.f90 
 ! Subroutine: 
 ! Parameters:subroutine imsr(n,x)
  integer :: n
  integer :: x(n)
  
  integer, allocatable :: iaux(:)
  
  integer :: iswap, iret, info, lp, k
  integer :: lswap

  if (n<0) then 
    write(0,*) 'Error: IMSR: N<0'
    return
  endif
  
  if (n<=1) return
  
  allocate(iaux(0:n+1),stat=info)
  if (info/=0) then 
    write(0,*) 'IMSR: memory allocation failed',info
    return
  endif
  

  call mrgsrt(n,x(1),iaux(1),iret)
  
  if (iret == 0) then 
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
      lswap    = iaux(lp)
      iaux(lp) = iaux(k)
      iaux(k)  = lswap
      lp = lswap 
      k  = k + 1
    enddo
  end if

  deallocate(iaux,stat=info)
  if (info/=0) then 
    write(0,*) 'IMSR: memory deallocation failed',info
  endif
  return
end subroutine imsr
