program pingpong
  use psb_base_mod
  implicit none
  integer(psb_ipk_) :: iam, np, icontxt, ip, jp, idummy
  integer(psb_ipk_), parameter  :: nmax=2**16
  integer(psb_ipk_) :: i,j,k,n
  real(psb_dpk_) :: v(nmax)
  real(psb_dpk_) :: t0, t1, t2, mbs, bytes


  call psb_init(icontxt)
  call psb_info(icontxt,iam,np)            
  !   have all processes check in 
  if ((iam >= 0).and.(iam < np)) then 
    if (iam == 0)  then 
      do ip = 1, np-1
        call psb_rcv(icontxt,idummy,ip)          
      enddo
      write(*,*) 'Hello, world: all ',np, &
           & ' processes checked in!'
    else
      ip = 0
      call psb_snd(icontxt,idummy,ip)
    endif
  end if

  n    = 1
  call psb_barrier(icontxt)
  if (iam == 0) then
    do i=1, 16
      ip = 1
      t0 = psb_wtime()
      call psb_snd(icontxt,v(1:n),ip)
      call psb_rcv(icontxt,v(1:n),ip)
      t1 = psb_wtime()
      bytes = done*n*psb_sizeof_dp
      mbs   = 2.d0*(bytes/(t1-t0))*1.d-6
      write(*,*) 'pingpong: ',n,bytes,mbs
      n = n * 2
    end do
  else if (iam == 1) then
    do i=1, 16
      ip = 0
      call psb_rcv(icontxt,v(1:n),ip)
      call psb_snd(icontxt,v(1:n),ip)
      n = n * 2
    end do
  end if

  call psb_barrier(icontxt)
  call psb_exit(icontxt)
end program pingpong
