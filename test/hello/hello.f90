program hello
  use iso_fortran_env
  implicit none
  integer me, np, icontxt, ip, jp, idummy
  integer :: snd_buf(4)[*]
  type(event_type), allocatable :: snd_copied(:)[:]


  me = this_image()
  np = num_images()
  
  write(*,*) 'Hello from ',me,' of:', np
  snd_buf(1:4) = me*(/1,2,3,4/)

  sync all
  if (me == 1) then
    do ip=1,np
      write(*,*) 'From ',ip,' :',snd_buf(:)[ip]
    end do
  end if
  
end program hello
