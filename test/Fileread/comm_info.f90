module comminfo
contains
  
  
  subroutine get_comminfo(icontxt,desc_a,comm_info)
    use typedesc
    implicit none

    type(desc_type) :: desc_a
    integer,pointer:: comm_info(:,:)
    integer        :: icontxt, nprow, npcol, myprow, mypcol,&
         & i,cnt,proc,n_elem_recv,n_elem_send
    integer,pointer:: sndbuf(:)



    call blacs_gridinfo(icontxt, nprow, npcol, myprow, mypcol)     
    
!    write(0,*)'inside comminfo',nprow,npcol,myprow,mypcol
    allocate(sndbuf(nprow))
    sndbuf(:)=0
    cnt=1

    do while(desc_a%halo_index(cnt).ne.-1)
       proc=desc_a%halo_index(cnt+proc_id_)
       n_elem_recv=desc_a%halo_index(cnt+n_elem_recv_)
       n_elem_send=desc_a%halo_index(cnt+n_elem_recv+n_elem_send_)
       cnt=cnt+n_elem_recv+n_elem_send+3
       sndbuf(proc+1)=n_elem_send
    end do

    
    if(myprow.eq.0) then
       comm_info(1,:)=sndbuf(:)
       deallocate(sndbuf)
       do i=1,nprow-1
          sndbuf=>comm_info(i+1,:)
!          call igerv2d( icontxt, 1, nprow, comm_info(i+1,:), nprow, i, 0)
          call igerv2d( icontxt, nprow,1, sndbuf, nprow, i, 0 )
!         write(0,'("Root has received from process n.",i3)'),i
!         write(0,*) comm_info(i+1,:)
!          write(0,*) sndbuf
       end do

    else

!      write(0,'("Process n.",i3," is sending to root")'),myprow
!      write(0,*) sndbuf
       call igesd2d( icontxt, nprow,1, sndbuf, nprow, 0, 0 )

    end if

  end subroutine get_comminfo


end module comminfo
