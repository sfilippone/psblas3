!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
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
