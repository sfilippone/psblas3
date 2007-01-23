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
! File: psb_cdprt.f90
!
! Subroutine: psb_cdprt
!    Prints the descriptor to an output file
! 
! Parameters: 
!    iout          - integer.                The output unit to print to.
!    desc_p        - type(<psb_desc_type>).  The communication descriptor to be printed.
!    glob          - logical(otpional).      Wheter to print out global or local data.
!    short         - logical(optional).      Used to choose a verbose output.
subroutine psb_cdprt(iout,desc_p,glob,short)
  use psb_const_mod
  use psb_descriptor_type
  implicit none 
  type(psb_desc_type), intent(in)    :: desc_p
  integer, intent(in)                :: iout
  logical, intent(in), optional      :: glob,short
  logical :: lshort, lglob

  integer :: m, n_row, n_col,counter,idx,n_elem_recv,n_elem_send,&
       & proc,i

  if (present(glob)) then 
    lglob = glob
  else
    lglob = .false.
  endif
  if (present(short)) then 
    lshort = short
  else
    lshort = .true.
  endif

  if (.not.lglob) then
    write(iout,*) 'Communication descriptor:',desc_p%matrix_data(1:10)
    m=desc_p%matrix_data(psb_m_)
    n_row=desc_p%matrix_data(psb_n_row_)
    n_col=desc_p%matrix_data(psb_n_col_)
    if (.not.lshort) &
         & write(iout,*) 'Loc_to_glob ',desc_p%loc_to_glob(1:n_row), ': ',&
         & desc_p%loc_to_glob(n_row+1:n_col)

!!$    if (.not.lshort) write(iout,*) 'glob_to_loc ',desc_p%glob_to_loc(1:m) 
    write(iout,*) 'Halo_index'
    counter      = 1
    Do
      proc=desc_p%halo_index(counter+psb_proc_id_)
      if (proc == -1) exit
      n_elem_recv=desc_p%halo_index(counter+psb_n_elem_recv_)
      n_elem_send=desc_p%halo_index(counter+n_elem_recv+psb_n_elem_send_)
      write(iout,*) 'Halo_index Receive',proc,n_elem_recv
      if (.not.lshort) write(iout,*) &
           & desc_p%halo_index(counter+psb_n_elem_recv_+1:counter+psb_n_elem_recv_+n_elem_recv)
      write(iout,*) 'Halo_index Send',proc,n_elem_send
      if (.not.lshort) write(iout,*) &
           & desc_p%halo_index(counter+n_elem_recv+psb_n_elem_send_+1: &
           &                   counter+n_elem_recv+psb_n_elem_send_+n_elem_send)

      counter   = counter+n_elem_recv+n_elem_send+3
    enddo

    write(iout,*) 'Ext_index'
    counter      = 1
    Do
      proc=desc_p%ext_index(counter+psb_proc_id_)
      if (proc == -1) exit
      n_elem_recv=desc_p%ext_index(counter+psb_n_elem_recv_)
      n_elem_send=desc_p%ext_index(counter+n_elem_recv+psb_n_elem_send_)
      write(iout,*) 'Ext_index Receive',proc,n_elem_recv
      if (.not.lshort) write(iout,*) &
           & desc_p%ext_index(counter+psb_n_elem_recv_+1:counter+psb_n_elem_recv_+n_elem_recv)
      write(iout,*) 'Ext_index Send',proc,n_elem_send
      if (.not.lshort) write(iout,*) &
           & desc_p%ext_index(counter+n_elem_recv+psb_n_elem_send_+1: &
           &                   counter+n_elem_recv+psb_n_elem_send_+n_elem_send)

      counter   = counter+n_elem_recv+n_elem_send+3
    enddo


    write(iout,*) 'Ovrlap_index'
    counter      = 1
    Do
      proc=desc_p%ovrlap_index(counter+psb_proc_id_)
      if (proc == -1) exit
      n_elem_recv=desc_p%ovrlap_index(counter+psb_n_elem_recv_)
      n_elem_send=desc_p%ovrlap_index(counter+n_elem_recv+psb_n_elem_send_)
      write(iout,*) 'Ovrlap_index Receive',proc,n_elem_recv
      if (.not.lshort) write(iout,*) &
           & desc_p%ovrlap_index(counter+psb_n_elem_recv_+1:&
           &    counter+psb_n_elem_recv_+n_elem_recv)
      write(iout,*) 'Ovrlap_index Send',proc,n_elem_send
      if (.not.lshort) write(iout,*) &
           & desc_p%ovrlap_index(counter+n_elem_recv+psb_n_elem_send_+1: &
           &                   counter+n_elem_recv+psb_n_elem_send_+n_elem_send)

      counter   = counter+n_elem_recv+n_elem_send+3
    enddo

    write(iout,*) 'Ovrlap_elem'
    counter      = 1
    Do
      idx=desc_p%ovrlap_elem(counter)
      if (idx == -1) exit
      n_elem_recv=desc_p%ovrlap_elem(counter+1)
      if (.not.lshort) write(iout,*) idx,n_elem_Recv
      counter   = counter+2
    enddo

  else if (lglob) then 

    write(iout,*) 'Communication descriptor:',desc_p%matrix_data(1:10)
    m=desc_p%matrix_data(psb_m_)
    n_row=desc_p%matrix_data(psb_n_row_)
    n_col=desc_p%matrix_data(psb_n_col_)
    if (.not.lshort) then 
      write(iout,*) 'Loc_to_glob '
      do i=1, n_row
        write(iout,*) i, desc_p%loc_to_glob(i)
      enddo
      write(iout,*) '........'
      do i=n_row+1,n_col
        write(iout,*) i, desc_p%loc_to_glob(i)
      enddo

!!$      write(iout,*) 'glob_to_loc '
!!$      do i=1,m
!!$        write(iout,*) i,desc_p%glob_to_loc(i)
!!$      enddo
    endif
    write(iout,*) 'Halo_index'
    counter      = 1
    Do
      proc=desc_p%halo_index(counter+psb_proc_id_)
      if (proc == -1) exit
      n_elem_recv=desc_p%halo_index(counter+psb_n_elem_recv_)
      n_elem_send=desc_p%halo_index(counter+n_elem_recv+psb_n_elem_send_)
      write(iout,*) 'Halo_index Receive',proc,n_elem_recv
      if (.not.lshort) then 
        do i=counter+psb_n_elem_recv_+1,counter+psb_n_elem_recv_+n_elem_recv
          write(iout,*) &
               & desc_p%loc_to_glob(desc_p%halo_index(i)),desc_p%halo_index(i)
        enddo
      endif
      write(iout,*) 'Halo_index Send',proc,n_elem_send
      if (.not.lshort) then 
        do i=counter+n_elem_recv+psb_n_elem_send_+1, &
             & counter+n_elem_recv+psb_n_elem_send_+n_elem_send
          write(iout,*) &
               & desc_p%loc_to_glob(desc_p%halo_index(i)),  desc_p%halo_index(i)
        enddo
      endif
      counter   = counter+n_elem_recv+n_elem_send+3
    enddo

    write(iout,*) 'Ext_index'
    counter      = 1
    Do
      proc=desc_p%ext_index(counter+psb_proc_id_)
      if (proc == -1) exit
      n_elem_recv=desc_p%ext_index(counter+psb_n_elem_recv_)
      n_elem_send=desc_p%ext_index(counter+n_elem_recv+psb_n_elem_send_)
      write(iout,*) 'Ext_index Receive',proc,n_elem_recv
      if (.not.lshort) then 
        do i=counter+psb_n_elem_recv_+1,counter+psb_n_elem_recv_+n_elem_recv
          write(iout,*) &
               & desc_p%loc_to_glob(desc_p%ext_index(i)),desc_p%ext_index(i)
        enddo
      endif
      write(iout,*) 'Ext_index Send',proc,n_elem_send
      if (.not.lshort) then 
        do i=counter+n_elem_recv+psb_n_elem_send_+1, &
             & counter+n_elem_recv+psb_n_elem_send_+n_elem_send
          write(iout,*) &
               & desc_p%loc_to_glob(desc_p%ext_index(i)),  desc_p%ext_index(i)
        enddo
      endif
      counter   = counter+n_elem_recv+n_elem_send+3
    enddo


    write(iout,*) 'Ovrlap_index'
    counter      = 1
    Do
      proc=desc_p%ovrlap_index(counter+psb_proc_id_)
      if (proc == -1) exit
      n_elem_recv=desc_p%ovrlap_index(counter+psb_n_elem_recv_)
      n_elem_send=desc_p%ovrlap_index(counter+n_elem_recv+psb_n_elem_send_)
      write(iout,*) 'Ovrlap_index Receive',proc,n_elem_recv
      if (.not.lshort) then 
        do i=counter+psb_n_elem_recv_+1,counter+psb_n_elem_recv_+n_elem_recv
          write(iout,*) desc_p%loc_to_glob(desc_p%ovrlap_index(i)),&
               & desc_p%ovrlap_index(i)
        enddo
      endif
      write(iout,*) 'Ovrlap_index Send',proc,n_elem_send
      if (.not.lshort) then 
        do i=counter+n_elem_recv+psb_n_elem_send_+1, &
             &                   counter+n_elem_recv+psb_n_elem_send_+n_elem_send
          write(iout,*) desc_p%loc_to_glob(desc_p%ovrlap_index(i)),&
               &   desc_p%ovrlap_index(i)
        enddo
      endif
      counter   = counter+n_elem_recv+n_elem_send+3
    enddo

    write(iout,*) 'Ovrlap_elem'
    counter      = 1
    if (.not.lshort) then 
      Do
        idx=desc_p%ovrlap_elem(counter)
        if (idx == -1) exit
        n_elem_recv=desc_p%ovrlap_elem(counter+1)
        write(iout,*) desc_p%loc_to_glob(idx),idx,n_elem_Recv
        counter   = counter+2
      enddo
    endif
  end if
end subroutine psb_cdprt
