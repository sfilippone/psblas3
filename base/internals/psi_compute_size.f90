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
subroutine psi_compute_size(desc_data, index_in, dl_lda, info)

  use psb_const_mod
  use psb_descriptor_type
  use psb_error_mod
  use psb_penv_mod
  implicit none

  !     ....scalars parameters....
  integer  :: info, dl_lda
  !     .....array parameters....
  integer  :: desc_data(:), index_in(:)
  !     ....local scalars....      
  integer  :: i,np,me,proc, max_index
  integer  :: ictxt, err_act
  !     ...local array...
  integer  :: int_err(5)
  integer, allocatable :: counter_recv(:), counter_dl(:)

  !     ...parameters
  logical, parameter :: debug=.false.
  character(len=20)  :: name

  name='psi_compute_size'
  call psb_get_erraction(err_act)

  info = 0
  ictxt = desc_data(psb_ctxt_)

  call psb_info(ictxt,me,np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(counter_dl(0:np-1),counter_recv(0:np-1),stat=info)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if

  !     ..initialize counters...
  do i=0,np-1
    counter_recv(i)=0
    counter_dl(i)=0
  enddo

  !     ....verify local correctness of halo_in....
  i=1
  do while (index_in(i).ne.-1)
    proc=index_in(i)
    if ((proc.gt.np-1).or.(proc.lt.0)) then
      info = 115
      int_err(1) = 11
      int_err(2) = proc
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    endif
    counter_dl(proc)=1

    !        ..update no of elements to receive from proc proc..         
    counter_recv(proc)=counter_recv(proc)+&
         & index_in(i+1)

    i=i+index_in(i+1)+2
  enddo

  !     ...computing max_halo: max halo points to be received from
  !                            same processor
  max_index=0
  dl_lda=0

  do i=0,np-1
    if (counter_recv(i).gt.max_index) max_index = counter_recv(i)
    if (counter_dl(i).eq.1) dl_lda = dl_lda+1
  enddo

  !     computing max global value of dl_lda
  call psb_amx(ictxt, dl_lda)

  if (debug) then 
    write(0,*) 'psi_compute_size: ',dl_lda
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psi_compute_size

         

