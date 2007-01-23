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
subroutine psi_crea_bnd_elem(bndel,desc_a,info)
  use psb_realloc_mod
  use psb_descriptor_type
  use psb_error_mod
  implicit none
  
  integer, allocatable :: bndel(:)
  type(psb_desc_type)  :: desc_a
  integer, intent(out) :: info

  integer, allocatable :: work(:)
  integer :: i, j, nr, ns, k, irv, err_act
  character(len=20)    :: name

  info = 0
  name='psi_crea_bnd_elem'
  call psb_erractionsave(err_act)

  allocate(work(size(desc_a%halo_index)),stat=info)
  if (info /= 0 ) then 
    info = 4000
    call psb_errpush(info,name)
    goto 9999
  end if

  i=0
  j=1
  do while(desc_a%halo_index(j) /= -1) 

    nr = desc_a%halo_index(j+1)
    ns = desc_a%halo_index(j+1+nr+1)
    do k=1, ns
      i = i + 1
      work(i) = desc_a%halo_index(j+1+nr+1+k)
    enddo
    j  = j + 1 + ns + 1 + nr + 1
  enddo

  if (i>0) then 
    call isr(i,work)
    j=1
    irv = work(1)
    do k=2, i
      if (work(k) /= irv) then
        irv = work(k)
        j = j + 1 
        work(j) = work(k)
      endif
    enddo
  else 
    j = 0
  endif


  if (.true.) then 
    if (j>0) then 
      call psb_realloc(j,bndel,info)
      if (info /= 0) then 
        call psb_errpush(4010,name,a_err='Allocate')
        goto 9999      
      end if
      bndel(1:j) = work(1:j)
    else
      if (allocated(bndel)) then 
        deallocate(bndel)
      end if
    end if
  else
    call psb_realloc(j+1,bndel,info)
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if
    bndel(1:j) = work(1:j)
    bndel(j+1) = -1
  endif

  deallocate(work)
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine psi_crea_bnd_elem
