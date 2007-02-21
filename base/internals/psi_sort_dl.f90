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
subroutine psi_sort_dl(dep_list,l_dep_list,np,info)
  !
  !     interface between former sort_dep_list subroutine
  !     and new srtlist
  !
  use psi_mod, psb_protect_name => psi_sort_dl
  use psb_const_mod
  use psb_error_mod
  implicit none

  integer :: np,dep_list(:,:), l_dep_list(:)
  integer :: idg, iupd, idgp, iedges, iidx, iich,ndgmx, isz, err_act
  integer :: i, info
  integer, allocatable   :: work(:)
  logical, parameter :: debug=.false.
  character(len=20)        :: name
  
  name='psi_sort_dl'
  if(psb_get_errstatus() /= 0) return 
  info=0
  call psb_erractionsave(err_act)
  
  info = 0
  ndgmx = 0
  do i=1,np
     ndgmx = ndgmx + l_dep_list(i)
     if (debug) write(0,*) i,l_dep_list(i)
  enddo
  idg = 1
  iupd = idg+np
  idgp = iupd+np
  iedges = idgp + ndgmx
  iidx = iedges + 2*ndgmx
  iich = iidx + ndgmx
  isz = iich + ndgmx
  if (debug)write(0,*) 'psi_sort_dl: ndgmx ',ndgmx,isz
  
  allocate(work(isz))
  ! call srtlist(dep_list, dl_lda, l_dep_list, np, info)
  call srtlist(dep_list,size(dep_list,1),l_dep_list,np,work(idg),&
       & work(idgp),work(iupd),work(iedges),work(iidx),work(iich),info)

  if (info  /=  0) then
     call psb_errpush(4010,name,a_err='srtlist')
     goto 9999
  endif
  
  deallocate(work)
  call psb_erractionrestore(err_act)
  return  
  
9999 continue
  call psb_erractionrestore(err_act)
  
  if (err_act == psb_act_abort_) then
     call psb_error()
     return
  end if
  return
end subroutine psi_sort_dl


      
