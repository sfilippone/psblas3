!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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

  integer(psb_ipk_) :: np,dep_list(:,:), l_dep_list(:)
  integer(psb_ipk_) :: idg, iupd, idgp, iedges, iidx, iich,ndgmx, isz, err_act
  integer(psb_ipk_) :: i, info
  integer(psb_ipk_), allocatable :: work(:)
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)    :: name
  
  name='psi_sort_dl'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  
  info = psb_success_
  ndgmx = 0
  do i=1,np
     ndgmx = ndgmx + l_dep_list(i)
     if (debug_level >= psb_debug_inner_)&
          & write(debug_unit,*) name,': ',i,l_dep_list(i)
  enddo
  idg = 1
  iupd = idg+np
  idgp = iupd+np
  iedges = idgp + ndgmx
  iidx = iedges + 2*ndgmx
  iich = iidx + ndgmx
  isz = iich + ndgmx
  if (debug_level >= psb_debug_inner_)&
       & write(debug_unit,*) name,': ndgmx ',ndgmx,isz

  allocate(work(isz))
  ! call srtlist(dep_list, dl_lda, l_dep_list, np, info)
  call srtlist(dep_list,size(dep_list,1,kind=psb_ipk_),l_dep_list,np,work(idg),&
       & work(idgp),work(iupd),work(iedges),work(iidx),work(iich),info)
  if (info  /=  psb_success_) then
     call psb_errpush(psb_err_from_subroutine_,name,a_err='srtlist')
     goto 9999
  endif
  
  deallocate(work)
  call psb_erractionrestore(err_act)
  return  
  
9999 call psb_error_handler(err_act)

  return

end subroutine psi_sort_dl


      
