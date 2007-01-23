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
subroutine psi_dl_check(dep_list,dl_lda,np,length_dl)

  use psb_const_mod
  use psb_descriptor_type
  implicit none

  integer  :: np,dl_lda,length_dl(0:np)
  integer  :: dep_list(dl_lda,0:np)
  ! locals
  integer  :: proc, proc2, i, j

  ! ...i must order communication in in halo

  ! ...if in dep_list of process i there is j
  !     and in dep_list of process j there isn't i,
  !     add to it process i...

  do proc=0,np-1
    i=1
    outer: do 
      if (i >length_dl(proc)) exit outer
      proc2=dep_list(i,proc)
      if (proc2.ne.psb_no_comm_) then
        ! ...search proc in proc2's dep_list....
        j=1
        p2loop:do 
          if (j > length_dl(proc2)) exit p2loop
          if (dep_list(j,proc2) == proc) exit p2loop
          j=j+1
        enddo p2loop

        if (j > length_dl(proc2)) then
          ! ...add proc to proc2 s dep_list.....',proc,proc2
          length_dl(proc2)     = length_dl(proc2)+1
          if (length_dl(proc2) > size(dep_list,1)) then
            write(0,*)'error in crea_halo', proc2,proc,&
                 & length_dl(proc2),'>',size(dep_list,1)
          endif
          dep_list(length_dl(proc2),proc2) = proc
        else if (dep_list(j,proc2) /= proc) then 
          write(0,*) 'PSI_DL_CHECK This should not happen!!! ',&
               & j,proc2,dep_list(j,proc2),proc
        endif
      endif
      i=i+1
    enddo outer
  enddo

end subroutine psi_dl_check
