!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!    
!
! File: psi_dl_check.f90
!
! Subroutine: psi_dl_check
!   Make sure a dependency list is symmetric, i.e. if process i depends on j
!   then process j should depend on i (even if the data to be sent in one of the
!   directions happens to be empty)
! 
! Arguments: 
!    dep_list(:,:) - integer             Initial dependency lists
!    dl_lda        - integer             Allocated size of dep_list
!    np            - integer             Total number of processes.
!    length_dl(:)  - integer             Items in dependency lists; updated on 
!                                        exit
! 
subroutine psi_dl_check(dep_list,dl_lda,np,length_dl)

  use psi_mod, psb_protect_name => psi_dl_check
  use psb_const_mod
  use psb_desc_mod
  implicit none

  integer(psb_ipk_) :: np,dl_lda,length_dl(0:np)
  integer(psb_ipk_) :: dep_list(dl_lda,0:np)
  ! locals
  integer(psb_ipk_) :: proc, proc2, i, j


  ! ...if j is in  dep_list of process i 
  !     and i is not in dep_list of process j 
  !     fix it.

  do proc=0,np-1
    i=1
    outer: do 
      if (i >length_dl(proc)) exit outer
      proc2=dep_list(i,proc)
      if ((proc2 /= -1).and.(proc2 /= proc)) then
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
            write(psb_err_unit,*)'error in dl_check', proc2,proc,&
                 & length_dl(proc2),'>',size(dep_list,1)
          endif
          dep_list(length_dl(proc2),proc2) = proc
        else if (dep_list(j,proc2) /= proc) then 
          write(psb_err_unit,*) 'PSI_DL_CHECK This should not happen!!! ',&
               & j,proc2,dep_list(j,proc2),proc
        endif
      endif
      i=i+1
    enddo outer
  enddo

end subroutine psi_dl_check
