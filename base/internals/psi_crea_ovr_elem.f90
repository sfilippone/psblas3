!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
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
!
! File: psi_crea_ovr_elem.f90
!
! Subroutine: psi_crea_ovr_elem
!   Creates the overlap_elem list: for each overlap index, store the index and 
!   the number of processes sharing it (minimum: 2). List is ended by -1.
!   See also description in base/modules/psb_desc_type.f90
! 
! Arguments: 
!    ovr_elem(:,:) - integer(psb_ipk_), allocatable  Array containing the output list              
!    desc_a   - type(psb_desc_type).       The communication descriptor.        
!    info     - integer.                   return code.
! 
subroutine psi_crea_ovr_elem(me,desc_overlap,ovr_elem,info)

  use psi_mod, psb_protect_name => psi_crea_ovr_elem
  use psb_realloc_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_serial_mod
  implicit none

  !     ...parameter arrays....      
  integer(psb_ipk_), intent(in)               :: me, desc_overlap(:)
  integer(psb_ipk_), allocatable, intent(out) :: ovr_elem(:,:)
  integer(psb_ipk_), intent(out)              :: info

  !     ...local scalars...
  integer(psb_ipk_) :: i,pnt_new_elem,ret,j
  integer(psb_ipk_) :: dim_ovr_elem
  integer(psb_ipk_) :: pairtree(2)

  !     ...external function...
  integer(psb_ipk_) :: psi_exist_ovr_elem
  external :: psi_exist_ovr_elem

  integer(psb_ipk_) :: nel, ip, ix, iel, insize, err_act, iproc
  integer(psb_ipk_), allocatable :: telem(:,:)

  character(len=20)    :: name


  info = psb_success_
  name='psi_crea_ovr_elem'


  if (allocated(ovr_elem)) then 
    dim_ovr_elem = size(ovr_elem,1)
  else
    dim_ovr_elem = 0
  endif


  insize = size(desc_overlap)
  insize = max(1,(insize+1)/2)
  allocate(telem(insize,3),stat=info)
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  endif
  i   = 1
  nel = 0
  do while (desc_overlap(i) /= -1)
    !        ...loop over all procs of desc_overlap list....
    iproc = desc_overlap(i)
    i     = i+1
    do j=1,desc_overlap(i)
      nel = nel + 1 
      telem(nel,1) = desc_overlap(i+j)
      telem(nel,2) = 1
      telem(nel,3) = iproc
    enddo
    i=i+2*desc_overlap(i)+2
  enddo

  if (nel > 0) then 
    call psb_msort(telem(1:nel,1),ix=telem(1:nel,3),flag=psb_sort_keep_idx_)

    iel        = telem(1,1)
    telem(1,2) = 2
    telem(1,3) = min(me,telem(1,3))
    ix = 1
    ip = 2
    do 
      if (ip > nel) exit
      if (telem(ip,1) == iel) then 
        telem(ix,2) = telem(ix,2) + 1
        telem(ix,3) = min(telem(ix,3),telem(ip,3))
      else
        ix = ix + 1
        telem(ix,1) = telem(ip,1)
        iel         = telem(ip,1)
        telem(ix,2) = 2
        telem(ix,3) = min(me,telem(ip,3))
      end if
      ip = ip + 1
    end do
  else
    ix = 0
  end if

  nel = ix

  call psb_realloc(nel,3,telem,info)
  call psb_move_alloc(telem,ovr_elem,info) 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psi_crea_ovr_elem
