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
subroutine psi_crea_ovr_elem(desc_overlap,ovr_elem,info)

  use psb_realloc_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  !     ...parameter arrays....      
  integer          :: desc_overlap(:)
  integer, allocatable, intent(inout)  :: ovr_elem(:)
  integer, intent(out) :: info

  !     ...local scalars...
  integer :: i,pnt_new_elem,ret,j,iret
  integer :: dim_ovr_elem
  integer :: pairtree(2)

  !     ...external function...
  integer  :: psi_exist_ovr_elem,dim
  external :: psi_exist_ovr_elem

  integer  :: nel, ip, ix, iel, insize, err_act
  integer, allocatable :: telem(:,:)

  logical, parameter :: usetree=.false.
  character(len=20)    :: name


  info = 0
  name='psi_crea_ovr_elem'
  

  if (allocated(ovr_elem)) then 
    dim_ovr_elem = size(ovr_elem)
  else
    dim_ovr_elem = 0
  endif


  if (usetree)  then 

    !
    ! While running through the column indices exchanged with other procs
    ! we have to record them in overlap_elem.  We do this by maintaining  
    ! an AVL balanced search tree: at each point counter_e is the next
    ! free index element. The search routine for gidx will return
    ! glx if gidx was already assigned a local index (glx<counter_e)
    ! but if gidx was a new index for this process, then it creates
    ! a new pair (gidx,counter_e), and glx==counter_e. In this case we
    ! need to record this for the overlap exchange. Otherwise it was 
    ! already there, so we need to record one more parnter in the exchange
    !

    i=1
    pnt_new_elem=1
    call initpairsearchtree(pairtree,info)
    do while (desc_overlap(i).ne.-1)
      !        ...loop over all procs of desc_overlap list....
      i=i+1
      do j=1,desc_overlap(i)
        !           ....loop over all overlap indices referred to act proc.....
        call searchinskeyval(pairtree,desc_overlap(i+j),pnt_new_elem,&
             & ret,info)
        if (ret == pnt_new_elem) ret=-1
        if (ret.eq.-1) then

          !            ...this point not exist in ovr_elem list:
          !               add to it.............................
          !              ...check if overflow element_d array......
          if ((pnt_new_elem +2) > dim_ovr_elem) then
            dim_ovr_elem=max(((3*dim_ovr_elem)/2+2),pnt_new_elem+100)
            call psb_realloc(dim_ovr_elem,ovr_elem,info)
            if (info /= 0) then 
              info = 4000
              call psb_errpush(info,name)
              goto 9999
            end if
          endif
          ovr_elem(pnt_new_elem)=desc_overlap(i+j)  
          ovr_elem(pnt_new_elem+1)=2             
          pnt_new_elem=pnt_new_elem+2              

        else
          !              ....this point already exist in ovr_elem list
          !                  its position is ret............................
          ovr_elem(ret+1)=ovr_elem(ret+1)+1
        endif
      enddo
      i=i+2*desc_overlap(i)+2
    enddo

    !  Add -1 at the end of output list. 
    !  And fix the size to the minimum necessary.
    dim_ovr_elem=pnt_new_elem
    call psb_realloc(dim_ovr_elem,ovr_elem,info)
    if (info /= 0) then 
      info = 4000
      call psb_errpush(info,name)
      goto 9999
    end if
    ovr_elem(pnt_new_elem)=-1
    call freepairsearchtree(pairtree)

  else

    insize = size(desc_overlap)
    insize = max(1,(insize+1)/2)
    allocate(telem(insize,2),stat=info)
    if (info /= 0) then 
      info = 4000
      call psb_errpush(info,name)
      goto 9999
    endif
    i   = 1
    nel = 0
    do while (desc_overlap(i).ne.-1)
      !        ...loop over all procs of desc_overlap list....

      i=i+1
      do j=1,desc_overlap(i)
        nel = nel + 1 
        telem(nel,1) = desc_overlap(i+j)
      enddo
      i=i+2*desc_overlap(i)+2
    enddo
    if (nel > 0) then 
      call imsr(nel,telem(:,1))
      iel        = telem(1,1)
      telem(1,2) = 2
      ix = 1
      ip = 2
      do 
        if (ip > nel) exit
        if (telem(ip,1) == iel) then 
          telem(ix,2) = telem(ix,2) + 1
        else
          ix = ix + 1
          telem(ix,1) = telem(ip,1)
          iel         = telem(ip,1)
          telem(ix,2) = 2
        end if
        ip = ip + 1
      end do
    else
      ix = 0
    end if
    dim_ovr_elem=2*ix+1
    call psb_realloc(dim_ovr_elem,ovr_elem,info)
    iel = 1
    do i=1, ix
      ovr_elem(iel)   = telem(i,1)
      ovr_elem(iel+1) = telem(i,2)
      iel = iel + 2
    end do
    ovr_elem(iel) = -1 
    deallocate(telem)
  endif
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
