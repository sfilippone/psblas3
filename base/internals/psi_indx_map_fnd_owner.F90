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
! File: psi_indx_map_fnd_owner.f90
!
! Subroutine: psi_indx_map_fnd_owner
!   Figure out who owns  global indices. 
! 
! Arguments: 
!    idx(:)   - integer                 Required indices on the calling process.
!                                       Note: the indices should be unique!
!    iprc(:)  - integer(psb_ipk_), allocatable    Output: process identifiers for the corresponding
!                                       indices
!    idxmap   - class(psb_indx_map).    The index map
!    info     - integer.                return code.
!
!  This is the default implementation of the FND_OWNER method.
!  If a particular index map class has additional information, it can override it
!  (see e.g. the GEN_BLOCK_MAP class).
!
!  1. Check if IDXM%PARTS is available, and use it; or
!  2. Check if TEMPVG(:) is allocated, and use it; or
!  3. Call the general method PSI_GRAPH_FND_OWNER. 
! 
subroutine psi_indx_map_fnd_owner(idx,iprc,idxmap,info)
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_realloc_mod
  use psb_indx_map_mod, psb_protect_name => psi_indx_map_fnd_owner
#ifdef MPI_MOD
  use mpi
#endif

  implicit none 
#ifdef MPI_H
  include 'mpif.h'
#endif
  integer(psb_lpk_), intent(in)      :: idx(:)
  integer(psb_ipk_), allocatable, intent(out) ::  iprc(:)
  class(psb_indx_map), intent(inout) :: idxmap
  integer(psb_ipk_), intent(out)     :: info


  integer(psb_ipk_), allocatable :: hhidx(:)
  integer(psb_mpk_) :: icomm, minfo, iictxt
  integer(psb_ipk_) :: i, err_act, hsize, nv
  integer(psb_lpk_) :: mglob
  integer(psb_ipk_) :: ictxt,np,me, nresp
  logical, parameter  :: gettime=.false.
  real(psb_dpk_)      :: t0, t1, t2, t3, t4, tamx, tidx
  character(len=20)   :: name

  info = psb_success_
  name = 'psb_indx_map_fnd_owner'
  call psb_erractionsave(err_act)

  ictxt   = idxmap%get_ctxt()
  icomm   = idxmap%get_mpic()
  mglob   = idxmap%get_gr()
  iictxt = ictxt 

  call psb_info(ictxt, me, np)

  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.(idxmap%is_valid())) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='invalid idxmap')
    goto 9999      
  end if

  if (gettime) then 
    t0 = psb_wtime()
  end if

  nv = size(idx)
  call psb_realloc(nv,iprc,info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_realloc')
    goto 9999      
  end if

  if (associated(idxmap%parts)) then 
    ! Use function shortcut
!!$    write(0,*) me,trim(name),' indxmap%parts shortcut'
    Allocate(hhidx(np), stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate') 
      goto 9999      
    end if
    do i=1, nv
      call idxmap%parts(idx(i),mglob,np,hhidx,nresp)
      if (nresp > 0) then
        iprc(i) = hhidx(1)
      else
        iprc(i) = -1 
      end if
    end do

  else if (allocated(idxmap%tempvg)) then 
!!$    write(0,*) me,trim(name),' indxmap%tempvg shortcut'
    ! Use temporary vector 
    do i=1, nv 
      iprc(i) = idxmap%tempvg(idx(i))
    end do

  else

    call psi_graph_fnd_owner(idx,iprc,idxmap,info)

  end if
  
  if (gettime) then 
    call psb_barrier(ictxt)
    t1 = psb_wtime()
    t1 = t1 -t0 - tamx - tidx   
    call psb_amx(ictxt,tamx)
    call psb_amx(ictxt,tidx)
    call psb_amx(ictxt,t1)
    if (me == psb_root_) then 
      write(psb_out_unit,'(" fnd_owner  idx time  : ",es10.4)') tidx
      write(psb_out_unit,'(" fnd_owner  amx time  : ",es10.4)') tamx
      write(psb_out_unit,'(" fnd_owner remainedr  : ",es10.4)') t1 
    endif
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psi_indx_map_fnd_owner
