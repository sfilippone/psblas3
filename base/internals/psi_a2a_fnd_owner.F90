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
! File: psi_a2a_fnd_owner.f90
!
! Subroutine: psi_a2a_fnd_owner
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
! This version does not assume any prior knowledge about the process topology,
! so it goes for an all-to-all by building an auxiliary neighbours list and
! reusing the neighbour version.
! 
subroutine psi_a2a_fnd_owner(idx,iprc,idxmap,info,samesize)
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_realloc_mod
  use psb_indx_map_mod, psb_protect_name => psi_a2a_fnd_owner
#ifdef MPI_MOD
  use mpi
#endif

  implicit none 
#ifdef MPI_H
  include 'mpif.h'
#endif
  integer(psb_lpk_), intent(in)   :: idx(:)
  integer(psb_ipk_), allocatable, intent(out) ::  iprc(:)
  class(psb_indx_map), intent(in) :: idxmap
  integer(psb_ipk_), intent(out)  :: info
  logical, intent(in), optional   :: samesize


  integer(psb_ipk_), allocatable :: tmpadj(:)
  integer(psb_lpk_), allocatable :: rmtidx(:)
  integer(psb_ipk_), allocatable :: tproc(:), lclidx(:)
  integer(psb_mpk_), allocatable :: hsz(:),hidx(:), sdidx(:), rvidx(:),&
       & sdsz(:), rvsz(:), sdhd(:), rvhd(:), p2pstat(:,:)
  integer(psb_mpk_) :: icomm, minfo, nv
  integer(psb_ipk_) :: i,n_row,n_col,err_act,gsz
  integer(psb_lpk_) :: mglob, ih
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np,me, nresp
  logical, parameter  :: use_psi_adj=.true.
  real(psb_dpk_)      :: t0, t1, t2, t3, t4, tamx, tidx
  character(len=20)   :: name
  logical             :: samesize_

  info = psb_success_
  name = 'psi_a2a_fnd_owner'
  call psb_erractionsave(err_act)

  ctxt   = idxmap%get_ctxt()
  icomm   = idxmap%get_mpic()
  mglob   = idxmap%get_gr()
  n_row   = idxmap%get_lr()
  n_col   = idxmap%get_lc()

  call psb_info(ctxt, me, np)

  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.(idxmap%is_valid())) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='invalid idxmap')
    goto 9999      
  end if

  if (present(samesize)) then
    samesize_ = samesize
  else
    samesize_ = .false.
  end if
  nv = size(idx)
  ! write(0,*) me,name,' :',use_psi_adj,samesize_,nv
  if (use_psi_adj) then 
    !
    ! Reuse the adjcncy version by tricking it with an adjcncy list
    ! that contains everybody but ME. 
    !
    call psb_realloc(np-1,tmpadj,info)
    tmpadj(1:me) = [(i,i=0,me-1)]
    tmpadj(me+1:np-1) = [(i,i=me+1,np-1)]
    call  psi_adjcncy_fnd_owner(idx,iprc,tmpadj,idxmap,info)

  else
    if (samesize_) then
      !
      ! Variant when IDX is guaranteed to have the same size on all
      ! processes. To be tested for performance: is it worth it?
      ! Probably yes. 
      !
      gsz = nv*np
      Allocate(rmtidx(gsz),lclidx(gsz),iprc(nv),stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
        goto 9999      
      end if
      call mpi_allgather(idx,nv,psb_mpi_lpk_,rmtidx,nv,psb_mpi_lpk_,icomm,minfo)
      call idxmap%g2l(rmtidx(1:gsz),lclidx(1:gsz),info,owned=.true.)
      !
      ! Reuse lclidx to encode owning process
      !
      do i=1, gsz
        if ((1<=lclidx(i)).and.(lclidx(i)<=n_row)) then
          lclidx(i) = me
        else
          lclidx(i) = -1
        end if
      end do
      call mpi_reduce_scatter_block(lclidx,iprc,nv,psb_mpi_ipk_,mpi_max,icomm,minfo)
      
    else
      !
      ! 1. allgetherv
      ! 2. local conversion
      ! 3. reduce_scatter
      !
      !
      ! The basic idea is very simple. 
      ! First we collect (to all) all the requests. 
      Allocate(hidx(np+1),hsz(np),stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate') 
        goto 9999      
      end if

      call mpi_allgather(nv,1,psb_mpi_mpk_,hsz,1,psb_mpi_mpk_,icomm,minfo)
      hidx(1)   = 0
      do i=1, np
        hidx(i+1) = hidx(i) + hsz(i)
      end do
      gsz = hidx(np+1)
      Allocate(rmtidx(gsz),lclidx(gsz),iprc(nv),stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
        goto 9999      
      end if

      call mpi_allgatherv(idx,hsz(me+1),psb_mpi_lpk_,&
           & rmtidx,hsz,hidx,psb_mpi_lpk_,&
           & icomm,minfo)

      call idxmap%g2l(rmtidx(1:gsz),lclidx(1:gsz),info,owned=.true.)
      !
      ! Reuse lclidx to encode owning process
      !
      do i=1, gsz
        if ((1<=lclidx(i)).and.(lclidx(i)<=n_row)) then
          lclidx(i) = me
        else
          lclidx(i) = -1
        end if
      end do
      call mpi_reduce_scatter(lclidx,iprc,hsz,psb_mpi_ipk_,mpi_max,icomm,minfo)
    end if
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psi_a2a_fnd_owner
