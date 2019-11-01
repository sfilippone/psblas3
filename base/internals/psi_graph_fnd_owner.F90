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
! File: psi_fnd_owner.f90
!
! Subroutine: psi_fnd_owner
!   Figure out who owns  global indices. 
! 
! Arguments: 
!    nv       - integer                 Number of indices required on  the calling
!                                       process 
!    idx(:)   - integer                 Required indices on the calling process.
!                                       Note: the indices should be unique!
!    iprc(:)  - integer(psb_ipk_), allocatable    Output: process identifiers for the corresponding
!                                       indices
!    desc_a   - type(psb_desc_type).    The communication descriptor.        
!    info     - integer.                return code.
! 
subroutine psi_graph_fnd_owner(idx,iprc,idxmap,info)
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_realloc_mod
  use psb_indx_map_mod, psb_protect_name => psi_graph_fnd_owner
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


  integer(psb_lpk_), allocatable :: answers(:,:), idxsrch(:,:), hproc(:), tidx(:)
  integer(psb_ipk_), allocatable :: helem(:), hhidx(:), tprc(:), tsmpl(:), ladj(:)  
  integer(psb_mpk_), allocatable :: hsz(:),hidx(:), &
       & sdsz(:),sdidx(:), rvsz(:), rvidx(:)
  integer(psb_mpk_) :: icomm, minfo, iictxt
  integer(psb_ipk_) :: i,n_row,n_col,err_act,hsize,ip,isz,j, k,&
       & last_ih, last_j, nv, n_answers, n_rest, n_samples, locr_max, nrest_max, nadj, maxspace
  integer(psb_lpk_) :: mglob, ih
  integer(psb_ipk_) :: ictxt,np,me, nresp
  logical, parameter  :: gettime=.false.
  real(psb_dpk_)      :: t0, t1, t2, t3, t4
  character(len=20)   :: name

  info = psb_success_
  name = 'psi_graph_fnd_owner'
  call psb_erractionsave(err_act)

  ictxt   = idxmap%get_ctxt()
  icomm   = idxmap%get_mpic()
  mglob   = idxmap%get_gr()
  n_row   = idxmap%get_lr()
  n_col   = idxmap%get_lc()
  iictxt  = ictxt 

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

  locr_max = n_row
  call psb_max(ictxt,locr_max)
  !
  ! Choice of maxspace should be adjusted to account for a default
  ! "sensible" size and/or a user-specified value 
  maxspace = 2*locr_max

  nv = size(idx)
  call psb_realloc(nv,iprc,info)
  if (info == psb_success_) call psb_realloc(nv,tidx,info)
  if (info == psb_success_) call psb_realloc(nv,tprc,info)
  if (info == psb_success_) call psb_realloc(nv,tsmpl,info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_realloc')
    goto 9999      
  end if
  iprc(:) = -1
  n_answers = 0
  if (.true.) then 
    !
    ! Start from the adjacncy list
    !
    ! Skip for the time being

    n_rest    = nv - n_answers
    nrest_max = n_rest
    call psb_max(ictxt,nrest_max)
    fnd_owner_loop: do while (nrest_max>0)
      !
      ! The basic idea of this loop is to alternate between
      ! searching through all processes and searching
      ! in the neighbourood.
      !
      ! 1. Select a sample to be sent to all processes; sample
      !    size is such that the total size is <= maxspace
      !    
      if (me == 0) write(0,*) 'Looping in graph_fnd_owner: ', nrest_max
      n_samples = min(n_rest,max(1,(maxspace+np-1)/np))
      !
      ! Choose a sample, should it be done in this simplistic way? 
      !
      write(0,*) me,' Into first sampling ',n_samples
      call psi_get_sample(idx,iprc,tidx,tsmpl,n_samples,k)      
      n_samples = min(k,n_samples)
      write(0,*) me,' From first sampling ',n_samples
      ! 
      ! 2. Do a search on all processes; this is supposed to find
      !    the owning process for all inputs;
      !    
      call psi_a2a_fnd_owner(tidx(1:n_samples),tprc,idxmap,info)
      call psi_cpy_out(iprc,tprc,tsmpl,n_samples,k)      
      if (k /= n_samples) then 
        write(0,*) me,'Warning: indices not found by a2a_fnd_owner ',k,n_samples
      end if
      n_answers = n_answers + k
      n_rest    = nv - n_answers
      !
      ! 3. Extract the resulting adjacency list and add it to the
      !    indxmap;
      !
      ladj = tprc(1:n_samples)
      call psb_msort_unique(ladj,nadj)
      !
      ! NOTE: should symmetrize the list...
      ! 
      call idxmap%xtnd_p_adjcncy(ladj(1:nadj)) 

      !
      ! 4. Extract a sample and do a neighbourhood search so that the total
      !    size is <= maxspace
      !    (will not be exact since nadj varies with process)
      ! 
      n_samples = min(n_rest,max(1,(maxspace+max(1,nadj)-1))/(max(1,nadj)))
      write(0,*) me,' Into second sampling ',n_samples
      call psi_get_sample(idx,iprc,tidx,tsmpl,n_samples,k)      
      n_samples = min(k,n_samples)
      write(0,*) me,' From second sampling ',n_samples
      call psi_adjcncy_fnd_owner(tidx(1:n_samples),tprc,ladj(1:nadj),idxmap,info)
      call psi_cpy_out(iprc,tprc,tsmpl,n_samples,k)      
      n_answers = n_answers + k
      n_rest    = nv - n_answers
      nrest_max = n_rest
      call psb_max(ictxt,nrest_max)
    end do fnd_owner_loop

  else
    call psi_a2a_fnd_owner(idx,iprc,idxmap,info)
  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

contains
  subroutine psi_get_sample(idx,iprc,tidx,tsmpl,ns_in,ns_out)      
    implicit none
    integer(psb_lpk_), intent(in)  :: idx(:)
    integer(psb_ipk_), intent(in)  :: ns_in, iprc(:)
    integer(psb_lpk_), intent(out) :: tidx(:)
    integer(psb_ipk_), intent(out) :: tsmpl(:), ns_out
    !
    integer(psb_ipk_) :: j, nv 

    nv = size(idx)
    !
    ! Choose a sample, should it be done in this simplistic way? 
    !
    ns_out = 0
    do j=1, nv
      if (iprc(j) == -1) then
        ns_out = ns_out + 1
        tsmpl(ns_out) = j
        tidx(ns_out)  = idx(j)
      end if
      if (ns_out >= ns_in) exit
    end do
  end subroutine psi_get_sample

  subroutine psi_cpy_out(iprc,tprc,tsmpl,ns_in,ns_out)      
    implicit none
    integer(psb_ipk_), intent(out) :: iprc(:)
    integer(psb_ipk_), intent(in)  :: ns_in
    integer(psb_ipk_), intent(in)  :: tprc(:), tsmpl(:)
    integer(psb_ipk_), intent(out) :: ns_out

    integer(psb_ipk_) :: j
    
    ns_out = 0 
    do j=1, ns_in
      if (tprc(j) /= -1) then
        ns_out = ns_out + 1
        iprc(tsmpl(j)) = tprc(j)
      end if
    end do
  end subroutine psi_cpy_out
  
end subroutine psi_graph_fnd_owner
