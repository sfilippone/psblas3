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
!
! File: psi_graph_fnd_owner.f90
!
! Subroutine: psi_graph_fnd_owner
!   Figure out who owns  global indices. 
! 
! Arguments: 
!    idx(:)   - integer                 Required indices on the calling process.
!                                       
!    iprc(:)  - integer(psb_ipk_), allocatable   Output: process identifiers
!                                                    for the corresponding indices
!    idxmap   - class(psb_indx_map).    The index map
!    info     - integer.                return code.
!
! This is the method to find out who owns a set of indices. 
! In principle we could do the following:
!   1. Do an allgatherv of IDX
!   2. For each of the collected indices figure out if current proces owns it
!   3. Scatter the results
!   4. Loop through the answers
! This method is guaranteed to find the owner, unless an input index has
! an invalid value, however it could easily require too much additional space
! because each block of indices is replicated to all processes.
! Therefore the current routine takes a different approach:
! -1. Figure out a maximum size for a buffer to collect the IDX; the buffer
!     should allow for at least one index from each process (i.e. min size NP); also
!     check if we have an adjacency list of processes on input; 
!  0. If the initial adjacency list is not empty, use psi_adj_fnd_sweep to go
!     through all indices and use  multiple calls to psi_adjcncy_fnd_owner
!     (up to the buffer size) to see if the owning processes are in the
!     initial neighbours list;
!  1. Extract a sample from IDX, up to the buffer size, and do a call
!     to psi_a2a_fnd_owner. This is guaranteed to find the owners of all indices
!     in the sample;
!  2. Build the list of processes that own the sample indices; these are
!     (a subset of) the topological neighbours, and store the list in IDXMAP
!  3. Use psi_adj_fnd_sweep to go through all remaining indices and use
!     multiple calls to psi_adjcncy_fnd_owner (up to the buffer size)
!     to see if the owning processes are in the current neighbours list;
!  4. If the input indices IDX have not been exhausted, cycle to 1.
!
!  Thus, we are alternating between asking all processes for a subset of indices,
!  and asking a subset of processes for all the indices,
!  thereby limiting the memory footprint to a predefined maximum
!  (that the user can force with psb_cd_set_maxspace()).
! 
subroutine psi_graph_fnd_owner(idx,iprc,idxmap,info)
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_realloc_mod
  use psb_timers_mod
  use psb_desc_mod, psb_protect_name => psi_graph_fnd_owner
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


  integer(psb_lpk_), allocatable :: tidx(:)
  integer(psb_ipk_), allocatable :: tprc(:), tsmpl(:), ladj(:)  
  integer(psb_mpk_) :: icomm, minfo
  integer(psb_ipk_) :: i,n_row,n_col,err_act,ip,j, nsampl_out,&
       & nv, n_answers, nqries, nsampl_in, locr_max, ist, iend,&
       & nqries_max, nadj, maxspace, nsampl, nlansw
  integer(psb_lpk_) :: mglob, ih
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np,me, nresp
  integer(psb_ipk_), parameter :: nt=4
  integer(psb_ipk_) :: tmpv(4)
  logical, parameter  :: do_timings=.false., trace=.false., debugsz=.false.
  integer(psb_ipk_), save  :: idx_sweep0=-1, idx_loop_a2a=-1, idx_loop_neigh=-1
  real(psb_dpk_)      :: t0, t1, t2, t3, t4
  character(len=20)   :: name

  info = psb_success_
  name = 'psi_graph_fnd_owner'
  call psb_erractionsave(err_act)

  ctxt   = idxmap%get_ctxt()
  icomm  = idxmap%get_mpic()
  mglob  = idxmap%get_gr()
  n_row  = idxmap%get_lr()
  n_col  = idxmap%get_lc()

  if ((do_timings).and.(idx_sweep0==-1))       &
       & idx_sweep0 = psb_get_timer_idx("GRPH_FND_OWN: Outer sweep")
  if ((do_timings).and.(idx_loop_a2a==-1))       &
       & idx_loop_a2a = psb_get_timer_idx("GRPH_FND_OWN: Loop a2a")
  if ((do_timings).and.(idx_loop_neigh==-1))       &
       & idx_loop_neigh = psb_get_timer_idx("GRPH_FND_OWN: Loop neigh")


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

  !
  !
  !
  nv = size(idx)
  call psb_realloc(nv,iprc,info)
  if (info == psb_success_) call psb_realloc(nv,tidx,info)
  if (info == psb_success_) call psb_realloc(nv,tprc,info)
  if (info == psb_success_) call psb_realloc(nv,tsmpl,info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_realloc')
    goto 9999      
  end if
  iprc(:)   = -1
  !
  ! Start from the adjacncy list
  !
  call psb_safe_ab_cpy(idxmap%p_adjcncy,ladj,info)
  nadj = psb_size(ladj)
  ! This makes ladj allocated with size 0 if needed, as opposed to unallocated
  call psb_realloc(nadj,ladj,info)
  !
  ! Throughout the subroutine, nqries is the number of local inquiries
  ! that have not been answered yet, n_answers is the number of queries
  ! that have been resolved. The queries/answers may be scattered
  ! through idx/iprc
  ! 
  n_answers  = 0
  nqries     = nv - n_answers
  nqries_max = nqries

  !
  ! Choice of maxspace should be adjusted to account for a default
  ! "sensible" size and/or a user-specified value
  ! Currently:
  ! 1. Use psb_cd_get_max_space()
  ! 2. If nt*locr_max > maxspace, use nt_loc_max
  ! 3. Should be at least NP 
  !
  tmpv(1) = nadj
  tmpv(2) = nqries_max
  tmpv(3) = n_row
  tmpv(4) = psb_cd_get_maxspace()
  call psb_max(ctxt,tmpv)
  nqries_max = tmpv(2)
  locr_max   = tmpv(3)
  maxspace   = nt*locr_max
  if (tmpv(4) > nt*locr_max)  maxspace = tmpv(4)
  maxspace = max(maxspace,np)
  if (trace.and.(me == 0)) write(0,*) ' Through graph_fnd_owner with maxspace:',&
       & maxspace, maxspace/np, nt*locr_max, psb_cd_get_maxspace()
  if (do_timings) call psb_tic(idx_sweep0)
  if ((tmpv(1) > 0).and.(tmpv(2) >0)) then
    !
    ! Do a preliminary run on the user-defined adjacency lists
    !
    if (debugsz) write(0,*) me,' Initial sweep on user-defined topology',nqries
    nsampl_in = min(nqries,max(1,(maxspace+max(1,nadj)-1))/(max(1,nadj)))
    if (trace.and.(me == 0)) write(0,*) ' Initial sweep on user-defined topology',&
         & nsampl_in
    call psi_adj_fnd_sweep(idx,iprc,ladj,idxmap,nsampl_in,n_answers)  
    call idxmap%xtnd_p_adjcncy(ladj) 
    nqries     = nv - n_answers
    nqries_max = nqries
    call psb_max(ctxt,nqries_max)
    if (trace.and.(me == 0)) write(0,*) ' After initial sweep:',nqries_max
    if (debugsz) write(0,*) me,' After sweep on user-defined topology',nqries_max
  end if
  if (do_timings) call psb_toc(idx_sweep0)
    
  fnd_owner_loop: do while (nqries_max>0)
    if (do_timings) call psb_tic(idx_loop_a2a)
    if (debugsz) write(0,*) me,' fnd_owner_loop',nqries_max
    !
    ! The basic idea of this loop is to alternate between
    ! searching through all processes and searching
    ! in the neighbourood.
    !
    ! 1. Select a sample such that the total size is <= maxspace
    !    sample query is then sent to all processes
    !    
    if (trace.and.(me == 0)) write(0,*) 'Looping in graph_fnd_owner: ', nqries_max
    nsampl_in = nqries
    nsampl_in = min(max(1,(maxspace+np-1)/np),nsampl_in)
    !
    ! Choose a sample, should it be done in this simplistic way?
    ! Note: nsampl_in is a hint, not an absolute, hence nsampl_out
    !
    call psi_get_sample(1,idx,iprc,tidx,tsmpl,iend,nsampl_in,nsampl_out)
    nsampl = min(nsampl_out,nsampl_in)
    if (debugsz) write(0,*) me,' From first sampling ',nsampl_in
    ! 
    ! 2. Do a search on all processes; this is supposed to find
    !    the owning process for all inputs;
    !    
    call psi_a2a_fnd_owner(tidx(1:nsampl),tprc,idxmap,info) 
    if (debugsz) write(0,*) me,' From a2a_fnd_owner ',info
        
    call psi_cpy_out(iprc,tprc,tsmpl,nsampl,nlansw)      
    if (nlansw < nsampl) then 
      write(0,*) me,'Warning: indices not found by a2a_fnd_owner ',nlansw,nsampl
    end if
    
    n_answers = n_answers + nlansw
    nqries    = nv - n_answers
    !
    ! 3. Extract the resulting adjacency list and add it to the
    !    indxmap;
    !
    ladj = tprc(1:nlansw)
    call psb_msort_unique(ladj,nadj)
    call psb_realloc(nadj,ladj,info)
    call idxmap%xtnd_p_adjcncy(ladj) 
    if (do_timings) call psb_toc(idx_loop_a2a)
    if (do_timings) call psb_tic(idx_loop_neigh)    
    !
    ! 4. Do a complete sweep over the queries using
    !    the adjacency list just computed.
    !    Rationale: 
    !    1. Only ask to the neighbours; any missing entries
    !       will eventually be found by the a2a step;
    !    2. Only use the adjacency list just recomputed: any
    !       currently open queries have already been tested
    !       on previous adjacency lists, no need to try them again.
    !
    if (trace.and.(me == 0)) write(0,*) ' Further sweep',nsampl_in
    call psi_adj_fnd_sweep(idx,iprc,ladj,&
         & idxmap,nsampl_in,n_answers)  

    nqries     = nv - n_answers
    nqries_max = nqries
    call psb_max(ctxt,nqries_max)
    if (trace.and.(me == 0)) write(0,*) ' fnd_owner_loop remaining:',nqries_max
    if (do_timings) call psb_toc(idx_loop_neigh)    
  end do fnd_owner_loop

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

contains

  !
  !  Get a sample.
  !  1. Start from  entry  istart;
  !  2. Collect unanswered queries
  !  3. Up to ns_in sample size;
  !  4. Record the actual sample size;
  !  5. Record where the sample stopped in case
  !     you need to complete a sweep through the data
  !  6. For each query, record where in the original vector
  !     it came from;
  !  7. There could be scattered answered/unanswered queries,
  !     so the code needs to skip existing answers; hence, the
  !     number of items sampled and the index where it stops
  !     differ. 
  !
  !
  subroutine psi_get_sample(istart,idx,iprc,tidx,tsmpl,iend,ns_in,ns_out)      
    implicit none
    
    integer(psb_lpk_), intent(in)    :: idx(:)
    integer(psb_ipk_), intent(in)    :: ns_in, iprc(:),istart
    integer(psb_lpk_), intent(out)   :: tidx(:)
    integer(psb_ipk_), intent(out)   :: tsmpl(:), ns_out,iend
    !
    integer(psb_ipk_) :: nv, ns, k, ipnt

    nv = size(idx)
    ns = ns_in
    !
    ! ns_in == 0 means that on the outside we figure there's
    ! nothing left, but we need to do something because the adj_sweep is
    ! sinchronized 
    ! Make sure we sweep through the entire vector immediately.
    ! But also make sure we do not overrun tsmpl
    !
    ns     = min(ns,size(tsmpl))
    ns_out = 0
    ipnt   = istart-1
    do while(ipnt<nv)
      ipnt = ipnt+1
      if (iprc(ipnt) == -1) then
        ns_out        = ns_out + 1
        tsmpl(ns_out) = ipnt
        tidx(ns_out)  = idx(ipnt)
      end if
      if (ns_out >= ns) exit
    end do
    iend = ipnt
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

  subroutine psi_adj_fnd_sweep(idx,iprc,adj,idxmap,n_samples,n_answers)      
    implicit none
    integer(psb_lpk_), intent(in)    :: idx(:)
    integer(psb_ipk_), intent(in)    :: n_samples
    integer(psb_ipk_), intent(inout) :: iprc(:), n_answers
    integer(psb_ipk_), intent(in)    :: adj(:)
    class(psb_indx_map), intent(inout) :: idxmap
    !
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: ipnt, ns_in, ns_out, n_rem, me, np, isw, n_reml,iend, nv
    integer(psb_lpk_), allocatable    :: tidx(:)
    integer(psb_ipk_), allocatable    :: tsmpl(:)

    ctxt   = idxmap%get_ctxt()
    call psb_info(ctxt,me,np)
    call psb_realloc(n_samples,tidx,info)
    call psb_realloc(n_samples,tsmpl,info)
    ipnt = 1
    isw  = 1
    nv   = size(idx)
    do
      ! Sweep through the vector, one section at a time,
      ! up to N_SAMPLES samples. The sections are unpredictable, because
      ! the queries are scattered; hence the need for get_sample
      ! to tell us where the current section ends
      !
      call psi_get_sample(ipnt,idx,iprc,tidx,tsmpl,iend,n_samples,ns_out)
      ns_in = min(n_samples,ns_out)
      !
      call psi_adjcncy_fnd_owner(tidx(1:ns_in),tprc,ladj,idxmap,info)
      call psi_cpy_out(iprc,tprc,tsmpl,ns_in,ns_out)
      !
      ! Update starting point of next sweep and number of remaining
      ! queries to check for end of loop.
      !
      n_answers = n_answers + ns_out
      ipnt      = iend + 1
      n_reml    = nv - ipnt + 1
      n_rem     = n_reml
      call psb_max(ctxt,n_rem)
      ! if (me == 0) write(0,*) me,' fnd_sweep Sweep ',isw,n_rem, ipnt, n_samples, n_reml
      isw = isw + 1 
      if (n_rem <= 0) exit
    end do
    !  if (me == 0) write(0,*)'adj_fnd_sweep: sweeps: ',isw

  end subroutine psi_adj_fnd_sweep

end subroutine psi_graph_fnd_owner
