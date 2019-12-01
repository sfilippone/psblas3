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
  integer(psb_mpk_), allocatable :: hsz(:),hidx(:)
  integer(psb_mpk_) :: icomm, minfo, iictxt
  integer(psb_ipk_) :: i,n_row,n_col,err_act,ip,j,ipnt, nsampl_out,&
       & nv, n_answers, n_rest, nsampl_in, locr_max, &
       & nrest_max, nadj, maxspace, mxnsin
  integer(psb_lpk_) :: mglob, ih
  integer(psb_ipk_) :: ictxt,np,me, nresp
  integer(psb_ipk_), parameter :: nt=4
  integer(psb_ipk_) :: tmpv(2)
  logical, parameter  :: gettime=.false., trace=.false.
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

  !
  ! Choice of maxspace should be adjusted to account for a default
  ! "sensible" size and/or a user-specified value
  !
  tmpv(1) = n_row
  tmpv(2) = psb_cd_get_maxspace()
  call psb_max(ictxt,tmpv)
  locr_max = tmpv(1)
  maxspace = min(nt*locr_max,tmpv(2))
  maxspace = max(maxspace,np)
  if (trace.and.(me == 0)) write(0,*) ' Through graph_fnd_owner with maxspace:',maxspace
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
  iprc(:) = -1
  n_answers = 0
  !
  ! Start from the adjacncy list
  !
  call psb_safe_ab_cpy(idxmap%p_adjcncy,ladj,info)
  nadj = psb_size(ladj)
  ! This makes ladj allocated with size 0 just in case 
  call psb_realloc(nadj,ladj,info)
  n_rest    = nv - n_answers
  nrest_max = n_rest

  tmpv(1) = nadj
  tmpv(2) = nrest_max
  call psb_max(ictxt,tmpv)
  if ((tmpv(1) > 0).and.(tmpv(2) >0)) then
    !
    ! Do a preliminary run on the user-defined adjacency lists
    !
    if (trace.and.(me == 0)) write(0,*) ' Initial sweep on user-defined topology'
    nsampl_in = min(n_rest,max(1,(maxspace+max(1,nadj)-1))/(max(1,nadj)))
    call psi_adj_fnd_sweep(idx,iprc,ladj,idxmap,nsampl_in,n_answers)  
    call idxmap%xtnd_p_adjcncy(ladj) 
    n_rest    = nv - n_answers
    nrest_max = n_rest
    call psb_max(ictxt,nrest_max)
    if (trace.and.(me == 0)) write(0,*) ' After initial sweep:',nrest_max
  end if

    
  fnd_owner_loop: do while (nrest_max>0)
    !
    ! The basic idea of this loop is to alternate between
    ! searching through all processes and searching
    ! in the neighbourood.
    !
    ! 1. Select a sample such that the total size is <= maxspace
    !    sample query is then sent to all processes
    !    
    ! if (trace.and.(me == 0)) write(0,*) 'Looping in graph_fnd_owner: ', nrest_max
    nsampl_in = min(n_rest,max(1,(maxspace+np-1)/np))
    !
    ! Choose a sample, should it be done in this simplistic way?
    ! Note: nsampl_in is a hint, not an absolute, hence nsampl_out
    !
    ipnt = 1
!!$      write(0,*) me,' Into first sampling ',nsampl_in
    call psi_get_sample(ipnt, idx,iprc,tidx,tsmpl,nsampl_in,nsampl_out)      
    nsampl_in = min(nsampl_out,nsampl_in)
!!$      write(0,*) me,' From first sampling ',nsampl_in
    ! 
    ! 2. Do a search on all processes; this is supposed to find
    !    the owning process for all inputs;
    !    
    call psi_a2a_fnd_owner(tidx(1:nsampl_in),tprc,idxmap,info)
    call psi_cpy_out(iprc,tprc,tsmpl,nsampl_in,nsampl_out)      
    if (nsampl_out /= nsampl_in) then 
      write(0,*) me,'Warning: indices not found by a2a_fnd_owner ',nsampl_out,nsampl_in
    end if
    n_answers = n_answers + nsampl_out
    n_rest    = nv - n_answers
    !
    ! 3. Extract the resulting adjacency list and add it to the
    !    indxmap;
    !
    ladj = tprc(1:nsampl_in)
    call psb_msort_unique(ladj,nadj)
    call psb_realloc(nadj,ladj,info)

    !
    ! 4. Extract again a sample and do a neighbourhood search
    !    so that the total size is <= maxspace
    !    (will not be exact since nadj varies with process)
    !    Need to set up a proper loop here to have a complete
    !    sweep over the input vector. Done inside adj_fnd_sweep. 
    !
!!$      write(0,*) me,' After a2a ',n_rest
    nsampl_in = min(n_rest,max(1,(maxspace+max(1,nadj)-1))/(max(1,nadj)))
    mxnsin = nsampl_in
    call psb_max(ictxt,mxnsin)
!!$      write(0,*) me, ' mxnsin ',mxnsin
    if (mxnsin>0) call psi_adj_fnd_sweep(idx,iprc,ladj,idxmap,nsampl_in,n_answers)  
    call idxmap%xtnd_p_adjcncy(ladj) 

    n_rest    = nv - n_answers
    nrest_max = n_rest
    call psb_max(ictxt,nrest_max)
    if (trace.and.(me == 0)) write(0,*) ' fnd_owner_loop remaining:',nrest_max
  end do fnd_owner_loop

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

contains

  subroutine psi_get_sample(ipntidx,idx,iprc,tidx,tsmpl,ns_in,ns_out)      
    implicit none
    integer(psb_ipk_), intent(inout) :: ipntidx
    integer(psb_lpk_), intent(in)    :: idx(:)
    integer(psb_ipk_), intent(in)    :: ns_in, iprc(:)
    integer(psb_lpk_), intent(out)   :: tidx(:)
    integer(psb_ipk_), intent(out)   :: tsmpl(:), ns_out
    !
    integer(psb_ipk_) :: nv, ns

    nv = size(idx)
    !
    ! Choose a sample, should it be done in this simplistic way? 
    !
    ns = ns_in
    !
    ! ns_in == 0 means that on the outside we figure there's
    ! nothing left, but we are here because we have to synchronize.
    ! Make sure we sweep through the entire vector immediately
    ! 
    if (ns == 0) ns = nv
    ns_out = 0

    do while (ipntidx<= nv)
      if (iprc(ipntidx) == -1) then
        ns_out        = ns_out + 1
        tsmpl(ns_out) = ipntidx
        tidx(ns_out)  = idx(ipntidx)
      end if
      ipntidx = ipntidx + 1
      if (ns_out >= ns) exit
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

  subroutine psi_adj_fnd_sweep(idx,iprc,adj,idxmap,n_samples,n_answers)      
    implicit none
    integer(psb_lpk_), intent(in)    :: idx(:)
    integer(psb_ipk_), intent(in)    :: n_samples
    integer(psb_ipk_), intent(inout) :: iprc(:), n_answers
    integer(psb_ipk_), intent(in)    :: adj(:)
    class(psb_indx_map), intent(inout) :: idxmap
    !
    integer(psb_ipk_) :: ipnt, ns_in, ns_out, n_rem, ictxt, me, np, isw
    integer(psb_lpk_), allocatable    :: tidx(:)
    integer(psb_ipk_), allocatable    :: tsmpl(:)

    ictxt   = idxmap%get_ctxt()
    call psb_info(ictxt,me,np)
    call psb_realloc(n_samples,tidx,info)
    call psb_realloc(n_samples,tsmpl,info)
    ipnt = 1
    isw  = 1
    do
      !write(0,*) me,' Into  sampling ',n_samples
      call psi_get_sample(ipnt, idx,iprc,tidx,tsmpl,n_samples,ns_out)      
      ns_in = min(n_samples,ns_out)
      !write(0,*) me,' From second sampling ',ns_out
      call psi_adjcncy_fnd_owner(tidx(1:ns_in),tprc,ladj,idxmap,info)
      call psi_cpy_out(iprc,tprc,tsmpl,ns_in,ns_out)
      !write(0,*) me,' Sweep ',isw,' answers:',ns_out
      n_answers = n_answers + ns_out
      n_rem = size(idx)-ipnt
      call psb_max(ictxt,n_rem)
      !write(0,*) me,' Sweep ',isw,n_rem, ipnt, n_samples
      if (n_rem <= 0) exit
      isw = isw + 1 
    end do


  end subroutine psi_adj_fnd_sweep

end subroutine psi_graph_fnd_owner
