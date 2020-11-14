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
! File:  psb_zovrl.f90
!
! Subroutine: psb_zovrl_vect
!   This subroutine performs the exchange of the overlap elements in a 
!    distributed dense vector between all the processes.
!
! Arguments:
!   x           -  type(psb_z_vect_type)  The local part of the vector
!   desc_a      -  type(psb_desc_type).    The communication descriptor.
!   info        -  integer.                  Return code.
!   jx          -  integer(optional).        The starting column of the global matrix
!   ik          -  integer(optional).        The number of columns to gather. 
!   work        -  complex(optional).           A work area.
!   update      -  integer(optional).        Type of update:
!                                            psb_none_   do nothing
!                                            psb_sum_    sum of overlaps
!                                            psb_avg_    average of overlaps
!   mode        -  integer(optional).        Choose the algorithm for data exchange: 
!                                       this is chosen through bit fields. 
!                                       - swap_mpi  = iand(flag,psb_swap_mpi_)  /= 0
!                                       - swap_sync = iand(flag,psb_swap_sync_) /= 0
!                                       - swap_send = iand(flag,psb_swap_send_) /= 0
!                                       - swap_recv = iand(flag,psb_swap_recv_) /= 0
!                                       - if (swap_mpi):  use underlying MPI_ALLTOALLV.
!                                       - if (swap_sync): use PSB_SND and PSB_RCV in 
!                                                       synchronized pairs
!                                       - if (swap_send .and. swap_recv): use mpi_irecv 
!                                                       and mpi_send
!                                       - if (swap_send): use psb_snd (but need another 
!                                                       call with swap_recv to complete)
!                                       - if (swap_recv): use psb_rcv (completing a 
!                                                       previous call with swap_send)
!
subroutine  psb_zovrl_vect(x,desc_a,info,work,update,mode)
  use psb_base_mod, psb_protect_name => psb_zovrl_vect
  use psi_mod
  implicit none

  type(psb_z_vect_type), intent(inout)   :: x
  type(psb_desc_type), intent(in)        :: desc_a
  integer(psb_ipk_), intent(out)                   :: info
  complex(psb_dpk_), optional, target, intent(inout) :: work(:)
  integer(psb_ipk_), intent(in), optional          :: update,mode

  ! locals
  type(psb_ctxt_type) :: ictxt
  integer(psb_ipk_) :: np, me, err_act, k, iix, jjx, &
       & nrow, ncol, ldx, liwork, data_, update_, mode_
  integer(psb_lpk_) :: m, n, ix, ijx
  complex(psb_dpk_),pointer :: iwork(:)
  logical                  :: do_swap
  character(len=20)        :: name, ch_err
  logical                  :: aliw

  name='psb_zovrlv'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ictxt=desc_a%get_context()

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(x%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = 1

  m = desc_a%get_global_rows()
  n = desc_a%get_global_cols()
  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  ldx = x%get_nrows()
  k = 1

  if (present(update)) then 
    update_ = update
  else
    update_ = psb_avg_
  endif

  if (present(mode)) then 
    mode_ = mode
  else
    mode_ = IOR(psb_swap_send_,psb_swap_recv_)
  endif
  do_swap = (mode_ /= 0)

  ! check vector correctness

  if ((info == 0).and.(ldx<ncol)) call x%reall(ncol,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_;     ch_err='psb_reall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  ! check for presence/size of a work area
  liwork=ncol
  if (present(work)) then
    if(size(work) >= liwork) then
      aliw=.false.
    else
      aliw=.true.
    end if
  else
    aliw=.true.
  end if
  if (aliw) then 
    allocate(iwork(liwork),stat=info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999
    end if
  else
    iwork => work    
  end if

  ! exchange overlap elements
  if (do_swap) then
    call psi_swapdata(mode_,zone,x%v,&
         & desc_a,iwork,info,data=psb_comm_ovr_)
  end if
  if (info == psb_success_) call psi_ovrl_upd(x%v,desc_a,update_,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Inner updates')
    goto 9999
  end if
  
  if (aliw) deallocate(iwork)
  nullify(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

    return
end subroutine psb_zovrl_vect

!
! Subroutine: psb_zovrl_multivect
!   This subroutine performs the exchange of the overlap elements in a 
!    distributed multivector between all the processes.
!
! Arguments:
!   x           -  type(psb_z_vect_type)  The local part of the multivector
!   desc_a      -  type(psb_desc_type).    The communication descriptor.
!   info        -  integer.                  Return code.
!   jx          -  integer(optional).        The starting column of the global matrix
!   ik          -  integer(optional).        The number of columns to gather. 
!   work        -  complex(optional).           A work area.
!   update      -  integer(optional).        Type of update:
!                                            psb_none_   do nothing
!                                            psb_sum_    sum of overlaps
!                                            psb_avg_    average of overlaps
!   mode        -  integer(optional).        Choose the algorithm for data exchange: 
!                                       this is chosen through bit fields. 
!                                       - swap_mpi  = iand(flag,psb_swap_mpi_)  /= 0
!                                       - swap_sync = iand(flag,psb_swap_sync_) /= 0
!                                       - swap_send = iand(flag,psb_swap_send_) /= 0
!                                       - swap_recv = iand(flag,psb_swap_recv_) /= 0
!                                       - if (swap_mpi):  use underlying MPI_ALLTOALLV.
!                                       - if (swap_sync): use PSB_SND and PSB_RCV in 
!                                                       synchronized pairs
!                                       - if (swap_send .and. swap_recv): use mpi_irecv 
!                                                       and mpi_send
!                                       - if (swap_send): use psb_snd (but need another 
!                                                       call with swap_recv to complete)
!                                       - if (swap_recv): use psb_rcv (completing a 
!                                                       previous call with swap_send)
!
subroutine  psb_zovrl_multivect(x,desc_a,info,work,update,mode)
  use psb_base_mod, psb_protect_name => psb_zovrl_multivect
  use psi_mod
  implicit none

  type(psb_z_multivect_type), intent(inout)   :: x
  type(psb_desc_type), intent(in)        :: desc_a
  integer(psb_ipk_), intent(out)                   :: info
  complex(psb_dpk_), optional, target, intent(inout) :: work(:)
  integer(psb_ipk_), intent(in), optional          :: update,mode

  ! locals
  type(psb_ctxt_type) :: ictxt
  integer(psb_ipk_) :: np, me, err_act, k, iix, jjx, &
       & nrow, ncol, ldx, liwork, data_, update_, mode_
  integer(psb_lpk_) :: m, n, ix, ijx
  complex(psb_dpk_),pointer :: iwork(:)
  logical                  :: do_swap
  character(len=20)        :: name, ch_err
  logical                  :: aliw

  name='psb_zovrlv'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ictxt=desc_a%get_context()

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(x%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = 1
  ijx = 1

  m = desc_a%get_global_rows()
  n = desc_a%get_global_cols()
  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  ldx = x%get_nrows()

  k = 1

  if (present(update)) then 
    update_ = update
  else
    update_ = psb_avg_
  endif

  if (present(mode)) then 
    mode_ = mode
  else
    mode_ = IOR(psb_swap_send_,psb_swap_recv_)
  endif
  do_swap = (mode_ /= 0)

  ! check vector correctness
  if (ldx < ncol) call x%reall(ncol,x%get_ncols(),info)
  
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_ ;    ch_err='psb_reall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  ! check for presence/size of a work area
  liwork=ncol
  if (present(work)) then
    if(size(work) >= liwork) then
      aliw=.false.
    else
      aliw=.true.
    end if
  else
    aliw=.true.
  end if
  if (aliw) then 
    allocate(iwork(liwork),stat=info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999
    end if
  else
    iwork => work    
  end if

  ! exchange overlap elements
  if (do_swap) then
    call psi_swapdata(mode_,zone,x%v,&
         & desc_a,iwork,info,data=psb_comm_ovr_)
  end if
  if (info == psb_success_) call psi_ovrl_upd(x%v,desc_a,update_,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Inner updates')
    goto 9999
  end if
  
  if (aliw) deallocate(iwork)
  nullify(iwork)
  
  call psb_erractionrestore(err_act)
  return  
  
9999 call psb_error_handler(ictxt,err_act)
  
  return
end subroutine psb_zovrl_multivect

