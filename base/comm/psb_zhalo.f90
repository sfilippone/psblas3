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
! File:  psb_zhalo.f90
!
! Subroutine: psb_zhalo_vect
!   This subroutine performs the exchange of the halo elements in a 
!    distributed vector between all the processes.
!
! Arguments:
!   x         -  type(psb_z_vect_type)    The local part of the vector
!   desc_a    -  type(psb_desc_type).        The communication descriptor.
!   info      -  integer.                      Return code
!   jx        -  integer(optional).            The starting column of the global matrix. 
!   ik        -  integer(optional).            The number of columns to gather. 
!   work      -  complex(optional).             Work  area.
!   tran      -  character(optional).          Transpose exchange.
!   mode      -  integer(optional).            Communication mode (see Swapdata)
!   data     - integer                 Which index list in desc_a should be used
!                                      to retrieve rows, default psb_comm_halo_
!                                       psb_comm_halo_    use halo_index
!                                       psb_comm_ext_     use ext_index 
!                                       psb_comm_ovrl_    use ovrl_index
!                                       psb_comm_mov_     use ovr_mst_idx
!
!
subroutine  psb_zhalo_vect(x,desc_a,info,work,tran,mode,data)
  use psb_base_mod, psb_protect_name => psb_zhalo_vect
  use psi_mod
  implicit none

  type(psb_z_vect_type), intent(inout)    :: x
  type(psb_desc_type), intent(in)         :: desc_a
  integer(psb_ipk_), intent(out)                    :: info
  complex(psb_dpk_), target, optional, intent(inout)  :: work(:)
  integer(psb_ipk_), intent(in), optional           :: mode,data
  character, intent(in), optional         :: tran

  ! locals
  type(psb_ctxt_type) :: ictxt
  integer(psb_ipk_) :: np, me, err_act, iix, jjx, &
       & nrow, ncol, lldx, imode, liwork,data_
  integer(psb_lpk_) :: m, n, ix, ijx
  complex(psb_dpk_),pointer :: iwork(:)
  character                 :: tran_
  character(len=20)         :: name, ch_err
  logical                   :: aliw

  name='psb_zhalov'
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
  lldx = x%get_nrows()

  if (present(tran)) then     
    tran_ = psb_toupper(tran)
  else
    tran_ = 'N'
  endif
  if (present(data)) then     
    data_ = data
  else
    data_ = psb_comm_halo_
  endif
  if (present(mode)) then 
    imode = mode
  else
    imode = IOR(psb_swap_send_,psb_swap_recv_)
  endif

  if ((info == 0).and.(lldx<ncol)) call x%reall(ncol,info)

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_ ;    ch_err='reall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  liwork=nrow
  if (present(work)) then
    if(size(work) >= liwork) then
      iwork => work
      aliw=.false.
    else
      aliw=.true.
      allocate(iwork(liwork),stat=info)
      if(info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='psb_realloc'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    end if
  else
    aliw=.true.
    allocate(iwork(liwork),stat=info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  end if

  ! exchange halo elements
  if(tran_ == 'N') then
    call psi_swapdata(imode,zzero,x%v,&
         & desc_a,iwork,info,data=data_)
  else if((tran_ == 'T').or.(tran_ == 'C')) then
    call psi_swaptran(imode,zone,x%v,&
         & desc_a,iwork,info)
  else
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='invalid tran')
    goto 9999      
  end if

  if (info /= psb_success_) then
    ch_err='PSI_swapdata'
    call psb_errpush(psb_err_from_subroutine_,name,a_err=ch_err)
    goto 9999
  end if

  if (aliw) deallocate(iwork)
  nullify(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

    return
end subroutine psb_zhalo_vect

!
! Subroutine: psb_zhalo_multivect
!   This subroutine performs the exchange of the halo elements in a 
!    distributed multivector between all the processes.
!
! Arguments:
!   x         -  type(psb_z_multivect_type)    The local part of the multivector
!   desc_a    -  type(psb_desc_type).        The communication descriptor.
!   info      -  integer.                      Return code
!   jx        -  integer(optional).            The starting column of the global matrix. 
!   ik        -  integer(optional).            The number of columns to gather. 
!   work      -  complex(optional).             Work  area.
!   tran      -  character(optional).          Transpose exchange.
!   mode      -  integer(optional).            Communication mode (see Swapdata)
!   data     - integer                 Which index list in desc_a should be used
!                                      to retrieve rows, default psb_comm_halo_
!                                       psb_comm_halo_    use halo_index
!                                       psb_comm_ext_     use ext_index 
!                                       psb_comm_ovrl_    use ovrl_index
!                                       psb_comm_mov_     use ovr_mst_idx
!
!
subroutine  psb_zhalo_multivect(x,desc_a,info,work,tran,mode,data)
  use psb_base_mod, psb_protect_name => psb_zhalo_multivect
  use psi_mod
  implicit none

  type(psb_z_multivect_type), intent(inout)    :: x
  type(psb_desc_type), intent(in)         :: desc_a
  integer(psb_ipk_), intent(out)                    :: info
  complex(psb_dpk_), target, optional, intent(inout)  :: work(:)
  integer(psb_ipk_), intent(in), optional           :: mode,data
  character, intent(in), optional         :: tran

  ! locals
  type(psb_ctxt_type) :: ictxt
  integer(psb_ipk_) :: np, me, err_act, iix, jjx, &
       & nrow, ncol, lldx, imode, liwork,data_
  integer(psb_lpk_) :: m, n, ix, ijx
  complex(psb_dpk_),pointer :: iwork(:)
  character                 :: tran_
  character(len=20)         :: name, ch_err
  logical                   :: aliw

  name='psb_zhalov'
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
  lldx = x%get_nrows()

  if (present(tran)) then     
    tran_ = psb_toupper(tran)
  else
    tran_ = 'N'
  endif
  if (present(data)) then     
    data_ = data
  else
    data_ = psb_comm_halo_
  endif
  if (present(mode)) then 
    imode = mode
  else
    imode = IOR(psb_swap_send_,psb_swap_recv_)
  endif

  if (lldx < ncol) call x%reall(ncol,x%get_ncols(),info)
  
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_;    ch_err='psb_reall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  liwork=nrow
  if (present(work)) then
    if(size(work) >= liwork) then
      iwork => work
      aliw=.false.
    else
      aliw=.true.
      allocate(iwork(liwork),stat=info)
      if(info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='psb_realloc'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    end if
  else
    aliw=.true.
    allocate(iwork(liwork),stat=info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  end if

  ! exchange halo elements
  if(tran_ == 'N') then
    call psi_swapdata(imode,zzero,x%v,&
         & desc_a,iwork,info,data=data_)
  else if((tran_ == 'T').or.(tran_ == 'C')) then
    call psi_swaptran(imode,zone,x%v,&
         & desc_a,iwork,info)
  else
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='invalid tran')
    goto 9999      
  end if

  if (info /= psb_success_) then
    ch_err='PSI_swapdata'
    call psb_errpush(psb_err_from_subroutine_,name,a_err=ch_err)
    goto 9999
  end if

  if (aliw) deallocate(iwork)
  nullify(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psb_zhalo_multivect

