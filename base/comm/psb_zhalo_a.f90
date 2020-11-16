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
! Subroutine: psb_zhalom
!   This subroutine performs the exchange of the halo elements in a 
!    distributed dense matrix between all the processes.
!
! Arguments:
!   x         -  complex,dimension(:,:).          The local part of the dense matrix.
!   desc_a    -  type(psb_desc_type).        The communication descriptor.
!   info      -  integer.                      Return code
!   jx        -  integer(optional).            The starting column of the global matrix. 
!   ik        -  integer(optional).            The number of columns to gather. 
!   work      -  complex(optional).            Work  area.
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
subroutine  psb_zhalom(x,desc_a,info,jx,ik,work,tran,mode,data)
  use psb_base_mod, psb_protect_name => psb_zhalom
  use psi_mod
  implicit none

  complex(psb_dpk_), intent(inout), target   :: x(:,:)
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)                      :: info
  complex(psb_dpk_), optional, target, intent(inout) :: work(:)
  integer(psb_ipk_), intent(in), optional             :: mode,jx,ik,data
  character, intent(in), optional           :: tran

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_mpk_) :: np, me
  integer(psb_ipk_) :: err_act, iix, jjx, k, maxk, nrow, imode, i,&
       & liwork,data_, ldx
  integer(psb_lpk_) :: m, n, ix, ijx
  complex(psb_dpk_),pointer :: iwork(:), xp(:,:)
  character                :: tran_
  character(len=20)        :: name, ch_err
  logical                  :: aliw

  name='psb_zhalom'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ctxt=desc_a%get_context()

  ! check on blacs grid 
  call psb_info(ctxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = 1
  if (present(jx)) then
    ijx = jx
  else
    ijx = 1
  endif

  m = desc_a%get_global_rows()
  n = desc_a%get_global_cols()
  nrow = desc_a%get_local_rows()

  maxk=size(x,2)-ijx+1

  if(present(ik)) then
    if(ik > maxk) then
      k=maxk
    else
      k=ik
    end if
  else
    k = maxk
  end if

  if (present(tran)) then     
    tran_ = psb_toupper(tran)
  else
    tran_ = 'N'
  endif
  if (present(mode)) then 
    imode = mode
  else
    imode = IOR(psb_swap_send_,psb_swap_recv_)
  endif

  if (present(data)) then     
    data_ = data
  else
    data_ = psb_comm_halo_
  endif
  ldx = size(x,1)
  ! check vector correctness
  call psb_chkvect(m,lone,ldx,ix,ijx,desc_a,info,iix,jjx,check_halo=.true.)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_ ;    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  liwork=nrow
  if (present(work)) then
    if(size(work) >= liwork) then
      aliw=.false.
      iwork => work
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
  xp => x(iix:size(x,1),jjx:jjx+k-1)
  if(tran_ == 'N') then
    call psi_swapdata(imode,k,zzero,xp,&
         & desc_a,iwork,info,data=data_)
  else if((tran_ == 'T').or.(tran_ == 'C')) then
    call psi_swaptran(imode,k,zone,xp,&
         &desc_a,iwork,info)
  else
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='invalid tran')
    goto 9999      
  end if

  if(info /= psb_success_) then
    ch_err='PSI_cswapdata'
    call psb_errpush(psb_err_from_subroutine_,name,a_err=ch_err)
    goto 9999
  end if

  if (aliw) deallocate(iwork)
  nullify(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ctxt,err_act)

    return
end subroutine psb_zhalom




!!$ 
!!$              Parallel Sparse BLAS  version 3.5
!!$    (C) Copyright 2006-2018
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari      
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
! Subroutine: psb_zhalov
!   This subroutine performs the exchange of the halo elements in a 
!    distributed dense vector between all the processes.
!
! Arguments:
!   x         -  real,dimension(:).            The local part of the dense vector.
!   desc_a    -  type(psb_desc_type).        The communication descriptor.
!   info      -  integer.                      Return code
!   jx        -  integer(optional).            The starting column of the global matrix. 
!   ik        -  integer(optional).            The number of columns to gather. 
!   work      -  complex(optional).            Work  area.
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
subroutine  psb_zhalov(x,desc_a,info,work,tran,mode,data)
  use psb_base_mod, psb_protect_name => psb_zhalov
  use psi_mod
  implicit none

  complex(psb_dpk_), intent(inout)          :: x(:)
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)                      :: info
  complex(psb_dpk_), target, optional, intent(inout) :: work(:)
  integer(psb_ipk_), intent(in), optional             :: mode,data
  character, intent(in), optional           :: tran

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_mpk_) :: np, me
  integer(psb_ipk_) :: err_act, ldx, iix, jjx, nrow, imode, liwork,data_
  integer(psb_lpk_) :: m, n, ix, ijx
  complex(psb_dpk_),pointer :: iwork(:)
  character                :: tran_
  character(len=20)        :: name, ch_err
  logical                  :: aliw

  name='psb_zhalov'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ctxt=desc_a%get_context()

  ! check on blacs grid 
  call psb_info(ctxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = 1
  ijx = 1

  m = desc_a%get_global_rows()
  n = desc_a%get_global_cols()
  nrow = desc_a%get_local_rows()

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
  ldx = size(x,1)
  ! check vector correctness
  call psb_chkvect(m,lone,ldx,ix,ijx,desc_a,info,iix,jjx,check_halo=.true.)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_ ;    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  liwork=nrow
  if (present(work)) then
    if(size(work) >= liwork) then
      aliw=.false.
      iwork => work
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
    call psi_swapdata(imode,zzero,x(iix:size(x)),&
         & desc_a,iwork,info,data=data_)
  else if((tran_ == 'T').or.(tran_ == 'C')) then
    call psi_swaptran(imode,zone,x(iix:size(x)),&
         & desc_a,iwork,info)
  else
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='invalid tran')
    goto 9999      
  end if

  if(info /= psb_success_) then
    ch_err='PSI_swapdata'
    call psb_errpush(psb_err_from_subroutine_,name,a_err=ch_err)
    goto 9999
  end if

  if (aliw) deallocate(iwork)
  nullify(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ctxt,err_act)

    return
end subroutine psb_zhalov

