!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
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
! File: psi_bld_tmphalo.f90
!
! Subroutine: psi_bld_tmphalo
!   Build initial versions of data exchange lists.
!   When the descriptor is for a large index space, we cannot build 
!   the data exchange lists "on-the-fly", but we also want to keep using the 
!   same format conversion routines we use in the small index space case, 
!   hence this adapter routine.
!   
! 
! Arguments: 
!    desc     - type(psb_desc_type).    The communication descriptor.        
!    info     - integer.                  return code.
!
subroutine psi_bld_tmphalo(desc,info)
  use psb_desc_mod
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_realloc_mod
  use psi_mod, psb_protect_name => psi_bld_tmphalo
  implicit none
  type(psb_desc_type), intent(inout) :: desc
  integer(psb_ipk_), intent(out) :: info

  integer(psb_ipk_),allocatable :: helem(:),hproc(:)
  integer(psb_ipk_),allocatable :: tmphl(:)

  integer(psb_ipk_) ::  i,j,np,me,lhalo,nhalo,&
       & n_col, err_act,  key, ih, nh, idx, nk,icomm
  integer(psb_ipk_) :: ictxt,n_row
  character(len=20)   :: name,ch_err

  info = psb_success_
  name = 'psi_bld_tmphalo'
  call psb_erractionsave(err_act)

  ictxt = desc%get_context()
  icomm = desc%get_mpic()
  n_row = desc%get_local_rows()
  n_col = desc%get_local_cols()

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.(desc%is_bld())) then 
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  ! Here we do not know yet who owns what, so we have 
  ! to call fnd_owner.
  nh = (n_col-n_row)
  Allocate(helem(max(1,nh)),stat=info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if

  do i=1, nh
    helem(i) = n_row+i ! desc%loc_to_glob(n_row+i)
  end do

  call desc%indxmap%l2gip(helem(1:nh),info)
  call desc%indxmap%fnd_owner(helem(1:nh),hproc,info)
      
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='fnd_owner')
    goto 9999      
  endif
  if (nh > size(hproc)) then 
    info=psb_err_from_subroutine_
    call psb_errpush(psb_err_from_subroutine_,name,a_err='nh > size(hproc)')
    goto 9999      
  end if

  allocate(tmphl((3*((n_col-n_row)+1)+1)),stat=info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if
  j  = 1
  do i=1,nh
    tmphl(j+0) = hproc(i)
    if (tmphl(j+0)<0) then 
      write(psb_err_unit,*) me,'Unrecoverable error: missing proc from asb',&
           & i, nh, n_row+i,helem(i),hproc(i)      
      info = psb_err_invalid_cd_state_
      call psb_errpush(info,name)
      goto 9999
    end if
    tmphl(j+1) = 1
    tmphl(j+2) = n_row+i
    j          = j + 3
  end do
  tmphl(j) = -1
  lhalo = j
  nhalo = (lhalo-1)/3

  call psb_move_alloc(tmphl,desc%halo_index,info)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error(ictxt)
  end if
  return


end subroutine psi_bld_tmphalo
