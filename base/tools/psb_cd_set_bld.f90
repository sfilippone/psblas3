!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
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
subroutine psb_cd_set_ovl_bld(desc,info)
  !
  ! Change state of a descriptor into ovl_build. 
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psi_mod
  use psb_descriptor_type
  use psb_tools_mod, psb_protect_name => psb_cd_set_ovl_bld
  implicit none
  type(psb_desc_type), intent(inout) :: desc
  integer                            :: info

  call psb_cd_set_bld(desc,info) 
  if (info == 0) desc%matrix_data(psb_dec_type_) = psb_cd_ovl_bld_ 

end subroutine psb_cd_set_ovl_bld

subroutine psb_cd_set_bld(desc,info)
  !
  ! Change state of a descriptor into BUILD. 
  ! If the descriptor is LARGE, check the  AVL search tree
  ! and initialize it if necessary.
  !
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psi_mod
  use psb_descriptor_type
  use psb_tools_mod, psb_protect_name => psb_cd_set_bld
  implicit none
  type(psb_desc_type), intent(inout) :: desc
  integer                            :: info
  !locals
  integer             :: np,me,ictxt, err_act,idx,gidx,lidx,nc
  logical, parameter  :: debug=.false.,debugprt=.false.
  character(len=20)   :: name
  if (debug) write(0,*) me,'Entered CDCPY'
  if (psb_get_errstatus() /= 0) return 
  info = 0
  call psb_erractionsave(err_act)
  name = 'psb_cd_set_bld'

  ictxt = psb_cd_get_context(desc)

  if (debug) write(0,*)'Entered CDSETBLD',ictxt
  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (debug) write(0,*) me,'Entered CDSETBLD'
  if (psb_is_asb_desc(desc)) then 
  end if

  desc%matrix_data(psb_dec_type_) = psb_desc_bld_ 
  
  if (psb_is_large_desc(desc)) then 
    !
    ! The idea: first build glb_lc with the info on
    ! rows we already have, then leave space in
    ! hash for newcomers (halo indices).
    ! The policy is to allocate for as many entries
    ! as there are rows; if we ever fill them up, we can
    ! try and enlarge again, but by the time the hash
    ! fills up it means we have as many halo as internals,
    ! therefore there are much worse problems ahead than
    ! the hash occupancy.
    !
    nc = psb_cd_get_local_cols(desc)
    if (.not.associated(desc%hash)) allocate(desc%hash,stat=info)
    if (info == 0)&
         & call psb_hash_init(nc,desc%hash,info)
    if (info == 0) call psi_bld_hash(desc,info)
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='hashInit')
      goto 9999      
    end if

  end if

  if (debug) write(0,*) me,'SET_BLD: done'
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
end subroutine psb_cd_set_bld
