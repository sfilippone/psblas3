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

Subroutine psb_cd_reinit(desc,info)

  use psb_base_mod, psb_protect_name => psb_cd_reinit
  use psi_mod

  Implicit None

  !     .. Array Arguments ..
  Type(psb_desc_type), Intent(inout) :: desc
  integer(psb_ipk_), intent(out)               :: info


  !     .. Local Scalars ..
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   ::  np, me, err_act
  integer(psb_mpk_)   :: icomm
  integer(psb_ipk_), allocatable :: tmp_halo(:),tmp_ext(:), tmp_ovr(:)
  integer(psb_ipk_)   :: debug_level, debug_unit
  character(len=20)   :: name, ch_err

  name='psb_cd_reinit'
  info  = psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt = desc%get_context()
  icomm = desc%get_mpic()
  Call psb_info(ctxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),': start'
  if (desc%is_asb()) then 
    call psb_cd_get_recv_idx(tmp_ovr,desc,psb_comm_ovr_,info)
    call psb_cd_get_recv_idx(tmp_halo,desc,psb_comm_halo_,info)    
    call psb_cd_get_recv_idx(tmp_ext,desc,psb_comm_ext_,info)        
    
    call psb_move_alloc(tmp_ovr,desc%ovrlap_index,info)
    call psb_move_alloc(tmp_halo,desc%halo_index,info)
    call psb_move_alloc(tmp_ext,desc%ext_index,info)
    call desc%indxmap%reinit(info)
!!$    if (me == 0) write(0,*) 'On cdreinit status :',&
!!$         & allocated(desc%indxmap%p_adjcncy),allocated(desc%indxmap%halo_owner), &
!!$         & desc%get_fmt()
    !  call psb_cd_set_bld(desc,info)
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),': end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

End Subroutine psb_cd_reinit
