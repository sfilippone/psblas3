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

subroutine psi_i_renum_index(iperm,idx,info)
  use psi_mod, psi_protect_name =>  psi_i_renum_index
  use psb_serial_mod
  implicit none

  integer(psb_ipk_), intent(out)   :: info
  integer(psb_ipk_), intent(in)    :: iperm(:)
  integer(psb_ipk_), intent(inout) :: idx(:)

  integer(psb_ipk_) :: i,j,k,nh

  i=1
  k=idx(i)
  do while (k /= -1)
    i = i+1
    nh = idx(i)
    do j = i+1, i+nh
      idx(j) = iperm(idx(j))
    enddo
    i  = i + nh + 1
    nh = idx(i)
    do j = i+1, i+nh
      idx(j) = iperm(idx(j))
    enddo
    i = i + nh + 1
    k = idx(i)
  enddo

end subroutine psi_i_renum_index

subroutine psi_i_cnv_dsc(halo_in,ovrlap_in,ext_in,cdesc, info, mold)

  use psi_mod, psi_protect_name =>  psi_i_cnv_dsc
  use psb_timers_mod
  use psb_realloc_mod
  implicit none

  !     ....scalars parameters....
  integer(psb_ipk_), intent(in)      :: halo_in(:), ovrlap_in(:),ext_in(:)
  type(psb_desc_type), intent(inout) :: cdesc
  integer(psb_ipk_), intent(out)     :: info
  class(psb_i_base_vect_type), optional, intent(in) :: mold

  !     ....local scalars....
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np,me
  integer(psb_ipk_)   :: err_act,nxch,nsnd,nrcv,j,k
  !     ...local array...
  integer(psb_ipk_), allocatable  :: idx_out(:), tmp_mst_idx(:)

  !     ...parameters
  integer(psb_ipk_) :: debug_level, debug_unit
  logical, parameter :: debug=.false.
  character(len=20)  :: name
  logical, parameter  :: do_timings=.false.
  integer(psb_ipk_), save  :: idx_phase1=-1, idx_phase2=-1, idx_phase3=-1
  integer(psb_ipk_), save  :: idx_phase11=-1, idx_phase12=-1, idx_phase13=-1

  name='psi_cnv_desc'
  call psb_get_erraction(err_act)
  debug_level = psb_get_debug_level()
  debug_unit  = psb_get_debug_unit()

  info = psb_success_
  ctxt = cdesc%get_context()

  call psb_info(ctxt,me,np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif
  if ((do_timings).and.(idx_phase1==-1))       &
       & idx_phase1 = psb_get_timer_idx("PSI_CNV_DSC: phase1 ")
  if ((do_timings).and.(idx_phase2==-1))       &
       & idx_phase2 = psb_get_timer_idx("PSI_CNV_DSC: phase2")
  if ((do_timings).and.(idx_phase3==-1))       &
       & idx_phase3 = psb_get_timer_idx("PSI_CNV_DSC: phase3")
  if ((do_timings).and.(idx_phase11==-1))       &
       & idx_phase11 = psb_get_timer_idx("PSI_CNV_DSC: phase11 ")
  if ((do_timings).and.(idx_phase12==-1))       &
       & idx_phase12 = psb_get_timer_idx("PSI_CNV_DSC: phase12")
  if ((do_timings).and.(idx_phase13==-1))       &
       & idx_phase13 = psb_get_timer_idx("PSI_CNV_DSC: phase13")


  if (do_timings) call psb_tic(idx_phase1)
  if (do_timings) call psb_tic(idx_phase11)

  ! first the halo index
  if (debug_level>0) write(debug_unit,*) me,'Calling crea_index on halo',&
       & size(halo_in)
  call psi_crea_index(cdesc,halo_in, idx_out,nxch,nsnd,nrcv,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_crea_index')
    goto 9999
  end if
  call psb_move_alloc(idx_out,cdesc%halo_index,info)

  if (debug_level>0) write(debug_unit,*) me,'Done crea_index on halo'
  if (debug_level>0) write(debug_unit,*) me,'Calling crea_index on ext'
  if (do_timings) call psb_toc(idx_phase11)
  if (do_timings) call psb_tic(idx_phase12)


  ! then ext index
  if (debug_level>0) write(debug_unit,*) me,'Calling crea_index on ext'
  call psi_crea_index(cdesc,ext_in, idx_out,nxch,nsnd,nrcv,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_crea_index')
    goto 9999
  end if
  call psb_move_alloc(idx_out,cdesc%ext_index,info)

  if (debug_level>0) write(debug_unit,*) me,'Done crea_index on ext'
  if (debug_level>0) write(debug_unit,*) me,'Calling crea_index on ovrlap'
  if (do_timings) call psb_toc(idx_phase12)
  if (do_timings) call psb_tic(idx_phase13)

  ! then the overlap index
  call psi_crea_index(cdesc,ovrlap_in, idx_out,nxch,nsnd,nrcv,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_crea_index')
    goto 9999
  end if
  call psb_move_alloc(idx_out,cdesc%ovrlap_index,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_move_alloc')
    goto 9999
  end if
  if (do_timings) call psb_toc(idx_phase13)
  if (do_timings) call psb_toc(idx_phase1)
  if (do_timings) call psb_tic(idx_phase2)


  ! next  ovrlap_elem
  if (debug_level>0) write(debug_unit,*) me,'Calling crea_ovr_elem'
  call psi_crea_ovr_elem(me,cdesc%ovrlap_index,cdesc%ovrlap_elem,info)
  if (debug_level>0) write(debug_unit,*) me,'Done crea_ovr_elem'
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_crea_ovr_elem')
    goto 9999
  end if
  ! Extract ovr_mst_idx from ovrlap_elem
  if (debug_level>0) write(debug_unit,*) me,'Calling bld_ovr_mst'
  call psi_bld_ovr_mst(me,cdesc%ovrlap_elem,tmp_mst_idx,info)
  if (info == psb_success_) call psi_crea_index(cdesc,&
       & tmp_mst_idx,idx_out,nxch,nsnd,nrcv,info)
  if (debug_level>0) write(debug_unit,*) me,'Done crea_indx'
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_bld_ovr_mst')
    goto 9999
  end if
  call psb_move_alloc(idx_out,cdesc%ovr_mst_idx,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_move_alloc')
    goto 9999
  end if
  if (do_timings) call psb_toc(idx_phase2)
  if (do_timings) call psb_tic(idx_phase3)

  ! finally bnd_elem
  call psi_crea_bnd_elem(idx_out,cdesc,info)
  if (info == psb_success_) call psb_move_alloc(idx_out,cdesc%bnd_elem,info)

  call cdesc%v_halo_index%bld(cdesc%halo_index,mold=mold)
  call cdesc%v_ext_index%bld(cdesc%ext_index,mold=mold)
  call cdesc%v_ovrlap_index%bld(cdesc%ovrlap_index,mold=mold)
  call cdesc%v_ovr_mst_idx%bld(cdesc%ovr_mst_idx,mold=mold)


  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_crea_bnd_elem')
    goto 9999
  end if
  if (debug_level>0) write(debug_unit,*) me,'Done crea_bnd_elem'
  if (do_timings) call psb_toc(idx_phase3)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psi_i_cnv_dsc

subroutine psi_i_bld_ovr_mst(me,ovrlap_elem,mst_idx,info)
  use psi_mod, psi_protect_name =>  psi_i_bld_ovr_mst

  use psb_realloc_mod
  implicit none

  !     ....scalars parameters....
  integer(psb_ipk_), intent(in)               :: me, ovrlap_elem(:,:)
  integer(psb_ipk_), allocatable, intent(out) :: mst_idx(:)
  integer(psb_ipk_), intent(out)              :: info

  integer(psb_ipk_) :: i, j, proc, nov,isz, ip, err_act, idx
  character(len=20)  :: name

  name='psi_bld_ovr_mst'
  call psb_get_erraction(err_act)

  nov = size(ovrlap_elem,1)
  isz = 3*nov+1
  call psb_realloc(isz,mst_idx,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='reallocate')
    goto 9999
  end if
  mst_idx = -1
  j = 1
  do i=1, nov
    proc = ovrlap_elem(i,3)
    if (me /= proc) then
      idx = ovrlap_elem(i,1)
      mst_idx(j+0) = proc
      mst_idx(j+1) = 1
      mst_idx(j+2) = idx
      j = j + 3
    end if
  end do
  mst_idx(j) = -1

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psi_i_bld_ovr_mst
