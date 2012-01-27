!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
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
Subroutine psb_cd_lstext(desc_a,in_list,desc_ov,info, mask,extype)

  use psb_base_mod, psb_protect_name => psb_cd_lstext
!!$  use psi_mod

  Implicit None

  !     .. Array Arguments ..
  Type(psb_desc_type), Intent(in), target :: desc_a
  integer(psb_ipk_), intent(in)                     :: in_list(:)
  Type(psb_desc_type), Intent(out)        :: desc_ov
  integer(psb_ipk_), intent(out)                    :: info
  logical, intent(in), optional, target   :: mask(:)
  integer(psb_ipk_), intent(in),optional            :: extype

  !     .. Local Scalars ..
  integer(psb_ipk_) ::  i, j, np, me,m,nnzero,&
       &  ictxt, lovr, lworks,lworkr, n_row,n_col, int_err(5),&
       &  index_dim,elem_dim, l_tmp_ovr_idx,l_tmp_halo, nztot,nhalo
  integer(psb_ipk_) :: counter,counter_h, counter_o, counter_e,idx,gidx,proc,n_elem_recv,&
       & n_elem_send,tot_recv,tot_elem,cntov_o,&
       & counter_t,n_elem,i_ovr,jj,proc_id,isz, nl, &
       & idxr, idxs, iszr, iszs, nxch, nsnd, nrcv,lidx, extype_
  integer(psb_ipk_) :: icomm, err_act

  integer(psb_ipk_), allocatable  :: tmp_halo(:),tmp_ovr_idx(:), orig_ovr(:)
  integer(psb_ipk_),allocatable   :: halo(:),works(:),workr(:),t_halo_in(:),&
       & t_halo_out(:),temp(:),maskr(:)
  integer(psb_ipk_),allocatable  :: brvindx(:),rvsz(:), bsdindx(:),sdsz(:)
  logical, allocatable, target :: lmask(:)
  logical, pointer     :: mask_(:)
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)    :: name, ch_err

  name='psb_cd_lstext'
  info  = psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_a%get_context()
  icomm = desc_a%get_mpic()
  Call psb_info(ictxt, me, np)

  If (debug_level >= psb_debug_outer_) &
       & Write(debug_unit,*) me,' ',trim(name),': start',size(in_list)



  m      = desc_a%get_local_rows()
  n_row  = desc_a%get_local_rows()
  n_col  = desc_a%get_local_cols()
  nhalo  = n_col-n_row

  nl = size(in_list)

  if (present(mask)) then 
    if (size(mask) < nl) then 
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='size of mask')
      goto 9999
    end if
    mask_ => mask
  else
    allocate(lmask(nl),stat=info) 
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='Allocat lmask')
      goto 9999
    end if
    lmask = .true.
    mask_ => lmask
  end if

  if (present(extype)) then
    extype_ = extype
  else
    extype_ = psb_ovt_xhal_  
  endif


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),':Calling desccpy'
  call psb_cdcpy(desc_a,desc_ov,info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_cdcpy'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),':From desccpy'


  call psb_cd_reinit(desc_ov,info)

  if (info == psb_success_) call psb_cdins(nl,in_list,desc_ov,info,mask=mask_)

  ! At this point we have added to the halo the indices in 
  ! in_list. Just call icdasb forcing to use 
  ! the halo_index provided. This is the same routine as gets 
  ! called inside CDASB.
  !

  if (debug_level >= psb_debug_outer_) then
    write(debug_unit,*) me,' ',trim(name),': converting indexes'
    call psb_barrier(ictxt)
  end if

  call psb_icdasb(desc_ov,info,ext_hv=.true.)

  call psb_cd_set_ovl_asb(desc_ov,info)

  if (info /= psb_success_) then
    ch_err='sp_free'
    call psb_errpush(psb_err_from_subroutine_ai_,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
    goto 9999
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),': end'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  Return

End Subroutine psb_cd_lstext
