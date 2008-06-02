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
function psb_c_inter_desc(map_kind,desc1, desc2, map_fw, map_bk, idx_fw, idx_bk)

  use psb_base_mod, psb_protect_name => psb_c_inter_desc

  implicit none 
  type(psb_inter_desc_type)         :: psb_c_inter_desc    
  type(psb_desc_type), target       :: desc1, desc2
  type(psb_cspmat_type), intent(in) :: map_fw, map_bk
  integer, intent(in)               :: map_kind,idx_fw(:), idx_bk(:)
  !
  type(psb_inter_desc_type)    :: this
  integer                      :: info
  character(len=20), parameter :: name='psb_inter_desc'

  info = 0 
  if (psb_is_ok_desc(desc1)) then 
    this%desc_1=>desc1
  else
    info = 2
  endif
  if (psb_is_ok_desc(desc2)) then 
    this%desc_2=>desc2
  else
    info = 3
  endif

  if (info == 0) call psb_sp_clone(map_fw,this%cmap%map_fw,info)
  if (info == 0) call psb_sp_clone(map_bk,this%cmap%map_bk,info)
  if (info == 0) call psb_safe_cpy(idx_fw,this%exch_fw_idx,info)
  if (info == 0) call psb_safe_cpy(idx_bk,this%exch_bk_idx,info)
  if (info == 0) call psb_realloc(psb_itd_data_size_,this%itd_data,info) 
  if (info == 0) then
    call psb_cd_set_map_kind(map_kind, this)
    call psb_cd_set_map_data(psb_map_complex_, this)
  end if
  if (info /= 0) then
    write(0,*) trim(name),' Invalid descriptor input'
    return
  end if

  psb_c_inter_desc = this

end function psb_c_inter_desc


function psb_c_inter_desc_noidx(map_kind,desc1, desc2, map_fw, map_bk)

  use psb_base_mod, psb_protect_name => psb_c_inter_desc_noidx

  implicit none 
  type(psb_inter_desc_type)         :: psb_c_inter_desc_noidx    
  type(psb_desc_type), target       :: desc1, desc2
  type(psb_cspmat_type), intent(in) :: map_fw, map_bk
  integer, intent(in)               :: map_kind
  !
  type(psb_inter_desc_type)    :: this
  integer                      :: info
  character(len=20), parameter :: name='psb_inter_desc'

  info = 0 
  select case(map_kind) 
  case (psb_map_aggr_)
    ! OK
  case default
    write(0,*) 'Bad map kind into psb_inter_desc ',map_kind
    info = 1
  end select

  if (psb_is_ok_desc(desc1)) then 
    this%desc_1=>desc1
  else
    info = 2
  endif
  if (psb_is_ok_desc(desc2)) then 
    this%desc_2=>desc2
  else
    info = 3
  endif

  if (info == 0) call psb_sp_clone(map_fw,this%cmap%map_fw,info)
  if (info == 0) call psb_sp_clone(map_bk,this%cmap%map_bk,info)
  if (info == 0) call psb_realloc(psb_itd_data_size_,this%itd_data,info) 
  if (info == 0) then
    call psb_cd_set_map_kind(map_kind, this)
    call psb_cd_set_map_data(psb_map_complex_, this)
  end if
  if (info /= 0) then
    write(0,*) trim(name),' Invalid descriptor input'
    return
  end if

  psb_c_inter_desc_noidx = this

end function psb_c_inter_desc_noidx

function psb_d_inter_desc(map_kind,desc1,desc2,map_fw,map_bk,idx_fw,idx_bk)

  use psb_base_mod, psb_protect_name => psb_d_inter_desc

  implicit none 
  type(psb_inter_desc_type)         :: psb_d_inter_desc    
  type(psb_desc_type), target       :: desc1, desc2
  type(psb_dspmat_type), intent(in) :: map_fw, map_bk
  integer, intent(in)               :: map_kind,idx_fw(:), idx_bk(:)
  !
  type(psb_inter_desc_type)    :: this
  integer                      :: info
  character(len=20), parameter :: name='psb_inter_desc'

  info = 0 
  if (psb_is_ok_desc(desc1)) then 
    this%desc_1=>desc1
  else
    info = 2
  endif
  if (psb_is_ok_desc(desc2)) then 
    this%desc_2=>desc2
  else
    info = 3
  endif

  if (info == 0) call psb_sp_clone(map_fw,this%dmap%map_fw,info)
  if (info == 0) call psb_sp_clone(map_bk,this%dmap%map_bk,info)
  if (info == 0) call psb_safe_cpy(idx_fw,this%exch_fw_idx,info)
  if (info == 0) call psb_safe_cpy(idx_bk,this%exch_bk_idx,info)
  if (info == 0) call psb_realloc(psb_itd_data_size_,this%itd_data,info) 
  if (info == 0) then
    call psb_cd_set_map_kind(map_kind, this)
    call psb_cd_set_map_data(psb_map_double_, this)
  end if
  if (info /= 0) then
    write(0,*) trim(name),' Invalid descriptor input'
    return
  end if

  psb_d_inter_desc = this

end function psb_d_inter_desc

function psb_d_inter_desc_noidx(map_kind,desc1, desc2, map_fw, map_bk)

  use psb_base_mod, psb_protect_name => psb_d_inter_desc_noidx

  implicit none 
  type(psb_inter_desc_type)         :: psb_d_inter_desc_noidx    
  type(psb_desc_type), target       :: desc1, desc2
  type(psb_dspmat_type), intent(in) :: map_fw, map_bk
  integer, intent(in)               :: map_kind
  !
  type(psb_inter_desc_type)    :: this
  integer                      :: info
  character(len=20), parameter :: name='psb_inter_desc'

  info = 0 
  select case(map_kind) 
  case (psb_map_aggr_)
    ! OK
  case default
    write(0,*) 'Bad map kind into psb_inter_desc ',map_kind
    info = 1
  end select

  if (psb_is_ok_desc(desc1)) then 
    this%desc_1=>desc1
  else
    info = 2
  endif
  if (psb_is_ok_desc(desc2)) then 
    this%desc_2=>desc2
  else
    info = 3
  endif

  if (info == 0) call psb_sp_clone(map_fw,this%dmap%map_fw,info)
  if (info == 0) call psb_sp_clone(map_bk,this%dmap%map_bk,info)
  if (info == 0) call psb_realloc(psb_itd_data_size_,this%itd_data,info) 
  if (info == 0) then
    call psb_cd_set_map_kind(map_kind, this)
    call psb_cd_set_map_data(psb_map_double_, this)
  end if
  if (info /= 0) then
    write(0,*) trim(name),' Invalid descriptor input'
    return
  end if

  psb_d_inter_desc_noidx = this

end function psb_d_inter_desc_noidx


function psb_s_inter_desc(map_kind,desc1,desc2,map_fw,map_bk,idx_fw,idx_bk)

  use psb_base_mod, psb_protect_name => psb_s_inter_desc

  implicit none 
  type(psb_inter_desc_type)         :: psb_s_inter_desc    
  type(psb_desc_type), target       :: desc1, desc2
  type(psb_sspmat_type), intent(in) :: map_fw, map_bk
  integer, intent(in)               :: map_kind,idx_fw(:), idx_bk(:)
  !
  type(psb_inter_desc_type)    :: this
  integer                      :: info
  character(len=20), parameter :: name='psb_inter_desc'

  info = 0 
  if (psb_is_ok_desc(desc1)) then 
    this%desc_1=>desc1
  else
    info = 2
  endif
  if (psb_is_ok_desc(desc2)) then 
    this%desc_2=>desc2
  else
    info = 3
  endif

  if (info == 0) call psb_sp_clone(map_fw,this%smap%map_fw,info)
  if (info == 0) call psb_sp_clone(map_bk,this%smap%map_bk,info)
  if (info == 0) call psb_safe_cpy(idx_fw,this%exch_fw_idx,info)
  if (info == 0) call psb_safe_cpy(idx_bk,this%exch_bk_idx,info)
  if (info == 0) call psb_realloc(psb_itd_data_size_,this%itd_data,info) 
  if (info == 0) then
    call psb_cd_set_map_kind(map_kind, this)
    call psb_cd_set_map_data(psb_map_single_, this)
  end if
  if (info /= 0) then
    write(0,*) trim(name),' Invalid descriptor input'
    return
  end if

  psb_s_inter_desc = this

end function psb_s_inter_desc

function psb_s_inter_desc_noidx(map_kind,desc1, desc2, map_fw, map_bk)

  use psb_base_mod, psb_protect_name => psb_s_inter_desc_noidx

  implicit none 
  type(psb_inter_desc_type)         :: psb_s_inter_desc_noidx    
  type(psb_desc_type), target       :: desc1, desc2
  type(psb_sspmat_type), intent(in) :: map_fw, map_bk
  integer, intent(in)               :: map_kind
  !
  type(psb_inter_desc_type)    :: this
  integer                      :: info
  character(len=20), parameter :: name='psb_inter_desc'

  info = 0 
  select case(map_kind) 
  case (psb_map_aggr_)
    ! OK
  case default
    write(0,*) 'Bad map kind into psb_inter_desc ',map_kind
    info = 1
  end select

  if (psb_is_ok_desc(desc1)) then 
    this%desc_1=>desc1
  else
    info = 2
  endif
  if (psb_is_ok_desc(desc2)) then 
    this%desc_2=>desc2
  else
    info = 3
  endif

  if (info == 0) call psb_sp_clone(map_fw,this%smap%map_fw,info)
  if (info == 0) call psb_sp_clone(map_bk,this%smap%map_bk,info)
  if (info == 0) call psb_realloc(psb_itd_data_size_,this%itd_data,info) 
  if (info == 0) then
    call psb_cd_set_map_kind(map_kind, this)
    call psb_cd_set_map_data(psb_map_single_, this)
  end if
  if (info /= 0) then
    write(0,*) trim(name),' Invalid descriptor input'
    return
  end if

  psb_s_inter_desc_noidx = this

end function psb_s_inter_desc_noidx

function psb_z_inter_desc(map_kind,desc1, desc2, map_fw, map_bk, idx_fw, idx_bk)

  use psb_base_mod, psb_protect_name => psb_z_inter_desc

  implicit none 
  type(psb_inter_desc_type)         :: psb_z_inter_desc    
  type(psb_desc_type), target       :: desc1, desc2
  type(psb_zspmat_type), intent(in) :: map_fw, map_bk
  integer, intent(in)               :: map_kind,idx_fw(:), idx_bk(:)
  !
  type(psb_inter_desc_type)    :: this
  integer                      :: info
  character(len=20), parameter :: name='psb_inter_desc'

  info = 0 
  if (psb_is_ok_desc(desc1)) then 
    this%desc_1=>desc1
  else
    info = 2
  endif
  if (psb_is_ok_desc(desc2)) then 
    this%desc_2=>desc2
  else
    info = 3
  endif

  if (info == 0) call psb_sp_clone(map_fw,this%zmap%map_fw,info)
  if (info == 0) call psb_sp_clone(map_bk,this%zmap%map_bk,info)
  if (info == 0) call psb_safe_cpy(idx_fw,this%exch_fw_idx,info)
  if (info == 0) call psb_safe_cpy(idx_bk,this%exch_bk_idx,info)
  if (info == 0) call psb_realloc(psb_itd_data_size_,this%itd_data,info) 
  if (info == 0) then
    call psb_cd_set_map_kind(map_kind, this)
    call psb_cd_set_map_data(psb_map_double_complex_, this)
  end if
  if (info /= 0) then
    write(0,*) trim(name),' Invalid descriptor input'
    return
  end if

  psb_z_inter_desc = this

end function psb_z_inter_desc

function psb_z_inter_desc_noidx(map_kind,desc1, desc2, map_fw, map_bk)

  use psb_base_mod, psb_protect_name => psb_z_inter_desc_noidx

  implicit none 
  type(psb_inter_desc_type)         :: psb_z_inter_desc_noidx    
  type(psb_desc_type), target       :: desc1, desc2
  type(psb_zspmat_type), intent(in) :: map_fw, map_bk
  integer, intent(in)               :: map_kind
  !
  type(psb_inter_desc_type)    :: this
  integer                      :: info
  character(len=20), parameter :: name='psb_inter_desc'

  info = 0 
  select case(map_kind) 
  case (psb_map_aggr_)
    ! OK
  case default
    write(0,*) 'Bad map kind into psb_inter_desc ',map_kind
    info = 1
  end select

  if (psb_is_ok_desc(desc1)) then 
    this%desc_1=>desc1
  else
    info = 2
  endif
  if (psb_is_ok_desc(desc2)) then 
    this%desc_2=>desc2
  else
    info = 3
  endif

  if (info == 0) call psb_sp_clone(map_fw,this%zmap%map_fw,info)
  if (info == 0) call psb_sp_clone(map_bk,this%zmap%map_bk,info)
  if (info == 0) call psb_realloc(psb_itd_data_size_,this%itd_data,info) 
  if (info == 0) then
    call psb_cd_set_map_kind(map_kind, this)
    call psb_cd_set_map_data(psb_map_double_complex_, this)
  end if
  if (info /= 0) then
    write(0,*) trim(name),' Invalid descriptor input'
    return
  end if

  psb_z_inter_desc_noidx = this

end function psb_z_inter_desc_noidx
