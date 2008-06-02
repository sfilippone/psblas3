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
!!$
!
! Takes a vector X from space desc%desc_1 and maps it onto
! desc%desc_2 under desc%map_fw possibly with communication
! due to exch_fw_idx
!
subroutine psb_s_forward_map(alpha,x,beta,y,desc,info,work)
  use psb_base_mod, psb_protect_name => psb_s_forward_map
  implicit none 
  type(psb_inter_desc_type), intent(in) :: desc
  real(psb_spk_), intent(in)     :: alpha,beta
  real(psb_spk_), intent(inout)  :: x(:)
  real(psb_spk_), intent(out)    :: y(:)
  integer, intent(out)             :: info 
  real(psb_spk_), optional       :: work(:)

  !
  real(psb_spk_), allocatable :: xt(:)
  integer                       :: itsz, i, j,totxch,totsnd,totrcv,&
       &  map_kind, map_data
  character(len=20), parameter  :: name='psb_forward_map'

  info = 0
  if (.not.psb_is_asb_desc(desc)) then 
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end if

  itsz     = psb_cd_get_fw_tmp_sz(desc)
  map_kind = psb_cd_get_map_kind(desc)
  map_data = psb_cd_get_map_data(desc)
  if (map_data /= psb_map_single_) then 
    write(0,*) trim(name),' Invalid descriptor inupt: map_data', &
         & map_data,psb_map_single_
    info = 1
    return 
  endif

  select case(map_kind)
  case(psb_map_aggr_)
    ! Ok, we just need to call a halo update on the base desc
    ! and a matrix-vector product. 
    call psb_halo(x,desc%desc_1,info,work=work)
    if (info == 0) call psb_csmm(alpha,desc%smap%map_fw,x,beta,y,info)

    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if

  case(psb_map_gen_linear_)

    call psb_linmap(alpha,x,beta,y,desc%smap%map_fw,&
         & desc%desc_fw,desc%desc_1,desc%desc_2)

    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if


  case default
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end select

end subroutine psb_s_forward_map


!
! Takes a vector X from space desc%desc_2 and maps it onto
! desc%desc_1 under desc%map_bk possibly with communication
! due to exch_bk_idx
!
subroutine psb_s_backward_map(alpha,x,beta,y,desc,info,work)
  use psb_base_mod, psb_protect_name => psb_s_backward_map

  implicit none 
  type(psb_inter_desc_type), intent(in) :: desc
  real(psb_spk_), intent(in)     :: alpha,beta
  real(psb_spk_), intent(inout)  :: x(:)
  real(psb_spk_), intent(out)    :: y(:)
  integer, intent(out)             :: info 
  real(psb_spk_), optional       :: work(:)

  !
  real(psb_spk_), allocatable :: xt(:)
  integer                       :: itsz, i, j,totxch,totsnd,totrcv,&
       & map_kind, map_data
  character(len=20), parameter  :: name='psb_backward_map'

  info = 0
  if (.not.psb_is_asb_desc(desc)) then 
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end if

  itsz     = psb_cd_get_bk_tmp_sz(desc)
  map_kind = psb_cd_get_map_kind(desc)
  map_data = psb_cd_get_map_data(desc)
  if (map_data /= psb_map_single_) then 
    write(0,*) trim(name),' Invalid descriptor inupt: map_data',&
         & map_data,psb_map_single_
    info = 1
    return 
  endif

  select case(map_kind)
  case(psb_map_aggr_)
    ! Ok, we just need to call a halo update and a matrix-vector product. 
    call psb_halo(x,desc%desc_2,info,work=work)
    if (info == 0) call psb_csmm(alpha,desc%smap%map_bk,x,beta,y,info)

    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if


  case(psb_map_gen_linear_)
    call psb_linmap(alpha,x,beta,y,desc%smap%map_bk,&
         & desc%desc_bk,desc%desc_2,desc%desc_1)
    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if

  case default
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end select

end subroutine psb_s_backward_map


!
! Takes a vector X from space desc%desc_1 and maps it onto
! desc%desc_2 under desc%map_fw possibly with communication
! due to exch_fw_idx
!
subroutine psb_d_forward_map(alpha,x,beta,y,desc,info,work)
  use psb_base_mod, psb_protect_name => psb_d_forward_map
  implicit none 
  type(psb_inter_desc_type), intent(in) :: desc
  real(psb_dpk_), intent(in)     :: alpha,beta
  real(psb_dpk_), intent(inout)  :: x(:)
  real(psb_dpk_), intent(out)    :: y(:)
  integer, intent(out)             :: info 
  real(psb_dpk_), optional       :: work(:)

  !
  real(psb_dpk_), allocatable :: xt(:)
  integer                       :: itsz, i, j,totxch,totsnd,totrcv,&
       &  map_kind, map_data
  character(len=20), parameter  :: name='psb_forward_map'

  info = 0
  if (.not.psb_is_asb_desc(desc)) then 
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end if

  itsz     = psb_cd_get_fw_tmp_sz(desc)
  map_kind = psb_cd_get_map_kind(desc)
  map_data = psb_cd_get_map_data(desc)
  if (map_data /= psb_map_double_) then 
    write(0,*) trim(name),' Invalid descriptor inupt: map_data', &
         & map_data,psb_map_double_
    info = 1
    return 
  endif

  select case(map_kind)
  case(psb_map_aggr_)
    ! Ok, we just need to call a halo update on the base desc
    ! and a matrix-vector product. 
    call psb_halo(x,desc%desc_1,info,work=work)
    if (info == 0) call psb_csmm(alpha,desc%dmap%map_fw,x,beta,y,info)

    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if

  case(psb_map_gen_linear_)

    call psb_linmap(alpha,x,beta,y,desc%dmap%map_fw,&
         & desc%desc_fw,desc%desc_1,desc%desc_2)

    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if


  case default
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end select

end subroutine psb_d_forward_map


!
! Takes a vector X from space desc%desc_2 and maps it onto
! desc%desc_1 under desc%map_bk possibly with communication
! due to exch_bk_idx
!
subroutine psb_d_backward_map(alpha,x,beta,y,desc,info,work)
  use psb_base_mod, psb_protect_name => psb_d_backward_map

  implicit none 
  type(psb_inter_desc_type), intent(in) :: desc
  real(psb_dpk_), intent(in)     :: alpha,beta
  real(psb_dpk_), intent(inout)  :: x(:)
  real(psb_dpk_), intent(out)    :: y(:)
  integer, intent(out)             :: info 
  real(psb_dpk_), optional       :: work(:)

  !
  real(psb_dpk_), allocatable :: xt(:)
  integer                       :: itsz, i, j,totxch,totsnd,totrcv,&
       & map_kind, map_data
  character(len=20), parameter  :: name='psb_backward_map'

  info = 0
  if (.not.psb_is_asb_desc(desc)) then 
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end if

  itsz     = psb_cd_get_bk_tmp_sz(desc)
  map_kind = psb_cd_get_map_kind(desc)
  map_data = psb_cd_get_map_data(desc)
  if (map_data /= psb_map_double_) then 
    write(0,*) trim(name),' Invalid descriptor inupt: map_data',&
         & map_data,psb_map_double_
    info = 1
    return 
  endif

  select case(map_kind)
  case(psb_map_aggr_)
    ! Ok, we just need to call a halo update and a matrix-vector product. 
    call psb_halo(x,desc%desc_2,info,work=work)
    if (info == 0) call psb_csmm(alpha,desc%dmap%map_bk,x,beta,y,info)

    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if


  case(psb_map_gen_linear_)
    call psb_linmap(alpha,x,beta,y,desc%dmap%map_bk,&
         & desc%desc_bk,desc%desc_2,desc%desc_1)
    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if

  case default
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end select

end subroutine psb_d_backward_map


!
! Takes a vector X from space desc%desc_1 and maps it onto
! desc%desc_2 under desc%map_fw possibly with communication
! due to exch_fw_idx
!
subroutine psb_c_forward_map(alpha,x,beta,y,desc,info,work)
  use psb_base_mod, psb_protect_name => psb_c_forward_map

  implicit none 
  type(psb_inter_desc_type), intent(in) :: desc
  complex(psb_spk_), intent(in)         :: alpha,beta
  complex(psb_spk_), intent(inout)      :: x(:)
  complex(psb_spk_), intent(out)        :: y(:)
  integer, intent(out)                  :: info 
  complex(psb_spk_), optional           :: work(:)

  !
  complex(psb_spk_), allocatable :: xt(:)
  integer                       :: itsz, i, j,totxch,totsnd,totrcv,&
       & map_kind, map_data
  character(len=20), parameter  :: name='psb_forward_map'

  info = 0
  if (.not.psb_is_asb_desc(desc)) then 
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end if

  itsz     = psb_cd_get_fw_tmp_sz(desc)
  map_kind = psb_cd_get_map_kind(desc)
  map_data = psb_cd_get_map_data(desc)
  if (map_data /= psb_map_complex_) then 
    write(0,*) trim(name),' Invalid descriptor inupt: map_data',&
         & map_data,psb_map_complex_
    info = 1
    return 
  endif

  select case(map_kind)
  case(psb_map_aggr_)
    ! Ok, we just need to call a halo update and a matrix-vector product. 
    call psb_halo(x,desc%desc_1,info,work=work)
    if (info == 0) call psb_csmm(alpha,desc%cmap%map_fw,x,beta,y,info)

    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if

  case(psb_map_gen_linear_)
    call psb_linmap(alpha,x,beta,y,desc%cmap%map_fw,&
         & desc%desc_fw,desc%desc_1,desc%desc_2)

    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if

  case default
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end select

end subroutine psb_c_forward_map


!
! Takes a vector X from space desc%desc_2 and maps it onto
! desc%desc_1 under desc%map_bk possibly with communication
! due to exch_bk_idx
!
subroutine psb_c_backward_map(alpha,x,beta,y,desc,info,work)
  use psb_base_mod, psb_protect_name => psb_c_backward_map

  implicit none 
  type(psb_inter_desc_type), intent(in) :: desc
  complex(psb_spk_), intent(in)       :: alpha,beta
  complex(psb_spk_), intent(inout)    :: x(:)
  complex(psb_spk_), intent(out)      :: y(:)
  integer, intent(out)                  :: info 
  complex(psb_spk_), optional         :: work(:)

  !
  complex(psb_spk_), allocatable :: xt(:)
  integer                       :: itsz, i, j,totxch,totsnd,totrcv,&
       & map_kind, map_data
  character(len=20), parameter  :: name='psb_backward_map'

  info = 0
  if (.not.psb_is_asb_desc(desc)) then 
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end if

  itsz     = psb_cd_get_bk_tmp_sz(desc)
  map_kind = psb_cd_get_map_kind(desc)
  map_data = psb_cd_get_map_data(desc)
  if (map_data /= psb_map_complex_) then 
    write(0,*) trim(name),' Invalid descriptor inupt: map_data',&
         & map_data,psb_map_complex_
    info = 1
    return 
  endif

  select case(map_kind)
  case(psb_map_aggr_)
    ! Ok, we just need to call a halo update and a matrix-vector product. 
    call psb_halo(x,desc%desc_2,info,work=work)
    if (info == 0) call psb_csmm(alpha,desc%cmap%map_bk,x,beta,y,info)

    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if


  case(psb_map_gen_linear_)
    call psb_linmap(alpha,x,beta,y,desc%cmap%map_bk,&
         & desc%desc_bk,desc%desc_1,desc%desc_2)

    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if

  case default
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end select

end subroutine psb_c_backward_map


!
! Takes a vector X from space desc%desc_1 and maps it onto
! desc%desc_2 under desc%map_fw possibly with communication
! due to exch_fw_idx
!
subroutine psb_z_forward_map(alpha,x,beta,y,desc,info,work)
  use psb_base_mod, psb_protect_name => psb_z_forward_map

  implicit none 
  type(psb_inter_desc_type), intent(in) :: desc
  complex(psb_dpk_), intent(in)       :: alpha,beta
  complex(psb_dpk_), intent(inout)    :: x(:)
  complex(psb_dpk_), intent(out)      :: y(:)
  integer, intent(out)                  :: info 
  complex(psb_dpk_), optional         :: work(:)

  !
  complex(psb_dpk_), allocatable :: xt(:)
  integer                       :: itsz, i, j,totxch,totsnd,totrcv,&
       & map_kind, map_data
  character(len=20), parameter  :: name='psb_forward_map'

  info = 0
  if (.not.psb_is_asb_desc(desc)) then 
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end if

  itsz     = psb_cd_get_fw_tmp_sz(desc)
  map_kind = psb_cd_get_map_kind(desc)
  map_data = psb_cd_get_map_data(desc)
  if (map_data /= psb_map_double_complex_) then 
    write(0,*) trim(name),' Invalid descriptor inupt: map_data',&
         & map_data,psb_map_double_complex_
    info = 1
    return 
  endif

  select case(map_kind)
  case(psb_map_aggr_)
    ! Ok, we just need to call a halo update and a matrix-vector product. 
    call psb_halo(x,desc%desc_1,info,work=work)
    if (info == 0) call psb_csmm(alpha,desc%zmap%map_fw,x,beta,y,info)

    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if

  case(psb_map_gen_linear_)
    call psb_linmap(alpha,x,beta,y,desc%zmap%map_fw,&
         & desc%desc_fw,desc%desc_1,desc%desc_2)

    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if

  case default
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end select

end subroutine psb_z_forward_map


!
! Takes a vector X from space desc%desc_2 and maps it onto
! desc%desc_1 under desc%map_bk possibly with communication
! due to exch_bk_idx
!
subroutine psb_z_backward_map(alpha,x,beta,y,desc,info,work)
  use psb_base_mod, psb_protect_name => psb_z_backward_map

  implicit none 
  type(psb_inter_desc_type), intent(in) :: desc
  complex(psb_dpk_), intent(in)       :: alpha,beta
  complex(psb_dpk_), intent(inout)    :: x(:)
  complex(psb_dpk_), intent(out)      :: y(:)
  integer, intent(out)                  :: info 
  complex(psb_dpk_), optional         :: work(:)

  !
  complex(psb_dpk_), allocatable :: xt(:)
  integer                       :: itsz, i, j,totxch,totsnd,totrcv,&
       & map_kind, map_data
  character(len=20), parameter  :: name='psb_backward_map'

  info = 0
  if (.not.psb_is_asb_desc(desc)) then 
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end if

  itsz     = psb_cd_get_bk_tmp_sz(desc)
  map_kind = psb_cd_get_map_kind(desc)
  map_data = psb_cd_get_map_data(desc)
  if (map_data /= psb_map_double_complex_) then 
    write(0,*) trim(name),' Invalid descriptor inupt: map_data',&
         & map_data,psb_map_double_complex_
    info = 1
    return 
  endif

  select case(map_kind)
  case(psb_map_aggr_)
    ! Ok, we just need to call a halo update and a matrix-vector product. 
    call psb_halo(x,desc%desc_2,info,work=work)
    if (info == 0) call psb_csmm(alpha,desc%zmap%map_bk,x,beta,y,info)

    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if


  case(psb_map_gen_linear_)
    call psb_linmap(alpha,x,beta,y,desc%zmap%map_bk,&
         & desc%desc_bk,desc%desc_1,desc%desc_2)

    if (info /= 0) then 
      write(0,*) trim(name),' Error from inner routines',info
      info = -1
    end if

  case default
    write(0,*) trim(name),' Invalid descriptor inupt'
    info = 1
    return 
  end select

end subroutine psb_z_backward_map



subroutine psb_s_apply_linmap(alpha,x,beta,y,a_map,cd_xt,descin,descout)
  use psb_base_mod, psb_protect_name => psb_s_apply_linmap

  implicit none 
  real(psb_spk_), intent(in)      :: alpha,beta
  real(psb_spk_), intent(inout)   :: x(:),y(:)
  type(psb_sspmat_type), intent(in) :: a_map
  type(psb_desc_type), intent(in)   :: cd_xt,descin, descout 

  integer :: nrt, nct, info
  real(psb_spk_), allocatable :: tmp(:)

  nrt = psb_cd_get_local_rows(cd_xt)
  nct = psb_cd_get_local_cols(cd_xt)
  allocate(tmp(nct),stat=info)
  if (info == 0) tmp(1:nrt) = x(1:nrt)
  if (info == 0) call psb_halo(tmp,cd_xt,info) 
  if (info == 0) call psb_csmm(alpha,a_map,tmp,beta,y,info)
  if (info /= 0) then 
    write(0,*) 'Error in apply_map'
  endif

end subroutine psb_s_apply_linmap


subroutine psb_d_apply_linmap(alpha,x,beta,y,a_map,cd_xt,descin,descout)
  use psb_base_mod, psb_protect_name => psb_d_apply_linmap

  implicit none 
  real(psb_dpk_), intent(in)      :: alpha,beta
  real(psb_dpk_), intent(inout)   :: x(:),y(:)
  type(psb_dspmat_type), intent(in) :: a_map
  type(psb_desc_type), intent(in)   :: cd_xt,descin, descout 

  integer :: nrt, nct, info
  real(psb_dpk_), allocatable :: tmp(:)

  nrt = psb_cd_get_local_rows(cd_xt)
  nct = psb_cd_get_local_cols(cd_xt)
  allocate(tmp(nct),stat=info)
  if (info == 0) tmp(1:nrt) = x(1:nrt)
  if (info == 0) call psb_halo(tmp,cd_xt,info) 
  if (info == 0) call psb_csmm(alpha,a_map,tmp,beta,y,info)
  if (info /= 0) then 
    write(0,*) 'Error in apply_map'
  endif

end subroutine psb_d_apply_linmap


subroutine psb_c_apply_linmap(alpha,x,beta,y,a_map,cd_xt,descin,descout)
  use psb_base_mod, psb_protect_name => psb_c_apply_linmap

  implicit none 
  complex(psb_spk_), intent(in)      :: alpha,beta
  complex(psb_spk_), intent(inout)   :: x(:),y(:)
  type(psb_cspmat_type), intent(in) :: a_map
  type(psb_desc_type), intent(in)   :: cd_xt,descin, descout 

  integer :: nrt, nct, info
  complex(psb_spk_), allocatable :: tmp(:)

  nrt = psb_cd_get_local_rows(cd_xt)
  nct = psb_cd_get_local_cols(cd_xt)
  allocate(tmp(nct),stat=info)
  if (info == 0) tmp(1:nrt) = x(1:nrt)
  if (info == 0) call psb_halo(tmp,cd_xt,info) 
  if (info == 0) call psb_csmm(alpha,a_map,tmp,beta,y,info)
  if (info /= 0) then 
    write(0,*) 'Error in apply_map'
  endif

end subroutine psb_c_apply_linmap

subroutine psb_z_apply_linmap(alpha,x,beta,y,a_map,cd_xt,descin,descout)

  use psb_base_mod, psb_protect_name => psb_z_apply_linmap

  implicit none 
  complex(psb_dpk_), intent(in)      :: alpha,beta
  complex(psb_dpk_), intent(inout)   :: x(:),y(:)
  type(psb_zspmat_type), intent(in) :: a_map
  type(psb_desc_type), intent(in)   :: cd_xt,descin, descout 

  integer :: nrt, nct, info
  complex(psb_dpk_), allocatable :: tmp(:)

  nrt = psb_cd_get_local_rows(cd_xt)
  nct = psb_cd_get_local_cols(cd_xt)
  allocate(tmp(nct),stat=info)
  if (info == 0) tmp(1:nrt) = x(1:nrt)
  if (info == 0) call psb_halo(tmp,cd_xt,info) 
  if (info == 0) call psb_csmm(alpha,a_map,tmp,beta,y,info)
  if (info /= 0) then 
    write(0,*) 'Error in apply_map'
  endif

end subroutine psb_z_apply_linmap
