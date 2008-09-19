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
!
!
! package: psb_inter_descriptor_type
!    Defines facilities for mapping between vectors belonging
!    to different spaces.
!
module psb_inter_desc_mod

  use psb_inter_descriptor_type
  use psb_descriptor_type

  interface psb_forward_map
    subroutine psb_s_forward_map(alpha,x,beta,y,desc,info,work)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type), intent(in) :: desc
      real(psb_spk_), intent(in)     :: alpha,beta
      real(psb_spk_), intent(inout)  :: x(:)
      real(psb_spk_), intent(out)    :: y(:)
      integer, intent(out)           :: info 
      real(psb_spk_), optional       :: work(:)
    end subroutine psb_s_forward_map
    subroutine psb_d_forward_map(alpha,x,beta,y,desc,info,work)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type), intent(in) :: desc
      real(psb_dpk_), intent(in)     :: alpha,beta
      real(psb_dpk_), intent(inout)  :: x(:)
      real(psb_dpk_), intent(out)    :: y(:)
      integer, intent(out)           :: info 
      real(psb_dpk_), optional       :: work(:)
    end subroutine psb_d_forward_map
    subroutine psb_c_forward_map(alpha,x,beta,y,desc,info,work)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type), intent(in) :: desc
      complex(psb_spk_), intent(in)         :: alpha,beta
      complex(psb_spk_), intent(inout)      :: x(:)
      complex(psb_spk_), intent(out)        :: y(:)
      integer, intent(out)                  :: info 
      complex(psb_spk_), optional           :: work(:)
    end subroutine psb_c_forward_map
    subroutine psb_z_forward_map(alpha,x,beta,y,desc,info,work)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type), intent(in) :: desc
      complex(psb_dpk_), intent(in)         :: alpha,beta
      complex(psb_dpk_), intent(inout)      :: x(:)
      complex(psb_dpk_), intent(out)        :: y(:)
      integer, intent(out)                  :: info 
      complex(psb_dpk_), optional           :: work(:)
    end subroutine psb_z_forward_map    
  end interface
  
  interface psb_backward_map
    subroutine psb_s_backward_map(alpha,x,beta,y,desc,info,work)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type), intent(in) :: desc
      real(psb_spk_), intent(in)     :: alpha,beta
      real(psb_spk_), intent(inout)  :: x(:)
      real(psb_spk_), intent(out)    :: y(:)
      integer, intent(out)           :: info 
      real(psb_spk_), optional       :: work(:)
    end subroutine psb_s_backward_map
    subroutine psb_d_backward_map(alpha,x,beta,y,desc,info,work)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type), intent(in) :: desc
      real(psb_dpk_), intent(in)     :: alpha,beta
      real(psb_dpk_), intent(inout)  :: x(:)
      real(psb_dpk_), intent(out)    :: y(:)
      integer, intent(out)           :: info 
      real(psb_dpk_), optional       :: work(:)
    end subroutine psb_d_backward_map
    subroutine psb_c_backward_map(alpha,x,beta,y,desc,info,work)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type), intent(in) :: desc
      complex(psb_spk_), intent(in)         :: alpha,beta
      complex(psb_spk_), intent(inout)      :: x(:)
      complex(psb_spk_), intent(out)        :: y(:)
      integer, intent(out)                  :: info 
      complex(psb_spk_), optional           :: work(:)
    end subroutine psb_c_backward_map
    subroutine psb_z_backward_map(alpha,x,beta,y,desc,info,work)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type), intent(in) :: desc
      complex(psb_dpk_), intent(in)         :: alpha,beta
      complex(psb_dpk_), intent(inout)      :: x(:)
      complex(psb_dpk_), intent(out)        :: y(:)
      integer, intent(out)                  :: info 
      complex(psb_dpk_), optional           :: work(:)
    end subroutine psb_z_backward_map    
  end interface



  interface psb_is_ok_desc
    module procedure psb_is_ok_inter_desc
  end interface

  interface psb_is_asb_desc
    module procedure psb_is_asb_inter_desc
  end interface

  interface psb_inter_desc
    function psb_s_inter_desc(map_kind,desc1,desc2,map_fw,map_bk,idx_fw,idx_bk)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type)         :: psb_s_inter_desc    
      type(psb_desc_type), target       :: desc1, desc2
      type(psb_sspmat_type), intent(in) :: map_fw, map_bk
      integer, intent(in)               :: map_kind,idx_fw(:), idx_bk(:)
    end function psb_s_inter_desc
    function psb_s_inter_desc_noidx(map_kind,desc1, desc2, map_fw, map_bk)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type)         :: psb_s_inter_desc_noidx    
      type(psb_desc_type), target       :: desc1, desc2
      type(psb_sspmat_type), intent(in) :: map_fw, map_bk
      integer, intent(in)               :: map_kind
    end function psb_s_inter_desc_noidx
    function psb_d_inter_desc(map_kind,desc1,desc2,map_fw,map_bk,idx_fw,idx_bk)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type)         :: psb_d_inter_desc    
      type(psb_desc_type), target       :: desc1, desc2
      type(psb_dspmat_type), intent(in) :: map_fw, map_bk
      integer, intent(in)               :: map_kind,idx_fw(:), idx_bk(:)
    end function psb_d_inter_desc
    function psb_d_inter_desc_noidx(map_kind,desc1, desc2, map_fw, map_bk)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type)         :: psb_d_inter_desc_noidx    
      type(psb_desc_type), target       :: desc1, desc2
      type(psb_dspmat_type), intent(in) :: map_fw, map_bk
      integer, intent(in)               :: map_kind
    end function psb_d_inter_desc_noidx
    function psb_c_inter_desc(map_kind,desc1,desc2,map_fw,map_bk,idx_fw,idx_bk)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type)         :: psb_c_inter_desc    
      type(psb_desc_type), target       :: desc1, desc2
      type(psb_cspmat_type), intent(in) :: map_fw, map_bk
      integer, intent(in)               :: map_kind,idx_fw(:), idx_bk(:)
    end function psb_c_inter_desc
    function psb_c_inter_desc_noidx(map_kind,desc1, desc2, map_fw, map_bk)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type)         :: psb_c_inter_desc_noidx    
      type(psb_desc_type), target       :: desc1, desc2
      type(psb_cspmat_type), intent(in) :: map_fw, map_bk
      integer, intent(in)               :: map_kind
    end function psb_c_inter_desc_noidx
    function psb_z_inter_desc(map_kind,desc1,desc2,map_fw,map_bk,idx_fw,idx_bk)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type)         :: psb_z_inter_desc    
      type(psb_desc_type), target       :: desc1, desc2
      type(psb_zspmat_type), intent(in) :: map_fw, map_bk
      integer, intent(in)               :: map_kind,idx_fw(:), idx_bk(:)
    end function psb_z_inter_desc
    function psb_z_inter_desc_noidx(map_kind,desc1, desc2, map_fw, map_bk)
      use psb_inter_descriptor_type
      implicit none 
      type(psb_inter_desc_type)         :: psb_z_inter_desc_noidx    
      type(psb_desc_type), target       :: desc1, desc2
      type(psb_zspmat_type), intent(in) :: map_fw, map_bk
      integer, intent(in)               :: map_kind
    end function psb_z_inter_desc_noidx
  end interface

  interface psb_sizeof
    module procedure psb_itd_sizeof,&
         & psb_s_map_sizeof, psb_c_map_sizeof,&
         & psb_d_map_sizeof, psb_z_map_sizeof
  end interface

  interface psb_linmap
    subroutine psb_s_apply_linmap(alpha,x,beta,y,a_map,cd_xt,descin,descout)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type, &
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_, psb_sizeof      
      use psb_descriptor_type, only: psb_desc_type
      implicit none 
      real(psb_spk_), intent(in)      :: alpha,beta
      real(psb_spk_), intent(inout)   :: x(:),y(:)
      type(psb_sspmat_type), intent(in) :: a_map
      type(psb_desc_type), intent(in)   :: cd_xt,descin, descout 
    end subroutine psb_s_apply_linmap
    subroutine psb_d_apply_linmap(alpha,x,beta,y,a_map,cd_xt,descin,descout)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type, &
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_, psb_sizeof
      use psb_descriptor_type, only: psb_desc_type
      implicit none 
      real(psb_dpk_), intent(in)      :: alpha,beta
      real(psb_dpk_), intent(inout)   :: x(:),y(:)
      type(psb_dspmat_type), intent(in) :: a_map
      type(psb_desc_type), intent(in)   :: cd_xt,descin, descout 
    end subroutine psb_d_apply_linmap
    subroutine psb_c_apply_linmap(alpha,x,beta,y,a_map,cd_xt,descin,descout)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type, &
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_, psb_sizeof
      use psb_descriptor_type, only: psb_desc_type
      implicit none 
      complex(psb_spk_), intent(in)     :: alpha,beta
      complex(psb_spk_), intent(inout)  :: x(:),y(:)
      type(psb_cspmat_type), intent(in) :: a_map
      type(psb_desc_type), intent(in)   :: cd_xt,descin, descout 
    end subroutine psb_c_apply_linmap
    subroutine psb_z_apply_linmap(alpha,x,beta,y,a_map,cd_xt,descin,descout)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type, &
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_, psb_sizeof
      use psb_descriptor_type, only: psb_desc_type
      implicit none 
      complex(psb_dpk_), intent(in)     :: alpha,beta
      complex(psb_dpk_), intent(inout)  :: x(:),y(:)
      type(psb_zspmat_type), intent(in) :: a_map
      type(psb_desc_type), intent(in)   :: cd_xt,descin, descout 
    end subroutine psb_z_apply_linmap
  end interface

contains
 
  function psb_cd_get_map_kind(desc)    
    implicit none
    type(psb_inter_desc_type), intent(in) :: desc
    Integer                      :: psb_cd_get_map_kind
    if (psb_is_ok_desc(desc)) then 
      psb_cd_get_map_kind = desc%itd_data(psb_map_kind_) 
    else    
      psb_cd_get_map_kind = -1
    end if
  end function psb_cd_get_map_kind
 
  subroutine psb_cd_set_map_kind(map_kind,desc)    
    implicit none
    integer, intent(in)          :: map_kind
    type(psb_inter_desc_type), intent(inout) :: desc

    desc%itd_data(psb_map_kind_) = map_kind

  end subroutine psb_cd_set_map_kind
 
  function psb_cd_get_map_data(desc)    
    implicit none
    type(psb_inter_desc_type), intent(in) :: desc
    Integer                      :: psb_cd_get_map_data
    if (psb_is_ok_desc(desc)) then 
      psb_cd_get_map_data = desc%itd_data(psb_map_data_) 
    else    
      psb_cd_get_map_data = -1
    end if
  end function psb_cd_get_map_data
 
  subroutine psb_cd_set_map_data(map_data,desc)    
    implicit none
    integer, intent(in)          :: map_data
    type(psb_inter_desc_type), intent(inout) :: desc

    
    desc%itd_data(psb_map_data_) = map_data

  end subroutine psb_cd_set_map_data

 
  function psb_cd_get_fw_tmp_sz(desc)    
    implicit none
    type(psb_inter_desc_type), intent(in) :: desc
    Integer                      :: psb_cd_get_fw_tmp_sz
    
    psb_cd_get_fw_tmp_sz = desc%itd_data(psb_fw_tmp_sz_) 
  end function psb_cd_get_fw_tmp_sz

  function psb_cd_get_bk_tmp_sz(desc)    
    implicit none
    type(psb_inter_desc_type), intent(in) :: desc
    Integer                      :: psb_cd_get_bk_tmp_sz
    
    psb_cd_get_bk_tmp_sz = desc%itd_data(psb_bk_tmp_sz_) 
  end function psb_cd_get_bk_tmp_sz

  subroutine psb_cd_set_fw_tmp_sz(isz,desc)    
    implicit none
    type(psb_inter_desc_type), intent(inout) :: desc
    integer, intent(in)                      :: isz 
    
    desc%itd_data(psb_fw_tmp_sz_) =isz
  end subroutine psb_cd_set_fw_tmp_sz

  subroutine psb_cd_set_bk_tmp_sz(isz,desc)    
    implicit none
    type(psb_inter_desc_type), intent(inout) :: desc
    integer, intent(in)                      :: isz 
    
    desc%itd_data(psb_bk_tmp_sz_) =isz

  end subroutine psb_cd_set_bk_tmp_sz


  logical function psb_is_asb_inter_desc(desc)
    implicit none 
    type(psb_inter_desc_type), intent(in) :: desc

    psb_is_asb_inter_desc = .false.
    if (.not.allocated(desc%itd_data)) return
    if (.not.associated(desc%desc_1)) return
    if (.not.associated(desc%desc_2)) return
    psb_is_asb_inter_desc = &
         & psb_is_asb_desc(desc%desc_1).and.psb_is_asb_desc(desc%desc_2)    

  end function psb_is_asb_inter_desc

  logical function psb_is_ok_inter_desc(desc)
    implicit none 
    type(psb_inter_desc_type), intent(in) :: desc

    psb_is_ok_inter_desc = .false.
    if (.not.allocated(desc%itd_data)) return
    select case(desc%itd_data(psb_map_data_))
    case(psb_map_integer_, psb_map_single_, psb_map_complex_,&
         &  psb_map_double_, psb_map_double_complex_) 
      ! Ok go ahead
    case default
      ! Since it's false so far, simply return
      return
    end select
    if (.not.associated(desc%desc_1)) return
    if (.not.associated(desc%desc_2)) return
    psb_is_ok_inter_desc = &
         & psb_is_ok_desc(desc%desc_1).and.psb_is_ok_desc(desc%desc_2)    

  end function psb_is_ok_inter_desc


  function psb_s_map_sizeof(map) result(val)
    implicit none
    type(psb_s_map_type), intent(in) :: map
    integer(psb_long_int_k_) :: val

    val = 0
    val = val + psb_sizeof(map%map_fw)
    val = val + psb_sizeof(map%map_bk)

  end function psb_s_map_sizeof

  function psb_d_map_sizeof(map) result(val)
    implicit none
    type(psb_d_map_type), intent(in) :: map
    integer(psb_long_int_k_) :: val

    val = 0
    val = val + psb_sizeof(map%map_fw)
    val = val + psb_sizeof(map%map_bk)

  end function psb_d_map_sizeof

  function psb_c_map_sizeof(map) result(val)
    implicit none
    type(psb_c_map_type), intent(in) :: map
    integer(psb_long_int_k_) :: val

    val = 0
    val = val + psb_sizeof(map%map_fw)
    val = val + psb_sizeof(map%map_bk)

  end function psb_c_map_sizeof

  function psb_z_map_sizeof(map) result(val)
    implicit none
    type(psb_z_map_type), intent(in) :: map
    integer(psb_long_int_k_) :: val

    val = 0
    val = val + psb_sizeof(map%map_fw)
    val = val + psb_sizeof(map%map_bk)

  end function psb_z_map_sizeof

  function psb_itd_sizeof(desc) result(val)
    implicit none 
    type(psb_inter_desc_type), intent(in) :: desc
    integer(psb_long_int_k_) :: val

    val = 0
    if (allocated(desc%itd_data))    val = val + psb_sizeof_int*size(desc%itd_data)
    if (allocated(desc%exch_fw_idx)) val = val + psb_sizeof_int*size(desc%exch_fw_idx)
    if (allocated(desc%exch_bk_idx)) val = val + psb_sizeof_int*size(desc%exch_bk_idx)
    val = val + psb_sizeof(desc%desc_fw)
    val = val + psb_sizeof(desc%desc_bk)
    val = val + psb_sizeof(desc%dmap)
    val = val + psb_sizeof(desc%zmap)

  end function psb_itd_sizeof




end module psb_inter_desc_mod
