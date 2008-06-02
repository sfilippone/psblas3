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
module psb_inter_descriptor_type
  use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type, &
       & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_, psb_sizeof

  use psb_descriptor_type, only: psb_desc_type
  

  ! Inter-descriptor mapping data structures. 
  integer, parameter :: psb_map_kind_     = 1
  integer, parameter :: psb_map_data_     = 2
  integer, parameter :: psb_map_integer_  = 1
  integer, parameter :: psb_map_single_   = 2
  integer, parameter :: psb_map_double_   = 3 
  integer, parameter :: psb_map_complex_  = 4
  integer, parameter :: psb_map_double_complex_ = 5 
 
  integer, parameter :: psb_fw_tmp_kind_ = 5 
  integer, parameter :: psb_fw_tmp_sz_   = 6 
  integer, parameter :: psb_bk_tmp_kind_ = 7
  integer, parameter :: psb_bk_tmp_sz_   = 8
  integer, parameter :: psb_itd_data_size_=20


  type psb_s_map_type
    type(psb_sspmat_type) :: map_fw, map_bk
  end type psb_s_map_type

  type psb_d_map_type
    type(psb_dspmat_type) :: map_fw, map_bk
  end type psb_d_map_type

  type psb_c_map_type
    type(psb_cspmat_type) :: map_fw, map_bk
  end type psb_c_map_type

  type psb_z_map_type
    type(psb_zspmat_type) :: map_fw, map_bk
  end type psb_z_map_type
  
  type psb_inter_desc_type 
    integer, allocatable :: itd_data(:)
    type(psb_desc_type), pointer :: desc_1=>null(), desc_2=>null()
    integer, allocatable :: exch_fw_idx(:), exch_bk_idx(:)
    type(psb_desc_type)  :: desc_fw, desc_bk
    type(psb_s_map_type) :: smap
    type(psb_d_map_type) :: dmap
    type(psb_c_map_type) :: cmap
    type(psb_z_map_type) :: zmap
  end type psb_inter_desc_type

end module psb_inter_descriptor_type

