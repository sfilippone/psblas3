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
! package: psb_linmap_type_mod
!    Defines data types for mapping between vectors belonging
!    to different spaces.
!
module psb_linmap_type_mod
  use psb_spmat_type, only : psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_, psb_sizeof

  use psb_mat_mod, only: psb_d_sparse_mat, psb_s_sparse_mat
  use psb_descriptor_type, only: psb_desc_type
  

  ! Inter-descriptor mapping data structures. 
  integer, parameter :: psb_map_kind_     = 1
  integer, parameter :: psb_map_data_     = 2
  integer, parameter :: psb_map_integer_  = 1
  integer, parameter :: psb_map_single_   = 2
  integer, parameter :: psb_map_double_   = 3 
  integer, parameter :: psb_map_complex_  = 4
  integer, parameter :: psb_map_double_complex_ = 5 
 
  integer, parameter :: psb_itd_data_size_=20


  type psb_slinmap_type 
    integer, allocatable   :: itd_data(:), iaggr(:), naggr(:)
    type(psb_desc_type), pointer :: p_desc_X=>null(), p_desc_Y=>null()
    type(psb_desc_type)    :: desc_X, desc_Y
    type(psb_s_sparse_mat) :: map_X2Y, map_Y2X
  end type psb_slinmap_type

  type psb_dlinmap_type 
    integer, allocatable   :: itd_data(:), iaggr(:), naggr(:)
    type(psb_desc_type), pointer :: p_desc_X=>null(), p_desc_Y=>null()
    type(psb_desc_type)    :: desc_X, desc_Y
    type(psb_d_sparse_mat) :: map_X2Y, map_Y2X
  end type psb_dlinmap_type

  type psb_clinmap_type 
    integer, allocatable  :: itd_data(:), iaggr(:), naggr(:)
    type(psb_desc_type), pointer :: p_desc_X=>null(), p_desc_Y=>null()
    type(psb_desc_type)   :: desc_X, desc_Y
    type(psb_cspmat_type) :: map_X2Y, map_Y2X
  end type psb_clinmap_type
  
  type psb_zlinmap_type 
    integer, allocatable  :: itd_data(:), iaggr(:), naggr(:)
    type(psb_desc_type), pointer :: p_desc_X=>null(), p_desc_Y=>null()
    type(psb_desc_type)   :: desc_X, desc_Y
    type(psb_zspmat_type) :: map_X2Y, map_Y2X
  end type psb_zlinmap_type

end module psb_linmap_type_mod

