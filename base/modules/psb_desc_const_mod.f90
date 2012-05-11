!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012
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
!
! package: psb_desc_const_mod
!    Auxiliary module for descriptor: constant values. 
!
module psb_desc_const_mod
  use psb_const_mod, only : psb_ipk_, psb_mpik_
  !
  !     Communication, prolongation & restriction
  !
  integer(psb_ipk_), parameter :: psb_nohalo_=0,  psb_halo_=1
  ! For overlap update. 
  integer(psb_ipk_), parameter :: psb_none_=0,  psb_sum_=1
  integer(psb_ipk_), parameter :: psb_avg_=2,  psb_square_root_=3
  integer(psb_ipk_), parameter :: psb_setzero_=4

  ! The following are bit fields. 
  integer(psb_ipk_), parameter :: psb_swap_send_=1, psb_swap_recv_=2
  integer(psb_ipk_), parameter :: psb_swap_sync_=4, psb_swap_mpi_=8
  ! Choice among lists on which to base data exchange
  integer(psb_ipk_), parameter :: psb_no_comm_=-1
  integer(psb_ipk_), parameter :: psb_comm_halo_=1, psb_comm_ovr_=2
  integer(psb_ipk_), parameter :: psb_comm_ext_=3,  psb_comm_mov_=4
  ! Types of mapping between descriptors.
  integer(psb_ipk_), parameter :: psb_map_xhal_        = 123
  integer(psb_ipk_), parameter :: psb_map_asov_        = psb_map_xhal_+1
  integer(psb_ipk_), parameter :: psb_map_aggr_        = psb_map_asov_+1 
  integer(psb_ipk_), parameter :: psb_map_gen_linear_  = psb_map_aggr_+1 

  integer(psb_ipk_), parameter :: psb_ovt_xhal_ = psb_map_xhal_, psb_ovt_asov_=psb_map_asov_
  !
  ! Entries and values in desc%matrix_data
  !
  integer(psb_ipk_), parameter :: psb_dec_type_  =  1
  integer(psb_ipk_), parameter :: psb_m_         =  2
  integer(psb_ipk_), parameter :: psb_n_         =  3
  integer(psb_ipk_), parameter :: psb_n_row_     =  4
  integer(psb_ipk_), parameter :: psb_n_col_     =  5
  integer(psb_ipk_), parameter :: psb_ctxt_      =  6
  integer(psb_ipk_), parameter :: psb_desc_size_ =  7
  integer(psb_ipk_), parameter :: psb_mpi_c_     =  9
  integer(psb_ipk_), parameter :: psb_pnt_h_     = 10
  integer(psb_ipk_), parameter :: psb_thal_xch_  = 11
  integer(psb_ipk_), parameter :: psb_thal_snd_  = 12
  integer(psb_ipk_), parameter :: psb_thal_rcv_  = 13
  integer(psb_ipk_), parameter :: psb_tovr_xch_  = 14
  integer(psb_ipk_), parameter :: psb_tovr_snd_  = 15
  integer(psb_ipk_), parameter :: psb_tovr_rcv_  = 16
  integer(psb_ipk_), parameter :: psb_text_xch_  = 17
  integer(psb_ipk_), parameter :: psb_text_snd_  = 18
  integer(psb_ipk_), parameter :: psb_text_rcv_  = 19
  integer(psb_ipk_), parameter :: psb_tmov_xch_  = 20
  integer(psb_ipk_), parameter :: psb_tmov_snd_  = 21
  integer(psb_ipk_), parameter :: psb_tmov_rcv_  = 22
  integer(psb_ipk_), parameter :: psb_mdata_size_= 24
  integer(psb_ipk_), parameter :: psb_desc_invalid_=-1
  integer(psb_ipk_), parameter :: psb_desc_null_=-1
  integer(psb_ipk_), parameter :: psb_desc_asb_=3099
  integer(psb_ipk_), parameter :: psb_desc_bld_=psb_desc_asb_+1
  integer(psb_ipk_), parameter :: psb_desc_upd_=psb_desc_bld_+1
  integer(psb_ipk_), parameter :: psb_desc_repl_=3199
  integer(psb_ipk_), parameter :: psb_desc_ovl_bld_=3399
  integer(psb_ipk_), parameter :: psb_desc_ovl_asb_=psb_desc_ovl_bld_+1
  ! these two are reserved for descriptors which are
  ! "overlap-extensions" of base descriptors. 
  integer(psb_ipk_), parameter :: psb_cd_ovl_bld_=psb_desc_ovl_bld_
  integer(psb_ipk_), parameter :: psb_cd_ovl_asb_=psb_desc_ovl_asb_
  integer(psb_ipk_), parameter :: psb_desc_normal_=3299
  integer(psb_ipk_), parameter :: psb_desc_large_=psb_desc_normal_+1
  !
  ! Constants for hashing into desc%hashv(:) and desc%glb_lc(:,:)
  !
  integer(psb_ipk_), parameter :: psb_hash_bits=16
  integer(psb_ipk_), parameter :: psb_max_hash_bits=22
  integer(psb_ipk_), parameter :: psb_hash_size=2**psb_hash_bits, psb_hash_mask=psb_hash_size-1
  integer(psb_ipk_), parameter :: psb_default_large_threshold=1*1024*1024   
  integer(psb_ipk_), parameter :: psb_hpnt_nentries_=7

  !
  !     Constants for desc_a handling
  !

  integer(psb_ipk_), parameter :: psb_upd_glbnum_=998
  integer(psb_ipk_), parameter :: psb_upd_locnum_=997
  integer(psb_ipk_), parameter :: psb_proc_id_=0, psb_n_elem_recv_=1
  integer(psb_ipk_), parameter :: psb_elem_recv_=2, psb_n_elem_send_=2
  integer(psb_ipk_), parameter :: psb_elem_send_=3, psb_n_ovrlp_elem_=1
  integer(psb_ipk_), parameter :: psb_ovrlp_elem_to_=2, psb_ovrlp_elem_=0
  integer(psb_ipk_), parameter :: psb_n_dom_ovr_=1

  interface 
    subroutine psb_parts(glob_index,nrow,np,pv,nv)
      import :: psb_ipk_
      integer(psb_ipk_), intent (in)  :: glob_index,nrow, np
      integer(psb_ipk_), intent (out) :: nv, pv(*)
    end subroutine psb_parts
  end interface
  
end module psb_desc_const_mod
