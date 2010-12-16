module psb_desc_const_mod
  !
  !     Communication, prolongation & restriction
  !
  integer, parameter :: psb_nohalo_=0,  psb_halo_=1
  ! For overlap update. 
  integer, parameter :: psb_none_=0,  psb_sum_=1
  integer, parameter :: psb_avg_=2,  psb_square_root_=3
  integer, parameter :: psb_setzero_=4

  ! The following are bit fields. 
  integer, parameter :: psb_swap_send_=1, psb_swap_recv_=2
  integer, parameter :: psb_swap_sync_=4, psb_swap_mpi_=8
  ! Choice among lists on which to base data exchange
  integer, parameter :: psb_no_comm_=-1
  integer, parameter :: psb_comm_halo_=1, psb_comm_ovr_=2
  integer, parameter :: psb_comm_ext_=3,  psb_comm_mov_=4
  ! Types of mapping between descriptors.
  integer, parameter :: psb_map_xhal_        = 123
  integer, parameter :: psb_map_asov_        = psb_map_xhal_+1
  integer, parameter :: psb_map_aggr_        = psb_map_asov_+1 
  integer, parameter :: psb_map_gen_linear_  = psb_map_aggr_+1 

  integer, parameter :: psb_ovt_xhal_ = psb_map_xhal_, psb_ovt_asov_=psb_map_asov_
  !
  ! Entries and values in desc%matrix_data
  !
  integer, parameter :: psb_dec_type_  =  1
  integer, parameter :: psb_m_         =  2
  integer, parameter :: psb_n_         =  3
  integer, parameter :: psb_n_row_     =  4
  integer, parameter :: psb_n_col_     =  5
  integer, parameter :: psb_ctxt_      =  6
  integer, parameter :: psb_desc_size_ =  7
  integer, parameter :: psb_mpi_c_     =  9
  integer, parameter :: psb_pnt_h_     = 10
  integer, parameter :: psb_thal_xch_  = 11
  integer, parameter :: psb_thal_snd_  = 12
  integer, parameter :: psb_thal_rcv_  = 13
  integer, parameter :: psb_tovr_xch_  = 14
  integer, parameter :: psb_tovr_snd_  = 15
  integer, parameter :: psb_tovr_rcv_  = 16
  integer, parameter :: psb_text_xch_  = 17
  integer, parameter :: psb_text_snd_  = 18
  integer, parameter :: psb_text_rcv_  = 19
  integer, parameter :: psb_tmov_xch_  = 20
  integer, parameter :: psb_tmov_snd_  = 21
  integer, parameter :: psb_tmov_rcv_  = 22
  integer, parameter :: psb_mdata_size_= 24
  integer, parameter :: psb_desc_invalid_=-1
  integer, parameter :: psb_desc_null_=-1
  integer, parameter :: psb_desc_asb_=3099
  integer, parameter :: psb_desc_bld_=psb_desc_asb_+1
  integer, parameter :: psb_desc_upd_=psb_desc_bld_+1
  integer, parameter :: psb_desc_repl_=3199
  integer, parameter :: psb_desc_ovl_bld_=3399
  integer, parameter :: psb_desc_ovl_asb_=psb_desc_ovl_bld_+1
  ! these two are reserved for descriptors which are
  ! "overlap-extensions" of base descriptors. 
  integer, parameter :: psb_cd_ovl_bld_=psb_desc_ovl_bld_
  integer, parameter :: psb_cd_ovl_asb_=psb_desc_ovl_asb_
  integer, parameter :: psb_desc_normal_=3299
  integer, parameter :: psb_desc_large_=psb_desc_normal_+1
  !
  ! Constants for hashing into desc%hashv(:) and desc%glb_lc(:,:)
  !
  integer, parameter :: psb_hash_bits=16
  integer, parameter :: psb_max_hash_bits=22
  integer, parameter :: psb_hash_size=2**psb_hash_bits, psb_hash_mask=psb_hash_size-1
  integer, parameter :: psb_default_large_threshold=1*1024*1024   
  integer, parameter :: psb_hpnt_nentries_=7

  !
  !     Constants for desc_a handling
  !

  integer, parameter :: psb_upd_glbnum_=998
  integer, parameter :: psb_upd_locnum_=997
  integer, parameter :: psb_proc_id_=0, psb_n_elem_recv_=1
  integer, parameter :: psb_elem_recv_=2, psb_n_elem_send_=2
  integer, parameter :: psb_elem_send_=3, psb_n_ovrlp_elem_=1
  integer, parameter :: psb_ovrlp_elem_to_=2, psb_ovrlp_elem_=0
  integer, parameter :: psb_n_dom_ovr_=1

end module psb_desc_const_mod
