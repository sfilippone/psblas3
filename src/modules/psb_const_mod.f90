!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
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

module psb_const_mod

  !
  !     Communication, prolongation & restriction
  !
  integer, parameter :: psb_nohalo_=0,  psb_halo_=4
  integer, parameter :: psb_none_=0,  psb_sum_=1
  integer, parameter :: psb_avg_=2,  psb_square_root_=3
  integer, parameter :: psb_swap_send_=1, psb_swap_recv_=2
  integer, parameter :: psb_swap_sync_=4, psb_swap_mpi_=8

  !
  !     Data checks
  !
  integer, parameter :: psb_deadlock_check_=0
  integer, parameter :: psb_local_mtrx_check_=1
  integer, parameter :: psb_local_comm_check_=2
  integer, parameter :: psb_consistency_check_=3
  integer, parameter :: psb_global_check_=4
  integer, parameter :: psb_order_communication_=5
  integer, parameter :: psb_change_represent_=6
  integer, parameter :: psb_loc_to_glob_check_=7
  integer, parameter :: psb_convert_halo_=1, psb_convert_ovrlap_=2
  integer, parameter :: psb_act_ret_=0, psb_act_abort_=1, no_err_=0
  !
  !     Entries and values in desc%matrix_data
  !
  integer, parameter :: psb_dec_type_=1, psb_m_=2,psb_n_=3
  integer, parameter :: psb_n_row_=4,  psb_n_col_=5,psb_ctxt_=6
  integer, parameter :: psb_loc_to_glob_=7
  integer, parameter :: psb_thal_xch_=11
  integer, parameter :: psb_thal_snd_=12
  integer, parameter :: psb_thal_rcv_=13
  integer, parameter :: psb_tovr_xch_=14
  integer, parameter :: psb_tovr_snd_=15
  integer, parameter :: psb_tovr_rcv_=16
  integer, parameter :: psb_mpi_c_=9,psb_mdata_size_=20
  integer, parameter :: psb_desc_asb_=3099
  integer, parameter :: psb_desc_bld_=psb_desc_asb_+1
  integer, parameter :: psb_desc_repl_=3199
  integer, parameter :: psb_desc_upd_=psb_desc_bld_+1
  integer, parameter :: psb_desc_upd_asb_=psb_desc_upd_+1
  integer, parameter :: psb_desc_large_asb_=psb_desc_upd_asb_+1
  integer, parameter :: psb_desc_large_bld_=psb_desc_large_asb_+1
  integer, parameter :: nbits=14
  integer, parameter :: hashsize=2**nbits, hashmask=hashsize-1
  integer, parameter :: psb_default_large_threshold=4*1024*1024   ! to be reviewed
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
  integer, parameter :: psb_no_comm_=-1
  integer, parameter :: psb_comm_halo_=0, psb_comm_ovr_=1

  !
  !     Queries into spmat%info
  !     
  integer, parameter :: psb_nztotreq_=1, psb_nzrowreq_=2
  integer, parameter :: psb_nzsizereq_=3
  !
  !     Entries and values for  spmat%info
  !     

  integer, parameter :: psb_nnz_=1
  integer, parameter :: psb_del_bnd_=7, psb_srtd_=8
  integer, parameter :: psb_state_=9
  integer, parameter :: psb_upd_pnt_=10
  integer, parameter :: psb_dupl_=11,  psb_upd_=12
  integer, parameter :: psb_ifasize_=16
  integer, parameter :: psb_spmat_null_=0, psb_spmat_bld_=1
  integer, parameter :: psb_spmat_asb_=2, psb_spmat_upd_=4
  integer, parameter :: psb_ireg_flgs_=10, psb_ip2_=0
  integer, parameter :: psb_iflag_=2, psb_ichk_=3
  integer, parameter :: psb_nnzt_=4, psb_zero_=5,psb_ipc_=6
  integer, parameter :: psb_dupl_ovwrt_ = 0
  integer, parameter :: psb_dupl_add_   = 1
  integer, parameter :: psb_dupl_err_   = 2
  integer, parameter :: psb_dupl_def_   = psb_dupl_ovwrt_
  integer, parameter :: psb_upd_dflt_   = 0
  integer, parameter :: psb_upd_srch_   = 98764
  integer, parameter :: psb_upd_perm_   = 98765
  integer, parameter :: psb_isrtdcoo_   = 98761
  integer, parameter :: psb_maxjdrows_=8, psb_minjdrows_=4
  integer, parameter :: psb_dbleint_=2
  !
  !     Error handling 
  !
  integer, parameter :: act_ret=0, act_abort=1, no_err=0

  !
  !     Handy & miscellaneous constants
  !
  integer, parameter :: ione=1, izero=0
  integer, parameter :: itwo=2, ithree=3,mone=-1, psb_root_=0
  real(kind(1.d0)), parameter :: psb_colrow_=0.33, psb_percent_=0.7
  real(kind(1.d0)), parameter :: dzero=0.d0, done=1.d0
  complex(kind(1.d0)), parameter :: zzero=(0.d0,0.0d0)
  complex(kind(1.d0)), parameter :: zone=(1.d0,0.0d0)
  real(kind(1.d0)), parameter :: epstol=1.d-32

  character, parameter :: psb_all_='A',  psb_topdef_=' '
  character(len=5)     :: psb_fidef_='CSR'

end module psb_const_mod
