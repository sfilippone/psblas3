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

module psb_const_mod
  ! This is the default integer
#if defined(LONG_INTEGERS)
  integer, parameter  :: ndig=12
#else  
  integer, parameter  :: ndig=8
#endif
  integer, parameter  :: psb_int_k_ = selected_int_kind(ndig)
  ! This is an 8-byte  integer, and normally different from default integer. 
  integer, parameter  :: longndig=12
  integer, parameter  :: psb_long_int_k_ = selected_int_kind(longndig)
  !
  ! These must be the kind parameter corresponding to MPI_DOUBLE_PRECISION
  ! and MPI_REAL
  !
  integer, parameter  :: psb_dpk_ = kind(1.d0)
  integer, parameter  :: psb_spk_ = kind(1.e0)
  integer, save       :: psb_sizeof_dp, psb_sizeof_sp
  integer, save       :: psb_sizeof_int, psb_sizeof_long_int
  integer, save       :: psb_mpi_integer

  !
  !     Handy & miscellaneous constants
  !
  integer, parameter             :: izero=0, ione=1
  integer, parameter             :: itwo=2, ithree=3,mone=-1, psb_root_=0
  real(psb_spk_), parameter      :: szero=0.e0, sone=1.e0
  real(psb_dpk_), parameter      :: dzero=0.d0, done=1.d0
  complex(psb_spk_), parameter   :: czero=(0.e0,0.0e0)
  complex(psb_spk_), parameter   :: cone=(1.e0,0.0e0)
  complex(psb_dpk_), parameter   :: zzero=(0.d0,0.0d0)
  complex(psb_dpk_), parameter   :: zone=(1.d0,0.0d0)
  real(psb_dpk_), parameter      :: d_epstol=1.1d-16 ! Unit roundoff.  
  real(psb_spk_), parameter      :: s_epstol=5.e-8   ! Is this right?
  character, parameter           :: psb_all_='A',  psb_topdef_=' '
  
  !
  ! Sparse matrix constants
  !

  !
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
  !  
  integer, parameter :: psb_invalid_ = -1 
  integer, parameter :: psb_spmat_null_=0, psb_spmat_bld_=1
  integer, parameter :: psb_spmat_asb_=2, psb_spmat_upd_=4

  integer, parameter :: psb_ireg_flgs_=10, psb_ip2_=0
  integer, parameter :: psb_iflag_=2, psb_ichk_=3
  integer, parameter :: psb_nnzt_=4, psb_zero_=5,psb_ipc_=6
  ! Duplicate coefficients handling
  ! These are usually set while calling spcnv as one of its
  ! optional arugments.
  integer, parameter :: psb_dupl_ovwrt_ = 0
  integer, parameter :: psb_dupl_add_   = 1
  integer, parameter :: psb_dupl_err_   = 2
  integer, parameter :: psb_dupl_def_   = psb_dupl_ovwrt_
  ! Matrix update mode
  integer, parameter :: psb_upd_srch_   = 98764
  integer, parameter :: psb_upd_perm_   = 98765
  integer, parameter :: psb_upd_dflt_   = psb_upd_srch_
  ! Mark a COO matrix with sorted entries.
  integer, parameter :: psb_isrtdcoo_   = 98761
  integer, parameter :: psb_maxjdrows_=8, psb_minjdrows_=4
  integer, parameter :: psb_dbleint_=2
  character(len=5)   :: psb_fidef_='CSR'

  !
  ! 
  !     Error constants
  integer, parameter, public :: psb_success_=0
  integer, parameter, public :: psb_err_pivot_too_small_=2
  integer, parameter, public :: psb_err_invalid_ovr_num_=3
  integer, parameter, public :: psb_err_invalid_input_=5
  integer, parameter, public :: psb_err_iarg_neg_=10
  integer, parameter, public :: psb_err_iarg_pos_=20
  integer, parameter, public :: psb_err_input_value_invalid_i_=30
  integer, parameter, public :: psb_err_input_asize_invalid_i_=35
  integer, parameter, public :: psb_err_iarg_invalid_i_=40
  integer, parameter, public :: psb_err_iarg_not_gtia_ii_=50
  integer, parameter, public :: psb_err_iarg_not_gteia_ii_=60
  integer, parameter, public :: psb_err_iarg_invalid_value_=70
  integer, parameter, public :: psb_err_asb_nrc_error_=71
  integer, parameter, public :: psb_err_iarg2_neg_=80
  integer, parameter, public :: psb_err_ia2_not_increasing_=90
  integer, parameter, public :: psb_err_ia1_not_increasing_=91
  integer, parameter, public :: psb_err_ia1_badindices_=100
  integer, parameter, public :: psb_err_invalid_args_combination_=110
  integer, parameter, public :: psb_err_invalid_pid_arg_=115
  integer, parameter, public :: psb_err_iarg_n_mbgtian_=120
  integer, parameter, public :: psb_err_dupl_cd_vl=123
  integer, parameter, public :: psb_err_duplicate_coo=130
  integer, parameter, public :: psb_err_invalid_input_format_=134
  integer, parameter, public :: psb_err_unsupported_format_=135
  integer, parameter, public :: psb_err_format_unknown_=136
  integer, parameter, public :: psb_err_iarray_outside_bounds_=140
  integer, parameter, public :: psb_err_iarray_outside_process_=150
  integer, parameter, public :: psb_err_forgot_geall_=290
  integer, parameter, public :: psb_err_forgot_spall_=295
  integer, parameter, public :: psb_err_wrong_ins_=298
  integer, parameter, public :: psb_err_iarg_mbeeiarra_i_=300
  integer, parameter, public :: psb_err_mpi_error_=400
  integer, parameter, public :: psb_err_parm_differs_among_procs_=550
  integer, parameter, public :: psb_err_entry_out_of_bounds_=551
  integer, parameter, public :: psb_err_inconsistent_index_lists_=552
  integer, parameter, public :: psb_err_partfunc_toomuchprocs_=570
  integer, parameter, public :: psb_err_partfunc_toofewprocs_=575
  integer, parameter, public :: psb_err_partfunc_wrong_pid_=580
  integer, parameter, public :: psb_err_no_optional_arg_=581
  integer, parameter, public :: psb_err_arg_m_required_=582
  integer, parameter, public :: psb_err_many_optional_arg_=583
  integer, parameter, public :: psb_err_spmat_invalid_state_=600
  integer, parameter, public :: psb_err_missing_override_method_=700
  integer, parameter, public :: psb_err_invalid_cd_state_=1122
  integer, parameter, public :: psb_err_invalid_a_and_cd_state_=1123
  integer, parameter, public :: psb_err_context_error_=2010
  integer, parameter, public :: psb_err_initerror_neugh_procs_=2011
  integer, parameter, public :: psb_err_invalid_matrix_input_state_=2231
  integer, parameter, public :: psb_err_input_no_regen_=2232
  integer, parameter, public :: psb_err_lld_case_not_implemented_=3010
  integer, parameter, public :: psb_err_transpose_unsupported_=3015
  integer, parameter, public :: psb_err_transpose_c_unsupported_=3020
  integer, parameter, public :: psb_err_transpose_not_n_unsupported_=3021
  integer, parameter, public :: psb_err_only_unit_diag_=3022
  integer, parameter, public :: psb_err_ja_nix_ia_niy_unsupported_=3030
  integer, parameter, public :: psb_err_ix_n1_iy_n1_unsupported_=3040
  integer, parameter, public :: psb_err_input_matrix_unassembled_=3110
  integer, parameter, public :: psb_err_alloc_dealloc_=4000
  integer, parameter, public :: psb_err_internal_error_=4001
  integer, parameter, public :: psb_err_from_subroutine_=4010
  integer, parameter, public :: psb_err_from_subroutine_i_=4012
  integer, parameter, public :: psb_err_from_subroutine_ai_=4013
  integer, parameter, public :: psb_err_alloc_request_=4025
  integer, parameter, public :: psb_err_from_subroutine_non_=4011
  integer, parameter, public :: psb_err_invalid_istop_=5001

end module psb_const_mod
