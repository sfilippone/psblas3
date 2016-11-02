!!$ 
!!$              Parallel Sparse BLAS  version 3.4
!!$    (C) Copyright 2006, 2010, 2015
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

module psb_const_mod
#if defined(HAVE_ISO_FORTRAN_ENV)
  use iso_fortran_env
  ! This is the default PSBLAS integer, can be 4 or 8 bytes.
#if defined(LONG_INTEGERS)
  integer, parameter  :: psb_ipk_ = int64
#else  
  integer, parameter  :: psb_ipk_ = int32
#endif
  ! This is always an 8-byte  integer.
  integer, parameter  :: psb_long_int_k_ = int64
  ! This is always a 4-byte integer, for MPI-related stuff
  integer, parameter  :: psb_mpik_ = int32
  !
  ! These must be the kind parameter corresponding to psb_mpi_r_dpk_
  ! and psb_mpi_r_spk_
  !
  integer(psb_mpik_), parameter  :: psb_spk_   = real32
  integer(psb_mpik_), parameter  :: psb_dpk_   = real64
#else
  ! This is the default PSBLAS integer, can be 4 or 8 bytes.
#if defined(LONG_INTEGERS)
  integer, parameter  :: ndig=12
#else  
  integer, parameter  :: ndig=8
#endif
  integer, parameter  :: psb_ipk_ = selected_int_kind(ndig)
  ! This is always an 8-byte  integer.
  integer, parameter  :: longndig=12
  integer, parameter  :: psb_long_int_k_ = selected_int_kind(longndig)
  ! This is always a 4-byte integer, for MPI-related stuff
  integer, parameter  :: psb_mpik_ = kind(1)
  !
  ! These must be the kind parameter corresponding to psb_mpi_r_dpk_
  ! and psb_mpi_r_spk_
  !
  integer(psb_mpik_), parameter  :: psb_spk_p_ = 6
  integer(psb_mpik_), parameter  :: psb_spk_r_ = 37
  integer(psb_mpik_), parameter  :: psb_spk_   = selected_real_kind(psb_spk_p_,psb_spk_r_)
  integer(psb_mpik_), parameter  :: psb_dpk_p_ = 15
  integer(psb_mpik_), parameter  :: psb_dpk_r_ = 307
  integer(psb_mpik_), parameter  :: psb_dpk_   = selected_real_kind(psb_dpk_p_,psb_dpk_r_)
#endif

  integer(psb_ipk_), save        :: psb_sizeof_dp, psb_sizeof_sp
  integer(psb_ipk_), save        :: psb_sizeof_int, psb_sizeof_long_int
  !
  ! Integer type identifiers for MPI operations. 
  !
  integer(psb_mpik_), save      :: psb_mpi_ipk_integer
  integer(psb_mpik_), save      :: psb_mpi_def_integer
  integer(psb_mpik_), save      :: psb_mpi_lng_integer
  integer(psb_mpik_), save      :: psb_mpi_r_spk_
  integer(psb_mpik_), save      :: psb_mpi_r_dpk_
  integer(psb_mpik_), save      :: psb_mpi_c_spk_
  integer(psb_mpik_), save      :: psb_mpi_c_dpk_
  ! 
  ! Version
  !
  character(len=*), parameter    :: psb_version_string_ = "3.4.0"
  integer(psb_ipk_), parameter   :: psb_version_major_  = 3
  integer(psb_ipk_), parameter   :: psb_version_minor_  = 4
  integer(psb_ipk_), parameter   :: psb_patchlevel_     = 0

  !
  !     Handy & miscellaneous constants
  !
  integer(psb_ipk_), parameter   :: izero=0, ione=1
  integer(psb_ipk_), parameter   :: itwo=2, ithree=3,mone=-1
  integer(psb_ipk_), parameter   :: psb_root_=0
  real(psb_spk_), parameter      :: szero=0.0_psb_spk_, sone=1.0_psb_spk_
  real(psb_dpk_), parameter      :: dzero=0.0_psb_dpk_, done=1.0_psb_dpk_
  complex(psb_spk_), parameter   :: czero=(0.0_psb_spk_,0.0_psb_spk_)
  complex(psb_spk_), parameter   :: cone=(1.0_psb_spk_,0.0_psb_spk_)
  complex(psb_dpk_), parameter   :: zzero=(0.0_psb_dpk_,0.0_psb_dpk_)
  complex(psb_dpk_), parameter   :: zone=(1.0_psb_dpk_,0.0_psb_dpk_)
  real(psb_dpk_), parameter      :: d_epstol=1.1e-16_psb_dpk_ ! Unit roundoff.  
  real(psb_spk_), parameter      :: s_epstol=5.e-8_psb_spk_   ! Is this right?
  character, parameter           :: psb_all_='A',  psb_topdef_=' '
  logical, parameter             :: psb_i_is_complex_ = .false.
  logical, parameter             :: psb_s_is_complex_ = .false.
  logical, parameter             :: psb_d_is_complex_ = .false.
  logical, parameter             :: psb_c_is_complex_ = .true.
  logical, parameter             :: psb_z_is_complex_ = .true.

  !
  ! Sort routines constants
  !
  ! 
  !  The up/down constant are defined in pairs having 
  !  opposite values. We make use of this fact in the heapsort routine.
  !
  integer(psb_ipk_), parameter :: psb_sort_up_      = 1, psb_sort_down_     = -1
  integer(psb_ipk_), parameter :: psb_lsort_up_     = 2, psb_lsort_down_    = -2
  integer(psb_ipk_), parameter :: psb_asort_up_     = 3, psb_asort_down_    = -3
  integer(psb_ipk_), parameter :: psb_alsort_up_    = 4, psb_alsort_down_   = -4
  integer(psb_ipk_), parameter :: psb_sort_ovw_idx_ = 0, psb_sort_keep_idx_ =  1
  integer(psb_ipk_), parameter :: psb_heap_resize   = 200


  !
  ! Sparse matrix constants
  !

  !
  ! State of matrices.
  !
  integer(psb_ipk_), parameter :: psb_invalid_ = -1 
  integer(psb_ipk_), parameter :: psb_spmat_null_=0, psb_spmat_bld_=1
  integer(psb_ipk_), parameter :: psb_spmat_asb_=2, psb_spmat_upd_=4

  integer(psb_ipk_), parameter :: psb_ireg_flgs_=10, psb_ip2_=0
  integer(psb_ipk_), parameter :: psb_iflag_=2, psb_ichk_=3
  integer(psb_ipk_), parameter :: psb_nnzt_=4, psb_zero_=5,psb_ipc_=6

  integer(psb_ipk_), parameter :: psb_unsorted_  = 0 
  integer(psb_ipk_), parameter :: psb_row_major_ = 1 
  integer(psb_ipk_), parameter :: psb_col_major_ = 2
  
  ! Duplicate coefficients handling
  ! These are usually set while calling spcnv as one of its
  ! optional arugments.
  integer(psb_ipk_), parameter :: psb_dupl_ovwrt_ = 0
  integer(psb_ipk_), parameter :: psb_dupl_add_   = 1
  integer(psb_ipk_), parameter :: psb_dupl_err_   = 2
  integer(psb_ipk_), parameter :: psb_dupl_def_   = psb_dupl_ovwrt_
  ! Matrix update mode
  integer(psb_ipk_), parameter :: psb_upd_srch_   = 98764
  integer(psb_ipk_), parameter :: psb_upd_perm_   = 98765
  integer(psb_ipk_), parameter :: psb_upd_dflt_   = psb_upd_srch_

#if defined(HAVE_ISO_FORTRAN_ENV) 
  integer(psb_ipk_), save :: psb_err_unit = error_unit  
  integer(psb_ipk_), save :: psb_inp_unit = input_unit  
  integer(psb_ipk_), save :: psb_out_unit = output_unit 
#else 
  integer(psb_ipk_), save :: psb_err_unit = 0
  integer(psb_ipk_), save :: psb_inp_unit = 5
  integer(psb_ipk_), save :: psb_out_unit = 6
#endif
  
  !
  ! 
  !     Error constants
  integer(psb_ipk_), parameter, public :: psb_success_=0
  integer(psb_ipk_), parameter, public :: psb_err_pivot_too_small_=2
  integer(psb_ipk_), parameter, public :: psb_err_invalid_ovr_num_=3
  integer(psb_ipk_), parameter, public :: psb_err_invalid_input_=5
  integer(psb_ipk_), parameter, public :: psb_err_iarg_neg_=10
  integer(psb_ipk_), parameter, public :: psb_err_iarg_pos_=20
  integer(psb_ipk_), parameter, public :: psb_err_input_value_invalid_i_=30
  integer(psb_ipk_), parameter, public :: psb_err_input_asize_invalid_i_=35
  integer(psb_ipk_), parameter, public :: psb_err_input_asize_small_i_=36
  integer(psb_ipk_), parameter, public :: psb_err_iarg_invalid_i_=40
  integer(psb_ipk_), parameter, public :: psb_err_iarg_not_gtia_ii_=50
  integer(psb_ipk_), parameter, public :: psb_err_iarg_not_gteia_ii_=60
  integer(psb_ipk_), parameter, public :: psb_err_iarg_invalid_value_=70
  integer(psb_ipk_), parameter, public :: psb_err_asb_nrc_error_=71
  integer(psb_ipk_), parameter, public :: psb_err_iarg2_neg_=80
  integer(psb_ipk_), parameter, public :: psb_err_ia2_not_increasing_=90
  integer(psb_ipk_), parameter, public :: psb_err_ia1_not_increasing_=91
  integer(psb_ipk_), parameter, public :: psb_err_ia1_badindices_=100
  integer(psb_ipk_), parameter, public :: psb_err_invalid_args_combination_=110
  integer(psb_ipk_), parameter, public :: psb_err_invalid_pid_arg_=115
  integer(psb_ipk_), parameter, public :: psb_err_iarg_n_mbgtian_=120
  integer(psb_ipk_), parameter, public :: psb_err_dupl_cd_vl=123
  integer(psb_ipk_), parameter, public :: psb_err_duplicate_coo=130
  integer(psb_ipk_), parameter, public :: psb_err_invalid_input_format_=134
  integer(psb_ipk_), parameter, public :: psb_err_unsupported_format_=135
  integer(psb_ipk_), parameter, public :: psb_err_format_unknown_=136
  integer(psb_ipk_), parameter, public :: psb_err_iarray_outside_bounds_=140
  integer(psb_ipk_), parameter, public :: psb_err_iarray_outside_process_=150
  integer(psb_ipk_), parameter, public :: psb_err_forgot_geall_=290
  integer(psb_ipk_), parameter, public :: psb_err_forgot_spall_=295
  integer(psb_ipk_), parameter, public :: psb_err_wrong_ins_=298
  integer(psb_ipk_), parameter, public :: psb_err_iarg_mbeeiarra_i_=300
  integer(psb_ipk_), parameter, public :: psb_err_mpi_error_=400
  integer(psb_ipk_), parameter, public :: psb_err_parm_differs_among_procs_=550
  integer(psb_ipk_), parameter, public :: psb_err_entry_out_of_bounds_=551
  integer(psb_ipk_), parameter, public :: psb_err_inconsistent_index_lists_=552
  integer(psb_ipk_), parameter, public :: psb_err_partfunc_toomuchprocs_=570
  integer(psb_ipk_), parameter, public :: psb_err_partfunc_toofewprocs_=575
  integer(psb_ipk_), parameter, public :: psb_err_partfunc_wrong_pid_=580
  integer(psb_ipk_), parameter, public :: psb_err_no_optional_arg_=581
  integer(psb_ipk_), parameter, public :: psb_err_arg_m_required_=582
  integer(psb_ipk_), parameter, public :: psb_err_many_optional_arg_=583
  integer(psb_ipk_), parameter, public :: psb_err_optional_arg_pair_=584
  integer(psb_ipk_), parameter, public :: psb_err_missing_override_method_=700
  integer(psb_ipk_), parameter, public :: psb_err_invalid_dynamic_type_=701
  integer(psb_ipk_), parameter, public :: psb_err_invalid_matrix_sizes_=1119
  integer(psb_ipk_), parameter, public :: psb_err_rectangular_mat_unsupported_=1120
  integer(psb_ipk_), parameter, public :: psb_err_invalid_mat_state_=1121
  integer(psb_ipk_), parameter, public :: psb_err_invalid_cd_state_=1122
  integer(psb_ipk_), parameter, public :: psb_err_invalid_a_and_cd_state_=1123
  integer(psb_ipk_), parameter, public :: psb_err_invalid_vect_state_=1124
  integer(psb_ipk_), parameter, public :: psb_err_context_error_=2010
  integer(psb_ipk_), parameter, public :: psb_err_initerror_neugh_procs_=2011
  integer(psb_ipk_), parameter, public :: psb_err_invalid_matrix_input_state_=2231
  integer(psb_ipk_), parameter, public :: psb_err_input_no_regen_=2232
  integer(psb_ipk_), parameter, public :: psb_err_lld_case_not_implemented_=3010
  integer(psb_ipk_), parameter, public :: psb_err_transpose_unsupported_=3015
  integer(psb_ipk_), parameter, public :: psb_err_transpose_c_unsupported_=3020
  integer(psb_ipk_), parameter, public :: psb_err_transpose_not_n_unsupported_=3021
  integer(psb_ipk_), parameter, public :: psb_err_only_unit_diag_=3022
  integer(psb_ipk_), parameter, public :: psb_err_ja_nix_ia_niy_unsupported_=3030
  integer(psb_ipk_), parameter, public :: psb_err_ix_n1_iy_n1_unsupported_=3040
  integer(psb_ipk_), parameter, public :: psb_err_input_matrix_unassembled_=3110
  integer(psb_ipk_), parameter, public :: psb_err_missing_aux_lib_=3999
  integer(psb_ipk_), parameter, public :: psb_err_alloc_dealloc_=4000
  integer(psb_ipk_), parameter, public :: psb_err_internal_error_=4001
  integer(psb_ipk_), parameter, public :: psb_err_from_subroutine_=4010
  integer(psb_ipk_), parameter, public :: psb_err_from_subroutine_i_=4012
  integer(psb_ipk_), parameter, public :: psb_err_from_subroutine_ai_=4013
  integer(psb_ipk_), parameter, public :: psb_err_alloc_request_=4025
  integer(psb_ipk_), parameter, public :: psb_err_from_subroutine_non_=4011
  integer(psb_ipk_), parameter, public :: psb_err_invalid_istop_=5001

end module psb_const_mod
