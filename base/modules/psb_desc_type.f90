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
! package: psb_descriptor_type
!    Defines a communication descriptor
!

module psb_descriptor_type
  use psb_const_mod
  use psb_hash_mod 

  implicit none

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
  integer, parameter :: psb_desc_asb_=3099
  integer, parameter :: psb_desc_bld_=psb_desc_asb_+1
  integer, parameter :: psb_desc_repl_=3199
  integer, parameter :: psb_desc_upd_=psb_desc_bld_+1
  ! these two are reserved for descriptors which are
  ! "overlap-extensions" of base descriptors. 
  integer, parameter :: psb_cd_ovl_bld_=3399
  integer, parameter :: psb_cd_ovl_asb_=psb_cd_ovl_bld_+1
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


  !
  !  type: psb_desc_type
  !  
  !  Communication Descriptor data structure.
  !| type psb_idxmap_type
  !|   integer              :: state 
  !|   integer, allocatable :: loc_to_glob(:)
  !|   integer, allocatable :: glob_to_loc(:)
  !|   integer              :: hashvsize, hashvmask
  !|   integer, allocatable :: hashv(:), glb_lc(:,:)
  !|   type(psb_hash_type)  :: hash
  !| end type psb_idxmap_type
  !
  !
  !|  type psb_desc_type
  !|     integer, allocatable  :: matrix_data(:)
  !|     integer, allocatable  :: halo_index(:), ext_index(:)
  !|     integer, allocatable  :: bnd_elem(:)
  !|     integer, allocatable  :: ovrlap_index(:)
  !|     integer, allocatable  :: ovrlap_elem(:,:)
  !|     integer, allocatable  :: ovr_mst_idx(:)
  !|     type(psb_idxmap_type) :: idxmap
  !|     integer, allocatable  :: lprm(:)
  !|     integer, allocatable  :: idx_space(:)
  !|     type(psb_desc_type), pointer :: base_desc => null()
  !|  end type psb_desc_type
  !
  !  
  !  This is the most important data structure: it holds all the data 
  !  necessary to organize data exchange. The pattern of communication 
  !  among processes depends not only on the allocation of portions of 
  !  the index space to the various processes, but also on the underlying
  !  mesh discretization pattern. Thus building a communication descriptor is
  !  very much linked to building a sparse matrix (since the matrix sparsity 
  !  pattern embodies the topology of the discretization graph).
  !  
  !  Most general info about the descriptor is stored in the matrix_data
  !  component, including the STATE which can be  PSB_DESC_BLD_, 
  !  PSB_DESC_ASB_ or PSB_DESC_REPL_. 
  !  Upon allocation with PSB_CDALL the descriptor enters the BLD state;
  !  then the user can specify the discretization pattern with PSB_CDINS;
  !  the call to PSB_CDASB puts the descriptor in the PSB_ASB_ state. 
  !
  !  PSB_DESC_REPL_ is a special value that specifies a replicated index space,
  !  and is only entered by the psb_cdrep call. Currently it is only 
  !  used in the last level of some multilevel preconditioners. 
  ! 
  !  The LOC_TO_GLOB, GLOB_TO_LOC, GLB_LC, HASHV and HASH  data structures
  !  inside IDXMAP  implement the  mapping between local and global indices,
  !  according to the following   guidelines:
  !
  !  1. Each global index I is owned by at least one process;
  !
  !  2. On each process, indices from 1 to N_ROW (desc%matrix_dat(psb_n_row_))
  !     are locally owned; the value of N_ROW can be determined upon allocation 
  !     based on the index distribution (see also the interface to CDALL).
  !
  !  3. If a global index is owned by more than one process, we have an OVERLAP
  !     in which case the sum of all the N_ROW values is greater than the total 
  !     size of the index space; 
  !
  !  4. During the buildup of the descriptor, according to the user specified 
  !     stencil, we also take notice of indices that are not owned by the current
  !     process, but whose value is needed to proceed with the computation; these 
  !     form the HALO of the current process. Halo indices are assigned local indices
  !     from N_ROW+1 to N_COL (inclusive).
  !
  !  5. Regardless of the descriptor state, LOC_TO_GLOB(I), I=1:N_COL always 
  !     contains the global index corresponding to local index I; the upper bound 
  !     N_COL moves during the descriptor build process (see CDINS). 
  !
  !  6. The descriptor also contains the inverse global-to-local mapping. This 
  !     mapping can take two forms according to the value returned by 
  !     psb_cd_choose_large_state:
  !     i. If the global index space size is not too large, it is possible to store 
  !        a complete mapping  in GLOB_TO_LOC: each entry contains the corresponding 
  !        local index (if there is one), or an encoded value identifying the process 
  !        owning that index. This array is filled in at initialization time CDALL, 
  !        and thus it is available throughout the insertion phase. The local storage 
  !        will thus be N_COL + N_GLOB
  !   ii.  If the global index space is very large (larger than the threshold value
  !        which may be set by the user), then it is not advisable to have such an 
  !        array.
  !        In this case we only record the global indices that do have a 
  !        local counterpart, so that the local storage will be proportional to 
  !        N_COL.
  !        The idea is that  glb_lc(:,1) will hold sorted global indices, and
  !        glb_lc(:,2) the corresponding local indices, so that we may do a binary search.
  !        To cut down  the search time we partition glb_lc into a set of lists
  !        addressed by  hashv(:) based on the value of the lowest
  !        PSB_HASH_BITS bits of the  global index. 
  !        During the build phase glb_lc() will store the indices of the internal points,
  !        i.e. local indices 1:NROW, since those are known ad CDALL time.
  !        The halo indices that we encounter during the build phase are put in
  !        a PSB_HASH_TYPE data structure, which implements a very simple hash; this
  !        hash  will nonetheless be quite fast at low occupancy rates.
  !        At assembly time, we move everything into hashv(:) and glb_lc(:,:).
  !
  !  7. The data exchange is based on lists of local indices to be exchanged; all the 
  !     lists have the same format, as follows:
  !     the list is  composed of variable dimension blocks, one for each process to 
  !     communicate with; each block contains indices of local elements to be 
  !     exchanged. We do choose the order of communications:  do not change 
  !     the sequence of blocks unless you know what you're doing, or you'll 
  !     risk a deadlock. NOTE: This is the format when the state is PSB_ASB_.
  !     See below for BLD. The end-of-list is marked with a -1. 
  !
  !|  notation        stored in                   explanation
  !|  --------------- --------------------------- -----------------------------------
  !|  process_id      index_v(p+proc_id_)      identifier of process with which 
  !|                                                data is  exchanged.
  !|  n_elements_recv index_v(p+n_elem_recv_)  number of elements to receive.
  !|  elements_recv   index_v(p+elem_recv_+i)  indexes of local elements to
  !|                                              receive. these are stored in the
  !|                                              array from location p+elem_recv_ to
  !|                                              location p+elem_recv_+
  !|                                              index_v(p+n_elem_recv_)-1.
  !|  n_elements_send index_v(p+n_elem_send_)  number of elements to send.
  !|  elements_send   index_v(p+elem_send_+i)  indexes of local elements to
  !|                                              send. these are stored in the
  !|                                              array from location p+elem_send_ to
  !|                                              location p+elem_send_+
  !|                                              index_v(p+n_elem_send_)-1.
  !
  !     This organization is valid for both halo and overlap indices; overlap entries
  !     need to be updated to ensure that a variable at a given global index 
  !     (assigned to multiple processes) has the same value. The way to resolve the 
  !     issue is to exchange the data and then sum (or average) the values. See
  !     psb_ovrl subroutine. 
  !  
  !  8. When the descriptor is in the BLD state the INDEX vectors contains only 
  !     the indices to be received, organized as  a sequence 
  !     of entries of the form (proc,N,(lx1,lx2,...,lxn)) with owning process,
  !     number of indices (most often but not necessarily N=1), list of local indices.  
  !     This is because we only know the list of halo indices to be received 
  !     as we go about building the sparse matrix pattern, and we want the build 
  !     phase to be loosely synchronized. Thus we record the indices we have to ask 
  !     for, and at the time we call PSB_CDASB we match all the requests to figure 
  !     out who should be sending what to whom.
  !     However this implies that we know who owns the indices; if we are in the 
  !     LARGE case (as described above) this is actually only true for the OVERLAP list 
  !     that is filled in at CDALL time, and not for the HALO (remember: we do not have 
  !     the space to encode the owning process index in the GLOB_TO_LOC mapping); thus 
  !     the HALO list is rebuilt during the CDASB process 
  !     (in the psi_ldsc_pre_halo subroutine). 
  !  
  !  9. Yet another twist comes about when building an extended descriptor with 
  !     the psb_cdbldext subroutine. In this case we are reaching out 
  !     layer by layer, but we may use the results in two ways:
  !      i. To build a descriptor with the same "owned" indices, but with an 
  !         extended halo, with additional layers; in this case the requests 
  !         go into halo_index;
  !     ii. To build a descriptor suitable for overlapped Schwarz-type computations.
  !         In this case we end up with more "owned" indices than in the base 
  !         descriptor, so that what was a halo index in the base becomes an overlap
  !         index in the extended descriptor. In this case we build three lists: 
  !         ovrlap_index    the indices that overlap
  !         halo_index      the halo indices (of the extended descriptor)
  !         ext_index       the indices of elements that need to be gathered to 
  !                         map the original index space onto the new (overlapped) 
  !                         index space. 
  !     So, ext_index has the same format as the others, but is only used in the 
  !     context of Schwarz-type computations; otherwise it is empty (i.e. 
  !     it only contains the end-of-list marker -1). 
  !
  ! 10. ovrlap_elem contains a list of overlap indices together with their degree
  !     of overlap, i.e. the number of processes "owning" the, and the "master"
  !     process whose value has to be considered authoritative when the need arises.
  !     
  ! 11. ovr_mst_idx is a list defining a retrieve of a copy of the values for
  !     overlap entries from their respecitve "master" processes by means of
  !     an halo exchange call. This is used for those cases where there is
  !     an overlap in the base data distribution.
  !
  ! It is complex, but it does the following:
  !  1. Allows a purely local matrix/stencil buildup phase, requiring only 
  !     one synch point at the end (CDASB)
  !  2. Takes shortcuts when the problem size is not too large (the default threshold
  !     assumes that you are willing to spend up to 4 MB on each process for the 
  !     glob_to_loc mapping)
  !  3. Supports restriction/prolongation operators with the same routines 
  !     just choosing (in the swapdata/swaptran internals) on which index list 
  !     they should work. 
  !
  !
  !
  type psb_idxmap_type
    integer              :: state 
    integer, allocatable :: loc_to_glob(:)
    integer, allocatable :: glob_to_loc(:)
    integer              :: hashvsize, hashvmask
    integer, allocatable :: hashv(:), glb_lc(:,:)
    type(psb_hash_type)  :: hash
  end type psb_idxmap_type

  type psb_desc_type
    integer, allocatable  :: matrix_data(:)
    integer, allocatable  :: halo_index(:)
    integer, allocatable  :: ext_index(:)
    integer, allocatable  :: ovrlap_index(:)
    integer, allocatable  :: ovrlap_elem(:,:)
    integer, allocatable  :: ovr_mst_idx(:)
    integer, allocatable  :: bnd_elem(:)
    type(psb_idxmap_type) :: idxmap 
    integer, allocatable  :: lprm(:)
    type(psb_desc_type), pointer     :: base_desc => null()
    integer, allocatable :: idx_space(:)
  end type psb_desc_type

  interface psb_sizeof
    module procedure psb_cd_sizeof, psb_idxmap_sizeof
  end interface

  interface psb_is_ok_desc
    module procedure psb_is_ok_desc
  end interface

  interface psb_is_asb_desc
    module procedure psb_is_asb_desc
  end interface

  interface psb_is_upd_desc
    module procedure psb_is_upd_desc
  end interface

  interface psb_is_ovl_desc
    module procedure psb_is_ovl_desc
  end interface

  interface psb_is_bld_desc
    module procedure psb_is_bld_desc
  end interface

  interface psb_is_large_desc
    module procedure psb_is_large_desc
  end interface


  interface psb_move_alloc
    module procedure psb_cdtransfer, psb_idxmap_transfer
  end interface


  interface psb_free
    module procedure psb_cdfree, psb_idxmap_free
  end interface

  interface psb_map_l2g
    module procedure psb_map_l2g_s1, psb_map_l2g_s2,&
         & psb_map_l2g_v1, psb_map_l2g_v2
  end interface

  integer, private, save :: cd_large_threshold=psb_default_large_threshold 


contains 

  function psb_idxmap_sizeof(map)  result(val)
    implicit none
    !....Parameters...

    Type(psb_idxmap_type), intent(in) :: map
    integer(psb_long_int_k_) :: val

    val = 3*psb_sizeof_int
    if (allocated(map%loc_to_glob))  val = val + psb_sizeof_int*size(map%loc_to_glob)
    if (allocated(map%glob_to_loc))  val = val + psb_sizeof_int*size(map%glob_to_loc)
    if (allocated(map%hashv))        val = val + psb_sizeof_int*size(map%hashv)
    if (allocated(map%glb_lc))       val = val + psb_sizeof_int*size(map%glb_lc)
    val = val + psb_sizeof(map%hash) 

  end function psb_idxmap_sizeof
    

  function psb_cd_sizeof(desc)  result(val)
    implicit none
    !....Parameters...

    Type(psb_desc_type), intent(in) :: desc
    integer(psb_long_int_k_) :: val

    val = 0 
    if (allocated(desc%matrix_data))  val = val + psb_sizeof_int*size(desc%matrix_data)
    if (allocated(desc%halo_index))   val = val + psb_sizeof_int*size(desc%halo_index)
    if (allocated(desc%ext_index))    val = val + psb_sizeof_int*size(desc%ext_index)
    if (allocated(desc%bnd_elem))     val = val + psb_sizeof_int*size(desc%bnd_elem)
    if (allocated(desc%ovrlap_index)) val = val + psb_sizeof_int*size(desc%ovrlap_index)
    if (allocated(desc%ovrlap_elem))  val = val + psb_sizeof_int*size(desc%ovrlap_elem)
    if (allocated(desc%ovr_mst_idx))  val = val + psb_sizeof_int*size(desc%ovr_mst_idx)
    if (allocated(desc%lprm))         val = val + psb_sizeof_int*size(desc%lprm)
    if (allocated(desc%idx_space))    val = val + psb_sizeof_int*size(desc%idx_space)
    val = val + psb_sizeof(desc%idxmap)

  end function psb_cd_sizeof



  subroutine psb_cd_set_large_threshold(ith)
    integer, intent(in) :: ith
    if (ith > 0) then 
      cd_large_threshold = ith
    end if
  end subroutine psb_cd_set_large_threshold

  integer function  psb_cd_get_large_threshold()
    psb_cd_get_large_threshold = cd_large_threshold 
  end function psb_cd_get_large_threshold

  logical  function  psb_cd_choose_large_state(ictxt,m)
    use psb_penv_mod

    implicit none
    integer, intent(in) :: ictxt,m
    !locals
    integer             :: np,me

    call psb_info(ictxt, me, np)
    ! 
    ! Since the hashed lists take up (somewhat) more than 2*N_COL integers,
    ! it makes no sense to use them if you don't have at least 
    ! 3 processes, no matter what the size of the process. 
    !
    psb_cd_choose_large_state = &
         & (m > psb_cd_get_large_threshold()) .and. &
         & (np > 2)
  end function psb_cd_choose_large_state

  subroutine psb_nullify_desc(desc)
    type(psb_desc_type), intent(inout) :: desc
    ! We have nothing left to do here.
    ! Perhaps we should delete this subroutine? 
    nullify(desc%base_desc)

  end subroutine psb_nullify_desc

  logical function psb_is_ok_desc(desc)

    type(psb_desc_type), intent(in) :: desc

    psb_is_ok_desc = psb_is_ok_dec(psb_cd_get_dectype(desc))

  end function psb_is_ok_desc

  logical function psb_is_bld_desc(desc)
    type(psb_desc_type), intent(in) :: desc

    psb_is_bld_desc = psb_is_bld_dec(psb_cd_get_dectype(desc))

  end function psb_is_bld_desc

  logical function psb_is_large_desc(desc)
    type(psb_desc_type), intent(in) :: desc

    psb_is_large_desc =(psb_desc_large_ == psb_cd_get_size(desc))

  end function psb_is_large_desc

  logical function psb_is_upd_desc(desc)
    type(psb_desc_type), intent(in) :: desc

    psb_is_upd_desc = psb_is_upd_dec(psb_cd_get_dectype(desc))

  end function psb_is_upd_desc

  logical function psb_is_repl_desc(desc)
    type(psb_desc_type), intent(in) :: desc

    psb_is_repl_desc = psb_is_repl_dec(psb_cd_get_dectype(desc))

  end function psb_is_repl_desc

  logical function psb_is_ovl_desc(desc)
    type(psb_desc_type), intent(in) :: desc

    psb_is_ovl_desc = psb_is_ovl_dec(psb_cd_get_dectype(desc))

  end function psb_is_ovl_desc


  logical function psb_is_asb_desc(desc)
    type(psb_desc_type), intent(in) :: desc

    psb_is_asb_desc = psb_is_asb_dec(psb_cd_get_dectype(desc))

  end function psb_is_asb_desc


  logical function psb_is_ok_dec(dectype)
    integer :: dectype

    psb_is_ok_dec = ((dectype == psb_desc_asb_).or.(dectype == psb_desc_bld_).or.&
         &(dectype == psb_cd_ovl_asb_).or.(dectype == psb_cd_ovl_bld_).or.&
         &(dectype == psb_desc_upd_).or.&
         &(dectype == psb_desc_repl_))
  end function psb_is_ok_dec

  logical function psb_is_bld_dec(dectype)
    integer :: dectype

    psb_is_bld_dec = (dectype == psb_desc_bld_).or.(dectype == psb_cd_ovl_bld_)
  end function psb_is_bld_dec

  logical function psb_is_upd_dec(dectype)          
    integer :: dectype

    psb_is_upd_dec = (dectype == psb_desc_upd_)

  end function psb_is_upd_dec

  logical function psb_is_repl_dec(dectype)          
    integer :: dectype

    psb_is_repl_dec = (dectype == psb_desc_repl_)

  end function psb_is_repl_dec


  logical function psb_is_asb_dec(dectype)          
    integer :: dectype

    psb_is_asb_dec = (dectype == psb_desc_asb_).or.&
         & (dectype == psb_desc_repl_).or.(dectype == psb_cd_ovl_asb_)

  end function psb_is_asb_dec

  logical function psb_is_ovl_dec(dectype)          
    integer :: dectype

    psb_is_ovl_dec = (dectype == psb_cd_ovl_bld_).or.&
         & (dectype == psb_cd_ovl_asb_)

  end function psb_is_ovl_dec


  integer function psb_cd_get_local_rows(desc)
    type(psb_desc_type), intent(in) :: desc

    if (psb_is_ok_desc(desc)) then 
      psb_cd_get_local_rows = desc%matrix_data(psb_n_row_)
    else
      psb_cd_get_local_rows = -1
    endif
  end function psb_cd_get_local_rows

  integer function psb_cd_get_local_cols(desc)
    type(psb_desc_type), intent(in) :: desc

    if (psb_is_ok_desc(desc)) then 
      psb_cd_get_local_cols = desc%matrix_data(psb_n_col_)
    else
      psb_cd_get_local_cols = -1
    endif
  end function psb_cd_get_local_cols

  integer function psb_cd_get_global_rows(desc)
    type(psb_desc_type), intent(in) :: desc

    if (psb_is_ok_desc(desc)) then 
      psb_cd_get_global_rows = desc%matrix_data(psb_m_)
    else
      psb_cd_get_global_rows = -1
    endif

  end function psb_cd_get_global_rows

  integer function psb_cd_get_global_cols(desc)
    type(psb_desc_type), intent(in) :: desc

    if (psb_is_ok_desc(desc)) then 
      psb_cd_get_global_cols = desc%matrix_data(psb_n_)
    else
      psb_cd_get_global_cols = -1
    endif

  end function psb_cd_get_global_cols

  integer function psb_cd_get_context(desc)
    use psb_error_mod
    type(psb_desc_type), intent(in) :: desc
    if (allocated(desc%matrix_data)) then 
      psb_cd_get_context = desc%matrix_data(psb_ctxt_)
    else
      psb_cd_get_context = -1
      call psb_errpush(psb_err_invalid_cd_state_,'psb_cd_get_context')
      call psb_error()
    end if
  end function psb_cd_get_context

  integer function psb_cd_get_dectype(desc)
    use psb_error_mod
    type(psb_desc_type), intent(in) :: desc

    if (allocated(desc%matrix_data)) then 
      psb_cd_get_dectype = desc%matrix_data(psb_dec_type_)
    else
      psb_cd_get_dectype = -1
      call psb_errpush(psb_err_invalid_cd_state_,'psb_cd_get_dectype')
      call psb_error()
    end if
      
  end function psb_cd_get_dectype

  integer function psb_cd_get_size(desc)
    use psb_error_mod
    type(psb_desc_type), intent(in) :: desc

    if (allocated(desc%matrix_data)) then 
      psb_cd_get_size = desc%idxmap%state
    else
      psb_cd_get_size = -1
      call psb_errpush(psb_err_invalid_cd_state_,'psb_cd_get_size')
      call psb_error()
    end if

  end function psb_cd_get_size

  integer function psb_cd_get_mpic(desc)
    use psb_error_mod
    type(psb_desc_type), intent(in) :: desc

    if (allocated(desc%matrix_data)) then 
      psb_cd_get_mpic = desc%matrix_data(psb_mpi_c_)
    else
      psb_cd_get_mpic = -1
      call psb_errpush(psb_err_invalid_cd_state_,'psb_cd_get_mpic')
      call psb_error()
    end if

  end function psb_cd_get_mpic


  subroutine psb_cd_set_ovl_asb(desc,info)
    !
    ! Change state of a descriptor into ovl_build. 
    implicit none
    type(psb_desc_type), intent(inout) :: desc
    integer                            :: info


    if (psb_is_asb_desc(desc)) desc%matrix_data(psb_dec_type_) = psb_cd_ovl_asb_ 

  end subroutine psb_cd_set_ovl_asb


  subroutine psb_get_xch_idx(idx,totxch,totsnd,totrcv)
    implicit none 
    integer, intent(in)  :: idx(:)
    integer, intent(out) :: totxch,totsnd,totrcv

    integer :: ip, nerv, nesd
    character(len=20), parameter  :: name='psb_get_xch_idx'    

    totxch = 0
    totsnd = 0
    totrcv = 0
    ip     = 1

    do 
      if (ip > size(idx)) then 
        write(0,*) trim(name),': Warning: out of size of input vector '
        exit
      end if
      if (idx(ip) == -1) exit
      totxch = totxch+1
      nerv   = idx(ip+psb_n_elem_recv_)
      nesd   = idx(ip+nerv+psb_n_elem_send_)
      totsnd = totsnd + nesd
      totrcv = totrcv + nerv
      ip     = ip+nerv+nesd+3
    end do

  end subroutine psb_get_xch_idx



  subroutine psb_cd_get_list(data,desc,ipnt,totxch,idxr,idxs,info)
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod

    implicit none
    integer, intent(in)          :: data
    integer, pointer             :: ipnt(:)
    type(psb_desc_type), target  :: desc
    integer, intent(out)         :: totxch,idxr,idxs,info

    !locals
    integer             :: np,me,ictxt,err_act, debug_level,debug_unit
    logical, parameter  :: debug=.false.,debugprt=.false.
    character(len=20), parameter  :: name='psb_cd_get_list'

    info = psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    ictxt = psb_cd_get_context(desc)

    call psb_info(ictxt, me, np)

    select case(data) 
    case(psb_comm_halo_) 
      ipnt   => desc%halo_index
    case(psb_comm_ovr_) 
      ipnt   => desc%ovrlap_index
    case(psb_comm_ext_) 
      ipnt   => desc%ext_index
      if (debug_level >= psb_debug_ext_) then
        if (.not.associated(desc%base_desc)) then
          write(debug_unit,*) trim(name),&
               & ': Warning: trying to get ext_index on a descriptor ',&
               & 'which does not have a base_desc!'
        end if
        if (.not.psb_is_ovl_desc(desc)) then
          write(debug_unit,*) trim(name),&
               & ': Warning: trying to get ext_index on a descriptor ',&
               & 'which is not overlap-extended!'
        end if
      end if
    case(psb_comm_mov_) 
      ipnt   => desc%ovr_mst_idx
    case default
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='wrong Data selector')
      goto 9999
    end select
    call psb_get_xch_idx(ipnt,totxch,idxs,idxr)


    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_ret_) then
      return
    else
      call psb_error(ictxt)
    end if
    return
  end subroutine psb_cd_get_list

  subroutine psb_idxmap_free(map,info)
    !...free descriptor structure...
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none
    !....parameters...
    type(psb_idxmap_type), intent(inout) :: map
    integer, intent(out)               :: info
    !...locals....
    integer             :: ictxt,np,me, err_act
    character(len=*), parameter ::  name = 'psb_idxmap_free'

    if(psb_get_errstatus() /= 0) return 
    info=psb_success_
    call psb_erractionsave(err_act)

    if (allocated(map%loc_to_glob)) then 
      deallocate(map%loc_to_glob,stat=info) 
    end if
    if ((info == psb_success_).and.allocated(map%glob_to_loc)) then 
      deallocate(map%glob_to_loc,stat=info) 
    end if
    if ((info == psb_success_).and.allocated(map%hashv)) then 
      deallocate(map%hashv,stat=info) 
    end if
    if ((info == psb_success_).and.allocated(map%glb_lc)) then 
      deallocate(map%glb_lc,stat=info) 
    end if
    if (info /= psb_success_) call psb_free(map%hash, info) 
    if (info /= psb_success_) then 
      info=2052
      call psb_errpush(info,name)
      goto 9999
    end if
    

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_ret_) then
      return
    else
      call psb_error(ictxt)
    end if
    return

  end subroutine psb_idxmap_free


  !
  ! Subroutine: psb_cdfree
  !   Frees a descriptor data structure.
  ! 
  ! Arguments: 
  !    desc_a   - type(psb_desc_type).         The communication descriptor to be freed.
  !    info     - integer.                       return code.
  subroutine psb_cdfree(desc_a,info)
    !...free descriptor structure...
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none
    !....parameters...
    type(psb_desc_type), intent(inout) :: desc_a
    integer, intent(out)               :: info
    !...locals....
    integer             :: ictxt,np,me, err_act
    character(len=20)   :: name

    if(psb_get_errstatus() /= 0) return 
    info=psb_success_
    call psb_erractionsave(err_act)
    name = 'psb_cdfree'


    if (.not.allocated(desc_a%matrix_data)) then
      info=psb_err_forgot_spall_
      call psb_errpush(info,name)
      return
    end if

    ictxt=psb_cd_get_context(desc_a)

    call psb_info(ictxt, me, np)
    !     ....verify blacs grid correctness..
    if (np == -1) then
      info = psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    endif

    call psb_free(desc_a%idxmap,info)

    if (.not.allocated(desc_a%halo_index)) then
      info=298
      call psb_errpush(info,name)
      goto 9999
    end if

    !deallocate halo_index field
    deallocate(desc_a%halo_index,stat=info)
    if (info /= psb_success_) then
      info=2053
      call psb_errpush(info,name)
      goto 9999
    end if

    if (.not.allocated(desc_a%bnd_elem)) then
!!$    info=296
!!$    call psb_errpush(info,name)
!!$    goto 9999
!!$  end if
    else
      !deallocate halo_index field
      deallocate(desc_a%bnd_elem,stat=info)
      if (info /= psb_success_) then
        info=2054
        call psb_errpush(info,name)
        goto 9999
      end if
    end if

    if (.not.allocated(desc_a%ovrlap_index)) then
      info=299
      call psb_errpush(info,name)
      goto 9999
    end if

    !deallocate ovrlap_index  field
    deallocate(desc_a%ovrlap_index,stat=info)
    if (info /= psb_success_) then
      info=2055
      call psb_errpush(info,name)
      goto 9999
    end if

    !deallocate ovrlap_elem  field
    deallocate(desc_a%ovrlap_elem,stat=info)
    if (info /= psb_success_) then 
      info=2056
      call psb_errpush(info,name)
      goto 9999
    end if

    !deallocate ovrlap_index  field
    deallocate(desc_a%ovr_mst_idx,stat=info)
    if (info /= psb_success_) then
      info=2055
      call psb_errpush(info,name)
      goto 9999
    end if


    deallocate(desc_a%lprm,stat=info)
    if (info /= psb_success_) then 
      info=2057
      call psb_errpush(info,name)
      goto 9999
    end if

    if (allocated(desc_a%idx_space)) then 
      deallocate(desc_a%idx_space,stat=info)
      if (info /= psb_success_) then 
        info=2056
        call psb_errpush(info,name)
        goto 9999
      end if
    end if

    deallocate(desc_a%matrix_data)

    call psb_nullify_desc(desc_a)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_ret_) then
      return
    else
      call psb_error(ictxt)
    end if
    return

  end subroutine psb_cdfree
  !
  ! Subroutine: psb_cdtransfer
  !   Transfers data and allocation from in to out; behaves like MOVE_ALLOC, i.e.
  !   the IN arg is empty (and deallocated) upon exit. 
  !
  ! 
  ! Arguments: 
  !    desc_in  - type(psb_desc_type).         The communication descriptor to be 
  !                                               transferred.
  !    desc_out - type(psb_desc_type).         The output communication descriptor.
  !    info     - integer.                       Return code.
  subroutine psb_cdtransfer(desc_in, desc_out, info)

    use psb_realloc_mod
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod

    implicit none
    !....parameters...

    type(psb_desc_type), intent(inout)  :: desc_in
    type(psb_desc_type), intent(inout)  :: desc_out
    integer, intent(out)                :: info

    !locals
    integer             :: np,me,ictxt, err_act
    integer             :: debug_level, debug_unit
    character(len=20)   :: name

    if (psb_get_errstatus() /= 0) return 
    info = psb_success_
    call psb_erractionsave(err_act)
    name = 'psb_cdtransfer'
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    ! Should not require ictxt to be present: this
    ! function might be called even when desc_in is
    ! empty. 

    call psb_move_alloc( desc_in%matrix_data ,    desc_out%matrix_data  , info)
    if (info == psb_success_)  &
         & call psb_move_alloc( desc_in%halo_index  ,    desc_out%halo_index   , info)
    if (info == psb_success_)  &
         & call psb_move_alloc( desc_in%bnd_elem    ,    desc_out%bnd_elem     , info)
    if (info == psb_success_)  &
         & call psb_move_alloc( desc_in%ovrlap_elem ,    desc_out%ovrlap_elem  , info)
    if (info == psb_success_)  &
         & call psb_move_alloc( desc_in%ovrlap_index,    desc_out%ovrlap_index , info)
    if (info == psb_success_)  &
         & call psb_move_alloc( desc_in%ovr_mst_idx ,    desc_out%ovr_mst_idx  , info)
    if (info == psb_success_)  &
         & call psb_move_alloc( desc_in%ext_index   ,    desc_out%ext_index    , info)
    if (info == psb_success_)  &
         & call psb_move_alloc( desc_in%lprm        ,    desc_out%lprm         , info)
    if (info == psb_success_)  &
         & call psb_move_alloc( desc_in%idx_space   ,    desc_out%idx_space    , info)
    if (info == psb_success_) &
         & call psb_move_alloc(desc_in%idxmap, desc_out%idxmap,info)
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info,name)
      goto 9999
    endif
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),': end'

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  end subroutine psb_cdtransfer

  subroutine psb_idxmap_transfer(map_in, map_out, info)

    use psb_realloc_mod
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod

    implicit none
    !....parameters...

    type(psb_idxmap_type), intent(inout)  :: map_in
    type(psb_idxmap_type), intent(inout)  :: map_out
    integer, intent(out)                :: info

    !locals
    integer             :: np,me,ictxt, err_act
    integer              :: debug_level, debug_unit
    character(len=*), parameter  ::  name = 'psb_idxmap_transfer'

    if (psb_get_errstatus() /= 0) return 
    info = psb_success_
    call psb_erractionsave(err_act)

    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    map_out%state     = map_in%state
    map_out%hashvsize = map_in%hashvsize
    map_out%hashvmask = map_in%hashvmask

    if (info == psb_success_)  &
         & call psb_move_alloc( map_in%loc_to_glob ,    map_out%loc_to_glob  , info)
    if (info == psb_success_)  &
         & call psb_move_alloc( map_in%glob_to_loc ,    map_out%glob_to_loc  , info)
    if (info == psb_success_)  &
         & call psb_move_alloc( map_in%hashv       ,    map_out%hashv        , info)
    if (info == psb_success_)  &
         & call psb_move_alloc( map_in%glb_lc      ,    map_out%glb_lc       , info)
    if (info == psb_success_)  &
         & call psb_move_alloc( map_in%hash        ,    map_out%hash        , info)
    
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info,name)
      goto 9999
    endif
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),': end'

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  end subroutine psb_idxmap_transfer

  subroutine psb_idxmap_copy(map_in, map_out, info)

    use psb_realloc_mod
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod

    implicit none
    !....parameters...

    type(psb_idxmap_type), intent(in)    :: map_in
    type(psb_idxmap_type), intent(inout) :: map_out
    integer, intent(out)                 :: info

    !locals
    integer             :: np,me,ictxt, err_act
    integer              :: debug_level, debug_unit
    character(len=*), parameter  ::  name = 'psb_idxmap_transfer'

    if (psb_get_errstatus() /= 0) return 
    info = psb_success_
    call psb_erractionsave(err_act)

    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    map_out%state     = map_in%state
    map_out%hashvsize = map_in%hashvsize
    map_out%hashvmask = map_in%hashvmask

    call psb_safe_ab_cpy( map_in%loc_to_glob ,    map_out%loc_to_glob  , info)
    if (info == psb_success_)  &
         & call psb_safe_ab_cpy( map_in%glob_to_loc ,    map_out%glob_to_loc  , info)
    if (info == psb_success_)  &
         & call psb_safe_ab_cpy( map_in%hashv       ,    map_out%hashv        , info)
    if (info == psb_success_)  &
         & call psb_safe_ab_cpy( map_in%glb_lc      ,    map_out%glb_lc       , info)
    if (info == psb_success_)  &
         & call psb_hash_copy( map_in%hash        ,    map_out%hash        , info)
    
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info,name)
      goto 9999
    endif
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),': end'

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  end subroutine psb_idxmap_copy

  subroutine psb_map_l2g_s1(idx,map,info)
    implicit none 
    integer, intent(inout) :: idx
    integer, intent(out)   :: info 
    type(psb_idxmap_type)  :: map
    integer :: nc

    info = psb_success_
    if (.not.allocated(map%loc_to_glob)) then 
      info = psb_err_iarray_outside_bounds_
      idx = -1 
      return
    end if
    nc = size(map%loc_to_glob) 
    if ((idx < 1).or.(idx>nc)) then 
      info = psb_err_iarray_outside_bounds_
      idx = -1 
      return
    end if
    idx = map%loc_to_glob(idx) 
    
  end subroutine psb_map_l2g_s1

  subroutine psb_map_l2g_s2(idx,gidx,map,info)
    implicit none 
    integer, intent(in)   :: idx
    integer, intent(out)  :: gidx, info 
    type(psb_idxmap_type) :: map
    integer :: nc

    info = psb_success_
    if (.not.allocated(map%loc_to_glob)) then 
      info = psb_err_iarray_outside_bounds_
      gidx = -1 
      return
    end if
    nc = size(map%loc_to_glob) 
    if ((idx < 1).or.(idx>nc)) then 
      info = psb_err_iarray_outside_bounds_
      gidx = -1 
      return
    end if
    gidx = map%loc_to_glob(idx) 
    
  end subroutine psb_map_l2g_s2

  subroutine psb_map_l2g_v1(idx,map,info)
    implicit none 
    integer, intent(inout) :: idx(:)
    integer, intent(out)   :: info 
    type(psb_idxmap_type)  :: map
    integer :: nc, i, ix

    info = psb_success_
    if (size(idx) == 0) return
    if (.not.allocated(map%loc_to_glob)) then 
      info = psb_err_iarray_outside_bounds_
      idx = -1 
      return
    end if
    nc = size(map%loc_to_glob) 
    do i=1, size(idx) 
      ix = idx(i)
      if ((ix < 1).or.(ix>nc)) then 
        info = psb_err_iarray_outside_bounds_
        idx(i) = -1 
      else        
        idx(i) = map%loc_to_glob(ix) 
      end if
    end do
    
  end subroutine psb_map_l2g_v1

  subroutine psb_map_l2g_v2(idx,gidx,map,info)
    implicit none 
    integer, intent(in)    :: idx(:)
    integer, intent(out)   :: gidx(:),info 
    type(psb_idxmap_type)  :: map
    integer :: nc, i, ix

    info = psb_success_
    if (size(idx) == 0) return
    if ((.not.allocated(map%loc_to_glob)).or.&
         & (size(gidx)<size(idx))) then 
      info = psb_err_iarray_outside_bounds_
      gidx = -1 
      return
    end if
      
    nc = size(map%loc_to_glob) 
    do i=1, size(idx) 
      ix = idx(i)
      if ((ix < 1).or.(ix>nc)) then 
        info = psb_err_iarray_outside_bounds_
        gidx(i) = -1 
      else        
        gidx(i) = map%loc_to_glob(ix) 
      end if
    end do
    
  end subroutine psb_map_l2g_v2


  Subroutine psb_cd_get_recv_idx(tmp,desc,data,info,toglob)

    use psb_error_mod
    use psb_penv_mod
    use psb_realloc_mod
    Implicit None
    integer, allocatable, intent(out)       :: tmp(:)
    integer, intent(in)                     :: data
    Type(psb_desc_type), Intent(in), target :: desc
    integer, intent(out)                    :: info
    logical, intent(in)                     :: toglob

    !     .. Local Scalars ..
    Integer ::  incnt, outcnt, j, np, me, ictxt, l_tmp,&
         & idx, gidx, proc, n_elem_send, n_elem_recv
    Integer, pointer   :: idxlist(:) 
    integer              :: debug_level, debug_unit, err_act
    character(len=20)    :: name

    name  = 'psb_cd_get_recv_idx'
    info  = psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    ictxt = psb_cd_get_context(desc)
    call psb_info(ictxt, me, np)
    
    select case(data)
    case(psb_comm_halo_)
      idxlist => desc%halo_index
    case(psb_comm_ovr_)
      idxlist => desc%ovrlap_index
    case(psb_comm_ext_)
      idxlist => desc%ext_index
    case(psb_comm_mov_)
      idxlist => desc%ovr_mst_idx
      write(0,*) 'Warning: unusual request getidx on ovr_mst_idx'
    case default
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='wrong Data selector')
      goto 9999
    end select
    
    l_tmp = 3*size(idxlist)

    allocate(tmp(l_tmp),stat=info)
    if (info /= psb_success_) then 
      info = psb_err_from_subroutine_
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999      
    end if
      
    incnt  = 1
    outcnt = 1
    tmp(:) = -1
    Do While (idxlist(incnt) /= -1)
      proc        = idxlist(incnt+psb_proc_id_)
      n_elem_recv = idxlist(incnt+psb_n_elem_recv_)
      n_elem_send = idxlist(incnt+n_elem_recv+psb_n_elem_send_)

      Do j=0,n_elem_recv-1
        idx = idxlist(incnt+psb_elem_recv_+j)
        call psb_ensure_size((outcnt+3),tmp,info,pad=-1)
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='psb_ensure_size')
          goto 9999
        end if
        if (toglob) then
          call psb_map_l2g(idx,gidx,desc%idxmap,info)
          If (gidx < 0) then 
            info=-3
            call psb_errpush(info,name)
            goto 9999
          endif
          tmp(outcnt)   = proc
          tmp(outcnt+1) = 1
          tmp(outcnt+2) = gidx
          tmp(outcnt+3) = -1
        else
          tmp(outcnt)   = proc
          tmp(outcnt+1) = 1
          tmp(outcnt+2) = idx
          tmp(outcnt+3) = -1
        end if
        outcnt          = outcnt+3
      end Do
      incnt = incnt+n_elem_recv+n_elem_send+3
    end Do
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    Return

  end Subroutine psb_cd_get_recv_idx



end module psb_descriptor_type
