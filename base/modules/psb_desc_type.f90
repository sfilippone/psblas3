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
!
!	Module to   define desc_a,
!      structure for coomunications.
!
! Typedef: psb_desc_type
!    Defines a communication descriptor


module psb_descriptor_type
  use psb_const_mod

  implicit none

  !
  !     Communication, prolongation & restriction
  !
  integer, parameter :: psb_nohalo_=0,  psb_halo_=4
  integer, parameter :: psb_none_=0,  psb_sum_=1
  integer, parameter :: psb_avg_=2,  psb_square_root_=3
  integer, parameter :: psb_swap_send_=1, psb_swap_recv_=2
  integer, parameter :: psb_swap_sync_=4, psb_swap_mpi_=8

  integer, parameter :: psb_no_comm_=-1
  integer, parameter :: psb_comm_halo_=0, psb_comm_ovr_=1, psb_comm_ext_=2
  integer, parameter :: psb_ovt_xhal_ = 123, psb_ovt_asov_=psb_ovt_xhal_+1

  !
  !     Entries and values in desc%matrix_data
  !
  integer, parameter :: psb_dec_type_=1, psb_m_=2,psb_n_=3
  integer, parameter :: psb_n_row_=4,  psb_n_col_=5,psb_ctxt_=6
  integer, parameter :: psb_desc_size_=7
  integer, parameter :: psb_mpi_c_=9
  integer, parameter :: psb_thal_xch_=11
  integer, parameter :: psb_thal_snd_=12
  integer, parameter :: psb_thal_rcv_=13
  integer, parameter :: psb_tovr_xch_=14
  integer, parameter :: psb_tovr_snd_=15
  integer, parameter :: psb_tovr_rcv_=16
  integer, parameter :: psb_text_xch_=17
  integer, parameter :: psb_text_snd_=18
  integer, parameter :: psb_text_rcv_=19
  integer, parameter :: psb_mdata_size_=20
  integer, parameter :: psb_desc_asb_=3099
  integer, parameter :: psb_desc_bld_=psb_desc_asb_+1
  integer, parameter :: psb_desc_repl_=3199
  integer, parameter :: psb_desc_upd_=psb_desc_bld_+1
  integer, parameter :: psb_desc_normal_=3299
  integer, parameter :: psb_desc_large_=psb_desc_normal_+1
  integer, parameter :: psb_cd_ovl_bld_=3399
  integer, parameter :: psb_cd_ovl_asb_=psb_cd_ovl_bld_+1
  integer, parameter :: psb_hash_bits=14
  integer, parameter :: psb_hash_size=2**psb_hash_bits, psb_hash_mask=psb_hash_size-1
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
  
  
  !
  !  DESC data structure. 
  !  This is the most important data structure: it holds all the data 
  !  necessary to organize data exchange. The pattern of communication 
  !  among processes depends not only on the allocation of portions of 
  !  the index space to the various processes, but also on the undelying
  !  mesh discretization pattern. Thus building a communication descriptor is
  !  very much linked to building a sparse matrix (since the matrix sparsity 
  !  pattern embodies the topology of the discretization graph).
  !  
  !  Most general info about the descriptor is stored in the matrix_dataq
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
  !  The LOC_TO_GLOB, GLOB_TO_LOC, GLB_LC, HASHV and PTREE arrays implement the 
  !  mapping between local and global indices, according to the following 
  !  guidelines:
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
  !        array; thus we only record the global indices that do have a 
  !        local counterpart. Thus the local storage will be proportional to 
  !        N_COL. During the build phase we keep the known global indices in an 
  !        AVL tree data structure whose pointer is stored in ptree(:), so that we 
  !        can do both search and insertions in log time. At assembly time, we move 
  !        the information into hashv(:) and glb_lc(:,:). The idea is that 
  !        glb_lc(:,1) will hold sorted global indices, and glb_lc(:,2) the 
  !        corresponding local indices, so that we may do a binary search. To cut down
  !        the search time we partition glb_lc into a set of lists addressed by 
  !        hashv(:) based on the value of the lowest PSB_HASH_BITS bits of the 
  !        global index. 
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
  !  notation        stored in		          explanation
  !  --------------- --------------------------- -----------------------------------
  !  process_id      index_v(p+proc_id_)      identifier of process with which 
  !                                                data is  exchanged.
  !  n_elements_recv index_v(p+n_elem_recv_)  number of elements to receive.
  !  elements_recv   index_v(p+elem_recv_+i)  indexes of local elements to
  !					          receive. these are stored in the
  !					          array from location p+elem_recv_ to
  !					          location p+elem_recv_+
  !						  index_v(p+n_elem_recv_)-1.
  !  n_elements_send index_v(p+n_elem_send_)  number of elements to send.
  !  elements_send   index_v(p+elem_send_+i)  indexes of local elements to
  !					          send. these are stored in the
  !					          array from location p+elem_send_ to
  !					          location p+elem_send_+
  !						  index_v(p+n_elem_send_)-1.
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
  !     number of indices (most often N=1), list of local indices. 
  !     This is because we only know the list of halo indices to be received 
  !     as we go about building the sparse matrix pattern, and we want the build 
  !     phase to be loosely synchronized. Thus we record the indices we have to ask 
  !     for, and at the time we call PSB_CDASB we match all the requests to figure 
  !     out who should be sending what. 
  !     However this implies that we know who owns the indices; if we are in the 
  !     LARGE case (as described above) this is actually only true for the OVERLAP list 
  !     that is filled in at CDALL time, and not for the HALO; thus the HALO list 
  !     is rebuilt during the CDASB process (in the psi_ldsc_pre_halo subroutine). 
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
  !     of overlap, i.e. the number of processes "owning" them.
  ! 
  !
  !  Yes, it is complex, but it does the following:
  !  1. Allows a purely local matrix/stencil buildup phase, requiring only 
  !     one synch point at the end (CDASB)
  !  2. Takes shortcuts when the problem size is not too large (the default threshold
  !     assumes that you are willing to spend up to 16 MB on each process for the 
  !     glob_to_loc mapping)
  !  3. Supports restriction/prolongation operators with the same routines 
  !     just choosing (in the swapdata/swaptran internals) on which index list 
  !     they should work. 
  !
  !
  !

  type psb_desc_type
     integer, allocatable :: matrix_data(:)
     integer, allocatable :: halo_index(:), ext_index(:)
     integer, allocatable :: bnd_elem(:)
     integer, allocatable :: ovrlap_index(:)
     integer, allocatable :: ovrlap_elem(:)
     integer, allocatable :: loc_to_glob(:)
     integer, allocatable :: glob_to_loc (:)
     integer, allocatable :: hashv(:), glb_lc(:,:), ptree(:)
     integer, allocatable :: lprm(:)
     integer, allocatable :: idx_space(:)
  end type psb_desc_type


  integer, private, save :: cd_large_threshold=psb_default_large_threshold 


contains 

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
    integer             :: np,me, isz, err_act,idx,gidx,lidx

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

    psb_is_large_desc =(psb_desc_large_==psb_cd_get_size(desc))

  end function psb_is_large_desc

  logical function psb_is_upd_desc(desc)
    type(psb_desc_type), intent(in) :: desc

    psb_is_upd_desc = psb_is_upd_dec(psb_cd_get_dectype(desc))

  end function psb_is_upd_desc


  logical function psb_is_asb_desc(desc)
    type(psb_desc_type), intent(in) :: desc

    psb_is_asb_desc = psb_is_asb_dec(psb_cd_get_dectype(desc))

  end function psb_is_asb_desc


  logical function psb_is_ok_dec(dectype)
    integer :: dectype

    psb_is_ok_dec = ((dectype == psb_desc_asb_).or.(dectype == psb_desc_bld_).or.&
         &(dectype == psb_desc_upd_).or.&
         &(dectype== psb_desc_repl_))
  end function psb_is_ok_dec

  logical function psb_is_bld_dec(dectype)
    integer :: dectype

    psb_is_bld_dec = (dectype == psb_desc_bld_)
  end function psb_is_bld_dec

  logical function psb_is_upd_dec(dectype)          
    integer :: dectype

    psb_is_upd_dec = (dectype == psb_desc_upd_)

  end function psb_is_upd_dec


  logical function psb_is_asb_dec(dectype)          
    integer :: dectype

    psb_is_asb_dec = (dectype == psb_desc_asb_).or.&
         & (dectype== psb_desc_repl_)

  end function psb_is_asb_dec


  integer function psb_cd_get_local_rows(desc)
    type(psb_desc_type), intent(in) :: desc
    
    psb_cd_get_local_rows = desc%matrix_data(psb_n_row_)
  end function psb_cd_get_local_rows

  integer function psb_cd_get_local_cols(desc)
    type(psb_desc_type), intent(in) :: desc
    
    psb_cd_get_local_cols = desc%matrix_data(psb_n_col_)
  end function psb_cd_get_local_cols

  integer function psb_cd_get_global_rows(desc)
    type(psb_desc_type), intent(in) :: desc
    
    psb_cd_get_global_rows = desc%matrix_data(psb_m_)
  end function psb_cd_get_global_rows

  integer function psb_cd_get_global_cols(desc)
    type(psb_desc_type), intent(in) :: desc
    
    psb_cd_get_global_cols = desc%matrix_data(psb_n_)
  end function psb_cd_get_global_cols

  integer function psb_cd_get_context(desc)
    type(psb_desc_type), intent(in) :: desc
    
    psb_cd_get_context = desc%matrix_data(psb_ctxt_)
  end function psb_cd_get_context

  integer function psb_cd_get_dectype(desc)
    type(psb_desc_type), intent(in) :: desc
    
    psb_cd_get_dectype = desc%matrix_data(psb_dec_type_)
  end function psb_cd_get_dectype

  integer function psb_cd_get_size(desc)
    type(psb_desc_type), intent(in) :: desc
    
    psb_cd_get_size = desc%matrix_data(psb_desc_size_)
  end function psb_cd_get_size

  integer function psb_cd_get_mpic(desc)
    type(psb_desc_type), intent(in) :: desc
    
    psb_cd_get_mpic = desc%matrix_data(psb_mpi_c_)
  end function psb_cd_get_mpic


  subroutine psb_cd_set_bld(desc,info)
    !
    ! Change state of a descriptor into BUILD. 
    ! If the descriptor is LARGE, check the  AVL search tree
    ! and initialize it if necessary.
    !
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod

    implicit none
    type(psb_desc_type), intent(inout) :: desc
    integer                            :: info
    !locals
    integer             :: np,me,ictxt, isz, err_act,idx,gidx,lidx
    logical, parameter  :: debug=.false.,debugprt=.false.
    character(len=20)   :: name, char_err
    if (debug) write(0,*) me,'Entered CDCPY'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    name = 'psb_cd_set_bld'

    ictxt = psb_cd_get_context(desc)

    ! check on blacs grid 
    call psb_info(ictxt, me, np)
    if (debug) write(0,*) me,'Entered CDCPY'

    if (psb_is_large_desc(desc)) then 
      if (.not.allocated(desc%ptree)) then 
        allocate(desc%ptree(2),stat=info)
        if (info /= 0) then 
          info=4000
          goto 9999
        endif
        call InitPairSearchTree(desc%ptree,info)
        do idx=1, psb_cd_get_local_cols(desc)
          gidx = desc%loc_to_glob(idx)
          call SearchInsKeyVal(desc%ptree,gidx,idx,lidx,info)        
          if (lidx /= idx) then 
            write(0,*) 'Warning from cdset: mismatch in PTREE ',idx,lidx
          endif
        enddo
      end if
    end if
    desc%matrix_data(psb_dec_type_) = psb_desc_bld_ 

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
  end subroutine psb_cd_set_bld
    
end module psb_descriptor_type
