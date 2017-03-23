!   
!                Parallel Sparse BLAS  version 3.4
!      (C) Copyright 2006, 2010, 2015
!                         Salvatore Filippone    University of Rome Tor Vergata
!                         Alfredo Buttari        CNRS-IRIT, Toulouse
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!    
!
!
! package: psb_desc_mod
!    Defines a communication descriptor
!

module psb_desc_mod
  use psb_const_mod
  use psb_desc_const_mod
  use psb_indx_map_mod
  use psb_i_vect_mod

  implicit none

  !
  !  type: psb_desc_type
  !  
  !  Communication Descriptor data structure.
  !
  !|  type psb_desc_type
  !|     class(psb_indx_map), allocatable :: indxmap
  !|     integer(psb_ipk_), allocatable  :: halo_index(:), ext_index(:)
  !|     integer(psb_ipk_), allocatable  :: bnd_elem(:)
  !|     integer(psb_ipk_), allocatable  :: ovrlap_index(:)
  !|     integer(psb_ipk_), allocatable  :: ovrlap_elem(:,:)
  !|     integer(psb_ipk_), allocatable  :: ovr_mst_idx(:)
  !|     integer(psb_ipk_), allocatable  :: lprm(:)
  !|     integer(psb_ipk_), allocatable  :: idx_space(:)
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
  !  This is a two-level data structure: it combines an INDX_MAP with
  !  a set of auxiliary lists.
  !  For a complete description of INDX_MAP see its own file, but the
  !  idea here is the following: the INDX_MAP contains information about
  !  the index space and its allocation to the various processors.
  !  In particular, besides the communicator, it contains the data relevant
  !  to the following queries:
  !  1. How many global rows/columns?
  !  2. How many local  rows/columns?
  !  3. Convert between local and global indices
  !  4. Add to local indices.
  !  5. Find (one of) the owner(s) of a given index
  !  Checking for the existence of overlap is very expensive, thus
  !  it is done at build time (for extended-halo cases it can be inferred from
  !  the construction process).
  !  There are multiple ways to represent an INDX_MAP internally, hence it is
  !  a CLASS variable, which can take different forms, more or less memory hungry. 
  !
  !  Guidelines
  !
  !  1. Each global index I is owned by at least one process;
  !
  !  2. On each process, indices from 1 to N_ROW (desc%indxmap%get_lr())
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
  !  5. The upper bound  N_COL moves during the descriptor build process (see CDINS). 
  !
  !  6. The descriptor also contains the inverse global-to-local mapping.
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
  !     However this implies that we know who owns the indices; 
  !     this is actually only true for the OVERLAP list 
  !     that is filled in at CDALL time, and not for the HALO (remember: we do not
  !     necessarily have the space to encode the owning process index); thus 
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
  !  2. Takes shortcuts when the problem size is not too large
  !  
  !  3. Supports restriction/prolongation operators with the same routines 
  !     just choosing (in the swapdata/swaptran internals) on which index list 
  !     they should work. 
  !
  !
  !


  type psb_desc_type
    class(psb_indx_map), allocatable :: indxmap

    integer(psb_ipk_), allocatable   :: halo_index(:)
    integer(psb_ipk_), allocatable   :: ext_index(:)
    integer(psb_ipk_), allocatable   :: ovrlap_index(:)
    integer(psb_ipk_), allocatable   :: ovr_mst_idx(:)

    type(psb_i_vect_type)            :: v_halo_index
    type(psb_i_vect_type)            :: v_ext_index
    type(psb_i_vect_type)            :: v_ovrlap_index
    type(psb_i_vect_type)            :: v_ovr_mst_idx 

    integer(psb_ipk_), allocatable   :: ovrlap_elem(:,:)
    integer(psb_ipk_), allocatable   :: bnd_elem(:)
    integer(psb_ipk_), allocatable   :: lprm(:)
    type(psb_desc_type), pointer     :: base_desc => null()
    integer(psb_ipk_), allocatable   :: idx_space(:)
  contains
    procedure, pass(desc) :: is_ok           => psb_is_ok_desc
    procedure, pass(desc) :: is_valid        => psb_is_valid_desc
    procedure, pass(desc) :: is_upd          => psb_is_upd_desc
    procedure, pass(desc) :: is_bld          => psb_is_bld_desc
    procedure, pass(desc) :: is_asb          => psb_is_asb_desc
    procedure, pass(desc) :: is_ovl          => psb_is_ovl_desc
    procedure, pass(desc) :: is_repl         => psb_is_repl_desc
    procedure, pass(desc) :: get_mpic        => psb_cd_get_mpic
    procedure, pass(desc) :: get_dectype     => psb_cd_get_dectype
    procedure, pass(desc) :: get_context     => psb_cd_get_context    
    procedure, pass(desc) :: get_ctxt        => psb_cd_get_context    
    procedure, pass(desc) :: get_local_rows  => psb_cd_get_local_rows
    procedure, pass(desc) :: get_local_cols  => psb_cd_get_local_cols
    procedure, pass(desc) :: get_global_rows => psb_cd_get_global_rows
    procedure, pass(desc) :: get_global_cols => psb_cd_get_global_cols
    procedure, pass(desc) :: get_global_indices => psb_cd_get_global_indices
    procedure, pass(desc) :: a_get_list      => psb_cd_get_list
    procedure, pass(desc) :: v_get_list      => psb_cd_v_get_list
    generic, public       :: get_list => a_get_list, v_get_list
    procedure, pass(desc) :: sizeof          => psb_cd_sizeof
    procedure, pass(desc) :: clone           => psb_cd_clone
    procedure, pass(desc) :: cnv             => psb_cd_cnv
    procedure, pass(desc) :: free            => psb_cdfree
    procedure, pass(desc) :: destroy         => psb_cd_destroy
    procedure, pass(desc) :: nullify         => nullify_desc

    procedure, pass(desc) :: get_fmt         => cd_get_fmt
    procedure, pass(desc) :: fnd_owner       => cd_fnd_owner
    procedure, pass(desc) :: l2gs1           => cd_l2gs1
    procedure, pass(desc) :: l2gs2           => cd_l2gs2
    procedure, pass(desc) :: l2gv1           => cd_l2gv1
    procedure, pass(desc) :: l2gv2           => cd_l2gv2
    generic, public       :: l2g             => l2gs2, l2gv2
    generic, public       :: l2gip           => l2gs1, l2gv1                            

    procedure, pass(desc) :: g2ls1           => cd_g2ls1
    procedure, pass(desc) :: g2ls2           => cd_g2ls2
    procedure, pass(desc) :: g2lv1           => cd_g2lv1
    procedure, pass(desc) :: g2lv2           => cd_g2lv2
    generic, public       :: g2l             => g2ls2, g2lv2
    generic, public       :: g2lip           => g2ls1, g2lv1

    procedure, pass(desc) :: g2ls1_ins       => cd_g2ls1_ins
    procedure, pass(desc) :: g2ls2_ins       => cd_g2ls2_ins
    procedure, pass(desc) :: g2lv1_ins       => cd_g2lv1_ins
    procedure, pass(desc) :: g2lv2_ins       => cd_g2lv2_ins
    generic, public       :: g2l_ins         => g2ls2_ins, g2lv2_ins
    generic, public       :: g2lip_ins       => g2ls1_ins, g2lv1_ins
    

  end type psb_desc_type

      
  interface psb_sizeof
    module procedure psb_cd_sizeof
  end interface psb_sizeof

  interface psb_move_alloc
    module procedure psb_cdtransfer
  end interface psb_move_alloc

  interface psb_free
    module procedure psb_cdfree
  end interface psb_free


  private :: nullify_desc, cd_get_fmt,&
       & cd_l2gs1, cd_l2gs2, cd_l2gv1, cd_l2gv2, cd_g2ls1,&
       & cd_g2ls2, cd_g2lv1, cd_g2lv2, cd_g2ls1_ins,&
       & cd_g2ls2_ins, cd_g2lv1_ins, cd_g2lv2_ins, cd_fnd_owner


  integer(psb_ipk_), private, save :: cd_large_threshold=psb_default_large_threshold 


contains 

  function psb_cd_sizeof(desc)  result(val)
    implicit none
    !....Parameters...

    class(psb_desc_type), intent(in) :: desc
    integer(psb_long_int_k_) :: val

    val = 0 
    val = val + psb_sizeof_int*psb_size(desc%halo_index)
    val = val + psb_sizeof_int*psb_size(desc%ext_index)
    val = val + psb_sizeof_int*psb_size(desc%bnd_elem)
    val = val + psb_sizeof_int*psb_size(desc%ovrlap_index)
    val = val + psb_sizeof_int*psb_size(desc%ovrlap_elem)
    val = val + psb_sizeof_int*psb_size(desc%ovr_mst_idx)
    val = val + psb_sizeof_int*psb_size(desc%lprm)
    val = val + psb_sizeof_int*psb_size(desc%idx_space)
    if (allocated(desc%indxmap))  val = val + desc%indxmap%sizeof()
    val = val + desc%v_halo_index%sizeof()
    val = val + desc%v_ext_index%sizeof()
    val = val + desc%v_ovrlap_index%sizeof()
    val = val + desc%v_ovr_mst_idx%sizeof()

  end function psb_cd_sizeof



  subroutine psb_cd_set_large_threshold(ith)
    implicit none 
    integer(psb_ipk_), intent(in) :: ith
    if (ith > 0) then 
      cd_large_threshold = ith
    end if
  end subroutine psb_cd_set_large_threshold

  function  psb_cd_get_large_threshold() result(val)
    implicit none 
    integer(psb_ipk_) :: val
    val  = cd_large_threshold 
  end function psb_cd_get_large_threshold

  logical  function  psb_cd_choose_large_state(ictxt,m)
    use psb_penv_mod

    implicit none
    integer(psb_ipk_), intent(in) :: ictxt,m
    !locals
    integer(psb_ipk_) :: np,me

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
    implicit none 
    type(psb_desc_type), intent(inout) :: desc
    ! We have nothing left to do here.
    ! Perhaps we should delete this subroutine? 
    nullify(desc%base_desc)

  end subroutine psb_nullify_desc

  subroutine nullify_desc(desc)
    implicit none 
    class(psb_desc_type), intent(inout) :: desc
    ! We have nothing left to do here.
    ! Perhaps we should delete this subroutine? 
    nullify(desc%base_desc)

  end subroutine nullify_desc

  function psb_is_ok_desc(desc) result(val)
    implicit none 
    class(psb_desc_type), intent(in) :: desc
    logical                         :: val 
    
    val = .false.
    if (allocated(desc%indxmap)) &
         & val = desc%indxmap%is_valid()

  end function psb_is_ok_desc

  function psb_is_valid_desc(desc) result(val)
    implicit none 
    class(psb_desc_type), intent(in) :: desc
    logical                         :: val 
    
    val = .false.
    if (allocated(desc%indxmap)) &
         & val = desc%indxmap%is_valid()

  end function psb_is_valid_desc

  function psb_is_bld_desc(desc) result(val)
    implicit none 
    class(psb_desc_type), intent(in) :: desc
    logical                         :: val 

    val = .false.
    if (allocated(desc%indxmap)) &
         & val = desc%indxmap%is_bld()

  end function psb_is_bld_desc

  function psb_is_upd_desc(desc)  result(val)
    implicit none 
    class(psb_desc_type), intent(in) :: desc
    logical                         :: val 

    val = .false.
    if (allocated(desc%indxmap)) &
         & val = desc%indxmap%is_upd()

  end function psb_is_upd_desc

  function psb_is_repl_desc(desc) result(val)
    implicit none 
    class(psb_desc_type), intent(in) :: desc
    logical                         :: val 

    val = .false.
    if (allocated(desc%indxmap)) &
         & val = desc%indxmap%is_repl()

  end function psb_is_repl_desc

  function psb_is_ovl_desc(desc) result(val)
    implicit none 
    class(psb_desc_type), intent(in) :: desc
    logical                         :: val 

    val = .false.
    if (allocated(desc%indxmap)) &
         & val = desc%indxmap%is_ovl()

  end function psb_is_ovl_desc


  function psb_is_asb_desc(desc) result(val)
    implicit none 
    class(psb_desc_type), intent(in) :: desc
    logical                         :: val 

    val = .false.
    if (allocated(desc%indxmap)) &
         & val = desc%indxmap%is_asb()

  end function psb_is_asb_desc

  function psb_cd_get_local_rows(desc) result(val)
    implicit none 
    integer(psb_ipk_) :: val 
    class(psb_desc_type), intent(in) :: desc

    if (allocated(desc%indxmap)) then 
      val = desc%indxmap%get_lr()
    else
      val = -1
    endif
  end function psb_cd_get_local_rows

  function psb_cd_get_local_cols(desc) result(val)
    implicit none 
    integer(psb_ipk_) :: val 
    class(psb_desc_type), intent(in) :: desc

    if (allocated(desc%indxmap)) then 
      val = desc%indxmap%get_lc()
    else
      val = -1
    endif
  end function psb_cd_get_local_cols

  function psb_cd_get_global_rows(desc) result(val)
    implicit none 
    integer(psb_ipk_) :: val 
    class(psb_desc_type), intent(in) :: desc

    if (allocated(desc%indxmap)) then 
      val = desc%indxmap%get_gr()
    else
      val = -1
    endif

  end function psb_cd_get_global_rows

  function psb_cd_get_global_cols(desc) result(val)
    implicit none 
    integer(psb_ipk_) :: val 
    class(psb_desc_type), intent(in) :: desc

    if (allocated(desc%indxmap)) then 
      val = desc%indxmap%get_gc()
    else
      val = -1
    endif

  end function psb_cd_get_global_cols

  function psb_cd_get_global_indices(desc,owned) result(val)
    implicit none 
    integer(psb_ipk_), allocatable   :: val(:)
    class(psb_desc_type), intent(in) :: desc
    logical, intent(in), optional    :: owned

    logical :: owned_
    integer(psb_ipk_) :: i, nr, info
    
    owned_=.true.
    if (present(owned)) owned_=owned

    
    if (allocated(desc%indxmap)) then 
      if (owned_) then
        nr = desc%get_local_rows()
      else
        nr = desc%get_local_cols()
      end if
      nr = max(nr,0)
      allocate(val(nr))
      do i=1, nr
        val(i) = i
      end do
      call desc%l2gip(val,info,owned=owned_)
    else
      allocate(val(0))
    endif

  end function psb_cd_get_global_indices



  
  function cd_get_fmt(desc) result(val)
    implicit none 
    character(len=5) :: val 
    class(psb_desc_type), intent(in) :: desc

    if (allocated(desc%indxmap)) then 
      val = desc%indxmap%get_fmt()
    else
      val = 'NULL'
    endif

  end function cd_get_fmt

  function psb_cd_get_context(desc) result(val)
    use psb_error_mod
    implicit none 
    integer(psb_ipk_) :: val 
    class(psb_desc_type), intent(in) :: desc
    if (allocated(desc%indxmap)) then
      val = desc%indxmap%get_ctxt()    
    else
      val = -1
      call psb_errpush(psb_err_invalid_cd_state_,'psb_cd_get_context')
      call psb_error()
    end if
  end function psb_cd_get_context

  function psb_cd_get_dectype(desc) result(val)
    use psb_error_mod
    implicit none 
    integer(psb_ipk_) :: val 
    class(psb_desc_type), intent(in) :: desc

    if (allocated(desc%indxmap)) then
      val = desc%indxmap%get_state()    
    else
      val = -1
      call psb_errpush(psb_err_invalid_cd_state_,'psb_cd_get_dectype')
      call psb_error()
    end if

  end function psb_cd_get_dectype

  function psb_cd_get_mpic(desc) result(val)
    use psb_error_mod
    implicit none 
    integer(psb_ipk_) :: val 
    class(psb_desc_type), intent(in) :: desc

    if (allocated(desc%indxmap)) then
      val = desc%indxmap%get_mpic()    
    else
      val = -1
      call psb_errpush(psb_err_invalid_cd_state_,'psb_cd_get_mpic')
      call psb_error()
    end if

  end function psb_cd_get_mpic


  subroutine psb_cd_set_ovl_asb(desc,info)
    !
    ! Change state of a descriptor into ovl_build. 
    implicit none
    type(psb_desc_type), intent(inout) :: desc
    integer(psb_ipk_) :: info

    info = 0
    if (psb_is_asb_desc(desc)) &
         & call desc%indxmap%set_state(psb_desc_ovl_asb_)

  end subroutine psb_cd_set_ovl_asb


  subroutine psb_get_xch_idx(idx,totxch,totsnd,totrcv)
    implicit none 
    integer(psb_ipk_), intent(in)  :: idx(:)
    integer(psb_ipk_), intent(out) :: totxch,totsnd,totrcv

    integer(psb_ipk_) :: ip, nerv, nesd
    character(len=20), parameter  :: name='psb_get_xch_idx'    

    totxch = 0
    totsnd = 0
    totrcv = 0
    ip     = 1

    do 
      if (ip > size(idx)) then 
        write(psb_err_unit,*) trim(name),': Warning: out of size of input vector '
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


  subroutine psb_get_v_xch_idx(idx,totxch,totsnd,totrcv)
    implicit none 
    class(psb_i_base_vect_type), intent(in)  :: idx
    integer(psb_ipk_), intent(out) :: totxch,totsnd,totrcv

    integer(psb_ipk_) :: ip, nerv, nesd
    character(len=20), parameter  :: name='psb_get_v_xch_idx'    

    call psb_get_xch_idx(idx%v,totxch,totsnd,totrcv)

  end subroutine psb_get_v_xch_idx



  subroutine psb_cd_get_list(data,desc,ipnt,totxch,idxr,idxs,info)
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod

    implicit none
    integer(psb_ipk_), intent(in)          :: data
    integer(psb_ipk_), pointer             :: ipnt(:)
    class(psb_desc_type), target  :: desc
    integer(psb_ipk_), intent(out)         :: totxch,idxr,idxs,info

    !locals
    integer(psb_ipk_) :: np,me,ictxt,err_act, debug_level,debug_unit
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


9999 call psb_error_handler(err_act)

    return

  end subroutine psb_cd_get_list


  subroutine psb_cd_v_get_list(data,desc,ipnt,totxch,idxr,idxs,info)
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none
    integer(psb_ipk_), intent(in)          :: data
    class(psb_i_base_vect_type), pointer   :: ipnt
    class(psb_desc_type), target  :: desc
    integer(psb_ipk_), intent(out)         :: totxch,idxr,idxs,info

    !locals
    integer(psb_ipk_) :: np,me,ictxt,err_act, debug_level,debug_unit
    logical, parameter  :: debug=.false.,debugprt=.false.
    character(len=20), parameter  :: name='psb_cd_v_get_list'

    info = psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    ictxt = psb_cd_get_context(desc)

    call psb_info(ictxt, me, np)

    select case(data) 
    case(psb_comm_halo_) 
      ipnt   => desc%v_halo_index%v
      if (.not.allocated(desc%v_halo_index%v)) &
           & info = psb_err_inconsistent_index_lists_
    case(psb_comm_ovr_) 
      ipnt   => desc%v_ovrlap_index%v
      if (.not.allocated(desc%v_ovrlap_index%v)) &
           & info = psb_err_inconsistent_index_lists_
    case(psb_comm_ext_) 
      ipnt   => desc%v_ext_index%v
      if (.not.allocated(desc%v_ext_index%v)) &
           & info = psb_err_inconsistent_index_lists_
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
      ipnt   => desc%v_ovr_mst_idx%v
      if (.not.allocated(desc%v_ovr_mst_idx%v)) &
           & info = psb_err_inconsistent_index_lists_
      
    case default
      info=psb_err_from_subroutine_
    end select
    if (info /= psb_success_) then
      call psb_errpush(info,name,a_err='wrong Data selector')
      goto 9999
    end if
    
    call psb_get_v_xch_idx(ipnt,totxch,idxs,idxr)


    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_cd_v_get_list

  !
  ! Subroutine: psb_cdfree
  !   Frees a descriptor data structure.
  ! 
  ! Arguments: 
  !    desc_a   - type(psb_desc_type).         The communication descriptor to be freed.
  !    info     - integer.                       return code.
  subroutine psb_cdfree(desc,info)
    !...free descriptor structure...
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none
    !....parameters...
    class(psb_desc_type), intent(inout) :: desc
    integer(psb_ipk_), intent(out)      :: info
    !...locals....
    integer(psb_ipk_) :: ictxt,np,me, err_act
    character(len=20)   :: name

    info=psb_success_
    call psb_erractionsave(err_act)
    name = 'psb_cdfree'

    call desc%destroy()

    call psb_erractionrestore(err_act)
    return


9999 call psb_error_handler(err_act)
    return

  end subroutine psb_cdfree

  !
  ! Subroutine: psb_cdfree
  !   Frees a descriptor data structure.
  ! 
  ! Arguments: 
  !    desc_a   - type(psb_desc_type).         The communication descriptor to be freed.
  subroutine psb_cd_destroy(desc)
    !...free descriptor structure...
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none
    !....parameters...
    class(psb_desc_type), intent(inout) :: desc
    !...locals....
    integer(psb_ipk_) :: info


    if (allocated(desc%halo_index)) &
         &  deallocate(desc%halo_index,stat=info)

    if (allocated(desc%bnd_elem)) &
         &    deallocate(desc%bnd_elem,stat=info)

    if (allocated(desc%ovrlap_index)) &
         & deallocate(desc%ovrlap_index,stat=info)
    
    if (allocated(desc%ovrlap_elem)) &
         & deallocate(desc%ovrlap_elem,stat=info)
    if (allocated(desc%ovr_mst_idx)) &
         & deallocate(desc%ovr_mst_idx,stat=info)

    if (allocated(desc%lprm)) &
         & deallocate(desc%lprm,stat=info)
    if (allocated(desc%idx_space)) &
         & deallocate(desc%idx_space,stat=info)

    if (allocated(desc%indxmap)) then 
      call desc%indxmap%free()
      deallocate(desc%indxmap, stat=info)
    end if
    call desc%v_halo_index%free(info)
    call desc%v_ovrlap_index%free(info)
    call desc%v_ext_index%free(info)
    call desc%v_ovr_mst_idx%free(info)

    call desc%nullify()

    return

  end subroutine psb_cd_destroy
  !
  ! Subroutine: psb_cdtransfer
  !   Transfers data and allocation from in to out; behaves like MOVE_ALLOC, i.e.
  !   the IN arg is empty (and deallocated) upon exit. 
  !
  ! 
  ! Arguments: 
  !    desc  - type(psb_desc_type).         The communication descriptor to be 
  !                                               transferred.
  !    desc_out - type(psb_desc_type).         The output communication descriptor.
  !    info     - integer.                       Return code.
  subroutine psb_cdtransfer(desc, desc_out, info)

    use psb_realloc_mod
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod

    implicit none
    !....parameters...

    type(psb_desc_type), intent(inout)  :: desc
    type(psb_desc_type), intent(inout)  :: desc_out
    integer(psb_ipk_), intent(out)                :: info

    !locals
    integer(psb_ipk_) :: np,me,ictxt, err_act
    integer(psb_ipk_) :: debug_level, debug_unit
    character(len=20)   :: name

    if (psb_get_errstatus() /= 0) return 
    info = psb_success_
    call psb_erractionsave(err_act)
    name = 'psb_cdtransfer'
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    !
    ! Note: this  might be called even
    ! when desc is empty.
    ! 
    if (desc%is_valid()) then 
      ictxt = psb_cd_get_context(desc)
      call psb_info(ictxt,me,np)

      if (info == psb_success_)  &
           & call psb_move_alloc( desc%halo_index  ,    desc_out%halo_index   , info)
      if (info == psb_success_)  &
           & call psb_move_alloc( desc%bnd_elem    ,    desc_out%bnd_elem     , info)
      if (info == psb_success_)  &
           & call psb_move_alloc( desc%ovrlap_elem ,    desc_out%ovrlap_elem  , info)
      if (info == psb_success_)  &
           & call psb_move_alloc( desc%ovrlap_index,    desc_out%ovrlap_index , info)
      if (info == psb_success_)  &
           & call psb_move_alloc( desc%ovr_mst_idx ,    desc_out%ovr_mst_idx  , info)
      if (info == psb_success_)  &
           & call psb_move_alloc( desc%ext_index   ,    desc_out%ext_index    , info)
      if (info == psb_success_)  &
           & call psb_move_alloc( desc%lprm        ,    desc_out%lprm         , info)
      if (info == psb_success_)  &
           & call psb_move_alloc( desc%idx_space   ,    desc_out%idx_space    , info)
      if (info == psb_success_) &
           & call move_alloc(desc%indxmap, desc_out%indxmap)
      if (info == psb_success_) &
           & call desc%v_halo_index%clone(desc_out%v_halo_index,info)
      if (info == psb_success_) &
           & call desc%v_ext_index%clone(desc_out%v_ext_index,info)
      if (info == psb_success_) &
           & call desc%v_ovrlap_index%clone(desc_out%v_ovrlap_index,info)
      if (info == psb_success_) &
           & call desc%v_ovr_mst_idx%clone(desc_out%v_ovr_mst_idx,info)

      if (info /= psb_success_) then
        info = psb_err_from_subroutine_
        call psb_errpush(info,name)
        goto 9999
      endif
      if (debug_level >= psb_debug_ext_) &
           & write(debug_unit,*) me,' ',trim(name),': end'
    else
      call desc_out%free(info)
    end if
    call desc%free(info)

    call psb_erractionrestore(err_act)
    return
    

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_cdtransfer
  !
  ! Subroutine: psb_cd_clone
  !   Copies data and allocation from in to out.
  !
  ! 
  ! Arguments: 
  !    desc  - type(psb_desc_type).         The communication descriptor to be 
  !                                               transferred.
  !    desc_out - type(psb_desc_type).         The output communication descriptor.
  !    info     - integer.                       Return code.
  subroutine psb_cd_clone(desc, desc_out, info)

    use psb_error_mod
    use psb_penv_mod
    use psb_realloc_mod
    implicit none
    !....parameters...

    class(psb_desc_type), intent(inout), target :: desc
    class(psb_desc_type), intent(inout)         :: desc_out
    integer(psb_ipk_), intent(out)              :: info
    !locals
    integer(psb_ipk_) :: np,me,ictxt, err_act
    integer(psb_ipk_) :: debug_level, debug_unit
    character(len=20)   :: name

    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    if (psb_get_errstatus() /= 0) return 
    info = psb_success_
    call psb_erractionsave(err_act)
    name = 'psb_cdcpy'

    if (desc%is_valid()) then 
      ictxt = desc%get_context()

      ! check on blacs grid 
      call psb_info(ictxt, me, np)
      if (debug_level >= psb_debug_ext_) &
           & write(debug_unit,*) me,' ',trim(name),': Entered'
      if (np == -1) then
        info = psb_err_context_error_
        call psb_errpush(info,name)
        goto 9999
      endif

      desc_out%base_desc  => desc%base_desc
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%halo_index,desc_out%halo_index,info)
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%ext_index,desc_out%ext_index,info)
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%ovrlap_index,&
           & desc_out%ovrlap_index,info)
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%bnd_elem,desc_out%bnd_elem,info)
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%ovrlap_elem,desc_out%ovrlap_elem,info)
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%ovr_mst_idx,desc_out%ovr_mst_idx,info)
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%lprm,desc_out%lprm,info)
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%idx_space,desc_out%idx_space,info)
!!$      if ((info == psb_success_).and.(allocated(desc%indxmap))) &
!!$           & call desc%indxmap%clone(desc_out%indxmap,info)
!!$      associate(indxin => desc%indxmap) 
!!$        if ((info == psb_success_).and.(allocated(desc%indxmap))) &
!!$             & call indxin%clone(desc_out%indxmap,info)
!!$      end associate
      if ((info == psb_success_).and.(allocated(desc%indxmap))) &
           & allocate(desc_out%indxmap,source=desc%indxmap,stat=info)
      if (info == psb_success_) &
           & call desc%v_halo_index%clone(desc_out%v_halo_index,info)
      if (info == psb_success_) &
           & call desc%v_ext_index%clone(desc_out%v_ext_index,info)
      if (info == psb_success_) &
           & call desc%v_ovrlap_index%clone(desc_out%v_ovrlap_index,info)
      if (info == psb_success_) &
           & call desc%v_ovr_mst_idx%clone(desc_out%v_ovr_mst_idx,info)


    else
      call desc_out%free(info)
    end if
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info,name)
      goto 9999
    endif
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),': Done'

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ictxt,err_act)

    return

  end subroutine psb_cd_clone

  
  Subroutine psb_cd_get_recv_idx(tmp,desc,data,info,toglob)

    use psb_error_mod
    use psb_penv_mod
    use psb_realloc_mod
    Implicit None
    integer(psb_ipk_), allocatable, intent(out)       :: tmp(:)
    integer(psb_ipk_), intent(in)                     :: data
    Type(psb_desc_type), Intent(in), target :: desc
    integer(psb_ipk_), intent(out)                    :: info
    logical, intent(in)                     :: toglob

    !     .. Local Scalars ..
    integer(psb_ipk_) ::  incnt, outcnt, j, np, me, ictxt, l_tmp,&
         & idx, gidx, proc, n_elem_send, n_elem_recv
    integer(psb_ipk_), pointer   :: idxlist(:) 
    integer(psb_ipk_) :: debug_level, debug_unit, err_act
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
      write(psb_err_unit,*) 'Warning: unusual request getidx on ovr_mst_idx'
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
        call psb_ensure_size((outcnt+3),tmp,info,pad=-ione)
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='psb_ensure_size')
          goto 9999
        end if
        if (toglob) then
          call desc%indxmap%l2g(idx,gidx,info)
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

9999 call psb_error_handler(ictxt,err_act)

    return

  end Subroutine psb_cd_get_recv_idx

  subroutine psb_cd_cnv(desc, mold)
    class(psb_desc_type), intent(inout), target :: desc
    class(psb_i_base_vect_type), intent(in)  :: mold
    
    call desc%v_halo_index%cnv(mold)
    call desc%v_ext_index%cnv(mold)
    call desc%v_ovrlap_index%cnv(mold)
    call desc%v_ovr_mst_idx%cnv(mold)

  end subroutine psb_cd_cnv


  subroutine cd_l2gs1(idx,desc,info,mask,owned)
    use psb_error_mod 
    implicit none 
    class(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(inout) :: idx
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='cd_l2g'

    info  = psb_success_
    call psb_erractionsave(err_act)

    if (allocated(desc%indxmap)) then 
      call desc%indxmap%l2gs1(idx,info,mask=mask,owned=owned)
    else
      info = psb_err_invalid_cd_state_
    end if
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return


9999 call psb_error_handler(err_act)

    return
  
  end subroutine cd_l2gs1

  subroutine cd_l2gs2(idxin,idxout,desc,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(in)    :: idxin
    integer(psb_ipk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='cd_l2g'
    logical, parameter :: debug=.false.

    info  = psb_success_
    call psb_erractionsave(err_act)

    if (allocated(desc%indxmap)) then 
      call desc%indxmap%l2gs2(idxin,idxout,info,mask=mask,owned=owned)
    else
      info = psb_err_invalid_cd_state_
    end if
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine cd_l2gs2


  subroutine cd_l2gv1(idx,desc,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='cd_l2g'
    logical, parameter :: debug=.false.

    info  = psb_success_
    call psb_erractionsave(err_act)

    if (allocated(desc%indxmap)) then 
      call desc%indxmap%l2gv1(idx,info,mask=mask,owned=owned)
    else
      info = psb_err_invalid_cd_state_
    end if
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine cd_l2gv1

  subroutine cd_l2gv2(idxin,idxout,desc,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='cd_l2g'
    logical, parameter :: debug=.false.

    info  = psb_success_
    call psb_erractionsave(err_act)

    if (allocated(desc%indxmap)) then 
      call desc%indxmap%l2gv2(idxin,idxout,info,mask=mask,owned=owned)
    else
      info = psb_err_invalid_cd_state_
    end if
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine cd_l2gv2


  subroutine cd_g2ls1(idx,desc,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(inout) :: idx
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='cd_g2l'
    logical, parameter :: debug=.false.

    info  = psb_success_
    call psb_erractionsave(err_act)

    if (allocated(desc%indxmap)) then 
      call desc%indxmap%g2ls1(idx,info,mask=mask,owned=owned)
    else
      info = psb_err_invalid_cd_state_
    end if
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine cd_g2ls1

  subroutine cd_g2ls2(idxin,idxout,desc,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(in)    :: idxin
    integer(psb_ipk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='cd_g2l'
    logical, parameter :: debug=.false.

    info  = psb_success_
    call psb_erractionsave(err_act)

    if (allocated(desc%indxmap)) then 
      call desc%indxmap%g2ls2(idxin,idxout,info,mask=mask,owned=owned)
    else
      info = psb_err_invalid_cd_state_
    end if
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine cd_g2ls2


  subroutine cd_g2lv1(idx,desc,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='cd_g2l'
    logical, parameter :: debug=.false.

    info  = psb_success_
    call psb_erractionsave(err_act)

    if (allocated(desc%indxmap)) then 
      call desc%indxmap%g2lv1(idx,info,mask=mask,owned=owned)
    else
      info = psb_err_invalid_cd_state_
    end if
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine cd_g2lv1

  subroutine cd_g2lv2(idxin,idxout,desc,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='cd_g2l'
    logical, parameter :: debug=.false.


    info  = psb_success_
    call psb_erractionsave(err_act)

    if (allocated(desc%indxmap)) then 
      call desc%indxmap%g2lv2(idxin,idxout,info,mask=mask,owned=owned)
    else
      info = psb_err_invalid_cd_state_
    end if
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine cd_g2lv2



  subroutine cd_g2ls1_ins(idx,desc,info,mask, lidx)
    use psb_error_mod
    implicit none 
    class(psb_desc_type), intent(inout) :: desc
    integer(psb_ipk_), intent(inout) :: idx
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    integer(psb_ipk_), intent(in), optional :: lidx
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='cd_g2l_ins'
    logical, parameter :: debug=.false.

    info  = psb_success_
    call psb_erractionsave(err_act)

    if (allocated(desc%indxmap)) then 
      call desc%indxmap%g2ls1_ins(idx,info,mask=mask,lidx=lidx)
    else
      info = psb_err_invalid_cd_state_
    end if
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine cd_g2ls1_ins

  subroutine cd_g2ls2_ins(idxin,idxout,desc,info,mask, lidx)
    use psb_error_mod
    implicit none 
    class(psb_desc_type), intent(inout) :: desc
    integer(psb_ipk_), intent(in)    :: idxin
    integer(psb_ipk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    integer(psb_ipk_), intent(in), optional :: lidx

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='cd_g2l_ins'
    logical, parameter :: debug=.false.


    info  = psb_success_
    call psb_erractionsave(err_act)

    if (allocated(desc%indxmap)) then 
      call desc%indxmap%g2ls2_ins(idxin,idxout,info,mask=mask,lidx=lidx)
    else
      info = psb_err_invalid_cd_state_
    end if
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    
    return

  end subroutine cd_g2ls2_ins


  subroutine cd_g2lv1_ins(idx,desc,info,mask, lidx)
    use psb_error_mod
    implicit none 
    class(psb_desc_type), intent(inout) :: desc
    integer(psb_ipk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    integer(psb_ipk_), intent(in), optional :: lidx(:)

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='cd_g2l_ins'
    logical, parameter :: debug=.false.


    info  = psb_success_
    call psb_erractionsave(err_act)

    if (allocated(desc%indxmap)) then 
      call desc%indxmap%g2lv1_ins(idx,info,mask=mask,lidx=lidx)
    else
      info = psb_err_invalid_cd_state_
    end if
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine cd_g2lv1_ins

  subroutine cd_g2lv2_ins(idxin,idxout,desc,info,mask,lidx)
    use psb_error_mod
    implicit none 
    class(psb_desc_type), intent(inout) :: desc
    integer(psb_ipk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    integer(psb_ipk_), intent(in), optional :: lidx(:)

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='cd_g2l_ins'
    logical, parameter :: debug=.false.


    info  = psb_success_
    call psb_erractionsave(err_act)

    if (allocated(desc%indxmap)) then 
      call desc%indxmap%g2lv2_ins(idxin,idxout,info,mask=mask,lidx=lidx)
    else
      info = psb_err_invalid_cd_state_
    end if
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    
    return

  end subroutine cd_g2lv2_ins


  subroutine cd_fnd_owner(idx,iprc,desc,info)
    use psb_error_mod
    implicit none 
    integer(psb_ipk_), intent(in) :: idx(:)
    integer(psb_ipk_), allocatable, intent(out) ::  iprc(:)
    class(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='cd_fnd_owner'
    logical, parameter :: debug=.false.


    info  = psb_success_
    call psb_erractionsave(err_act)

    if (allocated(desc%indxmap)) then 
      call desc%indxmap%fnd_owner(idx,iprc,info)
    else
      info = psb_err_invalid_cd_state_
    end if
    if (info /= psb_success_) then 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine cd_fnd_owner


end module psb_desc_mod
