!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
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
! package: psb_descriptor_type
!    Defines a communication descriptor
!

module psb_descriptor_type
  use psb_const_mod
  use psb_hash_mod 
  use psb_desc_const_mod
  use psb_indx_map_mod

  implicit none

  !
  !  type: psb_desc_type
  !  
  !  Communication Descriptor data structure.
  !
  !|  type psb_desc_type
  !|     class(psb_indx_map), allocatable :: indxmap
  !|     integer, allocatable  :: halo_index(:), ext_index(:)
  !|     integer, allocatable  :: bnd_elem(:)
  !|     integer, allocatable  :: ovrlap_index(:)
  !|     integer, allocatable  :: ovrlap_elem(:,:)
  !|     integer, allocatable  :: ovr_mst_idx(:)
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
    integer, allocatable  :: halo_index(:)
    integer, allocatable  :: ext_index(:)
    integer, allocatable  :: ovrlap_index(:)
    integer, allocatable  :: ovrlap_elem(:,:)
    integer, allocatable  :: ovr_mst_idx(:)
    integer, allocatable  :: bnd_elem(:)
    class(psb_indx_map), allocatable :: indxmap
    integer, allocatable  :: lprm(:)
    type(psb_desc_type), pointer     :: base_desc => null()
    integer, allocatable  :: idx_space(:)
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
    procedure, pass(desc) :: get_local_rows  => psb_cd_get_local_rows
    procedure, pass(desc) :: get_local_cols  => psb_cd_get_local_cols
    procedure, pass(desc) :: get_global_rows => psb_cd_get_global_rows
    procedure, pass(desc) :: get_global_cols => psb_cd_get_global_cols
  end type psb_desc_type

  interface psb_sizeof
    module procedure psb_cd_sizeof
  end interface psb_sizeof

!!$  interface psb_is_ok_desc
!!$    module procedure psb_is_ok_desc
!!$  end interface psb_is_ok_desc
!!$
!!$  interface psb_is_valid_desc
!!$    module procedure psb_is_valid_desc
!!$  end interface psb_is_valid_desc
!!$
!!$  interface psb_is_asb_desc
!!$    module procedure psb_is_asb_desc
!!$  end interface psb_is_asb_desc
!!$
!!$  interface psb_is_upd_desc
!!$    module procedure psb_is_upd_desc
!!$  end interface psb_is_upd_desc
!!$
!!$  interface psb_is_ovl_desc
!!$    module procedure psb_is_ovl_desc
!!$  end interface psb_is_ovl_desc
!!$
!!$  interface psb_is_bld_desc
!!$    module procedure psb_is_bld_desc
!!$  end interface psb_is_bld_desc
!!$

  interface psb_move_alloc
    module procedure psb_cdtransfer
  end interface psb_move_alloc

  interface psb_free
    module procedure psb_cdfree
  end interface psb_free


  integer, private, save :: cd_large_threshold=psb_default_large_threshold 


contains 

  function psb_cd_sizeof(desc)  result(val)
    implicit none
    !....Parameters...

    Type(psb_desc_type), intent(in) :: desc
    integer(psb_long_int_k_) :: val

    val = 0 
    if (allocated(desc%halo_index))   val = val + psb_sizeof_int*size(desc%halo_index)
    if (allocated(desc%ext_index))    val = val + psb_sizeof_int*size(desc%ext_index)
    if (allocated(desc%bnd_elem))     val = val + psb_sizeof_int*size(desc%bnd_elem)
    if (allocated(desc%ovrlap_index)) val = val + psb_sizeof_int*size(desc%ovrlap_index)
    if (allocated(desc%ovrlap_elem))  val = val + psb_sizeof_int*size(desc%ovrlap_elem)
    if (allocated(desc%ovr_mst_idx))  val = val + psb_sizeof_int*size(desc%ovr_mst_idx)
    if (allocated(desc%lprm))         val = val + psb_sizeof_int*size(desc%lprm)
    if (allocated(desc%idx_space))    val = val + psb_sizeof_int*size(desc%idx_space)
    if (allocated(desc%indxmap))      val = val + desc%indxmap%sizeof()

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

  function psb_is_ok_desc(desc) result(val)

    class(psb_desc_type), intent(in) :: desc
    logical                         :: val 
    
    val = .false.
    if (allocated(desc%indxmap)) &
         & val = desc%indxmap%is_valid()

  end function psb_is_ok_desc

  function psb_is_valid_desc(desc) result(val)

    class(psb_desc_type), intent(in) :: desc
    logical                         :: val 
    
    val = .false.
    if (allocated(desc%indxmap)) &
         & val = desc%indxmap%is_valid()

  end function psb_is_valid_desc

  function psb_is_bld_desc(desc) result(val)
    class(psb_desc_type), intent(in) :: desc
    logical                         :: val 

    val = .false.
    if (allocated(desc%indxmap)) &
         & val = desc%indxmap%is_bld()

  end function psb_is_bld_desc

  function psb_is_upd_desc(desc)  result(val)
    class(psb_desc_type), intent(in) :: desc
    logical                         :: val 

    val = .false.
    if (allocated(desc%indxmap)) &
         & val = desc%indxmap%is_upd()

  end function psb_is_upd_desc

  function psb_is_repl_desc(desc) result(val)
    class(psb_desc_type), intent(in) :: desc
    logical                         :: val 

    val = .false.
    if (allocated(desc%indxmap)) &
         & val = desc%indxmap%is_repl()

  end function psb_is_repl_desc

  function psb_is_ovl_desc(desc) result(val)
    class(psb_desc_type), intent(in) :: desc
    logical                         :: val 

    val = .false.
    if (allocated(desc%indxmap)) &
         & val = desc%indxmap%is_ovl()

  end function psb_is_ovl_desc


  function psb_is_asb_desc(desc) result(val)
    class(psb_desc_type), intent(in) :: desc
    logical                         :: val 

    val = .false.
    if (allocated(desc%indxmap)) &
         & val = desc%indxmap%is_asb()

  end function psb_is_asb_desc

  integer function psb_cd_get_local_rows(desc)
    class(psb_desc_type), intent(in) :: desc

    if (psb_is_ok_desc(desc)) then 
      psb_cd_get_local_rows = desc%indxmap%get_lr()
    else
      psb_cd_get_local_rows = -1
    endif
  end function psb_cd_get_local_rows

  integer function psb_cd_get_local_cols(desc)
    class(psb_desc_type), intent(in) :: desc

    if (psb_is_ok_desc(desc)) then 
      psb_cd_get_local_cols = desc%indxmap%get_lc()
    else
      psb_cd_get_local_cols = -1
    endif
  end function psb_cd_get_local_cols

  integer function psb_cd_get_global_rows(desc)
    class(psb_desc_type), intent(in) :: desc
    
    if (psb_is_ok_desc(desc)) then 
      psb_cd_get_global_rows = desc%indxmap%get_gr()
    else
      psb_cd_get_global_rows = -1
    endif

  end function psb_cd_get_global_rows

  integer function psb_cd_get_global_cols(desc)
    class(psb_desc_type), intent(in) :: desc

    if (psb_is_ok_desc(desc)) then 
      psb_cd_get_global_cols = desc%indxmap%get_gc()
    else
      psb_cd_get_global_cols = -1
    endif

  end function psb_cd_get_global_cols

  integer function psb_cd_get_context(desc)
    use psb_error_mod
    class(psb_desc_type), intent(in) :: desc
    if (allocated(desc%indxmap)) then
      psb_cd_get_context = desc%indxmap%get_ctxt()    
    else
      psb_cd_get_context = -1
      call psb_errpush(psb_err_invalid_cd_state_,'psb_cd_get_context')
      call psb_error()
    end if
  end function psb_cd_get_context

  integer function psb_cd_get_dectype(desc)
    use psb_error_mod
    class(psb_desc_type), intent(in) :: desc

    if (allocated(desc%indxmap)) then
      psb_cd_get_dectype = desc%indxmap%get_state()    
    else
      psb_cd_get_dectype = -1
      call psb_errpush(psb_err_invalid_cd_state_,'psb_cd_get_dectype')
      call psb_error()
    end if

  end function psb_cd_get_dectype

  integer function psb_cd_get_mpic(desc)
    use psb_error_mod
    class(psb_desc_type), intent(in) :: desc

    if (allocated(desc%indxmap)) then
      psb_cd_get_mpic = desc%indxmap%get_mpic()    
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

    info = 0
    if (psb_is_asb_desc(desc)) &
         & call desc%indxmap%set_state(psb_desc_ovl_asb_)

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


    ictxt=psb_cd_get_context(desc_a)

    call psb_info(ictxt, me, np)
    !     ....verify blacs grid correctness..
    if (np == -1) then
      info = psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    endif

    
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


    if (allocated(desc_a%lprm)) &
         & deallocate(desc_a%lprm,stat=info)
    if (info /= psb_success_) then 
      info=2057
      call psb_errpush(info,name)
      goto 9999
    end if

    if (allocated(desc_a%indxmap)) then 
      call desc_a%indxmap%free()
      deallocate(desc_a%indxmap, stat=info)
    end if
    if (allocated(desc_a%idx_space)) then 
      deallocate(desc_a%idx_space,stat=info)
      if (info /= psb_success_) then 
        info=2056
        call psb_errpush(info,name)
        goto 9999
      end if
    end if

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
    ictxt = psb_cd_get_context(desc_in)
    call psb_info(ictxt,me,np)
    ! Should not require ictxt to be present: this
    ! function might be called even when desc_in is
    ! empty. 

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
         & call move_alloc(desc_in%indxmap, desc_out%indxmap)
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
        call psb_ensure_size((outcnt+3),tmp,info,pad=-1)
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

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    Return

  end Subroutine psb_cd_get_recv_idx



end module psb_descriptor_type
