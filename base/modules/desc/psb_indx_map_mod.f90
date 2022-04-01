
!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
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
! package: psb_indx_map_mod
!    Defines the PSB_INDX_MAP class.
!
!
!
module psb_indx_map_mod
  use psb_const_mod
  use psb_desc_const_mod
  use psi_penv_mod, only : psb_ctxt_type
  !
  !> \namespace  psb_base_mod  \class  psb_indx_map
  !! \brief Object to handle the mapping between global and local indices.
  !!
  !!
  !!  In particular, besides the communicator, it contains the data relevant
  !!  to the following queries:
  !!  
  !!  -  How many global rows/columns?
  !!  
  !!  -  How many local  rows/columns?
  !!  
  !!  -  Convert between local and global indices
  !!  
  !!  -  Add to local indices.
  !!  
  !!  -  Find (one of) the owner(s) of a given index
  !!  
  !!  -  Query the indx state.
  !!    
  !!  -  Does the dynamic class support extensions of the rows? I.e., can
  !!     it have overlap? for instance, the BLOCK cannot, it would run afoul
  !!     of the glob_to_loc translation. 
  !!  
  !!  Checking for the existence of overlap is very expensive, thus
  !!  it is done at build time (for extended-halo cases it can be inferred from
  !!  the construction process).
  !!
  !!  The object can be in the NULL, BUILD or ASSEMBLED state.
  !!  
  !!  Rules/constraints:
  !!  
  !!  1. Each global index I is owned by at least one process;
  !!
  !!  2. On each process, indices from 1 to N_ROW (desc%indxmap%get_lr())
  !!     are locally owned; the value of N_ROW can be determined upon allocation 
  !!     based on the index distribution (see also the interface to CDALL).
  !!
  !!  3. If a global index is owned by more than one process, we have an OVERLAP
  !!     in which case the sum of all the N_ROW values is greater than the total 
  !!     size of the index space; 
  !!
  !!  4. During the buildup of the descriptor, according to the user specified 
  !!     stencil, we also take notice of indices that are not owned by the current
  !!     process, but whose value is needed to proceed with the computation; these 
  !!     form the HALO of the current process. Halo indices are assigned local indices
  !!     from N_ROW+1 to N_COL (inclusive).
  !!
  !!  5. The upper bound  N_COL moves during the descriptor build process (see CDINS). 
  !!
  !!
  !!  This is the base version of the class; as such, it only contains
  !!  methods for getting/setting the common attributes, whereas
  !!  the index translation methods are only implemented at the derived
  !!  class level.
  !!  Note that the INIT method is defined in the derived methods, and
  !!  is specialized for them; a better solution would have to have
  !!  a constructor for each specific class, with the name of the class,
  !!  but this is not yet working on many compilers, most notably GNU. 
  !!
  !!  Note: the CLONE method was implemented as a workaround for a problem
  !!  with SOURCE= allocation on GNU. Might be removed later on. 
  !!
  type      :: psb_indx_map
    !> State of the map 
    integer(psb_ipk_)   :: state        = psb_desc_null_    
    !> Communication context
    type(psb_ctxt_type) :: ctxt        
    !> MPI communicator
    integer(psb_mpk_)   :: mpic         = -1
    !> Number of global rows
    integer(psb_lpk_)   :: global_rows  = -1
    !> Number of global columns
    integer(psb_lpk_)   :: global_cols  = -1
    !> Number of local rows
    integer(psb_ipk_)   :: local_rows   = -1
    !> Number of local columns
    integer(psb_ipk_)   :: local_cols   = -1
    !> A pointer to the user-defined parts subroutine
    procedure(psb_parts), nopass, pointer  :: parts => null()
    !> The global vector assigning indices to processes, temp copy
    integer(psb_ipk_), allocatable :: tempvg(:)
    !> Reserved for future use. 
    integer(psb_ipk_), allocatable :: oracle(:,:)
    !> Halo owners
    integer(psb_ipk_), allocatable :: halo_owner(:)
    !> Adjacency list for processes
    integer(psb_ipk_), allocatable :: p_adjcncy(:)
    
  contains

    procedure, pass(idxmap)  :: get_state => base_get_state
    procedure, pass(idxmap)  :: set_state => base_set_state
    procedure, pass(idxmap)  :: is_null   => base_is_null
    procedure, nopass        :: is_repl   => base_is_repl
    procedure, pass(idxmap)  :: is_bld    => base_is_bld
    procedure, pass(idxmap)  :: is_upd    => base_is_upd
    procedure, pass(idxmap)  :: is_asb    => base_is_asb
    procedure, pass(idxmap)  :: is_valid  => base_is_valid
    procedure, pass(idxmap)  :: is_ovl    => base_is_ovl
    procedure, pass(idxmap)  :: get_gr    => base_get_gr
    procedure, pass(idxmap)  :: get_gc    => base_get_gc
    procedure, pass(idxmap)  :: get_lr    => base_get_lr
    procedure, pass(idxmap)  :: get_lc    => base_get_lc

    procedure, pass(idxmap)  :: get_p_adjcncy  => base_get_p_adjcncy


    procedure, pass(idxmap)  :: set_gri   => base_set_gri
    procedure, pass(idxmap)  :: set_gci   => base_set_gci
    procedure, pass(idxmap)  :: set_grl   => base_set_grl
    procedure, pass(idxmap)  :: set_gcl   => base_set_gcl
    generic, public          :: set_gr => set_grl
    generic, public          :: set_gc => set_gcl
    
    procedure, pass(idxmap)  :: set_lr    => base_set_lr
    procedure, pass(idxmap)  :: set_lc    => base_set_lc

    procedure, pass(idxmap)  :: set_p_adjcncy   => base_set_p_adjcncy
    procedure, pass(idxmap)  :: xtnd_p_adjcncy  => base_xtnd_p_adjcncy
    
    procedure, pass(idxmap)  :: set_ctxt  => base_set_ctxt
    procedure, pass(idxmap)  :: set_mpic  => base_set_mpic
    procedure, pass(idxmap)  :: get_ctxt  => base_get_ctxt
    procedure, pass(idxmap)  :: get_mpic  => base_get_mpic
    procedure, pass(idxmap)  :: sizeof    => base_sizeof
    procedure, pass(idxmap)  :: set_null  => base_set_null
    procedure, nopass        :: row_extendable => base_row_extendable

    procedure, nopass        :: get_fmt   => base_get_fmt

    procedure, pass(idxmap)  :: asb   => base_asb
    procedure, pass(idxmap)  :: free  => base_free
    procedure, pass(idxmap)  :: clone => base_clone
    procedure, pass(idxmap)  :: cpy   => base_cpy
    procedure, pass(idxmap)  :: reinit => base_reinit

!!$    procedure, pass(idxmap)  :: l2gs1   => base_l2gs1
!!$    procedure, pass(idxmap)  :: l2gs2   => base_l2gs2
!!$    procedure, pass(idxmap)  :: l2gv1   => base_l2gv1
!!$    procedure, pass(idxmap)  :: l2gv2   => base_l2gv2
    procedure, pass(idxmap)  :: ll2gs1  => base_ll2gs1
    procedure, pass(idxmap)  :: ll2gs2  => base_ll2gs2
    procedure, pass(idxmap)  :: ll2gv1  => base_ll2gv1
    procedure, pass(idxmap)  :: ll2gv2  => base_ll2gv2
!!$    generic, public          :: l2g =>   l2gs2, l2gv2
!!$    generic, public          :: l2gip => l2gs1, l2gv1
    generic, public          :: l2g =>   ll2gs2, ll2gv2
    generic, public          :: l2gip => ll2gs1, ll2gv1

!!$    procedure, pass(idxmap)  :: g2ls1   => base_g2ls1
!!$    procedure, pass(idxmap)  :: g2ls2   => base_g2ls2
!!$    procedure, pass(idxmap)  :: g2lv1   => base_g2lv1
!!$    procedure, pass(idxmap)  :: g2lv2   => base_g2lv2
    procedure, pass(idxmap)  :: lg2ls1  => base_lg2ls1
    procedure, pass(idxmap)  :: lg2ls2  => base_lg2ls2
    procedure, pass(idxmap)  :: lg2lv1  => base_lg2lv1
    procedure, pass(idxmap)  :: lg2lv2  => base_lg2lv2
!!$    generic, public          :: g2l =>   g2ls2, g2lv2
!!$    generic, public          :: g2lip => g2ls1, g2lv1
    generic, public          :: g2l =>   lg2ls2, lg2lv2
    generic, public          :: g2lip => lg2ls1, lg2lv1

!!$    procedure, pass(idxmap)  :: g2ls1_ins   => base_g2ls1_ins
!!$    procedure, pass(idxmap)  :: g2ls2_ins   => base_g2ls2_ins
!!$    procedure, pass(idxmap)  :: g2lv1_ins   => base_g2lv1_ins
!!$    procedure, pass(idxmap)  :: g2lv2_ins   => base_g2lv2_ins
    procedure, pass(idxmap)  :: lg2ls1_ins  => base_lg2ls1_ins
    procedure, pass(idxmap)  :: lg2ls2_ins  => base_lg2ls2_ins
    procedure, pass(idxmap)  :: lg2lv1_ins  => base_lg2lv1_ins
    procedure, pass(idxmap)  :: lg2lv2_ins  => base_lg2lv2_ins
!!$    generic, public          :: g2l_ins =>   g2ls2_ins, g2lv2_ins
!!$    generic, public          :: g2lip_ins => g2ls1_ins, g2lv1_ins
    generic, public          :: g2l_ins =>   lg2ls2_ins, lg2lv2_ins
    generic, public          :: g2lip_ins => lg2ls1_ins, lg2lv1_ins

    procedure, pass(idxmap)  :: set_halo_owner  => base_set_halo_owner
    procedure, pass(idxmap)  :: get_halo_owner  => base_get_halo_owner
    procedure, pass(idxmap)  :: qry_halo_owner_s => base_qry_halo_owner_s
    procedure, pass(idxmap)  :: qry_halo_owner_v => base_qry_halo_owner_v
    generic, public          :: qry_halo_owner => qry_halo_owner_s, qry_halo_owner_v
    
    procedure, pass(idxmap)  :: fnd_owner => psi_indx_map_fnd_owner
    procedure, pass(idxmap)  :: init_null => base_init_null 
    procedure, pass(idxmap)  :: init_vl   => base_init_vl
    generic, public          :: init      => init_vl

  end type psb_indx_map

  private :: base_get_state, base_set_state, base_is_repl, base_is_bld,&
       & base_is_upd, base_is_asb, base_is_valid, base_is_ovl,&
       & base_get_gr, base_get_gc, base_get_lr, base_get_lc, base_get_ctxt,&
       & base_get_mpic, base_sizeof, base_set_null, &
       & base_set_grl, base_set_gcl, &
       & base_set_lr, base_set_lc, base_set_ctxt,&
       & base_set_mpic, base_get_fmt, base_asb, base_free,&
       & base_l2gs1, base_l2gs2, base_l2gv1, base_l2gv2,&
       & base_g2ls1, base_g2ls2, base_g2lv1, base_g2lv2,&
       & base_g2ls1_ins, base_g2ls2_ins, base_g2lv1_ins, base_g2lv2_ins, &
       & base_ll2gs1, base_ll2gs2, base_ll2gv1, base_ll2gv2,&
       & base_lg2ls1, base_lg2ls2, base_lg2lv1, base_lg2lv2,&
       & base_lg2ls1_ins, base_lg2ls2_ins, base_lg2lv1_ins,&
       & base_lg2lv2_ins, base_init_vl, base_is_null, base_init_null, &
       & base_row_extendable, base_clone, base_cpy, base_reinit, &
       & base_set_halo_owner, base_get_halo_owner, &
       & base_qry_halo_owner_s, base_qry_halo_owner_v,&
       & base_get_p_adjcncy, base_set_p_adjcncy, base_xtnd_p_adjcncy
  
  !> Function: psi_indx_map_fnd_owner
  !! \memberof psb_indx_map
  !! \brief  Find the process owning indices
  !!
  !!  Given a list of indices IDX, return the processes owning
  !!  them. The base class provides the default implementation,
  !!  which is simply aggregating the requests and converting on
  !!  each proces who then builds its part of the solution.
  !!  This implies that in general this routine is a
  !!  synchronization point; for some derived classes it is
  !!  possible to answer locally, but this should not be relied
  !!  upon. 
  !!
  !!  \param idx(:) The set of indices (local to each process)
  !!  \param iprc(:) The processes owning them
  !!  \param info    return code. 
  !!

  interface 
    subroutine psi_indx_map_fnd_owner(idx,iprc,idxmap,info,adj)
      import :: psb_indx_map, psb_ipk_, psb_lpk_
      implicit none 
      integer(psb_lpk_), intent(in)      :: idx(:)
      integer(psb_ipk_), allocatable, intent(out) ::  iprc(:)
      class(psb_indx_map), intent(in) :: idxmap
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_ipk_), optional, allocatable, intent(out) ::  adj(:)
    end subroutine psi_indx_map_fnd_owner
  end interface

  interface 
    subroutine psi_a2a_fnd_owner(idx,iprc,idxmap,info,samesize)
      import :: psb_indx_map, psb_ipk_, psb_lpk_
      implicit none 
      integer(psb_lpk_), intent(in)      :: idx(:)
      integer(psb_ipk_), allocatable, intent(out) ::  iprc(:)
      class(psb_indx_map), intent(in)    :: idxmap
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: samesize
    end subroutine psi_a2a_fnd_owner
  end interface

  interface 
    subroutine psi_adjcncy_fnd_owner(idx,iprc,adj,idxmap,info)
      import :: psb_indx_map, psb_ipk_, psb_lpk_
      implicit none 
      integer(psb_lpk_), intent(in)   :: idx(:)
      integer(psb_ipk_), allocatable, intent(out)   ::  iprc(:)
      integer(psb_ipk_), intent(inout) :: adj(:)
      class(psb_indx_map), intent(in) :: idxmap
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_adjcncy_fnd_owner
  end interface

  interface 
    subroutine psi_graph_fnd_owner(idx,iprc,ladj,idxmap,info)
      import :: psb_indx_map, psb_ipk_, psb_lpk_
      implicit none 
      integer(psb_lpk_), intent(in)      :: idx(:)
      integer(psb_ipk_), allocatable, intent(out) ::  iprc(:)
      integer(psb_ipk_), allocatable, intent(out) ::  ladj(:)
      class(psb_indx_map), intent(in) :: idxmap
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_graph_fnd_owner
  end interface

  interface psb_cd_set_maxspace
    module procedure  psb_cd_set_maxspace
  end interface psb_cd_set_maxspace

  interface psb_cd_get_maxspace
    module procedure  psb_cd_get_maxspace
  end interface psb_cd_get_maxspace

  interface psb_cd_set_samplesize
    module procedure  psb_cd_set_samplesize
  end interface psb_cd_set_samplesize

  interface psb_cd_get_samplesize
    module procedure  psb_cd_get_samplesize
  end interface psb_cd_get_samplesize

  integer(psb_ipk_), private, save :: cd_maxspace   = -1
  integer(psb_ipk_), private, save :: samplesize    = 32
  
  integer, parameter :: psi_symm_flag_norv_ = 0
  integer, parameter :: psi_symm_flag_inrv_ = 1
  interface psi_symm_dep_list
    subroutine psi_symm_dep_list_inrv(rvsz,adj,ctxt,info)
      import :: psb_indx_map, psb_ipk_, psb_lpk_, psb_mpk_, &
           & psb_ctxt_type
      implicit none 
      integer(psb_mpk_), intent(inout)   :: rvsz(0:)
      integer(psb_ipk_), allocatable, intent(inout) :: adj(:)
      type(psb_ctxt_type), intent(in)      :: ctxt
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_symm_dep_list_inrv
    subroutine psi_symm_dep_list_norv(adj,ctxt,info)
      import :: psb_indx_map, psb_ipk_, psb_lpk_, psb_mpk_, &
           & psb_ctxt_type
      implicit none 
      integer(psb_ipk_), allocatable, intent(inout) :: adj(:)
      type(psb_ctxt_type), intent(in)      :: ctxt
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_symm_dep_list_norv    
  end interface psi_symm_dep_list

  integer(psb_mpk_), parameter :: psi_adj_fnd_irecv_  = 0
  integer(psb_mpk_), parameter :: psi_adj_fnd_a2av_   = 1
  integer(psb_mpk_), parameter :: psi_adj_fnd_pbrcv_  = 2
  integer(psb_mpk_), parameter :: psi_adj_alg_max_  = psi_adj_fnd_pbrcv_  
  integer(psb_mpk_), save      :: psi_adj_alg = psi_adj_fnd_irecv_
  
contains

  subroutine psi_set_adj_alg(ialg)
    integer(psb_mpk_), intent(in) :: ialg
    if ((ialg >=0) .and. (ialg <= psi_adj_alg_max_))&
         & psi_adj_alg = ialg
  end subroutine psi_set_adj_alg

  function psi_get_adj_alg() result(val)
    integer(psb_mpk_) :: val
    val = psi_adj_alg
  end function psi_get_adj_alg

  function psi_get_adj_alg_fmt() result(val)
    character(len=20) :: val
    select case(psi_adj_alg)
    case(psi_adj_fnd_a2av_)
      val = 'MPI_A2AV'
    case(psi_adj_fnd_irecv_)
      val = 'MPI_ISEND/IRECV'
    case(psi_adj_fnd_pbrcv_)
      val = 'PSB_SND/RCV'
    case default
      val = 'Unknown ?'
    end select
  end function psi_get_adj_alg_fmt

  subroutine psb_cd_set_maxspace(ith)
    implicit none 
    integer(psb_ipk_), intent(in) :: ith
    if (ith > 0) then 
      cd_maxspace = ith
    end if
  end subroutine psb_cd_set_maxspace

  function  psb_cd_get_maxspace() result(val)
    implicit none 
    integer(psb_ipk_) :: val
    val  = cd_maxspace
  end function psb_cd_get_maxspace
 

  subroutine psb_cd_set_samplesize(ith)
    implicit none 
    integer(psb_ipk_), intent(in) :: ith
    if (ith > 0) then 
      samplesize = ith
    end if
  end subroutine psb_cd_set_samplesize

  function  psb_cd_get_samplesize() result(val)
    implicit none 
    integer(psb_ipk_) :: val
    val  = samplesize
  end function psb_cd_get_samplesize
 
  
  !> 
  !! \memberof psb_indx_map
  !! \brief  Print a descriptive name
  function base_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'NULL'
  end function base_get_fmt


  function base_get_state(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_) :: val

    val = idxmap%state

  end function base_get_state


  function base_get_gr(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_lpk_) :: val

    val = idxmap%global_rows

  end function base_get_gr


  function base_get_gc(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_lpk_) :: val

    val = idxmap%global_cols

  end function base_get_gc


  function base_get_lr(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_) :: val

    val = idxmap%local_rows

  end function base_get_lr


  function base_get_lc(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_) :: val

    val = idxmap%local_cols

  end function base_get_lc

  function base_get_p_adjcncy(idxmap) result(val)
    use psb_realloc_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_), allocatable  :: val(:)
    integer(psb_ipk_) :: info

    call psb_safe_ab_cpy(idxmap%p_adjcncy,val,info)

  end function base_get_p_adjcncy

  function base_get_ctxt(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    type(psb_ctxt_type) :: val

    val = idxmap%ctxt

  end function base_get_ctxt


  function base_get_mpic(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_mpk_) :: val

    val = idxmap%mpic

  end function base_get_mpic


  subroutine base_set_state(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(in)  :: val

    idxmap%state = val
  end subroutine base_set_state

  subroutine base_set_ctxt(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    type(psb_ctxt_type), intent(in)  :: val

    idxmap%ctxt = val
  end subroutine base_set_ctxt

  subroutine base_set_gri(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(in)  :: val

    idxmap%global_rows = val
  end subroutine base_set_gri

  subroutine base_set_gci(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(in)  :: val

    idxmap%global_cols = val
  end subroutine base_set_gci

  subroutine base_set_grl(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_lpk_), intent(in)  :: val

    idxmap%global_rows = val
  end subroutine base_set_grl

  subroutine base_set_gcl(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_lpk_), intent(in)  :: val

    idxmap%global_cols = val
  end subroutine base_set_gcl

  subroutine base_set_lr(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(in)  :: val

    idxmap%local_rows = val
  end subroutine base_set_lr

  subroutine base_set_lc(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(in)  :: val

    idxmap%local_cols = val
  end subroutine base_set_lc

  subroutine base_set_p_adjcncy(idxmap,val)
    use psb_realloc_mod
    use psb_sort_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(in)  :: val(:)

    call idxmap%xtnd_p_adjcncy(val)

  end subroutine base_set_p_adjcncy
  
  subroutine base_xtnd_p_adjcncy(idxmap,val)
    use psb_realloc_mod
    use psb_sort_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(in)  :: val(:)
    integer(psb_ipk_) :: info, nv, nx
    
    nv = size(val)
    nx = psb_size(idxmap%p_adjcncy)
    call psb_realloc(nv+nx,idxmap%p_adjcncy,info)
    idxmap%p_adjcncy(nx+1:nx+nv) = val(1:nv)
    nx = size(idxmap%p_adjcncy)
    call psb_msort_unique(idxmap%p_adjcncy,nx,dir=psb_sort_up_)
    call psb_realloc(nx,idxmap%p_adjcncy,info)

  end subroutine base_xtnd_p_adjcncy
 
  subroutine base_set_mpic(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_mpk_), intent(in)  :: val

    idxmap%mpic = val
  end subroutine base_set_mpic


  !> 
  !! \memberof psb_indx_map
  !! \brief  Is the class capable of having overlapped rows?
  function base_row_extendable() result(val)
    implicit none 
    logical :: val
    val = .false.
  end function base_row_extendable

  function base_is_repl() result(val)
    implicit none 
    logical :: val
    val = .false.
  end function base_is_repl

  function base_is_null(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    logical :: val
    val =  (idxmap%state == psb_desc_null_)
  end function base_is_null


  function base_is_bld(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    logical :: val
    val = (idxmap%state == psb_desc_bld_).or.&
         & (idxmap%state == psb_desc_ovl_bld_)
  end function base_is_bld

  function base_is_upd(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    logical :: val
    val = (idxmap%state == psb_desc_upd_)
  end function base_is_upd

  function base_is_asb(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    logical :: val
    val = (idxmap%state == psb_desc_asb_).or.&
         & (idxmap%state == psb_desc_ovl_asb_)
  end function base_is_asb

  function base_is_valid(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    logical :: val
    val = idxmap%is_bld().or.idxmap%is_upd().or.idxmap%is_asb()
  end function base_is_valid


  function base_is_ovl(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    logical :: val
    val = (idxmap%state == psb_desc_ovl_bld_).or.&
         & (idxmap%state == psb_desc_ovl_asb_)
  end function base_is_ovl

  function base_sizeof(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_epk_) :: val

    val = 8 * psb_sizeof_ip
  end function base_sizeof


  !> 
  !! \memberof psb_indx_map
  !! \brief  Local to global, scalar, in place
  subroutine base_l2gs1(idx,idxmap,info,mask,owned)
    use psb_error_mod 
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(inout) :: idx
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_l2g'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_l2gs1

  subroutine base_l2gs2(idxin,idxout,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin
    integer(psb_ipk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_l2g'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)

  end subroutine base_l2gs2


  subroutine base_l2gv1(idx,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_l2g'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return
  end subroutine base_l2gv1

  subroutine base_l2gv2(idxin,idxout,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_l2g'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_l2gv2

  !> 
  !! \memberof psb_indx_map
  !! \brief  Local to global, scalar, in place
  subroutine base_ll2gs1(idx,idxmap,info,mask,owned)
    use psb_error_mod 
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_ll2g'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return
    
  end subroutine base_ll2gs1

  subroutine base_ll2gs2(idxin,idxout,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin
    integer(psb_lpk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_ll2g'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)

  end subroutine base_ll2gs2


  subroutine base_ll2gv1(idx,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_ll2g'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return
  end subroutine base_ll2gv1

  subroutine base_ll2gv2(idxin,idxout,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin(:)
    integer(psb_lpk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_ll2g'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_ll2gv2


  subroutine base_g2ls1(idx,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(inout) :: idx
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_g2l'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_g2ls1

  subroutine base_g2ls2(idxin,idxout,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin
    integer(psb_ipk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_g2l'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_g2ls2


  subroutine base_g2lv1(idx,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_g2l'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_g2lv1

  subroutine base_g2lv2(idxin,idxout,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_g2l'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return


  end subroutine base_g2lv2

  subroutine base_lg2ls1(idx,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_lg2l'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_lg2ls1

  subroutine base_lg2ls2(idxin,idxout,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(in)    :: idxin
    integer(psb_ipk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_lg2l'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_lg2ls2


  subroutine base_lg2lv1(idx,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_lg2l'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_lg2lv1

  subroutine base_lg2lv2(idxin,idxout,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_lg2l'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_lg2lv2


  subroutine base_g2ls1_ins(idx,idxmap,info,mask, lidx)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(inout) :: idx
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    integer(psb_ipk_), intent(in), optional :: lidx
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_g2l_ins'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_g2ls1_ins

  subroutine base_g2ls2_ins(idxin,idxout,idxmap,info,mask, lidx)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin
    integer(psb_ipk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    integer(psb_ipk_), intent(in), optional :: lidx

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_g2l_ins'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_g2ls2_ins


  subroutine base_g2lv1_ins(idx,idxmap,info,mask, lidx)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    integer(psb_ipk_), intent(in), optional :: lidx(:)

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_g2l_ins'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_g2lv1_ins

  subroutine base_g2lv2_ins(idxin,idxout,idxmap,info,mask,lidx)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    integer(psb_ipk_), intent(in), optional :: lidx(:)

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_g2l_ins'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_g2lv2_ins


  subroutine base_lg2ls1_ins(idx,idxmap,info,mask, lidx)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    integer(psb_ipk_), intent(in), optional :: lidx
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_lg2l_ins'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_lg2ls1_ins

  subroutine base_lg2ls2_ins(idxin,idxout,idxmap,info,mask, lidx)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_lpk_), intent(in)    :: idxin
    integer(psb_ipk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    integer(psb_ipk_), intent(in), optional :: lidx

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_lg2l_ins'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_lg2ls2_ins


  subroutine base_lg2lv1_ins(idx,idxmap,info,mask, lidx)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    integer(psb_ipk_), intent(in), optional :: lidx(:)

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_lg2l_ins'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_lg2lv1_ins

  subroutine base_lg2lv2_ins(idxin,idxout,idxmap,info,mask,lidx)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_lpk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    integer(psb_ipk_), intent(in), optional :: lidx(:)

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_lg2l_ins'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_lg2lv2_ins

  subroutine base_asb(idxmap,info)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_asb'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return

  end subroutine base_asb

  subroutine base_free(idxmap)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap

    ! almost nothing to be done here
    idxmap%state          = -1 
    if (allocated(idxmap%ctxt%ctxt)) deallocate(idxmap%ctxt%ctxt)
    idxmap%mpic           = -1
    idxmap%global_rows    = -1
    idxmap%global_cols    = -1
    idxmap%local_rows     = -1
    idxmap%local_cols     = -1

    return

  end subroutine base_free

  subroutine base_set_null(idxmap)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap

    idxmap%state          = psb_desc_null_
    if (allocated(idxmap%ctxt%ctxt)) deallocate(idxmap%ctxt%ctxt)
    idxmap%mpic           = -1
    idxmap%global_rows    = -1
    idxmap%global_cols    = -1
    idxmap%local_rows     = -1
    idxmap%local_cols     = -1

  end subroutine base_set_null

  subroutine base_init_null(idxmap,ctxt,info)
    class(psb_indx_map), intent(inout) :: idxmap
    type(psb_ctxt_type), intent(in)  :: ctxt
    integer(psb_ipk_), intent(out) :: info

    call idxmap%set_null()
    idxmap%ctxt = ctxt
    info = 0
    return
  end subroutine base_init_null
    
  subroutine base_init_vl(idxmap,ctxt,vl,info)
    use psb_penv_mod
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    type(psb_ctxt_type), intent(in)  :: ctxt
    integer(psb_lpk_), intent(in)  :: vl(:)
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_init_vl'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

  call psb_error_handler(err_act)
    return
  end subroutine base_init_vl

  subroutine base_clone(idxmap,outmap,info)
    use psb_penv_mod
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout)    :: idxmap
    class(psb_indx_map), allocatable, intent(out) :: outmap
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_clone'
    logical, parameter :: debug=.false.

    info = psb_success_

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return
  end subroutine base_clone

  !
  ! This is a kludge because we defined the outmap
  ! in base_clone to be allocatable intent(out).
  ! Should be revisited. 
  !
  subroutine base_cpy(idxmap,outmap,info)
    use psb_penv_mod
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    class(psb_indx_map), intent(out)   :: outmap
    integer(psb_ipk_), intent(out)     :: info
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_clone'
    logical, parameter :: debug=.false.

    info = psb_success_

    call psb_get_erraction(err_act)

    outmap%state       = idxmap%state      
    outmap%ctxt        = idxmap%ctxt      
    outmap%mpic        = idxmap%mpic       
    outmap%global_rows = idxmap%global_rows
    outmap%global_cols = idxmap%global_cols
    outmap%local_rows  = idxmap%local_rows 
    outmap%local_cols  = idxmap%local_cols 
    outmap%parts       => idxmap%parts 
    if (info == psb_success_)&
         &  call psb_safe_ab_cpy(idxmap%tempvg,outmap%tempvg,info)
    if (info == psb_success_)&
         &  call psb_safe_ab_cpy(idxmap%oracle,outmap%oracle,info)
    if (info == psb_success_)&
         &  call psb_safe_ab_cpy(idxmap%halo_owner,outmap%halo_owner,info)
    if (info == psb_success_)&
         &  call psb_safe_ab_cpy(idxmap%p_adjcncy,outmap%p_adjcncy,info)
    
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return
9999 call psb_error_handler(err_act)

    return
  end subroutine base_cpy


  subroutine base_reinit(idxmap,info)
    use psb_penv_mod
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout)    :: idxmap
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='base_reinit'
    logical, parameter :: debug=.false.

    info = psb_success_

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    call psb_error_handler(err_act)
    return
  end subroutine base_reinit


  subroutine base_set_halo_owner(idxmap,v,info)
    use psb_penv_mod
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(in)      :: v(:)
    integer(psb_ipk_), intent(out)     :: info
    !
    integer(psb_ipk_) :: me, np
    integer(psb_ipk_) :: i, j, nr, nc, nh

    call psb_info(idxmap%ctxt,me,np)
    ! The idea here is to store only the halo part
    nr = idxmap%local_rows
    nc = idxmap%local_cols
    nh = nc-nr
    if (size(v) < nh) then
      write(0,*) 'Error: set_halo_owner small size ',size(v),nh
    end if
    idxmap%halo_owner = v(1:min(size(v),nh))
    
  end subroutine base_set_halo_owner

  subroutine base_get_halo_owner(idxmap,v,info)
    use psb_realloc_mod
    use psb_penv_mod
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer(psb_ipk_), allocatable, intent(out) :: v(:)
    integer(psb_ipk_), intent(out)     :: info

    integer(psb_ipk_)  :: nh
    nh = size(idxmap%halo_owner)
    !v = idxmap%halo_owner(1:nh)
    call psb_safe_ab_cpy(idxmap%halo_owner,v,info)
  end subroutine base_get_halo_owner

  subroutine base_qry_halo_owner_s(idxmap,xin,xout,info)
    use psb_penv_mod
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)  ::  xin
    integer(psb_ipk_), intent(out) ::  xout
    integer(psb_ipk_), intent(out)     :: info

    integer(psb_ipk_)  :: i, j, nr, nc, nh
    nr = idxmap%local_rows
    nc = idxmap%local_cols
    nc = min(idxmap%local_cols, (nr+psb_size(idxmap%halo_owner)))    
    xout = -1
    if (.not.allocated(idxmap%halo_owner)) then
      !write(0,*) 'Halo_owner not allocated!', nr, nc, xin
      return
    end if
    if ((nr<xin).and.(xin <= nc)) then
!!$      if (size(idxmap%halo_owner)<(xin-nr)) then
!!$        !write(0,*) 'Halo_owner bad size',xin,nr,xin-nr,size(idxmap%halo_owner)
!!$        return
!!$      end if
      xout = idxmap%halo_owner(xin-nr)
    end if
    
  end subroutine base_qry_halo_owner_s

  subroutine base_qry_halo_owner_v(idxmap,xin,xout,info)
    use psb_penv_mod
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)  ::  xin(:)
    integer(psb_ipk_), intent(out) ::  xout(:)
    integer(psb_ipk_), intent(out)     :: info

    integer(psb_ipk_)  :: i, j, nr, nc, nh, sz
    nr = idxmap%local_rows
    nc = min(idxmap%local_cols, (nr+psb_size(idxmap%halo_owner)))
    sz = min(size(xin),size(xout))
    if (.not.allocated(idxmap%halo_owner)) then
      xout = -1
      return
    end if
    
    do i = 1, sz
      xout(i) = -1          
      if ((nr<xin(i)).and.(xin(i) <= nc)) xout(i) = idxmap%halo_owner(xin(i)-nr)
    end do
    do i=sz+1,size(xout)
      xout(i) = -1
    end do
  end subroutine base_qry_halo_owner_v

end module psb_indx_map_mod
