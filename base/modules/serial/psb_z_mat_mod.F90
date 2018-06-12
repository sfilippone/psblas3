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
! package: psb_z_mat_mod
!
! This module contains the definition of the psb_z_sparse type which
! is a generic container for a sparse matrix and it is mostly meant to
! provide a mean of switching, at run-time, among different formats,
! potentially unknown at the library compile-time by adding a layer of
! indirection. This type encapsulates the psb_z_base_sparse_mat class
! inside another class which is the one visible to the user. 
! Most methods of the psb_z_mat_mod simply call the methods of the
! encapsulated class.
! The exceptions are mainly cscnv and cp_from/cp_to; these provide
! the functionalities to have the encapsulated class change its
! type dynamically, and to extract/input an inner object.
!
! A sparse matrix has a state corresponding to its progression
! through the application life.
! In particular, computational methods can only be invoked when
! the matrix is in the ASSEMBLED state, whereas the other states are
! dedicated to operations on the internal matrix data. 
! A sparse matrix can move between states according to the 
! following state transition table. Associated with these states are
! the possible dynamic types of the inner matrix object.
! Only COO matrices can ever be in the BUILD state, whereas
! the ASSEMBLED and UPDATE state can be entered by any class. 
! 
! In           Out        Method    
!| ----------------------------------
!| Null         Build      csall
!| Build        Build      csput
!| Build        Assembled  cscnv
!| Assembled    Assembled  cscnv
!| Assembled    Update     reinit
!| Update       Update     csput
!| Update       Assembled  cscnv
!| *            unchanged  reall 
!| Assembled    Null       free
!
!
!
! We are also introducing the type psb_lzspmat_type.
! The basic difference with psb_zspmat_type is in the type
! of the indices, which are PSB_LPK_ so that the entries
! are guaranteed to be able to contain global indices.
! This type only supports data handling and preprocessing, it is
! not supposed to be used for computations. 
!
module psb_z_mat_mod

  use psb_z_base_mat_mod
  use psb_z_csr_mat_mod,  only : psb_z_csr_sparse_mat, psb_lz_csr_sparse_mat
  use psb_z_csc_mat_mod,  only : psb_z_csc_sparse_mat

  type :: psb_zspmat_type

    class(psb_z_base_sparse_mat), allocatable  :: a 

  contains
    ! Getters
    procedure, pass(a) :: get_nrows   => psb_z_get_nrows
    procedure, pass(a) :: get_ncols   => psb_z_get_ncols
    procedure, pass(a) :: get_nzeros  => psb_z_get_nzeros
    procedure, pass(a) :: get_nz_row  => psb_z_get_nz_row
    procedure, pass(a) :: get_size    => psb_z_get_size
    procedure, pass(a) :: get_dupl    => psb_z_get_dupl
    procedure, pass(a) :: is_null     => psb_z_is_null
    procedure, pass(a) :: is_bld      => psb_z_is_bld
    procedure, pass(a) :: is_upd      => psb_z_is_upd
    procedure, pass(a) :: is_asb      => psb_z_is_asb
    procedure, pass(a) :: is_sorted   => psb_z_is_sorted
    procedure, pass(a) :: is_by_rows  => psb_z_is_by_rows
    procedure, pass(a) :: is_by_cols  => psb_z_is_by_cols
    procedure, pass(a) :: is_upper    => psb_z_is_upper
    procedure, pass(a) :: is_lower    => psb_z_is_lower
    procedure, pass(a) :: is_triangle => psb_z_is_triangle
    procedure, pass(a) :: is_unit     => psb_z_is_unit
    procedure, pass(a) :: is_repeatable_updates => psb_z_is_repeatable_updates
    procedure, pass(a) :: get_fmt     => psb_z_get_fmt
    procedure, pass(a) :: sizeof      => psb_z_sizeof

    ! Setters
    procedure, pass(a) :: set_nrows    => psb_z_set_nrows
    procedure, pass(a) :: set_ncols    => psb_z_set_ncols
    procedure, pass(a) :: set_dupl     => psb_z_set_dupl
    procedure, pass(a) :: set_null     => psb_z_set_null
    procedure, pass(a) :: set_bld      => psb_z_set_bld
    procedure, pass(a) :: set_upd      => psb_z_set_upd
    procedure, pass(a) :: set_asb      => psb_z_set_asb
    procedure, pass(a) :: set_sorted   => psb_z_set_sorted
    procedure, pass(a) :: set_upper    => psb_z_set_upper
    procedure, pass(a) :: set_lower    => psb_z_set_lower
    procedure, pass(a) :: set_triangle => psb_z_set_triangle
    procedure, pass(a) :: set_unit     => psb_z_set_unit
    procedure, pass(a) :: set_repeatable_updates => psb_z_set_repeatable_updates

    ! Memory/data management 
    procedure, pass(a) :: csall       => psb_z_csall
    procedure, pass(a) :: free        => psb_z_free
    procedure, pass(a) :: trim        => psb_z_trim
    procedure, pass(a) :: csput_a     => psb_z_csput_a
    procedure, pass(a) :: csput_v     => psb_z_csput_v 
    generic, public    :: csput       => csput_a,  csput_v
    procedure, pass(a) :: csgetptn    => psb_z_csgetptn
    procedure, pass(a) :: csgetrow    => psb_z_csgetrow
    procedure, pass(a) :: csgetblk    => psb_z_csgetblk
    generic, public    :: csget       => csgetptn, csgetrow, csgetblk
#if defined(IPK4) && defined(LPK8)
    procedure, pass(a) :: lcsgetptn    => psb_z_lcsgetptn
    procedure, pass(a) :: lcsgetrow    => psb_z_lcsgetrow
    generic, public    :: csget        => lcsgetptn, lcsgetrow
#endif    
    procedure, pass(a) :: tril        => psb_z_tril
    procedure, pass(a) :: triu        => psb_z_triu
    procedure, pass(a) :: m_csclip    => psb_z_csclip
    procedure, pass(a) :: b_csclip    => psb_z_b_csclip
    generic, public    :: csclip      => b_csclip, m_csclip
    procedure, pass(a) :: clean_zeros => psb_z_clean_zeros
    procedure, pass(a) :: reall       => psb_z_reallocate_nz
    procedure, pass(a) :: get_neigh   => psb_z_get_neigh
    procedure, pass(a) :: reinit      => psb_z_reinit
    procedure, pass(a) :: print_i     => psb_z_sparse_print
    procedure, pass(a) :: print_n     => psb_z_n_sparse_print
    generic, public    :: print       => print_i, print_n
    procedure, pass(a) :: mold        => psb_z_mold
    procedure, pass(a) :: asb         => psb_z_asb
    procedure, pass(a) :: transp_1mat => psb_z_transp_1mat
    procedure, pass(a) :: transp_2mat => psb_z_transp_2mat
    generic, public    :: transp      => transp_1mat, transp_2mat
    procedure, pass(a) :: transc_1mat => psb_z_transc_1mat
    procedure, pass(a) :: transc_2mat => psb_z_transc_2mat
    generic, public    :: transc      => transc_1mat, transc_2mat

    !
    ! Sync: centerpiece of handling of external storage.
    ! Any derived class having extra storage upon sync
    ! will guarantee that both fortran/host side and
    ! external side contain the same data. The base
    ! version is only a placeholder. 
    !
    procedure, pass(a) :: sync        => z_mat_sync
    procedure, pass(a) :: is_host     => z_mat_is_host
    procedure, pass(a) :: is_dev      => z_mat_is_dev
    procedure, pass(a) :: is_sync     => z_mat_is_sync
    procedure, pass(a) :: set_host    => z_mat_set_host
    procedure, pass(a) :: set_dev     => z_mat_set_dev
    procedure, pass(a) :: set_sync    => z_mat_set_sync


    ! These are specific to this level of encapsulation.
    procedure, pass(a) :: mv_from_b   => psb_z_mv_from
    generic, public    :: mv_from     => mv_from_b
    procedure, pass(a) :: mv_to_b     => psb_z_mv_to
    generic, public    :: mv_to       => mv_to_b
    procedure, pass(a) :: cp_from_b   => psb_z_cp_from
    generic, public    :: cp_from     => cp_from_b
    procedure, pass(a) :: cp_to_b     => psb_z_cp_to
    generic, public    :: cp_to       => cp_to_b
    procedure, pass(a) :: clip_d_ip   => psb_z_clip_d_ip
    procedure, pass(a) :: clip_d      => psb_z_clip_d
    generic, public    :: clip_diag   => clip_d_ip, clip_d
    procedure, pass(a) :: cscnv_np    => psb_z_cscnv
    procedure, pass(a) :: cscnv_ip    => psb_z_cscnv_ip
    procedure, pass(a) :: cscnv_base  => psb_z_cscnv_base
    generic, public    :: cscnv       => cscnv_np, cscnv_ip, cscnv_base
    procedure, pass(a) :: clone       => psb_zspmat_clone
    !
    ! To/from lz
    !
    procedure, pass(a) :: mv_from_lb  => psb_z_mv_from_lb
    procedure, pass(a) :: mv_to_lb    => psb_z_mv_to_lb
    procedure, pass(a) :: cp_from_lb  => psb_z_cp_from_lb
    procedure, pass(a) :: cp_to_lb    => psb_z_cp_to_lb
    procedure, pass(a) :: mv_from_l   => psb_z_mv_from_l 
    procedure, pass(a) :: mv_to_l     => psb_z_mv_to_l 
    procedure, pass(a) :: cp_from_l   => psb_z_cp_from_l 
    procedure, pass(a) :: cp_to_l     => psb_z_cp_to_l 
    generic, public    :: mv_from     => mv_from_lb, mv_from_l
    generic, public    :: mv_to       => mv_to_lb, mv_to_l
    generic, public    :: cp_from     => cp_from_lb, cp_from_l
    generic, public    :: cp_to       => cp_to_lb, cp_to_l
    
    ! Computational routines 
    procedure, pass(a) :: get_diag => psb_z_get_diag
    procedure, pass(a) :: maxval   => psb_z_maxval
    procedure, pass(a) :: spnmi    => psb_z_csnmi
    procedure, pass(a) :: spnm1    => psb_z_csnm1
    procedure, pass(a) :: rowsum   => psb_z_rowsum
    procedure, pass(a) :: arwsum   => psb_z_arwsum
    procedure, pass(a) :: colsum   => psb_z_colsum
    procedure, pass(a) :: aclsum   => psb_z_aclsum
    procedure, pass(a) :: csmv_v   => psb_z_csmv_vect
    procedure, pass(a) :: csmv     => psb_z_csmv
    procedure, pass(a) :: csmm     => psb_z_csmm
    generic, public    :: spmm     => csmm, csmv, csmv_v
    procedure, pass(a) :: scals    => psb_z_scals
    procedure, pass(a) :: scalv    => psb_z_scal
    generic, public    :: scal     => scals, scalv
    procedure, pass(a) :: cssv_v   => psb_z_cssv_vect
    procedure, pass(a) :: cssv     => psb_z_cssv
    procedure, pass(a) :: cssm     => psb_z_cssm
    generic, public    :: spsm     => cssm, cssv, cssv_v

  end type psb_zspmat_type

  private :: psb_z_get_nrows, psb_z_get_ncols, &
       & psb_z_get_nzeros, psb_z_get_size, &
       & psb_z_get_dupl, psb_z_is_null, psb_z_is_bld, &
       & psb_z_is_upd, psb_z_is_asb, psb_z_is_sorted, &
       & psb_z_is_by_rows, psb_z_is_by_cols, psb_z_is_upper, &
       & psb_z_is_lower, psb_z_is_triangle, psb_z_get_nz_row, &
       & z_mat_sync, z_mat_is_host, z_mat_is_dev, &
       & z_mat_is_sync, z_mat_set_host, z_mat_set_dev,&
       & z_mat_set_sync



  class(psb_z_base_sparse_mat), allocatable, target, &
       & save, private :: psb_z_base_mat_default

  interface psb_set_mat_default
    module procedure psb_z_set_mat_default
  end interface

  interface psb_get_mat_default
    module procedure psb_z_get_mat_default
  end interface

  interface psb_sizeof
    module procedure psb_z_sizeof
  end interface


  type :: psb_lzspmat_type

    class(psb_lz_base_sparse_mat), allocatable  :: a 

  contains
    ! Getters
    procedure, pass(a) :: get_nrows   => psb_lz_get_nrows
    procedure, pass(a) :: get_ncols   => psb_lz_get_ncols
    procedure, pass(a) :: get_nzeros  => psb_lz_get_nzeros
    procedure, pass(a) :: get_nz_row  => psb_lz_get_nz_row
    procedure, pass(a) :: get_size    => psb_lz_get_size
    procedure, pass(a) :: get_dupl    => psb_lz_get_dupl
    procedure, pass(a) :: is_null     => psb_lz_is_null
    procedure, pass(a) :: is_bld      => psb_lz_is_bld
    procedure, pass(a) :: is_upd      => psb_lz_is_upd
    procedure, pass(a) :: is_asb      => psb_lz_is_asb
    procedure, pass(a) :: is_sorted   => psb_lz_is_sorted
    procedure, pass(a) :: is_by_rows  => psb_lz_is_by_rows
    procedure, pass(a) :: is_by_cols  => psb_lz_is_by_cols
    procedure, pass(a) :: is_upper    => psb_lz_is_upper
    procedure, pass(a) :: is_lower    => psb_lz_is_lower
    procedure, pass(a) :: is_triangle => psb_lz_is_triangle
    procedure, pass(a) :: is_unit     => psb_lz_is_unit
    procedure, pass(a) :: is_repeatable_updates => psb_lz_is_repeatable_updates
    procedure, pass(a) :: get_fmt     => psb_lz_get_fmt
    procedure, pass(a) :: sizeof      => psb_lz_sizeof

    ! Setters
    procedure, pass(a) :: set_nrows    => psb_lz_set_nrows
    procedure, pass(a) :: set_ncols    => psb_lz_set_ncols
    procedure, pass(a) :: set_dupl     => psb_lz_set_dupl
    procedure, pass(a) :: set_null     => psb_lz_set_null
    procedure, pass(a) :: set_bld      => psb_lz_set_bld
    procedure, pass(a) :: set_upd      => psb_lz_set_upd
    procedure, pass(a) :: set_asb      => psb_lz_set_asb
    procedure, pass(a) :: set_sorted   => psb_lz_set_sorted
    procedure, pass(a) :: set_upper    => psb_lz_set_upper
    procedure, pass(a) :: set_lower    => psb_lz_set_lower
    procedure, pass(a) :: set_triangle => psb_lz_set_triangle
    procedure, pass(a) :: set_unit     => psb_lz_set_unit
    procedure, pass(a) :: set_repeatable_updates => psb_lz_set_repeatable_updates

    ! Memory/data management 
    procedure, pass(a) :: csall       => psb_lz_csall
    procedure, pass(a) :: free        => psb_lz_free
    procedure, pass(a) :: trim        => psb_lz_trim
    procedure, pass(a) :: csput_a     => psb_lz_csput_a
    procedure, pass(a) :: csput_v     => psb_lz_csput_v 
    generic, public    :: csput       => csput_a,  csput_v
    procedure, pass(a) :: csgetptn    => psb_lz_csgetptn
    procedure, pass(a) :: csgetrow    => psb_lz_csgetrow
    procedure, pass(a) :: csgetblk    => psb_lz_csgetblk
    generic, public    :: csget       => csgetptn, csgetrow, csgetblk
#if defined(IPK4) && defined(LPK8)
    procedure, pass(a) :: icsgetptn    => psb_lz_icsgetptn
    procedure, pass(a) :: icsgetrow    => psb_lz_icsgetrow
    generic, public    :: csget        => icsgetptn, icsgetrow
#endif    
    procedure, pass(a) :: tril        => psb_lz_tril
    procedure, pass(a) :: triu        => psb_lz_triu
    procedure, pass(a) :: m_csclip    => psb_lz_csclip
    procedure, pass(a) :: b_csclip    => psb_lz_b_csclip
    generic, public    :: csclip      => b_csclip, m_csclip
    procedure, pass(a) :: clean_zeros => psb_lz_clean_zeros
    procedure, pass(a) :: reall       => psb_lz_reallocate_nz
    procedure, pass(a) :: get_neigh   => psb_lz_get_neigh
    procedure, pass(a) :: reinit      => psb_lz_reinit
    procedure, pass(a) :: print_i     => psb_lz_sparse_print
    procedure, pass(a) :: print_n     => psb_lz_n_sparse_print
    generic, public    :: print       => print_i, print_n
    procedure, pass(a) :: mold        => psb_lz_mold
    procedure, pass(a) :: asb         => psb_lz_asb
    procedure, pass(a) :: transp_1mat => psb_lz_transp_1mat
    procedure, pass(a) :: transp_2mat => psb_lz_transp_2mat
    generic, public    :: transp      => transp_1mat, transp_2mat
    procedure, pass(a) :: transc_1mat => psb_lz_transc_1mat
    procedure, pass(a) :: transc_2mat => psb_lz_transc_2mat
    generic, public    :: transc      => transc_1mat, transc_2mat

    !
    ! Sync: centerpiece of handling of external storage.
    ! Any derived class having extra storage upon sync
    ! will guarantee that both fortran/host side and
    ! external side contain the same data. The base
    ! version is only a placeholder. 
    !
    procedure, pass(a) :: sync        => lz_mat_sync
    procedure, pass(a) :: is_host     => lz_mat_is_host
    procedure, pass(a) :: is_dev      => lz_mat_is_dev
    procedure, pass(a) :: is_sync     => lz_mat_is_sync
    procedure, pass(a) :: set_host    => lz_mat_set_host
    procedure, pass(a) :: set_dev     => lz_mat_set_dev
    procedure, pass(a) :: set_sync    => lz_mat_set_sync


    ! These are specific to this level of encapsulation.
    procedure, pass(a) :: mv_from_b   => psb_lz_mv_from
    generic, public    :: mv_from     => mv_from_b
    procedure, pass(a) :: mv_to_b     => psb_lz_mv_to
    generic, public    :: mv_to       => mv_to_b
    procedure, pass(a) :: cp_from_b   => psb_lz_cp_from
    generic, public    :: cp_from     => cp_from_b
    procedure, pass(a) :: cp_to_b     => psb_lz_cp_to
    generic, public    :: cp_to       => cp_to_b
    procedure, pass(a) :: cscnv_np    => psb_lz_cscnv
    procedure, pass(a) :: cscnv_ip    => psb_lz_cscnv_ip
    procedure, pass(a) :: cscnv_base  => psb_lz_cscnv_base
    generic, public    :: cscnv       => cscnv_np, cscnv_ip, cscnv_base
    procedure, pass(a) :: clone       => psb_lzspmat_clone
    !
    ! To/from z
    !
    procedure, pass(a) :: mv_from_ib  => psb_lz_mv_from_ib
    procedure, pass(a) :: mv_to_ib    => psb_lz_mv_to_ib
    procedure, pass(a) :: cp_from_ib  => psb_lz_cp_from_ib
    procedure, pass(a) :: cp_to_ib    => psb_lz_cp_to_ib
    procedure, pass(a) :: mv_from_i   => psb_lz_mv_from_i 
    procedure, pass(a) :: mv_to_i     => psb_lz_mv_to_i 
    procedure, pass(a) :: cp_from_i   => psb_lz_cp_from_i 
    procedure, pass(a) :: cp_to_i     => psb_lz_cp_to_i 
    generic, public    :: mv_from     => mv_from_ib, mv_from_i
    generic, public    :: mv_to       => mv_to_ib, mv_to_i
    generic, public    :: cp_from     => cp_from_ib, cp_from_i
    generic, public    :: cp_to       => cp_to_ib, cp_to_i

    ! Computational routines 
    procedure, pass(a) :: get_diag => psb_lz_get_diag
    procedure, pass(a) :: scals    => psb_lz_scals
    procedure, pass(a) :: scalv    => psb_lz_scal
    generic, public    :: scal     => scals, scalv

  end type psb_lzspmat_type

  private :: psb_lz_get_nrows, psb_lz_get_ncols, &
       & psb_lz_get_nzeros, psb_lz_get_size, &
       & psb_lz_get_dupl, psb_lz_is_null, psb_lz_is_bld, &
       & psb_lz_is_upd, psb_lz_is_asb, psb_lz_is_sorted, &
       & psb_lz_is_by_rows, psb_lz_is_by_cols, psb_lz_is_upper, &
       & psb_lz_is_lower, psb_lz_is_triangle, psb_lz_get_nz_row, &
       & lz_mat_sync, lz_mat_is_host, lz_mat_is_dev, &
       & lz_mat_is_sync, lz_mat_set_host, lz_mat_set_dev,&
       & lz_mat_set_sync



  class(psb_lz_base_sparse_mat), allocatable, target, &
       & save, private :: psb_lz_base_mat_default

  interface psb_set_mat_default
    module procedure psb_lz_set_mat_default
  end interface

  interface psb_get_mat_default
    module procedure psb_lz_get_mat_default
  end interface


  ! == ===================================
  !
  !
  !
  ! Setters 
  !
  !
  !
  !
  !
  !
  ! == ===================================


  interface 
    subroutine  psb_z_set_nrows(m,a) 
      import :: psb_ipk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in) :: m
    end subroutine psb_z_set_nrows
  end interface
  
  interface 
    subroutine psb_z_set_ncols(n,a) 
      import :: psb_ipk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in) :: n
    end subroutine psb_z_set_ncols
  end interface
  
  interface 
    subroutine  psb_z_set_dupl(n,a) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in) :: n
    end subroutine psb_z_set_dupl
  end interface
  
  interface 
    subroutine psb_z_set_null(a) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_set_null
  end interface
  
  interface 
    subroutine psb_z_set_bld(a) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_set_bld
  end interface
  
  interface 
    subroutine psb_z_set_upd(a) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_set_upd
  end interface
  
  interface 
    subroutine psb_z_set_asb(a) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_set_asb
  end interface
  
  interface 
    subroutine psb_z_set_sorted(a,val) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_z_set_sorted
  end interface
  
  interface 
    subroutine psb_z_set_triangle(a,val) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_z_set_triangle
  end interface
  
  interface 
    subroutine psb_z_set_unit(a,val) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_z_set_unit
  end interface
  
  interface 
    subroutine psb_z_set_lower(a,val) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_z_set_lower
  end interface
  
  interface 
    subroutine psb_z_set_upper(a,val) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_z_set_upper
  end interface
  
  interface 
    subroutine psb_z_sparse_print(iout,a,iv,head,ivr,ivc)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      integer(psb_ipk_), intent(in)               :: iout
      class(psb_zspmat_type), intent(in) :: a   
      integer(psb_ipk_), intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer(psb_ipk_), intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_z_sparse_print
  end interface

  interface 
    subroutine psb_z_n_sparse_print(fname,a,iv,head,ivr,ivc)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      character(len=*), intent(in)      :: fname
      class(psb_zspmat_type), intent(in) :: a   
      integer(psb_ipk_), intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer(psb_ipk_), intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_z_n_sparse_print
  end interface
  
  interface 
    subroutine psb_z_get_neigh(a,idx,neigh,n,info,lev)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(in) :: a   
      integer(psb_ipk_), intent(in)                :: idx 
      integer(psb_ipk_), intent(out)               :: n   
      integer(psb_ipk_), allocatable, intent(out)  :: neigh(:)
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_ipk_), optional, intent(in)      :: lev 
    end subroutine psb_z_get_neigh
  end interface
  
  interface 
    subroutine psb_z_csall(nr,nc,a,info,nz) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in)             :: nr,nc
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: nz
    end subroutine psb_z_csall
  end interface
  
  interface 
    subroutine psb_z_reallocate_nz(nz,a) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      integer(psb_ipk_), intent(in) :: nz
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_reallocate_nz
  end interface
  
  interface 
    subroutine psb_z_free(a) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_free
  end interface
  
  interface 
    subroutine psb_z_trim(a) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_trim
  end interface
  
  interface 
    subroutine psb_z_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: val(:)
      integer(psb_ipk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: gtl(:)
    end subroutine psb_z_csput_a
  end interface

  
  interface 
    subroutine psb_z_csput_v(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      use psb_z_vect_mod, only : psb_z_vect_type
      use psb_i_vect_mod, only : psb_i_vect_type
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      type(psb_z_vect_type), intent(inout)  :: val
      type(psb_i_vect_type), intent(inout)  :: ia, ja
      integer(psb_ipk_), intent(in)             :: nz, imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: gtl(:)
    end subroutine psb_z_csput_v
  end interface
 
  interface 
    subroutine psb_z_csgetptn(imin,imax,a,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_z_csgetptn
  end interface
  
  interface 
    subroutine psb_z_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_z_csgetrow
  end interface
  
  interface 
    subroutine psb_z_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in)    :: a
      class(psb_zspmat_type), intent(inout) :: b
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_z_csgetblk
  end interface
  
  interface 
    subroutine psb_z_tril(a,l,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,u)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in)      :: a
      class(psb_zspmat_type), intent(inout)   :: l
      integer(psb_ipk_),intent(out)           :: info
      integer(psb_ipk_), intent(in), optional :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional           :: rscale,cscale
      class(psb_zspmat_type), optional, intent(inout)   :: u
    end subroutine psb_z_tril
  end interface
  
  interface 
    subroutine psb_z_triu(a,u,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,l)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in)      :: a
      class(psb_zspmat_type), intent(inout)   :: u
      integer(psb_ipk_),intent(out)           :: info
      integer(psb_ipk_), intent(in), optional :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional           :: rscale,cscale
      class(psb_zspmat_type), optional, intent(inout)   :: l
    end subroutine psb_z_triu
  end interface


  interface 
    subroutine psb_z_csclip(a,b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      class(psb_zspmat_type), intent(inout) :: b
      integer(psb_ipk_),intent(out)                  :: info
      integer(psb_ipk_), intent(in), optional        :: imin,imax,jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_z_csclip
  end interface
  
  interface 
    subroutine psb_z_b_csclip(a,b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_z_coo_sparse_mat
      class(psb_zspmat_type), intent(in) :: a
      type(psb_z_coo_sparse_mat), intent(out) :: b
      integer(psb_ipk_),intent(out)                  :: info
      integer(psb_ipk_), intent(in), optional        :: imin,imax,jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_z_b_csclip
  end interface
  
  interface 
    subroutine psb_z_mold(a,b)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_z_base_sparse_mat
      class(psb_zspmat_type), intent(inout)     :: a
      class(psb_z_base_sparse_mat), allocatable, intent(out) :: b
    end subroutine psb_z_mold
  end interface
  
  interface 
    subroutine psb_z_asb(a,mold) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_z_base_sparse_mat
      class(psb_zspmat_type), intent(inout) :: a
      class(psb_z_base_sparse_mat), optional, intent(in) :: mold
    end subroutine psb_z_asb
  end interface
  
  interface 
    subroutine psb_z_transp_1mat(a)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_transp_1mat
  end interface
  
  interface 
    subroutine psb_z_transp_2mat(a,b)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(in)  :: a
      class(psb_zspmat_type), intent(inout) :: b
    end subroutine psb_z_transp_2mat
  end interface
  
  interface 
    subroutine psb_z_transc_1mat(a)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_transc_1mat
  end interface
  
  interface 
    subroutine psb_z_transc_2mat(a,b)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(in)  :: a
      class(psb_zspmat_type), intent(inout) :: b
    end subroutine psb_z_transc_2mat
  end interface
  
  interface 
    subroutine psb_z_reinit(a,clear)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_z_reinit
    
  end interface


  !
  ! These methods are specific to the outer SPMAT_TYPE level, since
  ! they tamper with the inner BASE_SPARSE_MAT object.
  !
  !

  !
  ! CSCNV: switches to a different internal derived type.
  !   3 versions: copying to target
  !               copying to a base_sparse_mat object.
  !               in place
  !                 
  !
  interface 
    subroutine psb_z_cscnv(a,b,info,type,mold,upd,dupl)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_z_base_sparse_mat
      class(psb_zspmat_type), intent(in)      :: a
      class(psb_zspmat_type), intent(inout)   :: b
      integer(psb_ipk_), intent(out)                   :: info
      integer(psb_ipk_),optional, intent(in)           :: dupl, upd
      character(len=*), optional, intent(in) :: type
      class(psb_z_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_z_cscnv
  end interface
  

  interface 
    subroutine psb_z_cscnv_ip(a,iinfo,type,mold,dupl)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_z_base_sparse_mat
      class(psb_zspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(out)                   :: iinfo
      integer(psb_ipk_),optional, intent(in)           :: dupl
      character(len=*), optional, intent(in) :: type
      class(psb_z_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_z_cscnv_ip
  end interface
  

  interface 
    subroutine psb_z_cscnv_base(a,b,info,dupl)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_z_base_sparse_mat
      class(psb_zspmat_type), intent(in)       :: a
      class(psb_z_base_sparse_mat), intent(out) :: b
      integer(psb_ipk_), intent(out)                   :: info
      integer(psb_ipk_),optional, intent(in)           :: dupl
    end subroutine psb_z_cscnv_base
  end interface
  
  !
  ! Produce a version of the matrix with diagonal cut
  ! out; passes through a COO buffer. 
  !
  interface 
    subroutine psb_z_clip_d(a,b,info)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(in)    :: a
      class(psb_zspmat_type), intent(inout) :: b
      integer(psb_ipk_),intent(out)                  :: info
    end subroutine psb_z_clip_d
  end interface
  
  interface 
    subroutine psb_z_clip_d_ip(a,info)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      integer(psb_ipk_),intent(out)                  :: info
    end subroutine psb_z_clip_d_ip
  end interface
  
  !
  ! These four interfaces cut through the
  ! encapsulation between spmat_type and base_sparse_mat.
  !
  interface 
    subroutine psb_z_mv_from(a,b)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_z_base_sparse_mat
      class(psb_zspmat_type), intent(inout) :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
    end subroutine psb_z_mv_from
  end interface
  
  interface 
    subroutine psb_z_cp_from(a,b)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_z_base_sparse_mat
      class(psb_zspmat_type), intent(out) :: a
      class(psb_z_base_sparse_mat), intent(in) :: b
    end subroutine psb_z_cp_from
  end interface
  
  interface 
    subroutine psb_z_mv_to(a,b)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_z_base_sparse_mat
      class(psb_zspmat_type), intent(inout) :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
    end subroutine psb_z_mv_to
  end interface
  
  interface 
    subroutine psb_z_cp_to(a,b)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_z_base_sparse_mat    
      class(psb_zspmat_type), intent(in) :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
    end subroutine psb_z_cp_to
  end interface
  !
  ! Mixed type conversions
  !
  interface 
    subroutine psb_z_mv_from_lb(a,b)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_lz_base_sparse_mat
      class(psb_zspmat_type), intent(inout) :: a
      class(psb_lz_base_sparse_mat), intent(inout) :: b
    end subroutine psb_z_mv_from_lb
  end interface
  
  interface 
    subroutine psb_z_cp_from_lb(a,b)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_lz_base_sparse_mat
      class(psb_zspmat_type), intent(out) :: a
      class(psb_lz_base_sparse_mat), intent(in) :: b
    end subroutine psb_z_cp_from_lb
  end interface
  
  interface 
    subroutine psb_z_mv_to_lb(a,b)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_lz_base_sparse_mat
      class(psb_zspmat_type), intent(inout) :: a
      class(psb_lz_base_sparse_mat), intent(inout) :: b
    end subroutine psb_z_mv_to_lb
  end interface
  
  interface 
    subroutine psb_z_cp_to_lb(a,b)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_lz_base_sparse_mat    
      class(psb_zspmat_type), intent(in) :: a
      class(psb_lz_base_sparse_mat), intent(inout) :: b
    end subroutine psb_z_cp_to_lb
  end interface

    interface 
    subroutine psb_z_mv_from_l(a,b)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_lzspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      class(psb_lzspmat_type), intent(inout) :: b
    end subroutine psb_z_mv_from_l
  end interface
  
  interface 
    subroutine psb_z_cp_from_l(a,b)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_lzspmat_type
      class(psb_zspmat_type), intent(out) :: a
      class(psb_lzspmat_type), intent(in) :: b
    end subroutine psb_z_cp_from_l
  end interface
  
  interface 
    subroutine psb_z_mv_to_l(a,b)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_lzspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      class(psb_lzspmat_type), intent(inout) :: b
    end subroutine psb_z_mv_to_l
  end interface
  
  interface 
    subroutine psb_z_cp_to_l(a,b)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_, psb_lzspmat_type
      class(psb_zspmat_type), intent(in) :: a
      class(psb_lzspmat_type), intent(inout) :: b
    end subroutine psb_z_cp_to_l
  end interface

  !
  ! Transfer the internal allocation to the target.
  !  
  interface psb_move_alloc 
    subroutine psb_zspmat_type_move(a,b,info)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      class(psb_zspmat_type), intent(inout)   :: b
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine psb_zspmat_type_move
  end interface
  
  interface 
    subroutine psb_zspmat_clone(a,b,info)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      class(psb_zspmat_type), intent(inout) :: b
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_zspmat_clone
  end interface




  ! == ===================================
  !
  !
  !
  ! Computational routines
  !
  !
  !
  !
  !
  !
  ! == ===================================

  interface psb_csmm
    subroutine psb_z_csmm(alpha,a,x,beta,y,info,trans) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      complex(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
      complex(psb_dpk_), intent(inout) :: y(:,:)
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_z_csmm
    subroutine psb_z_csmv(alpha,a,x,beta,y,info,trans) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      complex(psb_dpk_), intent(in)    :: alpha, beta, x(:)
      complex(psb_dpk_), intent(inout) :: y(:)
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_z_csmv
    subroutine psb_z_csmv_vect(alpha,a,x,beta,y,info,trans) 
      use psb_z_vect_mod, only : psb_z_vect_type
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in)   :: a
      complex(psb_dpk_), intent(in)        :: alpha, beta
      type(psb_z_vect_type), intent(inout) :: x
      type(psb_z_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)                 :: info
      character, optional, intent(in)      :: trans
    end subroutine psb_z_csmv_vect
  end interface
  
  interface psb_cssm
    subroutine psb_z_cssm(alpha,a,x,beta,y,info,trans,scale,d) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      complex(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
      complex(psb_dpk_), intent(inout) :: y(:,:)
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in) :: trans, scale
      complex(psb_dpk_), intent(in), optional :: d(:)
    end subroutine psb_z_cssm
    subroutine psb_z_cssv(alpha,a,x,beta,y,info,trans,scale,d) 
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      complex(psb_dpk_), intent(in)    :: alpha, beta, x(:)
      complex(psb_dpk_), intent(inout) :: y(:)
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in) :: trans, scale
      complex(psb_dpk_), intent(in), optional :: d(:)
    end subroutine psb_z_cssv
    subroutine psb_z_cssv_vect(alpha,a,x,beta,y,info,trans,scale,d) 
      use psb_z_vect_mod, only : psb_z_vect_type
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in)   :: a
      complex(psb_dpk_), intent(in)        :: alpha, beta
      type(psb_z_vect_type), intent(inout) :: x
      type(psb_z_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)                 :: info
      character, optional, intent(in)      :: trans, scale
      type(psb_z_vect_type), optional, intent(inout)   :: d
    end subroutine psb_z_cssv_vect
  end interface
  
  interface 
    function psb_z_maxval(a) result(res)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_z_maxval
  end interface
  
  interface 
    function psb_z_csnmi(a) result(res)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_z_csnmi
  end interface
  
  interface 
    function psb_z_csnm1(a) result(res)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_z_csnm1
  end interface

  interface 
    function psb_z_rowsum(a,info) result(d)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      complex(psb_dpk_), allocatable      :: d(:)
      integer(psb_ipk_), intent(out)               :: info
    end function psb_z_rowsum
  end interface

  interface 
    function psb_z_arwsum(a,info) result(d)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      real(psb_dpk_), allocatable        :: d(:)
      integer(psb_ipk_), intent(out)               :: info
    end function psb_z_arwsum
  end interface
  
  interface 
    function psb_z_colsum(a,info) result(d)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      complex(psb_dpk_), allocatable      :: d(:)
      integer(psb_ipk_), intent(out)               :: info
    end function psb_z_colsum
  end interface

  interface 
    function psb_z_aclsum(a,info)  result(d)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      real(psb_dpk_), allocatable        :: d(:)
      integer(psb_ipk_), intent(out)        :: info
    end function psb_z_aclsum
  end interface

  interface 
    function psb_z_get_diag(a,info) result(d)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      complex(psb_dpk_), allocatable         :: d(:)
      integer(psb_ipk_), intent(out)       :: info
    end function psb_z_get_diag
  end interface
  
  interface psb_scal
    subroutine psb_z_scal(d,a,info,side)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(inout) :: a
      complex(psb_dpk_), intent(in)             :: d(:)
      integer(psb_ipk_), intent(out)                    :: info
      character, intent(in), optional :: side
    end subroutine psb_z_scal
    subroutine psb_z_scals(d,a,info)
      import :: psb_ipk_, psb_lpk_, psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(inout) :: a
      complex(psb_dpk_), intent(in)             :: d
      integer(psb_ipk_), intent(out)                    :: info
    end subroutine psb_z_scals
  end interface


  ! == ===================================
  !
  !
  !
  ! Setters 
  !
  !
  !
  !
  !
  !
  ! == ===================================


  interface 
    subroutine  psb_lz_set_nrows(m,a) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
      integer(psb_lpk_), intent(in) :: m
    end subroutine psb_lz_set_nrows
  end interface
  
  interface 
    subroutine psb_lz_set_ncols(n,a) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
      integer(psb_lpk_), intent(in) :: n
    end subroutine psb_lz_set_ncols
  end interface
  
  interface 
    subroutine  psb_lz_set_dupl(n,a) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in) :: n
    end subroutine psb_lz_set_dupl
  end interface
  
  interface 
    subroutine psb_lz_set_null(a) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
    end subroutine psb_lz_set_null
  end interface
  
  interface 
    subroutine psb_lz_set_bld(a) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
    end subroutine psb_lz_set_bld
  end interface
  
  interface 
    subroutine psb_lz_set_upd(a) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
    end subroutine psb_lz_set_upd
  end interface
  
  interface 
    subroutine psb_lz_set_asb(a) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
    end subroutine psb_lz_set_asb
  end interface
  
  interface 
    subroutine psb_lz_set_sorted(a,val) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_lz_set_sorted
  end interface
  
  interface 
    subroutine psb_lz_set_triangle(a,val) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_lz_set_triangle
  end interface
  
  interface 
    subroutine psb_lz_set_unit(a,val) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_lz_set_unit
  end interface
  
  interface 
    subroutine psb_lz_set_lower(a,val) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_lz_set_lower
  end interface
  
  interface 
    subroutine psb_lz_set_upper(a,val) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_lz_set_upper
  end interface
  
  interface 
    subroutine psb_lz_sparse_print(iout,a,iv,head,ivr,ivc)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      integer(psb_ipk_), intent(in)               :: iout
      class(psb_lzspmat_type), intent(in) :: a   
      integer(psb_lpk_), intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer(psb_lpk_), intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_lz_sparse_print
  end interface

  interface 
    subroutine psb_lz_n_sparse_print(fname,a,iv,head,ivr,ivc)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      character(len=*), intent(in)      :: fname
      class(psb_lzspmat_type), intent(in) :: a   
      integer(psb_lpk_), intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer(psb_lpk_), intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_lz_n_sparse_print
  end interface
  
  interface 
    subroutine psb_lz_get_neigh(a,idx,neigh,n,info,lev)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(in) :: a   
      integer(psb_lpk_), intent(in)                :: idx 
      integer(psb_lpk_), intent(out)               :: n   
      integer(psb_lpk_), allocatable, intent(out)  :: neigh(:)
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_lpk_), optional, intent(in)      :: lev 
    end subroutine psb_lz_get_neigh
  end interface
  
  interface 
    subroutine psb_lz_csall(nr,nc,a,info,nz) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
      integer(psb_lpk_), intent(in)             :: nr,nc
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_lpk_), intent(in), optional   :: nz
    end subroutine psb_lz_csall
  end interface
  
  interface 
    subroutine psb_lz_reallocate_nz(nz,a) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      integer(psb_lpk_), intent(in) :: nz
      class(psb_lzspmat_type), intent(inout) :: a
    end subroutine psb_lz_reallocate_nz
  end interface
  
  interface 
    subroutine psb_lz_free(a) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
    end subroutine psb_lz_free
  end interface
  
  interface 
    subroutine psb_lz_trim(a) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
    end subroutine psb_lz_trim
  end interface
  
  interface 
    subroutine psb_lz_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_
      class(psb_lzspmat_type), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: val(:)
      integer(psb_lpk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_lpk_), intent(in), optional   :: gtl(:)
    end subroutine psb_lz_csput_a
  end interface

  
  interface 
    subroutine psb_lz_csput_v(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      use psb_z_vect_mod, only : psb_z_vect_type
      use psb_l_vect_mod, only : psb_l_vect_type
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
      type(psb_z_vect_type), intent(inout)  :: val
      type(psb_l_vect_type), intent(inout)  :: ia, ja
      integer(psb_lpk_), intent(in)             :: nz, imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_lpk_), intent(in), optional   :: gtl(:)
    end subroutine psb_lz_csput_v
  end interface
 
  interface 
    subroutine psb_lz_csgetptn(imin,imax,a,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_
      class(psb_lzspmat_type), intent(in) :: a
      integer(psb_lpk_), intent(in)                  :: imin,imax
      integer(psb_lpk_), intent(out)                 :: nz
      integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_lpk_), intent(in), optional        :: iren(:)
      integer(psb_lpk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_lz_csgetptn
  end interface
  
  interface 
    subroutine psb_lz_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_
      class(psb_lzspmat_type), intent(in) :: a
      integer(psb_lpk_), intent(in)                  :: imin,imax
      integer(psb_lpk_), intent(out)                 :: nz
      integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_lpk_), intent(in), optional        :: iren(:)
      integer(psb_lpk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_lz_csgetrow
  end interface
  
  interface 
    subroutine psb_lz_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_
      class(psb_lzspmat_type), intent(in)    :: a
      class(psb_lzspmat_type), intent(inout) :: b
      integer(psb_lpk_), intent(in)                  :: imin,imax
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_lpk_), intent(in), optional        :: iren(:)
      integer(psb_lpk_), intent(in), optional        :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_lz_csgetblk
  end interface
  
  interface 
    subroutine psb_lz_tril(a,l,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,u)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_
      class(psb_lzspmat_type), intent(in)      :: a
      class(psb_lzspmat_type), intent(inout)   :: l
      integer(psb_ipk_),intent(out)           :: info
      integer(psb_lpk_), intent(in), optional :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional           :: rscale,cscale
      class(psb_lzspmat_type), optional, intent(inout)   :: u
    end subroutine psb_lz_tril
  end interface
  
  interface 
    subroutine psb_lz_triu(a,u,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,l)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_
      class(psb_lzspmat_type), intent(in)      :: a
      class(psb_lzspmat_type), intent(inout)   :: u
      integer(psb_ipk_),intent(out)           :: info
      integer(psb_lpk_), intent(in), optional :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional           :: rscale,cscale
      class(psb_lzspmat_type), optional, intent(inout)   :: l
    end subroutine psb_lz_triu
  end interface


  interface 
    subroutine psb_lz_csclip(a,b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_
      class(psb_lzspmat_type), intent(in) :: a
      class(psb_lzspmat_type), intent(inout) :: b
      integer(psb_ipk_),intent(out)                  :: info
      integer(psb_lpk_), intent(in), optional        :: imin,imax,jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_lz_csclip
  end interface
  
  interface 
    subroutine psb_lz_b_csclip(a,b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_lz_coo_sparse_mat
      class(psb_lzspmat_type), intent(in) :: a
      type(psb_lz_coo_sparse_mat), intent(out) :: b
      integer(psb_ipk_),intent(out)                  :: info
      integer(psb_lpk_), intent(in), optional        :: imin,imax,jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_lz_b_csclip
  end interface
  
  interface 
    subroutine psb_lz_mold(a,b)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_lz_base_sparse_mat
      class(psb_lzspmat_type), intent(inout)     :: a
      class(psb_lz_base_sparse_mat), allocatable, intent(out) :: b
    end subroutine psb_lz_mold
  end interface
  
  interface 
    subroutine psb_lz_asb(a,mold) 
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_lz_base_sparse_mat
      class(psb_lzspmat_type), intent(inout) :: a
      class(psb_lz_base_sparse_mat), optional, intent(in) :: mold
    end subroutine psb_lz_asb
  end interface
  
  interface 
    subroutine psb_lz_transp_1mat(a)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
    end subroutine psb_lz_transp_1mat
  end interface
  
  interface 
    subroutine psb_lz_transp_2mat(a,b)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(in)  :: a
      class(psb_lzspmat_type), intent(inout) :: b
    end subroutine psb_lz_transp_2mat
  end interface
  
  interface 
    subroutine psb_lz_transc_1mat(a)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
    end subroutine psb_lz_transc_1mat
  end interface
  
  interface 
    subroutine psb_lz_transc_2mat(a,b)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(in)  :: a
      class(psb_lzspmat_type), intent(inout) :: b
    end subroutine psb_lz_transc_2mat
  end interface
  
  interface 
    subroutine psb_lz_reinit(a,clear)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_lz_reinit
    
  end interface


  !
  ! These methods are specific to the outer SPMAT_TYPE level, since
  ! they tamper with the inner BASE_SPARSE_MAT object.
  !
  !

  !
  ! CSCNV: switches to a different internal derived type.
  !   3 versions: copying to target
  !               copying to a base_sparse_mat object.
  !               in place
  !                 
  !
  interface 
    subroutine psb_lz_cscnv(a,b,info,type,mold,upd,dupl)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_lz_base_sparse_mat
      class(psb_lzspmat_type), intent(in)      :: a
      class(psb_lzspmat_type), intent(inout)   :: b
      integer(psb_ipk_), intent(out)                   :: info
      integer(psb_ipk_),optional, intent(in)           :: dupl, upd
      character(len=*), optional, intent(in) :: type
      class(psb_lz_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_lz_cscnv
  end interface
  

  interface 
    subroutine psb_lz_cscnv_ip(a,iinfo,type,mold,dupl)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_lz_base_sparse_mat
      class(psb_lzspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(out)                   :: iinfo
      integer(psb_ipk_),optional, intent(in)           :: dupl
      character(len=*), optional, intent(in) :: type
      class(psb_lz_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_lz_cscnv_ip
  end interface
  

  interface 
    subroutine psb_lz_cscnv_base(a,b,info,dupl)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_lz_base_sparse_mat
      class(psb_lzspmat_type), intent(in)       :: a
      class(psb_lz_base_sparse_mat), intent(out) :: b
      integer(psb_ipk_), intent(out)                   :: info
      integer(psb_ipk_),optional, intent(in)           :: dupl
    end subroutine psb_lz_cscnv_base
  end interface
  
  
  !
  ! These four interfaces cut through the
  ! encapsulation between spmat_type and base_sparse_mat.
  !
  interface 
    subroutine psb_lz_mv_from(a,b)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_lz_base_sparse_mat
      class(psb_lzspmat_type), intent(inout) :: a
      class(psb_lz_base_sparse_mat), intent(inout) :: b
    end subroutine psb_lz_mv_from
  end interface
  
  interface 
    subroutine psb_lz_cp_from(a,b)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_lz_base_sparse_mat
      class(psb_lzspmat_type), intent(out) :: a
      class(psb_lz_base_sparse_mat), intent(in) :: b
    end subroutine psb_lz_cp_from
  end interface
  
  interface 
    subroutine psb_lz_mv_to(a,b)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_lz_base_sparse_mat
      class(psb_lzspmat_type), intent(inout) :: a
      class(psb_lz_base_sparse_mat), intent(inout) :: b
    end subroutine psb_lz_mv_to
  end interface
  
  interface 
    subroutine psb_lz_cp_to(a,b)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_lz_base_sparse_mat    
      class(psb_lzspmat_type), intent(in) :: a
      class(psb_lz_base_sparse_mat), intent(inout) :: b
    end subroutine psb_lz_cp_to
  end interface
  !
  ! Mixed type conversions
  !
  interface 
    subroutine psb_lz_mv_from_ib(a,b)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_z_base_sparse_mat
      class(psb_lzspmat_type), intent(inout) :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
    end subroutine psb_lz_mv_from_ib
  end interface
  
  interface 
    subroutine psb_lz_cp_from_ib(a,b)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_z_base_sparse_mat
      class(psb_lzspmat_type), intent(out) :: a
      class(psb_z_base_sparse_mat), intent(in) :: b
    end subroutine psb_lz_cp_from_ib
  end interface
  
  interface 
    subroutine psb_lz_mv_to_ib(a,b)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_z_base_sparse_mat
      class(psb_lzspmat_type), intent(inout) :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
    end subroutine psb_lz_mv_to_ib
  end interface
  
  interface 
    subroutine psb_lz_cp_to_ib(a,b)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_z_base_sparse_mat    
      class(psb_lzspmat_type), intent(in) :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
    end subroutine psb_lz_cp_to_ib
  end interface

    interface 
    subroutine psb_lz_mv_from_i(a,b)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_zspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
      class(psb_zspmat_type), intent(inout) :: b
    end subroutine psb_lz_mv_from_i
  end interface
  
  interface 
    subroutine psb_lz_cp_from_i(a,b)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_zspmat_type
      class(psb_lzspmat_type), intent(out) :: a
      class(psb_zspmat_type), intent(in) :: b
    end subroutine psb_lz_cp_from_i
  end interface
  
  interface 
    subroutine psb_lz_mv_to_i(a,b)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_zspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
      class(psb_zspmat_type), intent(inout) :: b
    end subroutine psb_lz_mv_to_i
  end interface
  
  interface 
    subroutine psb_lz_cp_to_i(a,b)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_, psb_zspmat_type
      class(psb_lzspmat_type), intent(in) :: a
      class(psb_zspmat_type), intent(inout) :: b
    end subroutine psb_lz_cp_to_i
  end interface

  
  !
  ! Transfer the internal allocation to the target.
  !  
  interface psb_move_alloc 
    subroutine psb_lzspmat_type_move(a,b,info)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
      class(psb_lzspmat_type), intent(inout)   :: b
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine psb_lzspmat_type_move
  end interface
  
  interface 
    subroutine psb_lzspmat_clone(a,b,info)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type
      class(psb_lzspmat_type), intent(inout) :: a
      class(psb_lzspmat_type), intent(inout) :: b
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_lzspmat_clone
  end interface



  interface 
    function psb_lz_get_diag(a,info) result(d)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_
      class(psb_lzspmat_type), intent(in) :: a
      complex(psb_dpk_), allocatable         :: d(:)
      integer(psb_ipk_), intent(out)       :: info
    end function psb_lz_get_diag
  end interface
  
  interface psb_scal
    subroutine psb_lz_scal(d,a,info,side)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_
      class(psb_lzspmat_type), intent(inout) :: a
      complex(psb_dpk_), intent(in)             :: d(:)
      integer(psb_ipk_), intent(out)                    :: info
      character, intent(in), optional :: side
    end subroutine psb_lz_scal
    subroutine psb_lz_scals(d,a,info)
      import :: psb_ipk_, psb_lpk_, psb_lzspmat_type, psb_dpk_
      class(psb_lzspmat_type), intent(inout) :: a
      complex(psb_dpk_), intent(in)             :: d
      integer(psb_ipk_), intent(out)                    :: info
    end subroutine psb_lz_scals
  end interface




  
contains 


  
  subroutine  psb_z_set_mat_default(a) 
    implicit none 
    class(psb_z_base_sparse_mat), intent(in) :: a
    
    if (allocated(psb_z_base_mat_default)) then 
      deallocate(psb_z_base_mat_default)
    end if
    allocate(psb_z_base_mat_default, mold=a)

  end subroutine psb_z_set_mat_default
  
  function psb_z_get_mat_default(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    class(psb_z_base_sparse_mat), pointer :: res
    
    res => psb_z_get_base_mat_default()
    
  end function psb_z_get_mat_default

  
  function psb_z_get_base_mat_default() result(res)
    implicit none 
    class(psb_z_base_sparse_mat), pointer :: res
    
    if (.not.allocated(psb_z_base_mat_default)) then 
      allocate(psb_z_csr_sparse_mat :: psb_z_base_mat_default)
    end if

    res => psb_z_base_mat_default
    
  end function psb_z_get_base_mat_default




  ! == ===================================
  !
  !
  !
  ! Getters 
  !
  !
  !
  !
  !
  ! == ===================================

  
  function psb_z_sizeof(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    integer(psb_epk_) :: res
    
    res = 0
    if (allocated(a%a)) then 
      res = a%a%sizeof()
    end if
    
  end function psb_z_sizeof


  function psb_z_get_fmt(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    character(len=5) :: res

    if (allocated(a%a)) then 
      res = a%a%get_fmt()
    else
      res = 'NULL'
    end if

  end function psb_z_get_fmt


  function psb_z_get_dupl(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    integer(psb_ipk_) :: res

    if (allocated(a%a)) then 
      res = a%a%get_dupl()
    else
      res = psb_invalid_
    end if
  end function psb_z_get_dupl

  function psb_z_get_nrows(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    integer(psb_ipk_) :: res

    if (allocated(a%a)) then 
      res = a%a%get_nrows()
    else
      res = 0
    end if

  end function psb_z_get_nrows

  function psb_z_get_ncols(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    integer(psb_ipk_) :: res

    if (allocated(a%a)) then 
      res = a%a%get_ncols()
    else
      res = 0
    end if

  end function psb_z_get_ncols

  function psb_z_is_triangle(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_triangle()
    else
      res = .false.
    end if

  end function psb_z_is_triangle

  function psb_z_is_unit(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_unit()
    else
      res = .false.
    end if

  end function psb_z_is_unit

  function psb_z_is_upper(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upper()
    else
      res = .false.
    end if

  end function psb_z_is_upper

  function psb_z_is_lower(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = .not. a%a%is_upper()
    else
      res = .false.
    end if

  end function psb_z_is_lower

  function psb_z_is_null(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_null() 
    else
      res = .true.
    end if

  end function psb_z_is_null

  function psb_z_is_bld(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_bld()
    else
      res = .false.
    end if

  end function psb_z_is_bld

  function psb_z_is_upd(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upd()
    else
      res = .false.
    end if

  end function psb_z_is_upd

  function psb_z_is_asb(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_asb()
    else
      res = .false.
    end if

  end function psb_z_is_asb

  function psb_z_is_sorted(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_sorted()
    else
      res = .false.
    end if

  end function psb_z_is_sorted

  function psb_z_is_by_rows(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_by_rows()
    else
      res = .false.
    end if

  end function psb_z_is_by_rows

  function psb_z_is_by_cols(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_by_cols()
    else
      res = .false.
    end if

  end function psb_z_is_by_cols


  !
  subroutine z_mat_sync(a)
    implicit none 
    class(psb_zspmat_type), target, intent(in) :: a
    
    if (allocated(a%a))  call a%a%sync()

  end subroutine z_mat_sync

  !
  subroutine z_mat_set_host(a)
    implicit none 
    class(psb_zspmat_type), intent(inout) :: a

    if (allocated(a%a))  call a%a%set_host()
    
  end subroutine z_mat_set_host


  !
  subroutine z_mat_set_dev(a)
    implicit none 
    class(psb_zspmat_type), intent(inout) :: a

    if (allocated(a%a))  call a%a%set_dev()
    
  end subroutine z_mat_set_dev

  !
  subroutine z_mat_set_sync(a)
    implicit none 
    class(psb_zspmat_type), intent(inout) :: a

    if (allocated(a%a))  call a%a%set_sync()
    
  end subroutine z_mat_set_sync

  !
  function z_mat_is_dev(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    logical  :: res
  
    if (allocated(a%a)) then
      res = a%a%is_dev()
    else
      res = .false.
    end if
  end function z_mat_is_dev
  
  !
  function z_mat_is_host(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    logical  :: res

  
    if (allocated(a%a)) then
      res = a%a%is_host()
    else
      res = .true.
    end if
  end function z_mat_is_host

  !
  function z_mat_is_sync(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    logical  :: res

  
    if (allocated(a%a)) then
      res = a%a%is_sync()
    else
      res = .true.
    end if

  end function z_mat_is_sync


  function psb_z_is_repeatable_updates(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_repeatable_updates()
    else
      res = .false.
    end if

  end function psb_z_is_repeatable_updates

  subroutine psb_z_set_repeatable_updates(a,val) 
    implicit none 
    class(psb_zspmat_type), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (allocated(a%a)) then 
      call a%a%set_repeatable_updates(val)
    end if
    
  end subroutine psb_z_set_repeatable_updates


  function psb_z_get_nzeros(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    integer(psb_ipk_) :: res

    res = 0
    if (allocated(a%a)) then 
      res = a%a%get_nzeros()
    end if

  end function psb_z_get_nzeros

  function psb_z_get_size(a) result(res)

    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    integer(psb_ipk_) :: res


    res = 0
    if (allocated(a%a)) then 
      res = a%a%get_size()
    end if

  end function psb_z_get_size


  function psb_z_get_nz_row(idx,a) result(res)
    implicit none 
    integer(psb_ipk_), intent(in)               :: idx
    class(psb_zspmat_type), intent(in) :: a
    integer(psb_ipk_) :: res

    res = 0
    
    if (allocated(a%a)) res = a%a%get_nz_row(idx)

  end function psb_z_get_nz_row

  subroutine psb_z_clean_zeros(a,info)
    implicit none 
    integer(psb_ipk_), intent(out)        :: info
    class(psb_zspmat_type), intent(inout) :: a

    info = 0 
    if (allocated(a%a)) call a%a%clean_zeros(info)

  end subroutine psb_z_clean_zeros

#if defined(IPK4) && defined(LPK8)
  subroutine psb_z_lcsgetptn(imin,imax,a,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in)                  :: imin,imax
    integer(psb_ipk_), intent(out)                 :: nz
    integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
    integer(psb_ipk_),intent(out)                  :: info
    logical, intent(in), optional        :: append
    integer(psb_ipk_), intent(in), optional        :: iren(:)
    integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
    logical, intent(in), optional        :: rscale,cscale

    ! Local
    integer(psb_ipk_), allocatable :: lia(:), lja(:)

    info = psb_success_
    !
    ! Note: in principle we could use reallocate on assignment,
    ! but GCC bug 52162 forces us to take defensive programming.
    !
    if (allocated(ia)) then
      call psb_realloc(size(ia),lia,info)
      if (info == psb_success_) lia(:) = ia(:)
    end if
    if (allocated(ja)) then
      call psb_realloc(size(ja),lja,info)
      if (info == psb_success_) lja(:) = ja(:)
    end if
    call a%csget(imin,imax,nz,lia,lja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
   
    call psb_ensure_size(size(lia),ia,info)
    if (info == psb_success_) ia(:) = lia(:)
    call psb_ensure_size(size(lja),ja,info)
    if (info == psb_success_) ja(:) = lja(:)
    
  end subroutine psb_z_lcsgetptn
  
  subroutine psb_z_lcsgetrow(imin,imax,a,nz,ia,ja,val,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in)                  :: imin,imax
    integer(psb_ipk_), intent(out)                 :: nz
    integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
    complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
    integer(psb_ipk_),intent(out)                  :: info
    logical, intent(in), optional        :: append
    integer(psb_ipk_), intent(in), optional        :: iren(:)
    integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
    logical, intent(in), optional        :: rscale,cscale
    ! Local
    integer(psb_ipk_), allocatable :: lia(:), lja(:)

    !
    ! Note: in principle we could use reallocate on assignment,
    ! but GCC bug 52162 forces us to take defensive programming.
    !
    if (allocated(ia)) then
      call psb_realloc(size(ia),lia,info)
      if (info == psb_success_) lia(:) = ia(:)
    end if
    if (allocated(ja)) then
      call psb_realloc(size(ja),lja,info)
      if (info == psb_success_) lja(:) = ja(:)
    end if

    call a%csget(imin,imax,nz,lia,lja,val,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
   
    call psb_ensure_size(size(lia),ia,info)
    if (info == psb_success_) ia(:) = lia(:)
    call psb_ensure_size(size(lja),ja,info)
    if (info == psb_success_) ja(:) = lja(:)
        
  end subroutine psb_z_lcsgetrow
#endif

  !
  ! lz methods
  !

  
  subroutine  psb_lz_set_mat_default(a) 
    implicit none 
    class(psb_lz_base_sparse_mat), intent(in) :: a
    
    if (allocated(psb_lz_base_mat_default)) then 
      deallocate(psb_lz_base_mat_default)
    end if
    allocate(psb_lz_base_mat_default, mold=a)

  end subroutine psb_lz_set_mat_default
  
  function psb_lz_get_mat_default(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    class(psb_lz_base_sparse_mat), pointer :: res
    
    res => psb_lz_get_base_mat_default()
    
  end function psb_lz_get_mat_default

  
  function psb_lz_get_base_mat_default() result(res)
    implicit none 
    class(psb_lz_base_sparse_mat), pointer :: res
    
    if (.not.allocated(psb_lz_base_mat_default)) then 
      allocate(psb_lz_csr_sparse_mat :: psb_lz_base_mat_default)
    end if

    res => psb_lz_base_mat_default
    
  end function psb_lz_get_base_mat_default




  ! == ===================================
  !
  !
  !
  ! Getters 
  !
  !
  !
  !
  !
  ! == ===================================

  
  function psb_lz_sizeof(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    integer(psb_epk_) :: res
    
    res = 0
    if (allocated(a%a)) then 
      res = a%a%sizeof()
    end if
    
  end function psb_lz_sizeof


  function psb_lz_get_fmt(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    character(len=5) :: res

    if (allocated(a%a)) then 
      res = a%a%get_fmt()
    else
      res = 'NULL'
    end if

  end function psb_lz_get_fmt


  function psb_lz_get_dupl(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    integer(psb_ipk_) :: res

    if (allocated(a%a)) then 
      res = a%a%get_dupl()
    else
      res = psb_invalid_
    end if
  end function psb_lz_get_dupl

  function psb_lz_get_nrows(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    integer(psb_lpk_) :: res

    if (allocated(a%a)) then 
      res = a%a%get_nrows()
    else
      res = 0
    end if

  end function psb_lz_get_nrows

  function psb_lz_get_ncols(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    integer(psb_lpk_) :: res

    if (allocated(a%a)) then 
      res = a%a%get_ncols()
    else
      res = 0
    end if

  end function psb_lz_get_ncols

  function psb_lz_is_triangle(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_triangle()
    else
      res = .false.
    end if

  end function psb_lz_is_triangle

  function psb_lz_is_unit(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_unit()
    else
      res = .false.
    end if

  end function psb_lz_is_unit

  function psb_lz_is_upper(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upper()
    else
      res = .false.
    end if

  end function psb_lz_is_upper

  function psb_lz_is_lower(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = .not. a%a%is_upper()
    else
      res = .false.
    end if

  end function psb_lz_is_lower

  function psb_lz_is_null(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_null() 
    else
      res = .true.
    end if

  end function psb_lz_is_null

  function psb_lz_is_bld(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_bld()
    else
      res = .false.
    end if

  end function psb_lz_is_bld

  function psb_lz_is_upd(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upd()
    else
      res = .false.
    end if

  end function psb_lz_is_upd

  function psb_lz_is_asb(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_asb()
    else
      res = .false.
    end if

  end function psb_lz_is_asb

  function psb_lz_is_sorted(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_sorted()
    else
      res = .false.
    end if

  end function psb_lz_is_sorted

  function psb_lz_is_by_rows(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_by_rows()
    else
      res = .false.
    end if

  end function psb_lz_is_by_rows

  function psb_lz_is_by_cols(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_by_cols()
    else
      res = .false.
    end if

  end function psb_lz_is_by_cols


  !
  subroutine lz_mat_sync(a)
    implicit none 
    class(psb_lzspmat_type), target, intent(in) :: a
    
    if (allocated(a%a))  call a%a%sync()

  end subroutine lz_mat_sync

  !
  subroutine lz_mat_set_host(a)
    implicit none 
    class(psb_lzspmat_type), intent(inout) :: a

    if (allocated(a%a))  call a%a%set_host()
    
  end subroutine lz_mat_set_host


  !
  subroutine lz_mat_set_dev(a)
    implicit none 
    class(psb_lzspmat_type), intent(inout) :: a

    if (allocated(a%a))  call a%a%set_dev()
    
  end subroutine lz_mat_set_dev

  !
  subroutine lz_mat_set_sync(a)
    implicit none 
    class(psb_lzspmat_type), intent(inout) :: a

    if (allocated(a%a))  call a%a%set_sync()
    
  end subroutine lz_mat_set_sync

  !
  function lz_mat_is_dev(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    logical  :: res
  
    if (allocated(a%a)) then
      res = a%a%is_dev()
    else
      res = .false.
    end if
  end function lz_mat_is_dev
  
  !
  function lz_mat_is_host(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    logical  :: res

  
    if (allocated(a%a)) then
      res = a%a%is_host()
    else
      res = .true.
    end if
  end function lz_mat_is_host

  !
  function lz_mat_is_sync(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    logical  :: res

  
    if (allocated(a%a)) then
      res = a%a%is_sync()
    else
      res = .true.
    end if

  end function lz_mat_is_sync


  function psb_lz_is_repeatable_updates(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_repeatable_updates()
    else
      res = .false.
    end if

  end function psb_lz_is_repeatable_updates

  subroutine psb_lz_set_repeatable_updates(a,val) 
    implicit none 
    class(psb_lzspmat_type), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (allocated(a%a)) then 
      call a%a%set_repeatable_updates(val)
    end if
    
  end subroutine psb_lz_set_repeatable_updates


  function psb_lz_get_nzeros(a) result(res)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    integer(psb_lpk_) :: res

    res = 0
    if (allocated(a%a)) then 
      res = a%a%get_nzeros()
    end if

  end function psb_lz_get_nzeros

  function psb_lz_get_size(a) result(res)

    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    integer(psb_lpk_) :: res


    res = 0
    if (allocated(a%a)) then 
      res = a%a%get_size()
    end if

  end function psb_lz_get_size


  function psb_lz_get_nz_row(idx,a) result(res)
    implicit none 
    integer(psb_lpk_), intent(in)               :: idx
    class(psb_lzspmat_type), intent(in) :: a
    integer(psb_lpk_) :: res

    res = 0
    
    if (allocated(a%a)) res = a%a%get_nz_row(idx)

  end function psb_lz_get_nz_row

  subroutine psb_lz_clean_zeros(a,info)
    implicit none 
    integer(psb_ipk_), intent(out)        :: info
    class(psb_lzspmat_type), intent(inout) :: a

    info = 0 
    if (allocated(a%a)) call a%a%clean_zeros(info)

  end subroutine psb_lz_clean_zeros

#if defined(IPK4) && defined(LPK8)
  subroutine psb_lz_icsgetptn(imin,imax,a,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in)                  :: imin,imax
    integer(psb_ipk_), intent(out)                 :: nz
    integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
    integer(psb_ipk_),intent(out)                  :: info
    logical, intent(in), optional        :: append
    integer(psb_ipk_), intent(in), optional        :: iren(:)
    integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
    logical, intent(in), optional        :: rscale,cscale

    ! Local
    integer(psb_ipk_), allocatable :: lia(:), lja(:)

    info = psb_success_
    !
    ! Note: in principle we could use reallocate on assignment,
    ! but GCC bug 52162 forces us to take defensive programming.
    !
    if (allocated(ia)) then
      call psb_realloc(size(ia),lia,info)
      if (info == psb_success_) lia(:) = ia(:)
    end if
    if (allocated(ja)) then
      call psb_realloc(size(ja),lja,info)
      if (info == psb_success_) lja(:) = ja(:)
    end if
    call a%csget(imin,imax,nz,lia,lja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
   
    call psb_ensure_size(size(lia),ia,info)
    if (info == psb_success_) ia(:) = lia(:)
    call psb_ensure_size(size(lja),ja,info)
    if (info == psb_success_) ja(:) = lja(:)
    
  end subroutine psb_lz_icsgetptn
  
  subroutine psb_lz_icsgetrow(imin,imax,a,nz,ia,ja,val,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
    implicit none 
    class(psb_lzspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in)                  :: imin,imax
    integer(psb_ipk_), intent(out)                 :: nz
    integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
    complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
    integer(psb_ipk_),intent(out)                  :: info
    logical, intent(in), optional        :: append
    integer(psb_ipk_), intent(in), optional        :: iren(:)
    integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
    logical, intent(in), optional        :: rscale,cscale
    ! Local
    integer(psb_ipk_), allocatable :: lia(:), lja(:)

    !
    ! Note: in principle we could use reallocate on assignment,
    ! but GCC bug 52162 forces us to take defensive programming.
    !
    if (allocated(ia)) then
      call psb_realloc(size(ia),lia,info)
      if (info == psb_success_) lia(:) = ia(:)
    end if
    if (allocated(ja)) then
      call psb_realloc(size(ja),lja,info)
      if (info == psb_success_) lja(:) = ja(:)
    end if

    call a%csget(imin,imax,nz,lia,lja,val,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
   
    call psb_ensure_size(size(lia),ia,info)
    if (info == psb_success_) ia(:) = lia(:)
    call psb_ensure_size(size(lja),ja,info)
    if (info == psb_success_) ja(:) = lja(:)
        
  end subroutine psb_lz_icsgetrow
#endif

end module psb_z_mat_mod
