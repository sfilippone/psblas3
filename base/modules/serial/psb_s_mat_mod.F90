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
! package: psb_s_mat_mod
!
! This module contains the definition of the psb_s_sparse type which
! is a generic container for a sparse matrix and it is mostly meant to
! provide a mean of switching, at run-time, among different formats,
! potentially unknown at the library compile-time by adding a layer of
! indirection. This type encapsulates the psb_s_base_sparse_mat class
! inside another class which is the one visible to the user. 
! Most methods of the psb_s_mat_mod simply call the methods of the
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


module psb_s_mat_mod

  use psb_s_base_mat_mod
  use psb_s_csr_mat_mod, only : psb_s_csr_sparse_mat
  use psb_s_csc_mat_mod, only : psb_s_csc_sparse_mat

  type :: psb_sspmat_type

    class(psb_s_base_sparse_mat), allocatable  :: a 

  contains
    ! Getters
    procedure, pass(a) :: get_nrows   => psb_s_get_nrows
    procedure, pass(a) :: get_ncols   => psb_s_get_ncols
    procedure, pass(a) :: get_nzeros  => psb_s_get_nzeros
    procedure, pass(a) :: get_nz_row  => psb_s_get_nz_row
    procedure, pass(a) :: get_size    => psb_s_get_size
    procedure, pass(a) :: get_dupl    => psb_s_get_dupl
    procedure, pass(a) :: is_null     => psb_s_is_null
    procedure, pass(a) :: is_bld      => psb_s_is_bld
    procedure, pass(a) :: is_upd      => psb_s_is_upd
    procedure, pass(a) :: is_asb      => psb_s_is_asb
    procedure, pass(a) :: is_sorted   => psb_s_is_sorted
    procedure, pass(a) :: is_by_rows  => psb_s_is_by_rows
    procedure, pass(a) :: is_by_cols  => psb_s_is_by_cols
    procedure, pass(a) :: is_upper    => psb_s_is_upper
    procedure, pass(a) :: is_lower    => psb_s_is_lower
    procedure, pass(a) :: is_triangle => psb_s_is_triangle
    procedure, pass(a) :: is_unit     => psb_s_is_unit
    procedure, pass(a) :: is_repeatable_updates => psb_s_is_repeatable_updates
    procedure, pass(a) :: get_fmt     => psb_s_get_fmt
    procedure, pass(a) :: sizeof      => psb_s_sizeof

    ! Setters
    procedure, pass(a) :: set_nrows    => psb_s_set_nrows
    procedure, pass(a) :: set_ncols    => psb_s_set_ncols
    procedure, pass(a) :: set_dupl     => psb_s_set_dupl
    procedure, pass(a) :: set_null     => psb_s_set_null
    procedure, pass(a) :: set_bld      => psb_s_set_bld
    procedure, pass(a) :: set_upd      => psb_s_set_upd
    procedure, pass(a) :: set_asb      => psb_s_set_asb
    procedure, pass(a) :: set_sorted   => psb_s_set_sorted
    procedure, pass(a) :: set_upper    => psb_s_set_upper
    procedure, pass(a) :: set_lower    => psb_s_set_lower
    procedure, pass(a) :: set_triangle => psb_s_set_triangle
    procedure, pass(a) :: set_unit     => psb_s_set_unit
    procedure, pass(a) :: set_repeatable_updates => psb_s_set_repeatable_updates

    ! Memory/data management 
    procedure, pass(a) :: csall       => psb_s_csall
    procedure, pass(a) :: free        => psb_s_free
    procedure, pass(a) :: trim        => psb_s_trim
    procedure, pass(a) :: csput_a     => psb_s_csput_a
    procedure, pass(a) :: csput_v     => psb_s_csput_v 
    generic, public    :: csput       => csput_a,  csput_v
    procedure, pass(a) :: csgetptn    => psb_s_csgetptn
    procedure, pass(a) :: csgetrow    => psb_s_csgetrow
    procedure, pass(a) :: csgetblk    => psb_s_csgetblk
    generic, public    :: csget       => csgetptn, csgetrow, csgetblk
#if defined(INT_I4_L8)
    procedure, pass(a) :: lcsgetptn    => psb_s_lcsgetptn
    procedure, pass(a) :: lcsgetrow    => psb_s_lcsgetrow
    generic, public    :: csget        => lcsgetptn, lcsgetrow
#endif    
    procedure, pass(a) :: tril        => psb_s_tril
    procedure, pass(a) :: triu        => psb_s_triu
    procedure, pass(a) :: m_csclip    => psb_s_csclip
    procedure, pass(a) :: b_csclip    => psb_s_b_csclip
    generic, public    :: csclip      => b_csclip, m_csclip
    procedure, pass(a) :: clean_zeros => psb_s_clean_zeros
    procedure, pass(a) :: reall       => psb_s_reallocate_nz
    procedure, pass(a) :: get_neigh   => psb_s_get_neigh
    procedure, pass(a) :: reinit      => psb_s_reinit
    procedure, pass(a) :: print_i     => psb_s_sparse_print
    procedure, pass(a) :: print_n     => psb_s_n_sparse_print
    generic, public    :: print       => print_i, print_n
    procedure, pass(a) :: mold        => psb_s_mold
    procedure, pass(a) :: asb         => psb_s_asb
    procedure, pass(a) :: transp_1mat => psb_s_transp_1mat
    procedure, pass(a) :: transp_2mat => psb_s_transp_2mat
    generic, public    :: transp      => transp_1mat, transp_2mat
    procedure, pass(a) :: transc_1mat => psb_s_transc_1mat
    procedure, pass(a) :: transc_2mat => psb_s_transc_2mat
    generic, public    :: transc      => transc_1mat, transc_2mat

    !
    ! Sync: centerpiece of handling of external storage.
    ! Any derived class having extra storage upon sync
    ! will guarantee that both fortran/host side and
    ! external side contain the same data. The base
    ! version is only a placeholder. 
    !
    procedure, pass(a) :: sync        => s_mat_sync
    procedure, pass(a) :: is_host     => s_mat_is_host
    procedure, pass(a) :: is_dev      => s_mat_is_dev
    procedure, pass(a) :: is_sync     => s_mat_is_sync
    procedure, pass(a) :: set_host    => s_mat_set_host
    procedure, pass(a) :: set_dev     => s_mat_set_dev
    procedure, pass(a) :: set_sync    => s_mat_set_sync


    ! These are specific to this level of encapsulation.
    procedure, pass(a) :: mv_from_b   => psb_s_mv_from
    generic, public    :: mv_from     => mv_from_b
    procedure, pass(a) :: mv_to_b     => psb_s_mv_to
    generic, public    :: mv_to       => mv_to_b
    procedure, pass(a) :: cp_from_b   => psb_s_cp_from
    generic, public    :: cp_from     => cp_from_b
    procedure, pass(a) :: cp_to_b     => psb_s_cp_to
    generic, public    :: cp_to       => cp_to_b
    procedure, pass(a) :: clip_d_ip   => psb_s_clip_d_ip
    procedure, pass(a) :: clip_d      => psb_s_clip_d
    generic, public    :: clip_diag   => clip_d_ip, clip_d
    procedure, pass(a) :: cscnv_np    => psb_s_cscnv
    procedure, pass(a) :: cscnv_ip    => psb_s_cscnv_ip
    procedure, pass(a) :: cscnv_base  => psb_s_cscnv_base
    generic, public    :: cscnv       => cscnv_np, cscnv_ip, cscnv_base
    procedure, pass(a) :: clone       => psb_sspmat_clone

    ! Computational routines 
    procedure, pass(a) :: get_diag => psb_s_get_diag
    procedure, pass(a) :: maxval   => psb_s_maxval
    procedure, pass(a) :: spnmi    => psb_s_csnmi
    procedure, pass(a) :: spnm1    => psb_s_csnm1
    procedure, pass(a) :: rowsum   => psb_s_rowsum
    procedure, pass(a) :: arwsum   => psb_s_arwsum
    procedure, pass(a) :: colsum   => psb_s_colsum
    procedure, pass(a) :: aclsum   => psb_s_aclsum
    procedure, pass(a) :: csmv_v   => psb_s_csmv_vect
    procedure, pass(a) :: csmv     => psb_s_csmv
    procedure, pass(a) :: csmm     => psb_s_csmm
    generic, public    :: spmm     => csmm, csmv, csmv_v
    procedure, pass(a) :: scals    => psb_s_scals
    procedure, pass(a) :: scalv    => psb_s_scal
    generic, public    :: scal     => scals, scalv
    procedure, pass(a) :: cssv_v   => psb_s_cssv_vect
    procedure, pass(a) :: cssv     => psb_s_cssv
    procedure, pass(a) :: cssm     => psb_s_cssm
    generic, public    :: spsm     => cssm, cssv, cssv_v

  end type psb_sspmat_type

  private :: psb_s_get_nrows, psb_s_get_ncols, &
       & psb_s_get_nzeros, psb_s_get_size, &
       & psb_s_get_dupl, psb_s_is_null, psb_s_is_bld, &
       & psb_s_is_upd, psb_s_is_asb, psb_s_is_sorted, &
       & psb_s_is_by_rows, psb_s_is_by_cols, psb_s_is_upper, &
       & psb_s_is_lower, psb_s_is_triangle, psb_s_get_nz_row, &
       & s_mat_sync, s_mat_is_host, s_mat_is_dev, &
       & s_mat_is_sync, s_mat_set_host, s_mat_set_dev,&
       & s_mat_set_sync



  class(psb_s_base_sparse_mat), allocatable, target, &
       & save, private :: psb_s_base_mat_default

  interface psb_set_mat_default
    module procedure psb_s_set_mat_default
  end interface

  interface psb_get_mat_default
    module procedure psb_s_get_mat_default
  end interface

  interface psb_sizeof
    module procedure psb_s_sizeof
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
    subroutine  psb_s_set_nrows(m,a) 
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in) :: m
    end subroutine psb_s_set_nrows
  end interface
  
  interface 
    subroutine psb_s_set_ncols(n,a) 
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in) :: n
    end subroutine psb_s_set_ncols
  end interface
  
  interface 
    subroutine  psb_s_set_dupl(n,a) 
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in) :: n
    end subroutine psb_s_set_dupl
  end interface
  
  interface 
    subroutine psb_s_set_null(a) 
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
    end subroutine psb_s_set_null
  end interface
  
  interface 
    subroutine psb_s_set_bld(a) 
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
    end subroutine psb_s_set_bld
  end interface
  
  interface 
    subroutine psb_s_set_upd(a) 
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
    end subroutine psb_s_set_upd
  end interface
  
  interface 
    subroutine psb_s_set_asb(a) 
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
    end subroutine psb_s_set_asb
  end interface
  
  interface 
    subroutine psb_s_set_sorted(a,val) 
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_s_set_sorted
  end interface
  
  interface 
    subroutine psb_s_set_triangle(a,val) 
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_s_set_triangle
  end interface
  
  interface 
    subroutine psb_s_set_unit(a,val) 
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_s_set_unit
  end interface
  
  interface 
    subroutine psb_s_set_lower(a,val) 
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_s_set_lower
  end interface
  
  interface 
    subroutine psb_s_set_upper(a,val) 
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_s_set_upper
  end interface
  
  interface 
    subroutine psb_s_sparse_print(iout,a,iv,head,ivr,ivc)
      import :: psb_ipk_, psb_sspmat_type
      integer(psb_ipk_), intent(in)               :: iout
      class(psb_sspmat_type), intent(in) :: a   
      integer(psb_ipk_), intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer(psb_ipk_), intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_s_sparse_print
  end interface

  interface 
    subroutine psb_s_n_sparse_print(fname,a,iv,head,ivr,ivc)
      import :: psb_ipk_, psb_sspmat_type
      character(len=*), intent(in)      :: fname
      class(psb_sspmat_type), intent(in) :: a   
      integer(psb_ipk_), intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer(psb_ipk_), intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_s_n_sparse_print
  end interface
  
  interface 
    subroutine psb_s_get_neigh(a,idx,neigh,n,info,lev)
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(in) :: a   
      integer(psb_ipk_), intent(in)                :: idx 
      integer(psb_ipk_), intent(out)               :: n   
      integer(psb_ipk_), allocatable, intent(out)  :: neigh(:)
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_ipk_), optional, intent(in)      :: lev 
    end subroutine psb_s_get_neigh
  end interface
  
  interface 
    subroutine psb_s_csall(nr,nc,a,info,nz) 
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in)             :: nr,nc
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: nz
    end subroutine psb_s_csall
  end interface
  
  interface 
    subroutine psb_s_reallocate_nz(nz,a) 
      import :: psb_ipk_, psb_sspmat_type
      integer(psb_ipk_), intent(in) :: nz
      class(psb_sspmat_type), intent(inout) :: a
    end subroutine psb_s_reallocate_nz
  end interface
  
  interface 
    subroutine psb_s_free(a) 
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
    end subroutine psb_s_free
  end interface
  
  interface 
    subroutine psb_s_trim(a) 
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
    end subroutine psb_s_trim
  end interface
  
  interface 
    subroutine psb_s_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(inout) :: a
      real(psb_spk_), intent(in)      :: val(:)
      integer(psb_ipk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: gtl(:)
    end subroutine psb_s_csput_a
  end interface

  
  interface 
    subroutine psb_s_csput_v(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      use psb_s_vect_mod, only : psb_s_vect_type
      use psb_i_vect_mod, only : psb_i_vect_type
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
      type(psb_s_vect_type), intent(inout)  :: val
      type(psb_i_vect_type), intent(inout)  :: ia, ja
      integer(psb_ipk_), intent(in)             :: nz, imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: gtl(:)
    end subroutine psb_s_csput_v
  end interface
 
  interface 
    subroutine psb_s_csgetptn(imin,imax,a,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_s_csgetptn
  end interface
  
  interface 
    subroutine psb_s_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      real(psb_spk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_s_csgetrow
  end interface
  
  interface 
    subroutine psb_s_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in)    :: a
      class(psb_sspmat_type), intent(inout) :: b
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_s_csgetblk
  end interface
  
  interface 
    subroutine psb_s_tril(a,l,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,u)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in)      :: a
      class(psb_sspmat_type), intent(inout)   :: l
      integer(psb_ipk_),intent(out)           :: info
      integer(psb_ipk_), intent(in), optional :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional           :: rscale,cscale
      class(psb_sspmat_type), optional, intent(inout)   :: u
    end subroutine psb_s_tril
  end interface
  
  interface 
    subroutine psb_s_triu(a,u,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,l)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in)      :: a
      class(psb_sspmat_type), intent(inout)   :: u
      integer(psb_ipk_),intent(out)           :: info
      integer(psb_ipk_), intent(in), optional :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional           :: rscale,cscale
      class(psb_sspmat_type), optional, intent(inout)   :: l
    end subroutine psb_s_triu
  end interface


  interface 
    subroutine psb_s_csclip(a,b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in) :: a
      class(psb_sspmat_type), intent(inout) :: b
      integer(psb_ipk_),intent(out)                  :: info
      integer(psb_ipk_), intent(in), optional        :: imin,imax,jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_s_csclip
  end interface
  
  interface 
    subroutine psb_s_b_csclip(a,b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_, psb_s_coo_sparse_mat
      class(psb_sspmat_type), intent(in) :: a
      type(psb_s_coo_sparse_mat), intent(out) :: b
      integer(psb_ipk_),intent(out)                  :: info
      integer(psb_ipk_), intent(in), optional        :: imin,imax,jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_s_b_csclip
  end interface
  
  interface 
    subroutine psb_s_mold(a,b)
      import :: psb_ipk_, psb_sspmat_type, psb_s_base_sparse_mat
      class(psb_sspmat_type), intent(inout)     :: a
      class(psb_s_base_sparse_mat), allocatable, intent(out) :: b
    end subroutine psb_s_mold
  end interface
  
  interface 
    subroutine psb_s_asb(a,mold) 
      import :: psb_ipk_, psb_sspmat_type, psb_s_base_sparse_mat
      class(psb_sspmat_type), intent(inout) :: a
      class(psb_s_base_sparse_mat), optional, intent(in) :: mold
    end subroutine psb_s_asb
  end interface
  
  interface 
    subroutine psb_s_transp_1mat(a)
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
    end subroutine psb_s_transp_1mat
  end interface
  
  interface 
    subroutine psb_s_transp_2mat(a,b)
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(in)  :: a
      class(psb_sspmat_type), intent(inout) :: b
    end subroutine psb_s_transp_2mat
  end interface
  
  interface 
    subroutine psb_s_transc_1mat(a)
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
    end subroutine psb_s_transc_1mat
  end interface
  
  interface 
    subroutine psb_s_transc_2mat(a,b)
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(in)  :: a
      class(psb_sspmat_type), intent(inout) :: b
    end subroutine psb_s_transc_2mat
  end interface
  
  interface 
    subroutine psb_s_reinit(a,clear)
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_s_reinit
    
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
    subroutine psb_s_cscnv(a,b,info,type,mold,upd,dupl)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_, psb_s_base_sparse_mat
      class(psb_sspmat_type), intent(in)      :: a
      class(psb_sspmat_type), intent(inout)   :: b
      integer(psb_ipk_), intent(out)                   :: info
      integer(psb_ipk_),optional, intent(in)           :: dupl, upd
      character(len=*), optional, intent(in) :: type
      class(psb_s_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_s_cscnv
  end interface
  

  interface 
    subroutine psb_s_cscnv_ip(a,iinfo,type,mold,dupl)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_, psb_s_base_sparse_mat
      class(psb_sspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(out)                   :: iinfo
      integer(psb_ipk_),optional, intent(in)           :: dupl
      character(len=*), optional, intent(in) :: type
      class(psb_s_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_s_cscnv_ip
  end interface
  

  interface 
    subroutine psb_s_cscnv_base(a,b,info,dupl)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_, psb_s_base_sparse_mat
      class(psb_sspmat_type), intent(in)       :: a
      class(psb_s_base_sparse_mat), intent(out) :: b
      integer(psb_ipk_), intent(out)                   :: info
      integer(psb_ipk_),optional, intent(in)           :: dupl
    end subroutine psb_s_cscnv_base
  end interface
  
  !
  ! Produce a version of the matrix with diagonal cut
  ! out; passes through a COO buffer. 
  !
  interface 
    subroutine psb_s_clip_d(a,b,info)
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(in)    :: a
      class(psb_sspmat_type), intent(inout) :: b
      integer(psb_ipk_),intent(out)                  :: info
    end subroutine psb_s_clip_d
  end interface
  
  interface 
    subroutine psb_s_clip_d_ip(a,info)
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
      integer(psb_ipk_),intent(out)                  :: info
    end subroutine psb_s_clip_d_ip
  end interface
  
  !
  ! These four interfaces cut through the
  ! encapsulation between spmat_type and base_sparse_mat.
  !
  interface 
    subroutine psb_s_mv_from(a,b)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_, psb_s_base_sparse_mat
      class(psb_sspmat_type), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(inout) :: b
    end subroutine psb_s_mv_from
  end interface
  
  interface 
    subroutine psb_s_cp_from(a,b)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_, psb_s_base_sparse_mat
      class(psb_sspmat_type), intent(out) :: a
      class(psb_s_base_sparse_mat), intent(in) :: b
    end subroutine psb_s_cp_from
  end interface
  
  interface 
    subroutine psb_s_mv_to(a,b)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_, psb_s_base_sparse_mat
      class(psb_sspmat_type), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(inout) :: b
    end subroutine psb_s_mv_to
  end interface
  
  interface 
    subroutine psb_s_cp_to(a,b)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_, psb_s_base_sparse_mat    
      class(psb_sspmat_type), intent(in) :: a
      class(psb_s_base_sparse_mat), intent(inout) :: b
    end subroutine psb_s_cp_to
  end interface
  
  !
  ! Transfer the internal allocation to the target.
  !  
  interface psb_move_alloc 
    subroutine psb_sspmat_type_move(a,b,info)
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
      class(psb_sspmat_type), intent(inout)   :: b
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine psb_sspmat_type_move
  end interface
  
  interface 
    subroutine psb_sspmat_clone(a,b,info)
      import :: psb_ipk_, psb_sspmat_type
      class(psb_sspmat_type), intent(inout) :: a
      class(psb_sspmat_type), intent(inout) :: b
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_sspmat_clone
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
    subroutine psb_s_csmm(alpha,a,x,beta,y,info,trans) 
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in) :: a
      real(psb_spk_), intent(in)    :: alpha, beta, x(:,:)
      real(psb_spk_), intent(inout) :: y(:,:)
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_s_csmm
    subroutine psb_s_csmv(alpha,a,x,beta,y,info,trans) 
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in) :: a
      real(psb_spk_), intent(in)    :: alpha, beta, x(:)
      real(psb_spk_), intent(inout) :: y(:)
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_s_csmv
    subroutine psb_s_csmv_vect(alpha,a,x,beta,y,info,trans) 
      use psb_s_vect_mod, only : psb_s_vect_type
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in)   :: a
      real(psb_spk_), intent(in)        :: alpha, beta
      type(psb_s_vect_type), intent(inout) :: x
      type(psb_s_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)                 :: info
      character, optional, intent(in)      :: trans
    end subroutine psb_s_csmv_vect
  end interface
  
  interface psb_cssm
    subroutine psb_s_cssm(alpha,a,x,beta,y,info,trans,scale,d) 
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in) :: a
      real(psb_spk_), intent(in)    :: alpha, beta, x(:,:)
      real(psb_spk_), intent(inout) :: y(:,:)
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in) :: trans, scale
      real(psb_spk_), intent(in), optional :: d(:)
    end subroutine psb_s_cssm
    subroutine psb_s_cssv(alpha,a,x,beta,y,info,trans,scale,d) 
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in) :: a
      real(psb_spk_), intent(in)    :: alpha, beta, x(:)
      real(psb_spk_), intent(inout) :: y(:)
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in) :: trans, scale
      real(psb_spk_), intent(in), optional :: d(:)
    end subroutine psb_s_cssv
    subroutine psb_s_cssv_vect(alpha,a,x,beta,y,info,trans,scale,d) 
      use psb_s_vect_mod, only : psb_s_vect_type
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in)   :: a
      real(psb_spk_), intent(in)        :: alpha, beta
      type(psb_s_vect_type), intent(inout) :: x
      type(psb_s_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)                 :: info
      character, optional, intent(in)      :: trans, scale
      type(psb_s_vect_type), optional, intent(inout)   :: d
    end subroutine psb_s_cssv_vect
  end interface
  
  interface 
    function psb_s_maxval(a) result(res)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_s_maxval
  end interface
  
  interface 
    function psb_s_csnmi(a) result(res)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_s_csnmi
  end interface
  
  interface 
    function psb_s_csnm1(a) result(res)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_s_csnm1
  end interface

  interface 
    function psb_s_rowsum(a,info) result(d)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in) :: a
      real(psb_spk_), allocatable      :: d(:)
      integer(psb_ipk_), intent(out)               :: info
    end function psb_s_rowsum
  end interface

  interface 
    function psb_s_arwsum(a,info) result(d)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in) :: a
      real(psb_spk_), allocatable        :: d(:)
      integer(psb_ipk_), intent(out)               :: info
    end function psb_s_arwsum
  end interface
  
  interface 
    function psb_s_colsum(a,info) result(d)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in) :: a
      real(psb_spk_), allocatable      :: d(:)
      integer(psb_ipk_), intent(out)               :: info
    end function psb_s_colsum
  end interface

  interface 
    function psb_s_aclsum(a,info)  result(d)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in) :: a
      real(psb_spk_), allocatable        :: d(:)
      integer(psb_ipk_), intent(out)        :: info
    end function psb_s_aclsum
  end interface

  interface 
    function psb_s_get_diag(a,info) result(d)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(in) :: a
      real(psb_spk_), allocatable         :: d(:)
      integer(psb_ipk_), intent(out)       :: info
    end function psb_s_get_diag
  end interface
  
  interface psb_scal
    subroutine psb_s_scal(d,a,info,side)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(inout) :: a
      real(psb_spk_), intent(in)             :: d(:)
      integer(psb_ipk_), intent(out)                    :: info
      character, intent(in), optional :: side
    end subroutine psb_s_scal
    subroutine psb_s_scals(d,a,info)
      import :: psb_ipk_, psb_sspmat_type, psb_spk_
      class(psb_sspmat_type), intent(inout) :: a
      real(psb_spk_), intent(in)             :: d
      integer(psb_ipk_), intent(out)                    :: info
    end subroutine psb_s_scals
  end interface


contains 


  
  subroutine  psb_s_set_mat_default(a) 
    implicit none 
    class(psb_s_base_sparse_mat), intent(in) :: a
    
    if (allocated(psb_s_base_mat_default)) then 
      deallocate(psb_s_base_mat_default)
    end if
    allocate(psb_s_base_mat_default, mold=a)

  end subroutine psb_s_set_mat_default
  
  function psb_s_get_mat_default(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    class(psb_s_base_sparse_mat), pointer :: res
    
    res => psb_s_get_base_mat_default()
    
  end function psb_s_get_mat_default

  
  function psb_s_get_base_mat_default() result(res)
    implicit none 
    class(psb_s_base_sparse_mat), pointer :: res
    
    if (.not.allocated(psb_s_base_mat_default)) then 
      allocate(psb_s_csr_sparse_mat :: psb_s_base_mat_default)
    end if

    res => psb_s_base_mat_default
    
  end function psb_s_get_base_mat_default




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

  
  function psb_s_sizeof(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    integer(psb_epk_) :: res
    
    res = 0
    if (allocated(a%a)) then 
      res = a%a%sizeof()
    end if
    
  end function psb_s_sizeof


  function psb_s_get_fmt(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    character(len=5) :: res

    if (allocated(a%a)) then 
      res = a%a%get_fmt()
    else
      res = 'NULL'
    end if

  end function psb_s_get_fmt


  function psb_s_get_dupl(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    integer(psb_ipk_) :: res

    if (allocated(a%a)) then 
      res = a%a%get_dupl()
    else
      res = psb_invalid_
    end if
  end function psb_s_get_dupl

  function psb_s_get_nrows(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    integer(psb_ipk_) :: res

    if (allocated(a%a)) then 
      res = a%a%get_nrows()
    else
      res = 0
    end if

  end function psb_s_get_nrows

  function psb_s_get_ncols(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    integer(psb_ipk_) :: res

    if (allocated(a%a)) then 
      res = a%a%get_ncols()
    else
      res = 0
    end if

  end function psb_s_get_ncols

  function psb_s_is_triangle(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_triangle()
    else
      res = .false.
    end if

  end function psb_s_is_triangle

  function psb_s_is_unit(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_unit()
    else
      res = .false.
    end if

  end function psb_s_is_unit

  function psb_s_is_upper(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upper()
    else
      res = .false.
    end if

  end function psb_s_is_upper

  function psb_s_is_lower(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = .not. a%a%is_upper()
    else
      res = .false.
    end if

  end function psb_s_is_lower

  function psb_s_is_null(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_null() 
    else
      res = .true.
    end if

  end function psb_s_is_null

  function psb_s_is_bld(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_bld()
    else
      res = .false.
    end if

  end function psb_s_is_bld

  function psb_s_is_upd(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upd()
    else
      res = .false.
    end if

  end function psb_s_is_upd

  function psb_s_is_asb(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_asb()
    else
      res = .false.
    end if

  end function psb_s_is_asb

  function psb_s_is_sorted(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_sorted()
    else
      res = .false.
    end if

  end function psb_s_is_sorted

  function psb_s_is_by_rows(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_by_rows()
    else
      res = .false.
    end if

  end function psb_s_is_by_rows

  function psb_s_is_by_cols(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_by_cols()
    else
      res = .false.
    end if

  end function psb_s_is_by_cols


  !
  subroutine s_mat_sync(a)
    implicit none 
    class(psb_sspmat_type), target, intent(in) :: a
    
    if (allocated(a%a))  call a%a%sync()

  end subroutine s_mat_sync

  !
  subroutine s_mat_set_host(a)
    implicit none 
    class(psb_sspmat_type), intent(inout) :: a

    if (allocated(a%a))  call a%a%set_host()
    
  end subroutine s_mat_set_host


  !
  subroutine s_mat_set_dev(a)
    implicit none 
    class(psb_sspmat_type), intent(inout) :: a

    if (allocated(a%a))  call a%a%set_dev()
    
  end subroutine s_mat_set_dev

  !
  subroutine s_mat_set_sync(a)
    implicit none 
    class(psb_sspmat_type), intent(inout) :: a

    if (allocated(a%a))  call a%a%set_sync()
    
  end subroutine s_mat_set_sync

  !
  function s_mat_is_dev(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    logical  :: res
  
    if (allocated(a%a)) then
      res = a%a%is_dev()
    else
      res = .false.
    end if
  end function s_mat_is_dev
  
  !
  function s_mat_is_host(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    logical  :: res

  
    if (allocated(a%a)) then
      res = a%a%is_host()
    else
      res = .true.
    end if
  end function s_mat_is_host

  !
  function s_mat_is_sync(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    logical  :: res

  
    if (allocated(a%a)) then
      res = a%a%is_sync()
    else
      res = .true.
    end if

  end function s_mat_is_sync


  function psb_s_is_repeatable_updates(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_repeatable_updates()
    else
      res = .false.
    end if

  end function psb_s_is_repeatable_updates

  subroutine psb_s_set_repeatable_updates(a,val) 
    implicit none 
    class(psb_sspmat_type), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (allocated(a%a)) then 
      call a%a%set_repeatable_updates(val)
    end if
    
  end subroutine psb_s_set_repeatable_updates


  function psb_s_get_nzeros(a) result(res)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    integer(psb_ipk_) :: res

    res = 0
    if (allocated(a%a)) then 
      res = a%a%get_nzeros()
    end if

  end function psb_s_get_nzeros

  function psb_s_get_size(a) result(res)

    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    integer(psb_ipk_) :: res


    res = 0
    if (allocated(a%a)) then 
      res = a%a%get_size()
    end if

  end function psb_s_get_size


  function psb_s_get_nz_row(idx,a) result(res)
    implicit none 
    integer(psb_ipk_), intent(in)               :: idx
    class(psb_sspmat_type), intent(in) :: a
    integer(psb_ipk_) :: res

    res = 0
    
    if (allocated(a%a)) res = a%a%get_nz_row(idx)

  end function psb_s_get_nz_row

  subroutine psb_s_clean_zeros(a,info)
    implicit none 
    integer(psb_ipk_), intent(out)        :: info
    class(psb_sspmat_type), intent(inout) :: a

    info = 0 
    if (allocated(a%a)) call a%a%clean_zeros(info)

  end subroutine psb_s_clean_zeros

#if defined(INT_I4_L8)
  subroutine psb_s_lcsgetptn(imin,imax,a,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
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
    
  end subroutine psb_s_lcsgetptn
  
  subroutine psb_s_lcsgetrow(imin,imax,a,nz,ia,ja,val,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
    implicit none 
    class(psb_sspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in)                  :: imin,imax
    integer(psb_ipk_), intent(out)                 :: nz
    integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
    real(psb_spk_), allocatable,  intent(inout)    :: val(:)
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
        
  end subroutine psb_s_lcsgetrow
#endif
end module psb_s_mat_mod
