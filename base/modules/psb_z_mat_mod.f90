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
! package: psb_z_mat_mod
!
! This module contains the definition of the psb_z_sparse type which
! is a generic container for a sparse matrix and it is mostly meant to
! provide a mean of switching, at run-time, among different formats,
! potentially unknown at the library compile-time by adding a layer of
! indirection. This type encapsulates the psb_z_base_sparse_mat class
! inside another class which is the one visible to the user. All the
! methods of the psb_z_mat_mod simply call the methods of the
! encapsulated class.


module psb_z_mat_mod

  use psb_z_base_mat_mod
  use psb_z_csr_mat_mod, only : psb_z_csr_sparse_mat
  use psb_z_csc_mat_mod, only : psb_z_csc_sparse_mat

  type :: psb_zspmat_type

    class(psb_z_base_sparse_mat), allocatable  :: a 

  contains
    ! Getters
    procedure, pass(a) :: get_nrows   => psb_z_get_nrows
    procedure, pass(a) :: get_ncols   => psb_z_get_ncols
    procedure, pass(a) :: get_nzeros  => psb_z_get_nzeros
    procedure, pass(a) :: get_nz_row  => psb_z_get_nz_row
    procedure, pass(a) :: get_size    => psb_z_get_size
    procedure, pass(a) :: get_state   => psb_z_get_state
    procedure, pass(a) :: get_dupl    => psb_z_get_dupl
    procedure, pass(a) :: is_null     => psb_z_is_null
    procedure, pass(a) :: is_bld      => psb_z_is_bld
    procedure, pass(a) :: is_upd      => psb_z_is_upd
    procedure, pass(a) :: is_asb      => psb_z_is_asb
    procedure, pass(a) :: is_sorted   => psb_z_is_sorted
    procedure, pass(a) :: is_upper    => psb_z_is_upper
    procedure, pass(a) :: is_lower    => psb_z_is_lower
    procedure, pass(a) :: is_triangle => psb_z_is_triangle
    procedure, pass(a) :: is_unit     => psb_z_is_unit
    procedure, pass(a) :: get_fmt     => psb_z_get_fmt
    procedure, pass(a) :: sizeof      => psb_z_sizeof

    ! Setters
    procedure, pass(a) :: set_nrows    => psb_z_set_nrows
    procedure, pass(a) :: set_ncols    => psb_z_set_ncols
    procedure, pass(a) :: set_dupl     => psb_z_set_dupl
    procedure, pass(a) :: set_state    => psb_z_set_state
    procedure, pass(a) :: set_null     => psb_z_set_null
    procedure, pass(a) :: set_bld      => psb_z_set_bld
    procedure, pass(a) :: set_upd      => psb_z_set_upd
    procedure, pass(a) :: set_asb      => psb_z_set_asb
    procedure, pass(a) :: set_sorted   => psb_z_set_sorted
    procedure, pass(a) :: set_upper    => psb_z_set_upper
    procedure, pass(a) :: set_lower    => psb_z_set_lower
    procedure, pass(a) :: set_triangle => psb_z_set_triangle
    procedure, pass(a) :: set_unit     => psb_z_set_unit

    ! Memory/data management 
    procedure, pass(a) :: csall         => psb_z_csall
    procedure, pass(a) :: free          => psb_z_free
    procedure, pass(a) :: trim          => psb_z_trim
    procedure, pass(a) :: csput         => psb_z_csput 
    procedure, pass(a) :: z_csgetptn    => psb_z_csgetptn
    procedure, pass(a) :: z_csgetrow    => psb_z_csgetrow
    procedure, pass(a) :: z_csgetblk    => psb_z_csgetblk
    generic, public    :: csget         => z_csgetptn, z_csgetrow, z_csgetblk 
    procedure, pass(a) :: z_csclip      => psb_z_csclip
    procedure, pass(a) :: z_b_csclip    => psb_z_b_csclip
    generic, public    :: csclip        => z_b_csclip, z_csclip
    procedure, pass(a) :: z_clip_d_ip   => psb_z_clip_d_ip
    procedure, pass(a) :: z_clip_d      => psb_z_clip_d
    generic, public    :: clip_diag     => z_clip_d_ip, z_clip_d
    procedure, pass(a) :: reall         => psb_z_reallocate_nz
    procedure, pass(a) :: get_neigh     => psb_z_get_neigh
    procedure, pass(a) :: z_cscnv       => psb_z_cscnv
    procedure, pass(a) :: z_cscnv_ip    => psb_z_cscnv_ip
    procedure, pass(a) :: z_cscnv_base  => psb_z_cscnv_base
    generic, public    :: cscnv         => z_cscnv, z_cscnv_ip, z_cscnv_base
    procedure, pass(a) :: reinit        => psb_z_reinit
    procedure, pass(a) :: print_i       => psb_z_sparse_print
    procedure, pass(a) :: print_n       => psb_z_n_sparse_print
    generic, public    :: print         => print_i, print_n
    procedure, pass(a) :: z_mv_from     => psb_z_mv_from
    generic, public    :: mv_from       => z_mv_from
    procedure, pass(a) :: z_mv_to       => psb_z_mv_to
    generic, public    :: mv_to         => z_mv_to
    procedure, pass(a) :: z_cp_from     => psb_z_cp_from
    generic, public    :: cp_from       => z_cp_from
    procedure, pass(a) :: z_cp_to       => psb_z_cp_to
    generic, public    :: cp_to         => z_cp_to
    procedure, pass(a) :: z_transp_1mat => psb_z_transp_1mat
    procedure, pass(a) :: z_transp_2mat => psb_z_transp_2mat
    generic, public    :: transp        => z_transp_1mat, z_transp_2mat
    procedure, pass(a) :: z_transc_1mat => psb_z_transc_1mat
    procedure, pass(a) :: z_transc_2mat => psb_z_transc_2mat
    generic, public    :: transc        => z_transc_1mat, z_transc_2mat

    
    
    ! Computational routines 
    procedure, pass(a) :: get_diag => psb_z_get_diag
    procedure, pass(a) :: csnmi    => psb_z_csnmi
    procedure, pass(a) :: z_csmv   => psb_z_csmv
    procedure, pass(a) :: z_csmm   => psb_z_csmm
    generic, public    :: csmm     => z_csmm, z_csmv
    procedure, pass(a) :: z_scals  => psb_z_scals
    procedure, pass(a) :: z_scal   => psb_z_scal
    generic, public    :: scal     => z_scals, z_scal 
    procedure, pass(a) :: z_cssv   => psb_z_cssv
    procedure, pass(a) :: z_cssm   => psb_z_cssm
    generic, public    :: cssm     => z_cssm, z_cssv

  end type psb_zspmat_type

  private :: psb_z_get_nrows, psb_z_get_ncols, psb_z_get_nzeros, psb_z_get_size, &
       & psb_z_get_state, psb_z_get_dupl, psb_z_is_null, psb_z_is_bld, psb_z_is_upd, &
       & psb_z_is_asb, psb_z_is_sorted, psb_z_is_upper, psb_z_is_lower, psb_z_is_triangle,&
       & psb_z_get_nz_row

  interface psb_sizeof
    module procedure psb_z_sizeof
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
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      integer, intent(in) :: m
    end subroutine psb_z_set_nrows
  end interface
  
  interface 
    subroutine psb_z_set_ncols(n,a) 
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      integer, intent(in) :: n
    end subroutine psb_z_set_ncols
  end interface
  
  interface 
    subroutine  psb_z_set_state(n,a) 
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      integer, intent(in) :: n
    end subroutine psb_z_set_state
  end interface
  
  interface 
    subroutine  psb_z_set_dupl(n,a) 
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      integer, intent(in) :: n
    end subroutine psb_z_set_dupl
  end interface
  
  interface 
    subroutine psb_z_set_null(a) 
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_set_null
  end interface
  
  interface 
    subroutine psb_z_set_bld(a) 
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_set_bld
  end interface
  
  interface 
    subroutine psb_z_set_upd(a) 
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_set_upd
  end interface
  
  interface 
    subroutine psb_z_set_asb(a) 
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_set_asb
  end interface
  
  interface 
    subroutine psb_z_set_sorted(a,val) 
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_z_set_sorted
  end interface
  
  interface 
    subroutine psb_z_set_triangle(a,val) 
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_z_set_triangle
  end interface
  
  interface 
    subroutine psb_z_set_unit(a,val) 
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_z_set_unit
  end interface
  
  interface 
    subroutine psb_z_set_lower(a,val) 
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_z_set_lower
  end interface
  
  interface 
    subroutine psb_z_set_upper(a,val) 
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_z_set_upper
  end interface
  
  interface 
    subroutine psb_z_sparse_print(iout,a,iv,eirs,eics,head,ivr,ivc)
      import :: psb_zspmat_type
      integer, intent(in)               :: iout
      class(psb_zspmat_type), intent(in) :: a   
      integer, intent(in), optional     :: iv(:)
      integer, intent(in), optional     :: eirs,eics
      character(len=*), optional        :: head
      integer, intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_z_sparse_print
  end interface

  interface 
    subroutine psb_z_n_sparse_print(fname,a,iv,eirs,eics,head,ivr,ivc)
      import :: psb_zspmat_type
      character(len=*), intent(in)      :: fname
      class(psb_zspmat_type), intent(in) :: a   
      integer, intent(in), optional     :: iv(:)
      integer, intent(in), optional     :: eirs,eics
      character(len=*), optional        :: head
      integer, intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_z_n_sparse_print
  end interface
  
  interface 
    subroutine psb_z_get_neigh(a,idx,neigh,n,info,lev)
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(in) :: a   
      integer, intent(in)                :: idx 
      integer, intent(out)               :: n   
      integer, allocatable, intent(out)  :: neigh(:)
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: lev 
    end subroutine psb_z_get_neigh
  end interface
  
  interface 
    subroutine psb_z_csall(nr,nc,a,info,nz) 
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(out) :: a
      integer, intent(in)             :: nr,nc
      integer, intent(out)            :: info
      integer, intent(in), optional   :: nz
    end subroutine psb_z_csall
  end interface
  
  interface 
    subroutine psb_z_reallocate_nz(nz,a) 
      import :: psb_zspmat_type
      integer, intent(in) :: nz
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_reallocate_nz
  end interface
  
  interface 
    subroutine psb_z_free(a) 
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_free
  end interface
  
  interface 
    subroutine psb_z_trim(a) 
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_trim
  end interface
  
  interface 
    subroutine psb_z_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      import :: psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: val(:)
      integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
      integer, intent(out)            :: info
      integer, intent(in), optional   :: gtl(:)
    end subroutine psb_z_csput
  end interface
  
  interface 
    subroutine psb_z_csgetptn(imin,imax,a,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
      import :: psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_z_csgetptn
  end interface
  
  interface 
    subroutine psb_z_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import :: psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_z_csgetrow
  end interface
  
  interface 
    subroutine psb_z_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
      import :: psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      class(psb_zspmat_type), intent(out) :: b
      integer, intent(in)                  :: imin,imax
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_z_csgetblk
  end interface
  
  interface 
    subroutine psb_z_csclip(a,b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
      import :: psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      class(psb_zspmat_type), intent(out) :: b
      integer,intent(out)                  :: info
      integer, intent(in), optional        :: imin,imax,jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_z_csclip
  end interface
  
  interface 
    subroutine psb_z_b_csclip(a,b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
      import :: psb_zspmat_type, psb_dpk_, psb_z_coo_sparse_mat
      class(psb_zspmat_type), intent(in) :: a
      type(psb_z_coo_sparse_mat), intent(out) :: b
      integer,intent(out)                  :: info
      integer, intent(in), optional        :: imin,imax,jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_z_b_csclip
  end interface
  
  interface 
    subroutine psb_z_cscnv(a,b,info,type,mold,upd,dupl)
      import :: psb_zspmat_type, psb_dpk_, psb_z_base_sparse_mat
      class(psb_zspmat_type), intent(in)    :: a
      class(psb_zspmat_type), intent(out)   :: b
      integer, intent(out)                   :: info
      integer,optional, intent(in)           :: dupl, upd
      character(len=*), optional, intent(in) :: type
      class(psb_z_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_z_cscnv
  end interface
  

  interface 
    subroutine psb_z_cscnv_ip(a,iinfo,type,mold,dupl)
      import :: psb_zspmat_type, psb_dpk_, psb_z_base_sparse_mat
      class(psb_zspmat_type), intent(inout) :: a
      integer, intent(out)                   :: iinfo
      integer,optional, intent(in)           :: dupl
      character(len=*), optional, intent(in) :: type
      class(psb_z_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_z_cscnv_ip
  end interface
  

  interface 
    subroutine psb_z_cscnv_base(a,b,info,dupl)
      import :: psb_zspmat_type, psb_dpk_, psb_z_base_sparse_mat
      class(psb_zspmat_type), intent(in)       :: a
      class(psb_z_base_sparse_mat), intent(out) :: b
      integer, intent(out)                   :: info
      integer,optional, intent(in)           :: dupl
    end subroutine psb_z_cscnv_base
  end interface
  
  interface 
    subroutine psb_z_clip_d(a,b,info)
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(in) :: a
      class(psb_zspmat_type), intent(out) :: b
      integer,intent(out)                  :: info
    end subroutine psb_z_clip_d
  end interface
  
  interface 
    subroutine psb_z_clip_d_ip(a,info)
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      integer,intent(out)                  :: info
    end subroutine psb_z_clip_d_ip
  end interface
  
  interface 
    subroutine psb_z_mv_from(a,b)
      import :: psb_zspmat_type, psb_dpk_, psb_z_base_sparse_mat
      class(psb_zspmat_type), intent(out) :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
    end subroutine psb_z_mv_from
  end interface
  
  interface 
    subroutine psb_z_cp_from(a,b)
      import :: psb_zspmat_type, psb_dpk_, psb_z_base_sparse_mat
      class(psb_zspmat_type), intent(out) :: a
      class(psb_z_base_sparse_mat), intent(inout), allocatable :: b
    end subroutine psb_z_cp_from
  end interface
  
  interface 
    subroutine psb_z_mv_to(a,b)
      import :: psb_zspmat_type, psb_dpk_, psb_z_base_sparse_mat
      class(psb_zspmat_type), intent(inout) :: a
      class(psb_z_base_sparse_mat), intent(out) :: b
    end subroutine psb_z_mv_to
  end interface
  
  interface 
    subroutine psb_z_cp_to(a,b)
      import :: psb_zspmat_type, psb_dpk_, psb_z_base_sparse_mat    
      class(psb_zspmat_type), intent(in) :: a
      class(psb_z_base_sparse_mat), intent(out) :: b
    end subroutine psb_z_cp_to
  end interface
  
  interface psb_move_alloc 
    subroutine psb_zspmat_type_move(a,b,info)
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
      class(psb_zspmat_type), intent(out)   :: b
      integer, intent(out)                   :: info
    end subroutine psb_zspmat_type_move
  end interface
  

  interface psb_clone
    subroutine psb_zspmat_type_clone(a,b,info)
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(in)  :: a
      class(psb_zspmat_type), intent(out) :: b
      integer, intent(out)                 :: info
    end subroutine psb_zspmat_type_clone
  end interface
  
  interface 
    subroutine psb_z_transp_1mat(a)
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_transp_1mat
  end interface
  
  interface 
    subroutine psb_z_transp_2mat(a,b)
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(out) :: a
      class(psb_zspmat_type), intent(in)  :: b
    end subroutine psb_z_transp_2mat
  end interface
  
  interface 
    subroutine psb_z_transc_1mat(a)
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a
    end subroutine psb_z_transc_1mat
  end interface
  
  interface 
    subroutine psb_z_transc_2mat(a,b)
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(out) :: a
      class(psb_zspmat_type), intent(in)  :: b
    end subroutine psb_z_transc_2mat
  end interface
  
  interface 
    subroutine psb_z_reinit(a,clear)
      import :: psb_zspmat_type
      class(psb_zspmat_type), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_z_reinit
    
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
      import :: psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      complex(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
      complex(psb_dpk_), intent(inout) :: y(:,:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_z_csmm
    subroutine psb_z_csmv(alpha,a,x,beta,y,info,trans) 
      import :: psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      complex(psb_dpk_), intent(in)    :: alpha, beta, x(:)
      complex(psb_dpk_), intent(inout) :: y(:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_z_csmv
  end interface
  
  interface psb_cssm
    subroutine psb_z_cssm(alpha,a,x,beta,y,info,trans,scale,d) 
      import :: psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      complex(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
      complex(psb_dpk_), intent(inout) :: y(:,:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans, scale
      complex(psb_dpk_), intent(in), optional :: d(:)
    end subroutine psb_z_cssm
    subroutine psb_z_cssv(alpha,a,x,beta,y,info,trans,scale,d) 
      import :: psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      complex(psb_dpk_), intent(in)    :: alpha, beta, x(:)
      complex(psb_dpk_), intent(inout) :: y(:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans, scale
      complex(psb_dpk_), intent(in), optional :: d(:)
    end subroutine psb_z_cssv
  end interface
  
  interface 
    function psb_z_csnmi(a) result(res)
      import :: psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_z_csnmi
  end interface
  
  interface 
    subroutine psb_z_get_diag(a,d,info)
      import :: psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(in) :: a
      complex(psb_dpk_), intent(out)          :: d(:)
      integer, intent(out)                 :: info
    end subroutine psb_z_get_diag
  end interface
  
  interface psb_scal
    subroutine psb_z_scal(d,a,info)
      import :: psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(inout) :: a
      complex(psb_dpk_), intent(in)              :: d(:)
      integer, intent(out)                    :: info
    end subroutine psb_z_scal
    subroutine psb_z_scals(d,a,info)
      import :: psb_zspmat_type, psb_dpk_
      class(psb_zspmat_type), intent(inout) :: a
      complex(psb_dpk_), intent(in)              :: d
      integer, intent(out)                    :: info
    end subroutine psb_z_scals
  end interface




contains 


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
    integer(psb_long_int_k_) :: res
    
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
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_dupl()
    else
      res = psb_invalid_
    end if
  end function psb_z_get_dupl


  function psb_z_get_state(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_state()
    else
      res = psb_spmat_null_
    end if
  end function psb_z_get_state

  function psb_z_get_nrows(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_nrows()
    else
      res = 0
    end if

  end function psb_z_get_nrows

  function psb_z_get_ncols(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    integer :: res

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



  function psb_z_get_nzeros(a) result(res)
    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    integer :: res

    res = 0
    if (allocated(a%a)) then 
      res = a%a%get_nzeros()
    end if

  end function psb_z_get_nzeros

  function psb_z_get_size(a) result(res)

    implicit none 
    class(psb_zspmat_type), intent(in) :: a
    integer :: res


    res = 0
    if (allocated(a%a)) then 
      res = a%a%get_size()
    end if

  end function psb_z_get_size


  function psb_z_get_nz_row(idx,a) result(res)
    implicit none 
    integer, intent(in)               :: idx
    class(psb_zspmat_type), intent(in) :: a
    integer :: res

    res = 0
    
    if (allocated(a%a)) res = a%a%get_nz_row(idx)

  end function psb_z_get_nz_row


end module psb_z_mat_mod
