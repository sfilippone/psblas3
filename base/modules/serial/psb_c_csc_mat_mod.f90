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
! package: psb_c_csc_mat_mod
!
! This module contains the definition of the psb_c_csc_sparse_mat type
! which implements an actual storage format (the CSC in this case) for
! a sparse matrix as well as the related methods (those who are
! specific to the type and could not be defined higher in the
! hierarchy). We are at the bottom level of the inheritance chain.
!
! Please refere to psb_c_base_mat_mod for a detailed description
! of the various methods, and to psb_c_csc_impl for implementation details.
!
module psb_c_csc_mat_mod

  use psb_c_base_mat_mod

  !> \namespace  psb_base_mod  \class  psb_c_csc_sparse_mat
  !! \extends psb_c_base_mat_mod::psb_c_base_sparse_mat
  !!
  !! psb_c_csc_sparse_mat type and the related methods.
  !!
  type, extends(psb_c_base_sparse_mat) :: psb_c_csc_sparse_mat

    !> Pointers to beginning of cols in IA and VAL.
    integer(psb_ipk_), allocatable :: icp(:)
    !> Row indices.
    integer(psb_ipk_), allocatable :: ia(:)
    !> Coefficient values.
    complex(psb_spk_), allocatable :: val(:)

  contains
    procedure, pass(a) :: is_by_cols  => c_csc_is_by_cols
    procedure, pass(a) :: get_size    => c_csc_get_size
    procedure, pass(a) :: get_nzeros  => c_csc_get_nzeros
    procedure, nopass  :: get_fmt     => c_csc_get_fmt
    procedure, pass(a) :: sizeof      => c_csc_sizeof
    procedure, pass(a) :: csmm        => psb_c_csc_csmm
    procedure, pass(a) :: csmv        => psb_c_csc_csmv
    procedure, pass(a) :: inner_cssm  => psb_c_csc_cssm
    procedure, pass(a) :: inner_cssv  => psb_c_csc_cssv
    procedure, pass(a) :: scals       => psb_c_csc_scals
    procedure, pass(a) :: scalv       => psb_c_csc_scal
    procedure, pass(a) :: maxval      => psb_c_csc_maxval
    procedure, pass(a) :: spnm1       => psb_c_csc_csnm1
    procedure, pass(a) :: rowsum      => psb_c_csc_rowsum
    procedure, pass(a) :: arwsum      => psb_c_csc_arwsum
    procedure, pass(a) :: colsum      => psb_c_csc_colsum
    procedure, pass(a) :: aclsum      => psb_c_csc_aclsum
    procedure, pass(a) :: reallocate_nz => psb_c_csc_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_c_csc_allocate_mnnz
    procedure, pass(a) :: cp_to_coo   => psb_c_cp_csc_to_coo
    procedure, pass(a) :: cp_from_coo => psb_c_cp_csc_from_coo
    procedure, pass(a) :: cp_to_fmt   => psb_c_cp_csc_to_fmt
    procedure, pass(a) :: cp_from_fmt => psb_c_cp_csc_from_fmt
    procedure, pass(a) :: mv_to_coo   => psb_c_mv_csc_to_coo
    procedure, pass(a) :: mv_from_coo => psb_c_mv_csc_from_coo
    procedure, pass(a) :: mv_to_fmt   => psb_c_mv_csc_to_fmt
    procedure, pass(a) :: mv_from_fmt => psb_c_mv_csc_from_fmt
    procedure, pass(a) :: clean_zeros => psb_c_csc_clean_zeros
    procedure, pass(a) :: csput_a      => psb_c_csc_csput_a
    procedure, pass(a) :: get_diag    => psb_c_csc_get_diag
    procedure, pass(a) :: csgetptn    => psb_c_csc_csgetptn
    procedure, pass(a) :: csgetrow   => psb_c_csc_csgetrow
    procedure, pass(a) :: get_nz_col  => c_csc_get_nz_col
    procedure, pass(a) :: reinit      => psb_c_csc_reinit
    procedure, pass(a) :: trim        => psb_c_csc_trim
    procedure, pass(a) :: print       => psb_c_csc_print
    procedure, pass(a) :: free        => c_csc_free
    procedure, pass(a) :: mold        => psb_c_csc_mold

  end type psb_c_csc_sparse_mat

  private :: c_csc_get_nzeros, c_csc_free,  c_csc_get_fmt, &
       & c_csc_get_size, c_csc_sizeof, c_csc_get_nz_col


  !> \namespace  psb_base_mod  \class  psb_c_csc_sparse_mat
  !! \extends psb_c_base_mat_mod::psb_c_base_sparse_mat
  !!
  !! psb_c_csc_sparse_mat type and the related methods.
  !!
  type, extends(psb_lc_base_sparse_mat) :: psb_lc_csc_sparse_mat

    !> Pointers to beginning of cols in IA and VAL.
    integer(psb_lpk_), allocatable :: icp(:)
    !> Row indices.
    integer(psb_lpk_), allocatable :: ia(:)
    !> Coefficient values.
    complex(psb_spk_), allocatable :: val(:)

  contains
    procedure, pass(a) :: is_by_cols  => lc_csc_is_by_cols
    procedure, pass(a) :: get_size    => lc_csc_get_size
    procedure, pass(a) :: get_nzeros  => lc_csc_get_nzeros
    procedure, nopass  :: get_fmt     => lc_csc_get_fmt
    procedure, pass(a) :: sizeof      => lc_csc_sizeof
    procedure, pass(a) :: scals       => psb_lc_csc_scals
    procedure, pass(a) :: scalv       => psb_lc_csc_scal
    procedure, pass(a) :: maxval      => psb_lc_csc_maxval
    procedure, pass(a) :: spnm1       => psb_lc_csc_csnm1
    procedure, pass(a) :: rowsum      => psb_lc_csc_rowsum
    procedure, pass(a) :: arwsum      => psb_lc_csc_arwsum
    procedure, pass(a) :: colsum      => psb_lc_csc_colsum
    procedure, pass(a) :: aclsum      => psb_lc_csc_aclsum
    procedure, pass(a) :: reallocate_nz => psb_lc_csc_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_lc_csc_allocate_mnnz
    procedure, pass(a) :: cp_to_coo   => psb_lc_cp_csc_to_coo
    procedure, pass(a) :: cp_from_coo => psb_lc_cp_csc_from_coo
    procedure, pass(a) :: cp_to_fmt   => psb_lc_cp_csc_to_fmt
    procedure, pass(a) :: cp_from_fmt => psb_lc_cp_csc_from_fmt
    procedure, pass(a) :: mv_to_coo   => psb_lc_mv_csc_to_coo
    procedure, pass(a) :: mv_from_coo => psb_lc_mv_csc_from_coo
    procedure, pass(a) :: mv_to_fmt   => psb_lc_mv_csc_to_fmt
    procedure, pass(a) :: mv_from_fmt => psb_lc_mv_csc_from_fmt
    procedure, pass(a) :: clean_zeros => psb_lc_csc_clean_zeros
    procedure, pass(a) :: csput_a     => psb_lc_csc_csput_a
    procedure, pass(a) :: get_diag    => psb_lc_csc_get_diag
    procedure, pass(a) :: csgetptn    => psb_lc_csc_csgetptn
    procedure, pass(a) :: csgetrow    => psb_lc_csc_csgetrow
    procedure, pass(a) :: get_nz_col  => lc_csc_get_nz_col
    procedure, pass(a) :: reinit      => psb_lc_csc_reinit
    procedure, pass(a) :: trim        => psb_lc_csc_trim
    procedure, pass(a) :: print       => psb_lc_csc_print
    procedure, pass(a) :: free        => lc_csc_free
    procedure, pass(a) :: mold        => psb_lc_csc_mold

  end type psb_lc_csc_sparse_mat

  private :: lc_csc_get_nzeros, lc_csc_free,  lc_csc_get_fmt, &
       & lc_csc_get_size, lc_csc_sizeof, lc_csc_get_nz_col

  !> \memberof psb_c_csc_sparse_mat
  !| \see psb_base_mat_mod::psb_base_reallocate_nz
  interface
    subroutine  psb_c_csc_reallocate_nz(nz,a)
      import
      integer(psb_ipk_), intent(in) :: nz
      class(psb_c_csc_sparse_mat), intent(inout) :: a
    end subroutine psb_c_csc_reallocate_nz
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !| \see psb_base_mat_mod::psb_base_reinit
  interface
    subroutine psb_c_csc_reinit(a,clear)
      import
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      logical, intent(in), optional :: clear
    end subroutine psb_c_csc_reinit
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !| \see psb_base_mat_mod::psb_base_trim
  interface
    subroutine  psb_c_csc_trim(a)
      import
      class(psb_c_csc_sparse_mat), intent(inout) :: a
    end subroutine psb_c_csc_trim
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !| \see psb_base_mat_mod::psb_base_mold
  interface
    subroutine psb_c_csc_mold(a,b,info)
      import
      class(psb_c_csc_sparse_mat), intent(in)                  :: a
      class(psb_c_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_c_csc_mold
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !| \see psb_base_mat_mod::psb_base_allocate_mnnz
  interface
    subroutine  psb_c_csc_allocate_mnnz(m,n,a,nz)
      import
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in), optional :: nz
    end subroutine psb_c_csc_allocate_mnnz
  end interface


  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_print
  interface
    subroutine psb_c_csc_print(iout,a,iv,head,ivr,ivc)
      import
      integer(psb_ipk_), intent(in)               :: iout
      class(psb_c_csc_sparse_mat), intent(in) :: a
      integer(psb_lpk_), intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer(psb_lpk_), intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_c_csc_print
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_cp_to_coo
  interface
    subroutine psb_c_cp_csc_to_coo(a,b,info)
      import
      class(psb_c_csc_sparse_mat), intent(in) :: a
      class(psb_c_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_c_cp_csc_to_coo
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_cp_from_coo
  interface
    subroutine psb_c_cp_csc_from_coo(a,b,info)
      import
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_c_cp_csc_from_coo
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_cp_to_fmt
  interface
    subroutine psb_c_cp_csc_to_fmt(a,b,info)
      import
      class(psb_c_csc_sparse_mat), intent(in)   :: a
      class(psb_c_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                       :: info
    end subroutine psb_c_cp_csc_to_fmt
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_cp_from_fmt
  interface
    subroutine psb_c_cp_csc_from_fmt(a,b,info)
      import
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      class(psb_c_base_sparse_mat), intent(in)   :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_c_cp_csc_from_fmt
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_mv_to_coo
  interface
    subroutine psb_c_mv_csc_to_coo(a,b,info)
      import
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(inout)   :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_c_mv_csc_to_coo
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_mv_from_coo
  interface
    subroutine psb_c_mv_csc_from_coo(a,b,info)
      import
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_c_mv_csc_from_coo
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_mv_to_fmt
  interface
    subroutine psb_c_mv_csc_to_fmt(a,b,info)
      import
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      class(psb_c_base_sparse_mat), intent(inout)  :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_c_mv_csc_to_fmt
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_mv_from_fmt
  interface
    subroutine psb_c_mv_csc_from_fmt(a,b,info)
      import
      class(psb_c_csc_sparse_mat), intent(inout)  :: a
      class(psb_c_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                         :: info
    end subroutine psb_c_mv_csc_from_fmt
  end interface

  !
  !>
  !! \memberof  psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_clean_zeros
  !
  interface
    subroutine  psb_c_csc_clean_zeros(a, info)
      import
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_c_csc_clean_zeros
  end interface


  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_cp_from
  interface
    subroutine psb_c_csc_cp_from(a,b)
      import
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      type(psb_c_csc_sparse_mat), intent(in)   :: b
    end subroutine psb_c_csc_cp_from
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_mv_from
  interface
    subroutine psb_c_csc_mv_from(a,b)
      import
      class(psb_c_csc_sparse_mat), intent(inout)  :: a
      type(psb_c_csc_sparse_mat), intent(inout) :: b
    end subroutine psb_c_csc_mv_from
  end interface


  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_csput_a
  interface
    subroutine psb_c_csc_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info)
      import
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)      :: val(:)
      integer(psb_ipk_), intent(in)             :: nz,ia(:), ja(:),&
           &  imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_c_csc_csput_a
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_base_mat_mod::psb_base_csgetptn
  interface
    subroutine psb_c_csc_csgetptn(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import
      class(psb_c_csc_sparse_mat), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_c_csc_csgetptn
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_csgetrow
  interface
    subroutine psb_c_csc_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale,chksz)
      import :: psb_ipk_, psb_c_csc_sparse_mat, psb_spk_
      class(psb_c_csc_sparse_mat), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale,chksz
    end subroutine psb_c_csc_csgetrow
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_csgetblk
  interface
    subroutine psb_c_csc_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale,chksz)
      import
      class(psb_c_csc_sparse_mat), intent(in) :: a
      class(psb_c_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale,chksz
    end subroutine psb_c_csc_csgetblk
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_cssv
  interface
    subroutine psb_c_csc_cssv(alpha,a,x,beta,y,info,trans)
      import
      class(psb_c_csc_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:)
      complex(psb_spk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_c_csc_cssv
  end interface
  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_cssm
  interface
    subroutine psb_c_csc_cssm(alpha,a,x,beta,y,info,trans)
      import
      class(psb_c_csc_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
      complex(psb_spk_), intent(inout)       :: y(:,:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_c_csc_cssm
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_csmv
  interface
    subroutine psb_c_csc_csmv(alpha,a,x,beta,y,info,trans,ivshft)
      import
      class(psb_c_csc_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:)
      complex(psb_spk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
      integer(psb_ipk_), optional, intent(in) :: ivshft
    end subroutine psb_c_csc_csmv
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_csmm
  interface
    subroutine psb_c_csc_csmm(alpha,a,x,beta,y,info,trans,ivshft)
      import
      class(psb_c_csc_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
      complex(psb_spk_), intent(inout)       :: y(:,:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
      integer(psb_ipk_), optional, intent(in) :: ivshft
    end subroutine psb_c_csc_csmm
  end interface


  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_maxval
  interface
    function psb_c_csc_maxval(a) result(res)
      import
      class(psb_c_csc_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_c_csc_maxval
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_csnm1
  interface
    function psb_c_csc_csnm1(a) result(res)
      import
      class(psb_c_csc_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_c_csc_csnm1
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_rowsum
  interface
    subroutine psb_c_csc_rowsum(d,a)
      import
      class(psb_c_csc_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_c_csc_rowsum
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_arwsum
  interface
    subroutine psb_c_csc_arwsum(d,a)
      import
      class(psb_c_csc_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_c_csc_arwsum
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_colsum
  interface
    subroutine psb_c_csc_colsum(d,a)
      import
      class(psb_c_csc_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_c_csc_colsum
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_aclsum
  interface
    subroutine psb_c_csc_aclsum(d,a)
      import
      class(psb_c_csc_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_c_csc_aclsum
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_get_diag
  interface
    subroutine psb_c_csc_get_diag(a,d,info)
      import
      class(psb_c_csc_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(out)     :: d(:)
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_c_csc_get_diag
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_scal
  interface
    subroutine psb_c_csc_scal(d,a,info,side)
      import
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)      :: d(:)
      integer(psb_ipk_), intent(out)            :: info
      character, intent(in), optional :: side
    end subroutine psb_c_csc_scal
  end interface

  !> \memberof psb_c_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_c_base_scals
  interface
    subroutine psb_c_csc_scals(d,a,info)
      import
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)      :: d
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_c_csc_scals
  end interface

  !
  ! lc
  !
  !> \memberof psb_lc_csc_sparse_mat
  !| \see psb_base_mat_mod::psb_base_reallocate_nz
  interface
    subroutine  psb_lc_csc_reallocate_nz(nz,a)
      import
      integer(psb_lpk_), intent(in) :: nz
      class(psb_lc_csc_sparse_mat), intent(inout) :: a
    end subroutine psb_lc_csc_reallocate_nz
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !| \see psb_base_mat_mod::psb_base_reinit
  interface
    subroutine psb_lc_csc_reinit(a,clear)
      import
      class(psb_lc_csc_sparse_mat), intent(inout) :: a
      logical, intent(in), optional :: clear
    end subroutine psb_lc_csc_reinit
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !| \see psb_base_mat_mod::psb_base_trim
  interface
    subroutine  psb_lc_csc_trim(a)
      import
      class(psb_lc_csc_sparse_mat), intent(inout) :: a
    end subroutine psb_lc_csc_trim
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !| \see psb_base_mat_mod::psb_base_mold
  interface
    subroutine psb_lc_csc_mold(a,b,info)
      import
      class(psb_lc_csc_sparse_mat), intent(in)                  :: a
      class(psb_lc_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_lc_csc_mold
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !| \see psb_base_mat_mod::psb_base_allocate_mnnz
  interface
    subroutine  psb_lc_csc_allocate_mnnz(m,n,a,nz)
      import
      integer(psb_lpk_), intent(in) :: m,n
      class(psb_lc_csc_sparse_mat), intent(inout) :: a
      integer(psb_lpk_), intent(in), optional :: nz
    end subroutine psb_lc_csc_allocate_mnnz
  end interface


  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_print
  interface
    subroutine psb_lc_csc_print(iout,a,iv,head,ivr,ivc)
      import
      integer(psb_ipk_), intent(in)               :: iout
      class(psb_lc_csc_sparse_mat), intent(in) :: a
      integer(psb_lpk_), intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer(psb_lpk_), intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_lc_csc_print
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_cp_to_coo
  interface
    subroutine psb_lc_cp_csc_to_coo(a,b,info)
      import
      class(psb_lc_csc_sparse_mat), intent(in) :: a
      class(psb_lc_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lc_cp_csc_to_coo
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_cp_from_coo
  interface
    subroutine psb_lc_cp_csc_from_coo(a,b,info)
      import
      class(psb_lc_csc_sparse_mat), intent(inout) :: a
      class(psb_lc_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_lc_cp_csc_from_coo
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_cp_to_fmt
  interface
    subroutine psb_lc_cp_csc_to_fmt(a,b,info)
      import
      class(psb_lc_csc_sparse_mat), intent(in)   :: a
      class(psb_lc_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                       :: info
    end subroutine psb_lc_cp_csc_to_fmt
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_cp_from_fmt
  interface
    subroutine psb_lc_cp_csc_from_fmt(a,b,info)
      import
      class(psb_lc_csc_sparse_mat), intent(inout) :: a
      class(psb_lc_base_sparse_mat), intent(in)   :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_lc_cp_csc_from_fmt
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_mv_to_coo
  interface
    subroutine psb_lc_mv_csc_to_coo(a,b,info)
      import
      class(psb_lc_csc_sparse_mat), intent(inout) :: a
      class(psb_lc_coo_sparse_mat), intent(inout)   :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lc_mv_csc_to_coo
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_mv_from_coo
  interface
    subroutine psb_lc_mv_csc_from_coo(a,b,info)
      import
      class(psb_lc_csc_sparse_mat), intent(inout) :: a
      class(psb_lc_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_lc_mv_csc_from_coo
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_mv_to_fmt
  interface
    subroutine psb_lc_mv_csc_to_fmt(a,b,info)
      import
      class(psb_lc_csc_sparse_mat), intent(inout) :: a
      class(psb_lc_base_sparse_mat), intent(inout)  :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_lc_mv_csc_to_fmt
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_mv_from_fmt
  interface
    subroutine psb_lc_mv_csc_from_fmt(a,b,info)
      import
      class(psb_lc_csc_sparse_mat), intent(inout)  :: a
      class(psb_lc_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                         :: info
    end subroutine psb_lc_mv_csc_from_fmt
  end interface

  !
  !>
  !! \memberof  psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_clean_zeros
  !
  interface
    subroutine  psb_lc_csc_clean_zeros(a, info)
      import
      class(psb_lc_csc_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_lc_csc_clean_zeros
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_cp_from
  interface
    subroutine psb_lc_csc_cp_from(a,b)
      import
      class(psb_lc_csc_sparse_mat), intent(inout) :: a
      type(psb_lc_csc_sparse_mat), intent(in)   :: b
    end subroutine psb_lc_csc_cp_from
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_mv_from
  interface
    subroutine psb_lc_csc_mv_from(a,b)
      import
      class(psb_lc_csc_sparse_mat), intent(inout)  :: a
      type(psb_lc_csc_sparse_mat), intent(inout) :: b
    end subroutine psb_lc_csc_mv_from
  end interface


  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_csput_a
  interface
    subroutine psb_lc_csc_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info)
      import
      class(psb_lc_csc_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)      :: val(:)
      integer(psb_lpk_), intent(in)             :: nz,ia(:), ja(:),&
           &  imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lc_csc_csput_a
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_base_mat_mod::psb_base_csgetptn
  interface
    subroutine psb_lc_csc_csgetptn(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import
      class(psb_lc_csc_sparse_mat), intent(in) :: a
      integer(psb_lpk_), intent(in)                  :: imin,imax
      integer(psb_lpk_), intent(out)                 :: nz
      integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_lpk_), intent(in), optional        :: iren(:)
      integer(psb_lpk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_lc_csc_csgetptn
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_csgetrow
  interface
    subroutine psb_lc_csc_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import
      class(psb_lc_csc_sparse_mat), intent(in) :: a
      integer(psb_lpk_), intent(in)                  :: imin,imax
      integer(psb_lpk_), intent(out)                 :: nz
      integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_lpk_), intent(in), optional        :: iren(:)
      integer(psb_lpk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_lc_csc_csgetrow
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_csgetblk
  interface
    subroutine psb_lc_csc_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
      import
      class(psb_lc_csc_sparse_mat), intent(in) :: a
      class(psb_lc_coo_sparse_mat), intent(inout) :: b
      integer(psb_lpk_), intent(in)                  :: imin,imax
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_lpk_), intent(in), optional        :: iren(:)
      integer(psb_lpk_), intent(in), optional        :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_lc_csc_csgetblk
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_get_diag
  interface
    subroutine psb_lc_csc_get_diag(a,d,info)
      import
      class(psb_lc_csc_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(out)     :: d(:)
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lc_csc_get_diag
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_lc_base_maxval
  interface
    function psb_lc_csc_maxval(a) result(res)
      import
      class(psb_lc_csc_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_lc_csc_maxval
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_lc_base_csnm1
  interface
    function psb_lc_csc_csnm1(a) result(res)
      import
      class(psb_lc_csc_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_lc_csc_csnm1
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_lc_base_rowsum
  interface
    subroutine psb_lc_csc_rowsum(d,a)
      import
      class(psb_lc_csc_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_lc_csc_rowsum
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_lc_base_arwsum
  interface
    subroutine psb_lc_csc_arwsum(d,a)
      import
      class(psb_lc_csc_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_lc_csc_arwsum
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_lc_base_colsum
  interface
    subroutine psb_lc_csc_colsum(d,a)
      import
      class(psb_lc_csc_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_lc_csc_colsum
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_c_base_mat_mod::psb_lc_base_aclsum
  interface
    subroutine psb_lc_csc_aclsum(d,a)
      import
      class(psb_lc_csc_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_lc_csc_aclsum
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_scal
  interface
    subroutine psb_lc_csc_scal(d,a,info,side)
      import
      class(psb_lc_csc_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)      :: d(:)
      integer(psb_ipk_), intent(out)            :: info
      character, intent(in), optional :: side
    end subroutine psb_lc_csc_scal
  end interface

  !> \memberof psb_lc_csc_sparse_mat
  !! \see psb_lc_base_mat_mod::psb_lc_base_scals
  interface
    subroutine psb_lc_csc_scals(d,a,info)
      import
      class(psb_lc_csc_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)      :: d
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lc_csc_scals
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


  function c_csc_is_by_cols(a) result(res)
    implicit none
    class(psb_c_csc_sparse_mat), intent(in) :: a
    logical  :: res
    res = .true.

  end function c_csc_is_by_cols


  function c_csc_sizeof(a) result(res)
    implicit none
    class(psb_c_csc_sparse_mat), intent(in) :: a
    integer(psb_epk_) :: res
    res = 2 * psb_sizeof_ip
    res = res + (2*psb_sizeof_sp)  * psb_size(a%val)
    res = res + psb_sizeof_ip * psb_size(a%icp)
    res = res + psb_sizeof_ip * psb_size(a%ia)

  end function c_csc_sizeof

  function c_csc_get_fmt() result(res)
    implicit none
    character(len=5) :: res
    res = 'CSC'
  end function c_csc_get_fmt

  function c_csc_get_nzeros(a) result(res)
    implicit none
    class(psb_c_csc_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res
    res = a%icp(a%get_ncols()+1)-1
  end function c_csc_get_nzeros

  function c_csc_get_size(a) result(res)
    implicit none
    class(psb_c_csc_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res

    res = -1

    if (allocated(a%ia)) then
      res = size(a%ia)
    end if
    if (allocated(a%val)) then
      if (res >= 0) then
        res = min(res,size(a%val))
      else
        res = size(a%val)
      end if
    end if

  end function c_csc_get_size



  function  c_csc_get_nz_col(idx,a) result(res)
    use psb_const_mod
    implicit none

    class(psb_c_csc_sparse_mat), intent(in) :: a
    integer(psb_ipk_), intent(in)                  :: idx
    integer(psb_ipk_) :: res

    res = 0

    if ((1<=idx).and.(idx<=a%get_ncols())) then
      res = a%icp(idx+1)-a%icp(idx)
    end if

  end function c_csc_get_nz_col



  ! == ===================================
  !
  !
  !
  ! Data management
  !
  !
  !
  !
  !
  ! == ===================================


  subroutine  c_csc_free(a)
    implicit none

    class(psb_c_csc_sparse_mat), intent(inout) :: a

    if (allocated(a%icp)) deallocate(a%icp)
    if (allocated(a%ia)) deallocate(a%ia)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(0_psb_ipk_)
    call a%set_ncols(0_psb_ipk_)

    return

  end subroutine c_csc_free



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


  function lc_csc_is_by_cols(a) result(res)
    implicit none
    class(psb_lc_csc_sparse_mat), intent(in) :: a
    logical  :: res
    res = .true.

  end function lc_csc_is_by_cols

  !
  ! lc
  !

  function lc_csc_sizeof(a) result(res)
    implicit none
    class(psb_lc_csc_sparse_mat), intent(in) :: a
    integer(psb_epk_) :: res
    res = 2*psb_sizeof_lp
    res = res + (2*psb_sizeof_sp)  * psb_size(a%val)
    res = res + psb_sizeof_lp * psb_size(a%icp)
    res = res + psb_sizeof_lp * psb_size(a%ia)

  end function lc_csc_sizeof

  function lc_csc_get_fmt() result(res)
    implicit none
    character(len=5) :: res
    res = 'CSC'
  end function lc_csc_get_fmt

  function lc_csc_get_nzeros(a) result(res)
    implicit none
    class(psb_lc_csc_sparse_mat), intent(in) :: a
    integer(psb_lpk_) :: res
    res = a%icp(a%get_ncols()+1)-1
  end function lc_csc_get_nzeros

  function lc_csc_get_size(a) result(res)
    implicit none
    class(psb_lc_csc_sparse_mat), intent(in) :: a
    integer(psb_lpk_) :: res

    res = -1

    if (allocated(a%ia)) then
      res = size(a%ia)
    end if
    if (allocated(a%val)) then
      if (res >= 0) then
        res = min(res,size(a%val))
      else
        res = size(a%val)
      end if
    end if

  end function lc_csc_get_size



  function  lc_csc_get_nz_col(idx,a) result(res)
    use psb_const_mod
    implicit none

    class(psb_lc_csc_sparse_mat), intent(in) :: a
    integer(psb_lpk_), intent(in)                  :: idx
    integer(psb_lpk_) :: res

    res = 0

    if ((1<=idx).and.(idx<=a%get_ncols())) then
      res = a%icp(idx+1)-a%icp(idx)
    end if

  end function lc_csc_get_nz_col



  ! == ===================================
  !
  !
  !
  ! Data management
  !
  !
  !
  !
  !
  ! == ===================================


  subroutine  lc_csc_free(a)
    implicit none

    class(psb_lc_csc_sparse_mat), intent(inout) :: a

    if (allocated(a%icp)) deallocate(a%icp)
    if (allocated(a%ia)) deallocate(a%ia)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(0_psb_lpk_)
    call a%set_ncols(0_psb_lpk_)

    return

  end subroutine lc_csc_free



end module psb_c_csc_mat_mod
