
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
! package: psb_d_csr_mat_mod
!
! This module contains the definition of the psb_d_csr_sparse_mat type
! which implements an actual storage format (the CSR in this case) for
! a sparse matrix as well as the related methods (those who are
! specific to the type and could not be defined higher in the
! hierarchy). We are at the bottom level of the inheritance chain.
!
! Please refere to psb_d_base_mat_mod for a detailed description
! of the various methods, and to psb_d_csr_impl for implementation details.
!
module psb_d_csr_mat_mod

  use psb_d_base_mat_mod

  !> \namespace  psb_base_mod  \class  psb_d_csr_sparse_mat
  !! \extends psb_d_base_mat_mod::psb_d_base_sparse_mat
  !!
  !! psb_d_csr_sparse_mat type and the related methods.
  !! This is a very common storage type, and is the default for assembled
  !! matrices in our library
  type, extends(psb_d_base_sparse_mat) :: psb_d_csr_sparse_mat

    !> Pointers to beginning of rows in JA and VAL.
    integer(psb_ipk_), allocatable :: irp(:)
    !> Column indices.
    integer(psb_ipk_), allocatable :: ja(:)
    !> Coefficient values.
    real(psb_dpk_), allocatable :: val(:)

  contains
    procedure, pass(a) :: is_by_rows  => d_csr_is_by_rows
    procedure, pass(a) :: get_size    => d_csr_get_size
    procedure, pass(a) :: get_nzeros  => d_csr_get_nzeros
    procedure, nopass  :: get_fmt     => d_csr_get_fmt
    procedure, pass(a) :: sizeof      => d_csr_sizeof
    procedure, pass(a) :: csmm        => psb_d_csr_csmm
    procedure, pass(a) :: csmv        => psb_d_csr_csmv
    procedure, pass(a) :: inner_cssm  => psb_d_csr_cssm
    procedure, pass(a) :: inner_cssv  => psb_d_csr_cssv
    procedure, pass(a) :: scals       => psb_d_csr_scals
    procedure, pass(a) :: scalv       => psb_d_csr_scal
    procedure, pass(a) :: maxval      => psb_d_csr_maxval
    procedure, pass(a) :: spnmi       => psb_d_csr_csnmi
    procedure, pass(a) :: rowsum      => psb_d_csr_rowsum
    procedure, pass(a) :: arwsum      => psb_d_csr_arwsum
    procedure, pass(a) :: colsum      => psb_d_csr_colsum
    procedure, pass(a) :: aclsum      => psb_d_csr_aclsum
    procedure, pass(a) :: reallocate_nz => psb_d_csr_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_d_csr_allocate_mnnz
    procedure, pass(a) :: tril          => psb_d_csr_tril
    procedure, pass(a) :: triu          => psb_d_csr_triu
    procedure, pass(a) :: cp_to_coo   => psb_d_cp_csr_to_coo
    procedure, pass(a) :: cp_from_coo => psb_d_cp_csr_from_coo
    procedure, pass(a) :: cp_to_fmt   => psb_d_cp_csr_to_fmt
    procedure, pass(a) :: cp_from_fmt => psb_d_cp_csr_from_fmt
    procedure, pass(a) :: mv_to_coo   => psb_d_mv_csr_to_coo
    procedure, pass(a) :: mv_from_coo => psb_d_mv_csr_from_coo
    procedure, pass(a) :: mv_to_fmt   => psb_d_mv_csr_to_fmt
    procedure, pass(a) :: mv_from_fmt => psb_d_mv_csr_from_fmt
    procedure, pass(a) :: clean_zeros => psb_d_csr_clean_zeros
    procedure, pass(a) :: csput_a     => psb_d_csr_csput_a
    procedure, pass(a) :: get_diag    => psb_d_csr_get_diag
    procedure, pass(a) :: csgetptn    => psb_d_csr_csgetptn
    procedure, pass(a) :: csgetrow   => psb_d_csr_csgetrow
    procedure, pass(a) :: get_nz_row  => d_csr_get_nz_row
    procedure, pass(a) :: reinit      => psb_d_csr_reinit
    procedure, pass(a) :: trim        => psb_d_csr_trim
    procedure, pass(a) :: print       => psb_d_csr_print
    procedure, pass(a) :: free        => d_csr_free
    procedure, pass(a) :: mold        => psb_d_csr_mold

  end type psb_d_csr_sparse_mat

  private :: d_csr_get_nzeros, d_csr_free,  d_csr_get_fmt, &
       & d_csr_get_size, d_csr_sizeof, d_csr_get_nz_row, &
       & d_csr_is_by_rows

  !> \memberof psb_d_csr_sparse_mat
  !| \see psb_base_mat_mod::psb_base_reallocate_nz
  interface
    subroutine  psb_d_csr_reallocate_nz(nz,a)
      import
      integer(psb_ipk_), intent(in) :: nz
      class(psb_d_csr_sparse_mat), intent(inout) :: a
    end subroutine psb_d_csr_reallocate_nz
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !| \see psb_base_mat_mod::psb_base_reinit
  interface
    subroutine psb_d_csr_reinit(a,clear)
      import
      class(psb_d_csr_sparse_mat), intent(inout) :: a
      logical, intent(in), optional :: clear
    end subroutine psb_d_csr_reinit
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !| \see psb_base_mat_mod::psb_base_trim
  interface
    subroutine  psb_d_csr_trim(a)
      import
      class(psb_d_csr_sparse_mat), intent(inout) :: a
    end subroutine psb_d_csr_trim
  end interface


  !> \memberof psb_d_csr_sparse_mat
  !| \see psb_base_mat_mod::psb_base_mold
  interface
    subroutine psb_d_csr_mold(a,b,info)
      import
      class(psb_d_csr_sparse_mat), intent(in)                  :: a
      class(psb_d_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_d_csr_mold
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !| \see psb_base_mat_mod::psb_base_allocate_mnnz
  interface
    subroutine  psb_d_csr_allocate_mnnz(m,n,a,nz)
      import
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_d_csr_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in), optional :: nz
    end subroutine psb_d_csr_allocate_mnnz
  end interface


  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_print
  interface
    subroutine psb_d_csr_print(iout,a,iv,head,ivr,ivc)
      import
      integer(psb_ipk_), intent(in)               :: iout
      class(psb_d_csr_sparse_mat), intent(in) :: a
      integer(psb_lpk_), intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer(psb_lpk_), intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_d_csr_print
  end interface
  !
  !> Function  tril:
  !! \memberof  psb_d_base_sparse_mat
  !! \brief  Copy the lower triangle, i.e. all entries
  !!         A(I,J) such that J-I <= DIAG
  !!         default value is DIAG=0, i.e. lower triangle up to
  !!         the main diagonal.
  !!         DIAG=-1 means copy the strictly lower triangle
  !!         DIAG= 1 means copy the lower triangle plus the first diagonal
  !!                 of the upper triangle.
  !!         Moreover, apply a clipping by copying entries A(I,J) only if
  !!         IMIN<=I<=IMAX
  !!         JMIN<=J<=JMAX
  !!
  !!  \param l     the output (sub)matrix
  !!  \param info  return code
  !!  \param diag [0] the last diagonal (J-I) to be considered.
  !!  \param imin [1] the minimum row index we are interested in
  !!  \param imax [a\%get_nrows()] the minimum row index we are interested in
  !!  \param jmin [1] minimum col index
  !!  \param jmax [a\%get_ncols()] maximum col index
  !!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !!  \param u  [none]  copy of the complementary triangle
  !!
  !
  interface
    subroutine psb_d_csr_tril(a,l,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,u)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      class(psb_d_coo_sparse_mat), intent(out) :: l
      integer(psb_ipk_),intent(out)              :: info
      integer(psb_ipk_), intent(in), optional    :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional              :: rscale,cscale
      class(psb_d_coo_sparse_mat), optional, intent(out) :: u
    end subroutine psb_d_csr_tril
  end interface

  !
  !> Function  triu:
  !! \memberof  psb_d_csr_sparse_mat
  !! \brief  Copy the upper triangle, i.e. all entries
  !!         A(I,J) such that DIAG <= J-I
  !!         default value is DIAG=0, i.e. upper triangle from
  !!         the main diagonal up.
  !!         DIAG= 1 means copy the strictly upper triangle
  !!         DIAG=-1 means copy the upper triangle plus the first diagonal
  !!                 of the lower triangle.
  !!         Moreover, apply a clipping by copying entries A(I,J) only if
  !!         IMIN<=I<=IMAX
  !!         JMIN<=J<=JMAX
  !!         Optionally copies the lower triangle at the same time
  !!
  !!  \param u     the output (sub)matrix
  !!  \param info  return code
  !!  \param diag [0] the last diagonal (J-I) to be considered.
  !!  \param imin [1] the minimum row index we are interested in
  !!  \param imax [a\%get_nrows()] the minimum row index we are interested in
  !!  \param jmin [1] minimum col index
  !!  \param jmax [a\%get_ncols()] maximum col index
  !!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !!  \param l  [none]  copy of the complementary triangle
  !!
  !
  interface
    subroutine psb_d_csr_triu(a,u,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,l)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      class(psb_d_coo_sparse_mat), intent(out) :: u
      integer(psb_ipk_),intent(out)              :: info
      integer(psb_ipk_), intent(in), optional    :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional              :: rscale,cscale
      class(psb_d_coo_sparse_mat), optional, intent(out) :: l
    end subroutine psb_d_csr_triu
  end interface

  !
  !>
  !! \memberof  psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_clean_zeros
  !
  interface
    subroutine  psb_d_csr_clean_zeros(a, info)
      import
      class(psb_d_csr_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_d_csr_clean_zeros
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_cp_to_coo
  interface
    subroutine psb_d_cp_csr_to_coo(a,b,info)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      class(psb_d_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_d_cp_csr_to_coo
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_cp_from_coo
  interface
    subroutine psb_d_cp_csr_from_coo(a,b,info)
      import
      class(psb_d_csr_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_d_cp_csr_from_coo
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_cp_to_fmt
  interface
    subroutine psb_d_cp_csr_to_fmt(a,b,info)
      import
      class(psb_d_csr_sparse_mat), intent(in)   :: a
      class(psb_d_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                       :: info
    end subroutine psb_d_cp_csr_to_fmt
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_cp_from_fmt
  interface
    subroutine psb_d_cp_csr_from_fmt(a,b,info)
      import
      class(psb_d_csr_sparse_mat), intent(inout) :: a
      class(psb_d_base_sparse_mat), intent(in)   :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_d_cp_csr_from_fmt
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_mv_to_coo
  interface
    subroutine psb_d_mv_csr_to_coo(a,b,info)
      import
      class(psb_d_csr_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(inout)   :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_d_mv_csr_to_coo
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_mv_from_coo
  interface
    subroutine psb_d_mv_csr_from_coo(a,b,info)
      import
      class(psb_d_csr_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_d_mv_csr_from_coo
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_mv_to_fmt
  interface
    subroutine psb_d_mv_csr_to_fmt(a,b,info)
      import
      class(psb_d_csr_sparse_mat), intent(inout) :: a
      class(psb_d_base_sparse_mat), intent(inout)  :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_d_mv_csr_to_fmt
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_mv_from_fmt
  interface
    subroutine psb_d_mv_csr_from_fmt(a,b,info)
      import
      class(psb_d_csr_sparse_mat), intent(inout)  :: a
      class(psb_d_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                         :: info
    end subroutine psb_d_mv_csr_from_fmt
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_cp_from
  interface
    subroutine psb_d_csr_cp_from(a,b)
      import
      class(psb_d_csr_sparse_mat), intent(inout) :: a
      type(psb_d_csr_sparse_mat), intent(in)   :: b
    end subroutine psb_d_csr_cp_from
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_mv_from
  interface
    subroutine psb_d_csr_mv_from(a,b)
      import
      class(psb_d_csr_sparse_mat), intent(inout)  :: a
      type(psb_d_csr_sparse_mat), intent(inout) :: b
    end subroutine psb_d_csr_mv_from
  end interface


  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_csput_a
  interface
    subroutine psb_d_csr_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info)
      import
      class(psb_d_csr_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: val(:)
      integer(psb_ipk_), intent(in)             :: nz,ia(:), ja(:),&
           &  imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_d_csr_csput_a
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_base_mat_mod::psb_base_csgetptn
  interface
    subroutine psb_d_csr_csgetptn(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_d_csr_csgetptn
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_csgetrow
  interface
    subroutine psb_d_csr_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale,chksz)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale,chksz
    end subroutine psb_d_csr_csgetrow
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_cssv
  interface
    subroutine psb_d_csr_cssv(alpha,a,x,beta,y,info,trans)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      real(psb_dpk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_d_csr_cssv
  end interface
  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_cssm
  interface
    subroutine psb_d_csr_cssm(alpha,a,x,beta,y,info,trans)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_dpk_), intent(inout)       :: y(:,:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_d_csr_cssm
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_csmv
  interface
    subroutine psb_d_csr_csmv(alpha,a,x,beta,y,info,trans)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      real(psb_dpk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_d_csr_csmv
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_csmm
  interface
    subroutine psb_d_csr_csmm(alpha,a,x,beta,y,info,trans)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_dpk_), intent(inout)       :: y(:,:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_d_csr_csmm
  end interface


  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_maxval
  interface
    function psb_d_csr_maxval(a) result(res)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_d_csr_maxval
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_csnmi
  interface
    function psb_d_csr_csnmi(a) result(res)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_d_csr_csnmi
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_rowsum
  interface
    subroutine psb_d_csr_rowsum(d,a)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_d_csr_rowsum
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_arwsum
  interface
    subroutine psb_d_csr_arwsum(d,a)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_d_csr_arwsum
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_colsum
  interface
    subroutine psb_d_csr_colsum(d,a)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_d_csr_colsum
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_aclsum
  interface
    subroutine psb_d_csr_aclsum(d,a)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_d_csr_aclsum
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_get_diag
  interface
    subroutine psb_d_csr_get_diag(a,d,info)
      import
      class(psb_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)     :: d(:)
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_d_csr_get_diag
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_scal
  interface
    subroutine psb_d_csr_scal(d,a,info,side)
      import
      class(psb_d_csr_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: d(:)
      integer(psb_ipk_), intent(out)            :: info
      character, intent(in), optional :: side
    end subroutine psb_d_csr_scal
  end interface

  !> \memberof psb_d_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_scals
  interface
    subroutine psb_d_csr_scals(d,a,info)
      import
      class(psb_d_csr_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: d
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_d_csr_scals
  end interface

  
  type, extends(psb_d_csr_sparse_mat) :: psb_d_ecsr_sparse_mat
    
    !> Number of non-empty rows
    integer(psb_ipk_) :: nnerws
    !> Indices of non-empty rows
    integer(psb_ipk_), allocatable :: nerwp(:)

  contains
    procedure, nopass  :: get_fmt     => d_ecsr_get_fmt

    !    procedure, pass(a) :: csmm        => psb_d_ecsr_csmm
    procedure, pass(a) :: csmv        => psb_d_ecsr_csmv

    procedure, pass(a) :: cp_from_coo => psb_d_cp_ecsr_from_coo
    procedure, pass(a) :: cp_from_fmt => psb_d_cp_ecsr_from_fmt
    procedure, pass(a) :: mv_from_coo => psb_d_mv_ecsr_from_coo
    procedure, pass(a) :: mv_from_fmt => psb_d_mv_ecsr_from_fmt

    procedure, pass(a) :: cmp_nerwp   => psb_d_ecsr_cmp_nerwp
    procedure, pass(a) :: free        => d_ecsr_free
    procedure, pass(a) :: mold        => psb_d_ecsr_mold

  end type psb_d_ecsr_sparse_mat
  !> \memberof psb_d_ecsr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_csmv
  interface
    subroutine psb_d_ecsr_csmv(alpha,a,x,beta,y,info,trans)
      import
      class(psb_d_ecsr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      real(psb_dpk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_d_ecsr_csmv
  end interface

  !> \memberof psb_d_ecsr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_cp_from_coo
  interface
    subroutine psb_d_ecsr_cmp_nerwp(a,info)
      import
      class(psb_d_ecsr_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_d_ecsr_cmp_nerwp
  end interface

  !> \memberof psb_d_ecsr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_cp_from_coo
  interface
    subroutine psb_d_cp_ecsr_from_coo(a,b,info)
      import
      class(psb_d_ecsr_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_d_cp_ecsr_from_coo
  end interface

  !> \memberof psb_d_ecsr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_cp_from_fmt
  interface
    subroutine psb_d_cp_ecsr_from_fmt(a,b,info)
      import
      class(psb_d_ecsr_sparse_mat), intent(inout) :: a
      class(psb_d_base_sparse_mat), intent(in)   :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_d_cp_ecsr_from_fmt
  end interface

  !> \memberof psb_d_ecsr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_mv_from_coo
  interface
    subroutine psb_d_mv_ecsr_from_coo(a,b,info)
      import
      class(psb_d_ecsr_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_d_mv_ecsr_from_coo
  end interface

  !> \memberof psb_d_ecsr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_d_base_mv_from_fmt
  interface
    subroutine psb_d_mv_ecsr_from_fmt(a,b,info)
      import
      class(psb_d_ecsr_sparse_mat), intent(inout)  :: a
      class(psb_d_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                         :: info
    end subroutine psb_d_mv_ecsr_from_fmt
  end interface

  !> \memberof psb_d_ecsr_sparse_mat
  !| \see psb_base_mat_mod::psb_base_mold
  interface
    subroutine psb_d_ecsr_mold(a,b,info)
      import
      class(psb_d_ecsr_sparse_mat), intent(in)                  :: a
      class(psb_d_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_d_ecsr_mold
  end interface


  
  !> \namespace  psb_base_mod  \class  psb_ld_csr_sparse_mat
  !! \extends psb_ld_base_mat_mod::psb_ld_base_sparse_mat
  !!
  !! psb_ld_csr_sparse_mat type and the related methods.
  !! This is a very common storage type, and is the default for assembled
  !! matrices in our library
  type, extends(psb_ld_base_sparse_mat) :: psb_ld_csr_sparse_mat

    !> Pointers to beginning of rows in JA and VAL.
    integer(psb_lpk_), allocatable :: irp(:)
    !> Column indices.
    integer(psb_lpk_), allocatable :: ja(:)
    !> Coefficient values.
    real(psb_dpk_), allocatable :: val(:)

  contains
    procedure, pass(a) :: is_by_rows  => ld_csr_is_by_rows
    procedure, pass(a) :: get_size    => ld_csr_get_size
    procedure, pass(a) :: get_nzeros  => ld_csr_get_nzeros
    procedure, nopass  :: get_fmt     => ld_csr_get_fmt
    procedure, pass(a) :: sizeof      => ld_csr_sizeof
    procedure, pass(a) :: reallocate_nz => psb_ld_csr_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_ld_csr_allocate_mnnz
    procedure, pass(a) :: tril        => psb_ld_csr_tril
    procedure, pass(a) :: triu        => psb_ld_csr_triu
    procedure, pass(a) :: cp_to_coo   => psb_ld_cp_csr_to_coo
    procedure, pass(a) :: cp_from_coo => psb_ld_cp_csr_from_coo
    procedure, pass(a) :: cp_to_fmt   => psb_ld_cp_csr_to_fmt
    procedure, pass(a) :: cp_from_fmt => psb_ld_cp_csr_from_fmt
    procedure, pass(a) :: mv_to_coo   => psb_ld_mv_csr_to_coo
    procedure, pass(a) :: mv_from_coo => psb_ld_mv_csr_from_coo
    procedure, pass(a) :: mv_to_fmt   => psb_ld_mv_csr_to_fmt
    procedure, pass(a) :: mv_from_fmt => psb_ld_mv_csr_from_fmt
    procedure, pass(a) :: clean_zeros => psb_ld_csr_clean_zeros
    procedure, pass(a) :: csput_a     => psb_ld_csr_csput_a
    procedure, pass(a) :: get_diag    => psb_ld_csr_get_diag
    procedure, pass(a) :: csgetptn    => psb_ld_csr_csgetptn
    procedure, pass(a) :: csgetrow    => psb_ld_csr_csgetrow
    procedure, pass(a) :: get_nz_row  => ld_csr_get_nz_row
    procedure, pass(a) :: reinit      => psb_ld_csr_reinit
    procedure, pass(a) :: trim        => psb_ld_csr_trim
    procedure, pass(a) :: print       => psb_ld_csr_print
    procedure, pass(a) :: free        => ld_csr_free
    procedure, pass(a) :: mold        => psb_ld_csr_mold
    procedure, pass(a) :: scals       => psb_ld_csr_scals
    procedure, pass(a) :: scalv       => psb_ld_csr_scal
    procedure, pass(a) :: maxval      => psb_ld_csr_maxval
    procedure, pass(a) :: spnmi       => psb_ld_csr_csnmi
    procedure, pass(a) :: rowsum      => psb_ld_csr_rowsum
    procedure, pass(a) :: arwsum      => psb_ld_csr_arwsum
    procedure, pass(a) :: colsum      => psb_ld_csr_colsum
    procedure, pass(a) :: aclsum      => psb_ld_csr_aclsum

  end type psb_ld_csr_sparse_mat

  private :: ld_csr_get_nzeros, ld_csr_free,  ld_csr_get_fmt, &
       & ld_csr_get_size, ld_csr_sizeof, ld_csr_get_nz_row, &
       & ld_csr_is_by_rows

  !> \memberof psb_ld_csr_sparse_mat
  !| \see psb_base_mat_mod::psb_base_reallocate_nz
  interface
    subroutine  psb_ld_csr_reallocate_nz(nz,a)
      import
      integer(psb_lpk_), intent(in) :: nz
      class(psb_ld_csr_sparse_mat), intent(inout) :: a
    end subroutine psb_ld_csr_reallocate_nz
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !| \see psb_base_mat_mod::psb_base_reinit
  interface
    subroutine psb_ld_csr_reinit(a,clear)
      import
      class(psb_ld_csr_sparse_mat), intent(inout) :: a
      logical, intent(in), optional :: clear
    end subroutine psb_ld_csr_reinit
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !| \see psb_base_mat_mod::psb_base_trim
  interface
    subroutine  psb_ld_csr_trim(a)
      import
      class(psb_ld_csr_sparse_mat), intent(inout) :: a
    end subroutine psb_ld_csr_trim
  end interface


  !> \memberof psb_ld_csr_sparse_mat
  !| \see psb_base_mat_mod::psb_base_mold
  interface
    subroutine psb_ld_csr_mold(a,b,info)
      import
      class(psb_ld_csr_sparse_mat), intent(in)                  :: a
      class(psb_ld_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_ld_csr_mold
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !| \see psb_base_mat_mod::psb_base_allocate_mnnz
  interface
    subroutine  psb_ld_csr_allocate_mnnz(m,n,a,nz)
      import
      integer(psb_lpk_), intent(in) :: m,n
      class(psb_ld_csr_sparse_mat), intent(inout) :: a
      integer(psb_lpk_), intent(in), optional :: nz
    end subroutine psb_ld_csr_allocate_mnnz
  end interface


  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_print
  interface
    subroutine psb_ld_csr_print(iout,a,iv,head,ivr,ivc)
      import
      integer(psb_ipk_), intent(in)               :: iout
      class(psb_ld_csr_sparse_mat), intent(in) :: a
      integer(psb_lpk_), intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer(psb_lpk_), intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_ld_csr_print
  end interface
  !
  !> Function  tril:
  !! \memberof  psb_d_base_sparse_mat
  !! \brief  Copy the lower triangle, i.e. all entries
  !!         A(I,J) such that J-I <= DIAG
  !!         default value is DIAG=0, i.e. lower triangle up to
  !!         the main diagonal.
  !!         DIAG=-1 means copy the strictly lower triangle
  !!         DIAG= 1 means copy the lower triangle plus the first diagonal
  !!                 of the upper triangle.
  !!         Moreover, apply a clipping by copying entries A(I,J) only if
  !!         IMIN<=I<=IMAX
  !!         JMIN<=J<=JMAX
  !!
  !!  \param l     the output (sub)matrix
  !!  \param info  return code
  !!  \param diag [0] the last diagonal (J-I) to be considered.
  !!  \param imin [1] the minimum row index we are interested in
  !!  \param imax [a\%get_nrows()] the minimum row index we are interested in
  !!  \param jmin [1] minimum col index
  !!  \param jmax [a\%get_ncols()] maximum col index
  !!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !!  \param u  [none]  copy of the complementary triangle
  !!
  !
  interface
    subroutine psb_ld_csr_tril(a,l,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,u)
      import
      class(psb_ld_csr_sparse_mat), intent(in) :: a
      class(psb_ld_coo_sparse_mat), intent(out) :: l
      integer(psb_ipk_),intent(out)              :: info
      integer(psb_lpk_), intent(in), optional    :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional              :: rscale,cscale
      class(psb_ld_coo_sparse_mat), optional, intent(out) :: u
    end subroutine psb_ld_csr_tril
  end interface

  !
  !> Function  triu:
  !! \memberof  psb_d_csr_sparse_mat
  !! \brief  Copy the upper triangle, i.e. all entries
  !!         A(I,J) such that DIAG <= J-I
  !!         default value is DIAG=0, i.e. upper triangle from
  !!         the main diagonal up.
  !!         DIAG= 1 means copy the strictly upper triangle
  !!         DIAG=-1 means copy the upper triangle plus the first diagonal
  !!                 of the lower triangle.
  !!         Moreover, apply a clipping by copying entries A(I,J) only if
  !!         IMIN<=I<=IMAX
  !!         JMIN<=J<=JMAX
  !!         Optionally copies the lower triangle at the same time
  !!
  !!  \param u     the output (sub)matrix
  !!  \param info  return code
  !!  \param diag [0] the last diagonal (J-I) to be considered.
  !!  \param imin [1] the minimum row index we are interested in
  !!  \param imax [a\%get_nrows()] the minimum row index we are interested in
  !!  \param jmin [1] minimum col index
  !!  \param jmax [a\%get_ncols()] maximum col index
  !!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !!  \param l  [none]  copy of the complementary triangle
  !!
  !
  interface
    subroutine psb_ld_csr_triu(a,u,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,l)
      import
      class(psb_ld_csr_sparse_mat), intent(in) :: a
      class(psb_ld_coo_sparse_mat), intent(out) :: u
      integer(psb_ipk_),intent(out)              :: info
      integer(psb_lpk_), intent(in), optional    :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional              :: rscale,cscale
      class(psb_ld_coo_sparse_mat), optional, intent(out) :: l
    end subroutine psb_ld_csr_triu
  end interface

    !
  !>
  !! \memberof  psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_clean_zeros
  !
  interface
    subroutine  psb_ld_csr_clean_zeros(a, info)
      import
      class(psb_ld_csr_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_ld_csr_clean_zeros
  end interface



  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_cp_to_coo
  interface
    subroutine psb_ld_cp_csr_to_coo(a,b,info)
      import
      class(psb_ld_csr_sparse_mat), intent(in) :: a
      class(psb_ld_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_ld_cp_csr_to_coo
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_cp_from_coo
  interface
    subroutine psb_ld_cp_csr_from_coo(a,b,info)
      import
      class(psb_ld_csr_sparse_mat), intent(inout) :: a
      class(psb_ld_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_ld_cp_csr_from_coo
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_cp_to_fmt
  interface
    subroutine psb_ld_cp_csr_to_fmt(a,b,info)
      import
      class(psb_ld_csr_sparse_mat), intent(in)   :: a
      class(psb_ld_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                       :: info
    end subroutine psb_ld_cp_csr_to_fmt
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_cp_from_fmt
  interface
    subroutine psb_ld_cp_csr_from_fmt(a,b,info)
      import
      class(psb_ld_csr_sparse_mat), intent(inout) :: a
      class(psb_ld_base_sparse_mat), intent(in)   :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_ld_cp_csr_from_fmt
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_mv_to_coo
  interface
    subroutine psb_ld_mv_csr_to_coo(a,b,info)
      import
      class(psb_ld_csr_sparse_mat), intent(inout) :: a
      class(psb_ld_coo_sparse_mat), intent(inout)   :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_ld_mv_csr_to_coo
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_mv_from_coo
  interface
    subroutine psb_ld_mv_csr_from_coo(a,b,info)
      import
      class(psb_ld_csr_sparse_mat), intent(inout) :: a
      class(psb_ld_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_ld_mv_csr_from_coo
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_mv_to_fmt
  interface
    subroutine psb_ld_mv_csr_to_fmt(a,b,info)
      import
      class(psb_ld_csr_sparse_mat), intent(inout) :: a
      class(psb_ld_base_sparse_mat), intent(inout)  :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_ld_mv_csr_to_fmt
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_mv_from_fmt
  interface
    subroutine psb_ld_mv_csr_from_fmt(a,b,info)
      import
      class(psb_ld_csr_sparse_mat), intent(inout)  :: a
      class(psb_ld_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                         :: info
    end subroutine psb_ld_mv_csr_from_fmt
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_cp_from
  interface
    subroutine psb_ld_csr_cp_from(a,b)
      import
      class(psb_ld_csr_sparse_mat), intent(inout) :: a
      type(psb_ld_csr_sparse_mat), intent(in)   :: b
    end subroutine psb_ld_csr_cp_from
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_mv_from
  interface
    subroutine psb_ld_csr_mv_from(a,b)
      import
      class(psb_ld_csr_sparse_mat), intent(inout)  :: a
      type(psb_ld_csr_sparse_mat), intent(inout) :: b
    end subroutine psb_ld_csr_mv_from
  end interface


  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_csput_a
  interface
    subroutine psb_ld_csr_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info)
      import
      class(psb_ld_csr_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: val(:)
      integer(psb_lpk_), intent(in)             :: nz,ia(:), ja(:),&
           &  imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_ld_csr_csput_a
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_base_mat_mod::psb_base_csgetptn
  interface
    subroutine psb_ld_csr_csgetptn(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import
      class(psb_ld_csr_sparse_mat), intent(in) :: a
      integer(psb_lpk_), intent(in)                  :: imin,imax
      integer(psb_lpk_), intent(out)                 :: nz
      integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_lpk_), intent(in), optional        :: iren(:)
      integer(psb_lpk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_ld_csr_csgetptn
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_csgetrow
  interface
    subroutine psb_ld_csr_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import
      class(psb_ld_csr_sparse_mat), intent(in) :: a
      integer(psb_lpk_), intent(in)                  :: imin,imax
      integer(psb_lpk_), intent(out)                 :: nz
      integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
      real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_lpk_), intent(in), optional        :: iren(:)
      integer(psb_lpk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_ld_csr_csgetrow
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_get_diag
  interface
    subroutine psb_ld_csr_get_diag(a,d,info)
      import
      class(psb_ld_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)     :: d(:)
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_ld_csr_get_diag
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_scal
  interface
    subroutine psb_ld_csr_scal(d,a,info,side)
      import
      class(psb_ld_csr_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: d(:)
      integer(psb_ipk_), intent(out)            :: info
      character, intent(in), optional :: side
    end subroutine psb_ld_csr_scal
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_ld_base_mat_mod::psb_ld_base_scals
  interface
    subroutine psb_ld_csr_scals(d,a,info)
      import
      class(psb_ld_csr_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: d
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_ld_csr_scals
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_ld_base_maxval
  interface
    function psb_ld_csr_maxval(a) result(res)
      import
      class(psb_ld_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_ld_csr_maxval
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_ld_base_csnmi
  interface
    function psb_ld_csr_csnmi(a) result(res)
      import
      class(psb_ld_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_ld_csr_csnmi
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_ld_base_rowsum
  interface
    subroutine psb_ld_csr_rowsum(d,a)
      import
      class(psb_ld_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_ld_csr_rowsum
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_ld_base_arwsum
  interface
    subroutine psb_ld_csr_arwsum(d,a)
      import
      class(psb_ld_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_ld_csr_arwsum
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_ld_base_colsum
  interface
    subroutine psb_ld_csr_colsum(d,a)
      import
      class(psb_ld_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_ld_csr_colsum
  end interface

  !> \memberof psb_ld_csr_sparse_mat
  !! \see psb_d_base_mat_mod::psb_ld_base_aclsum
  interface
    subroutine psb_ld_csr_aclsum(d,a)
      import
      class(psb_ld_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_ld_csr_aclsum
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



  function d_csr_is_by_rows(a) result(res)
    implicit none
    class(psb_d_csr_sparse_mat), intent(in) :: a
    logical  :: res
    res = .true.

  end function d_csr_is_by_rows


  function d_csr_sizeof(a) result(res)
    implicit none
    class(psb_d_csr_sparse_mat), intent(in) :: a
    integer(psb_epk_) :: res
    res = 2 * psb_sizeof_ip
    res = res + psb_sizeof_dp  * psb_size(a%val)
    res = res + psb_sizeof_ip * psb_size(a%irp)
    res = res + psb_sizeof_ip * psb_size(a%ja)

  end function d_csr_sizeof

  function d_csr_get_fmt() result(res)
    implicit none
    character(len=5) :: res
    res = 'CSR'
  end function d_csr_get_fmt

  function d_csr_get_nzeros(a) result(res)
    implicit none
    class(psb_d_csr_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res
    res = a%irp(a%get_nrows()+1)-1
  end function d_csr_get_nzeros

  function d_csr_get_size(a) result(res)
    implicit none
    class(psb_d_csr_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res

    res = -1

    if (allocated(a%ja)) then
      res = size(a%ja)
    end if
    if (allocated(a%val)) then
      if (res >= 0) then
        res = min(res,size(a%val))
      else
        res = size(a%val)
      end if
    end if

  end function d_csr_get_size



  function  d_csr_get_nz_row(idx,a) result(res)

    implicit none

    class(psb_d_csr_sparse_mat), intent(in) :: a
    integer(psb_ipk_), intent(in)                  :: idx
    integer(psb_ipk_) :: res

    res = 0

    if ((1<=idx).and.(idx<=a%get_nrows())) then
      res = a%irp(idx+1)-a%irp(idx)
    end if

  end function d_csr_get_nz_row



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

  subroutine  d_csr_free(a)
    implicit none

    class(psb_d_csr_sparse_mat), intent(inout) :: a

    if (allocated(a%irp)) deallocate(a%irp)
    if (allocated(a%ja)) deallocate(a%ja)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(0_psb_ipk_)
    call a%set_ncols(0_psb_ipk_)

    return

  end subroutine d_csr_free



  function d_ecsr_get_fmt() result(res)
    implicit none
    character(len=5) :: res
    res = 'ECSR'
  end function d_ecsr_get_fmt

  subroutine  d_ecsr_free(a)
    implicit none

    class(psb_d_ecsr_sparse_mat), intent(inout) :: a

    
    if (allocated(a%nerwp)) deallocate(a%nerwp)
    a%nnerws = 0
    call a%psb_d_csr_sparse_mat%free()

    return
  end subroutine d_ecsr_free


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



  function ld_csr_is_by_rows(a) result(res)
    implicit none
    class(psb_ld_csr_sparse_mat), intent(in) :: a
    logical  :: res
    res = .true.

  end function ld_csr_is_by_rows


  function ld_csr_sizeof(a) result(res)
    implicit none
    class(psb_ld_csr_sparse_mat), intent(in) :: a
    integer(psb_epk_) :: res
    res = 2 * psb_sizeof_lp
    res = res + psb_sizeof_dp  * psb_size(a%val)
    res = res + psb_sizeof_lp * psb_size(a%irp)
    res = res + psb_sizeof_lp * psb_size(a%ja)

  end function ld_csr_sizeof

  function ld_csr_get_fmt() result(res)
    implicit none
    character(len=5) :: res
    res = 'CSR'
  end function ld_csr_get_fmt

  function ld_csr_get_nzeros(a) result(res)
    implicit none
    class(psb_ld_csr_sparse_mat), intent(in) :: a
    integer(psb_lpk_) :: res
    res = a%irp(a%get_nrows()+1)-1
  end function ld_csr_get_nzeros

  function ld_csr_get_size(a) result(res)
    implicit none
    class(psb_ld_csr_sparse_mat), intent(in) :: a
    integer(psb_lpk_) :: res

    res = -1

    if (allocated(a%ja)) then
      res = size(a%ja)
    end if
    if (allocated(a%val)) then
      if (res >= 0) then
        res = min(res,size(a%val))
      else
        res = size(a%val)
      end if
    end if

  end function ld_csr_get_size



  function  ld_csr_get_nz_row(idx,a) result(res)

    implicit none

    class(psb_ld_csr_sparse_mat), intent(in) :: a
    integer(psb_lpk_), intent(in)                  :: idx
    integer(psb_lpk_) :: res

    res = 0

    if ((1<=idx).and.(idx<=a%get_nrows())) then
      res = a%irp(idx+1)-a%irp(idx)
    end if

  end function ld_csr_get_nz_row



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

  subroutine  ld_csr_free(a)
    implicit none

    class(psb_ld_csr_sparse_mat), intent(inout) :: a

    if (allocated(a%irp)) deallocate(a%irp)
    if (allocated(a%ja)) deallocate(a%ja)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(0_psb_lpk_)
    call a%set_ncols(0_psb_lpk_)

    return

  end subroutine ld_csr_free

end module psb_d_csr_mat_mod
