
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
! package: psb_z_csrli_mat_mod
!
! This module contains the definition of the psb_z_csr_sparse_mat type
! which implements an actual storage format (the CSR in this case) for
! a sparse matrix as well as the related methods (those who are
! specific to the type and could not be defined higher in the
! hierarchy). We are at the bottom level of the inheritance chain.
!
! Please refere to psb_z_base_mat_mod for a detailed description
! of the various methods, and to psb_z_csr_impl for implementation details.
!
module psb_z_csrli_mat_mod

  use psb_z_csr_mat_mod

  type, extends(psb_z_csr_sparse_mat) :: psb_z_csrli_sparse_mat

    complex(psb_dpk_) :: lambda=zzero
    
  contains
    procedure, nopass  :: get_fmt     => z_csrli_get_fmt
    procedure, pass(a) :: csmm        => psb_z_csrli_csmm
    procedure, pass(a) :: csmv        => psb_z_csrli_csmv
    procedure, pass(a) :: inner_cssm  => psb_z_csrli_cssm
    procedure, pass(a) :: inner_cssv  => psb_z_csrli_cssv
    procedure, pass(a) :: scals       => psb_z_csrli_scals
    procedure, pass(a) :: scalv       => psb_z_csrli_scal
    procedure, pass(a) :: maxval      => psb_z_csrli_maxval
    procedure, pass(a) :: spnmi       => psb_z_csrli_csnmi
    procedure, pass(a) :: rowsum      => psb_z_csrli_rowsum
    procedure, pass(a) :: arwsum      => psb_z_csrli_arwsum
    procedure, pass(a) :: colsum      => psb_z_csrli_colsum
    procedure, pass(a) :: aclsum      => psb_z_csrli_aclsum
    procedure, pass(a) :: tril          => psb_z_csrli_tril
    procedure, pass(a) :: triu          => psb_z_csrli_triu
    procedure, pass(a) :: cp_to_coo   => psb_z_cp_csrli_to_coo
    procedure, pass(a) :: cp_from_coo => psb_z_cp_csrli_from_coo
    !procedure, pass(a) :: cp_to_fmt   => psb_z_cp_csrli_to_fmt
    !procedure, pass(a) :: cp_from_fmt => psb_z_cp_csrli_from_fmt
    procedure, pass(a) :: mv_to_coo   => psb_z_mv_csrli_to_coo
    procedure, pass(a) :: mv_from_coo => psb_z_mv_csrli_from_coo
    !procedure, pass(a) :: mv_to_fmt   => psb_z_mv_csrli_to_fmt
    !procedure, pass(a) :: mv_from_fmt => psb_z_mv_csrli_from_fmt
    procedure, pass(a) :: get_diag    => psb_z_csrli_get_diag
    procedure, pass(a) :: csgetrow   => psb_z_csrli_csgetrow
    procedure, pass(a) :: reinit      => psb_z_csrli_reinit
    procedure, pass(a) :: print       => psb_z_csrli_print
    procedure, pass(a) :: free        => z_csrli_free
    procedure, pass(a) :: mold        => psb_z_csrli_mold

    procedure, pass(a) :: set_lambda  => z_csrli_set_lambda
    procedure, pass(a) :: get_lambda  => z_csrli_get_lambda
    
  end type psb_z_csrli_sparse_mat

  private :: z_csrli_get_nzeros, z_csrli_free,  z_csrli_get_fmt, &
       & z_csrli_get_size, z_csrli_sizeof, z_csrli_get_nz_row, &
       & z_csrli_is_by_rows

  !> \memberof psb_z_csrli_sparse_mat
  !| \see psb_base_mat_mod::psb_base_reallocate_nz
  interface
    subroutine  psb_z_csrli_reallocate_nz(nz,a)
      import
      integer(psb_ipk_), intent(in) :: nz
      class(psb_z_csrli_sparse_mat), intent(inout) :: a
    end subroutine psb_z_csrli_reallocate_nz
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !| \see psb_base_mat_mod::psb_base_reinit
  interface
    subroutine psb_z_csrli_reinit(a,clear)
      import
      class(psb_z_csrli_sparse_mat), intent(inout) :: a
      logical, intent(in), optional :: clear
    end subroutine psb_z_csrli_reinit
  end interface

!!$  !> \memberof psb_z_csrli_sparse_mat
!!$  !| \see psb_base_mat_mod::psb_base_trim
!!$  interface
!!$    subroutine  psb_z_csrli_trim(a)
!!$      import
!!$      class(psb_z_csrli_sparse_mat), intent(inout) :: a
!!$    end subroutine psb_z_csrli_trim
!!$  end interface


  !> \memberof psb_z_csrli_sparse_mat
  !| \see psb_base_mat_mod::psb_base_mold
  interface
    subroutine psb_z_csrli_mold(a,b,info)
      import
      class(psb_z_csrli_sparse_mat), intent(in)                  :: a
      class(psb_z_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_z_csrli_mold
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !| \see psb_base_mat_mod::psb_base_allocate_mnnz
  interface
    subroutine  psb_z_csrli_allocate_mnnz(m,n,a,nz)
      import
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_z_csrli_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in), optional :: nz
    end subroutine psb_z_csrli_allocate_mnnz
  end interface


  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_print
  interface
    subroutine psb_z_csrli_print(iout,a,iv,head,ivr,ivc)
      import
      integer(psb_ipk_), intent(in)               :: iout
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      integer(psb_lpk_), intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer(psb_lpk_), intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_z_csrli_print
  end interface
  !
  !> Function  tril:
  !! \memberof  psb_z_base_sparse_mat
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
    subroutine psb_z_csrli_tril(a,l,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,u)
      import
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      class(psb_z_coo_sparse_mat), intent(out) :: l
      integer(psb_ipk_),intent(out)              :: info
      integer(psb_ipk_), intent(in), optional    :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional              :: rscale,cscale
      class(psb_z_coo_sparse_mat), optional, intent(out) :: u
    end subroutine psb_z_csrli_tril
  end interface

  !
  !> Function  triu:
  !! \memberof  psb_z_csrli_sparse_mat
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
    subroutine psb_z_csrli_triu(a,u,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,l)
      import
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      class(psb_z_coo_sparse_mat), intent(out) :: u
      integer(psb_ipk_),intent(out)              :: info
      integer(psb_ipk_), intent(in), optional    :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional              :: rscale,cscale
      class(psb_z_coo_sparse_mat), optional, intent(out) :: l
    end subroutine psb_z_csrli_triu
  end interface

  !
  !>
  !! \memberof  psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_clean_zeros
  !
!!$  interface
!!$    subroutine  psb_z_csrli_clean_zeros(a, info)
!!$      import
!!$      class(psb_z_csrli_sparse_mat), intent(inout) :: a
!!$      integer(psb_ipk_), intent(out)              :: info
!!$    end subroutine psb_z_csrli_clean_zeros
!!$  end interface
!!$
  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_cp_to_coo
  interface
    subroutine psb_z_cp_csrli_to_coo(a,b,info)
      import
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      class(psb_z_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_cp_csrli_to_coo
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_cp_from_coo
  interface
    subroutine psb_z_cp_csrli_from_coo(a,b,info)
      import
      class(psb_z_csrli_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_z_cp_csrli_from_coo
  end interface

!!$  !> \memberof psb_z_csrli_sparse_mat
!!$  !! \see psb_z_base_mat_mod::psb_z_base_cp_to_fmt
!!$  interface
!!$    subroutine psb_z_cp_csrli_to_fmt(a,b,info)
!!$      import
!!$      class(psb_z_csrli_sparse_mat), intent(in)   :: a
!!$      class(psb_z_base_sparse_mat), intent(inout) :: b
!!$      integer(psb_ipk_), intent(out)                       :: info
!!$    end subroutine psb_z_cp_csrli_to_fmt
!!$  end interface
!!$
!!$  !> \memberof psb_z_csrli_sparse_mat
!!$  !! \see psb_z_base_mat_mod::psb_z_base_cp_from_fmt
!!$  interface
!!$    subroutine psb_z_cp_csrli_from_fmt(a,b,info)
!!$      import
!!$      class(psb_z_csrli_sparse_mat), intent(inout) :: a
!!$      class(psb_z_base_sparse_mat), intent(in)   :: b
!!$      integer(psb_ipk_), intent(out)                        :: info
!!$    end subroutine psb_z_cp_csrli_from_fmt
!!$  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_mv_to_coo
  interface
    subroutine psb_z_mv_csrli_to_coo(a,b,info)
      import
      class(psb_z_csrli_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(inout)   :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_mv_csrli_to_coo
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_mv_from_coo
  interface
    subroutine psb_z_mv_csrli_from_coo(a,b,info)
      import
      class(psb_z_csrli_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_z_mv_csrli_from_coo
  end interface

!!$  !> \memberof psb_z_csrli_sparse_mat
!!$  !! \see psb_z_base_mat_mod::psb_z_base_mv_to_fmt
!!$  interface
!!$    subroutine psb_z_mv_csrli_to_fmt(a,b,info)
!!$      import
!!$      class(psb_z_csrli_sparse_mat), intent(inout) :: a
!!$      class(psb_z_base_sparse_mat), intent(inout)  :: b
!!$      integer(psb_ipk_), intent(out)                        :: info
!!$    end subroutine psb_z_mv_csrli_to_fmt
!!$  end interface
!!$
!!$  !> \memberof psb_z_csrli_sparse_mat
!!$  !! \see psb_z_base_mat_mod::psb_z_base_mv_from_fmt
!!$  interface
!!$    subroutine psb_z_mv_csrli_from_fmt(a,b,info)
!!$      import
!!$      class(psb_z_csrli_sparse_mat), intent(inout)  :: a
!!$      class(psb_z_base_sparse_mat), intent(inout) :: b
!!$      integer(psb_ipk_), intent(out)                         :: info
!!$    end subroutine psb_z_mv_csrli_from_fmt
!!$  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_cp_from
  interface
    subroutine psb_z_csrli_cp_from(a,b)
      import
      class(psb_z_csrli_sparse_mat), intent(inout) :: a
      type(psb_z_csrli_sparse_mat), intent(in)   :: b
    end subroutine psb_z_csrli_cp_from
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_mv_from
  interface
    subroutine psb_z_csrli_mv_from(a,b)
      import
      class(psb_z_csrli_sparse_mat), intent(inout)  :: a
       type(psb_z_csrli_sparse_mat), intent(inout) :: b
    end subroutine psb_z_csrli_mv_from
  end interface


!!$  !> \memberof psb_z_csrli_sparse_mat
!!$  !! \see psb_z_base_mat_mod::psb_z_base_csput_a
!!$  interface
!!$    subroutine psb_z_csrli_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info)
!!$      import
!!$      class(psb_z_csrli_sparse_mat), intent(inout) :: a
!!$      complex(psb_dpk_), intent(in)      :: val(:)
!!$      integer(psb_ipk_), intent(in)             :: nz,ia(:), ja(:),&
!!$           &  imin,imax,jmin,jmax
!!$      integer(psb_ipk_), intent(out)            :: info
!!$    end subroutine psb_z_csrli_csput_a
!!$  end interface
!!$
!!$  !> \memberof psb_z_csrli_sparse_mat
!!$  !! \see psb_base_mat_mod::psb_base_csgetptn
!!$  interface
!!$    subroutine psb_z_csrli_csgetptn(imin,imax,a,nz,ia,ja,info,&
!!$         & jmin,jmax,iren,append,nzin,rscale,cscale)
!!$      import
!!$      class(psb_z_csrli_sparse_mat), intent(in) :: a
!!$      integer(psb_ipk_), intent(in)                  :: imin,imax
!!$      integer(psb_ipk_), intent(out)                 :: nz
!!$      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
!!$      integer(psb_ipk_),intent(out)                  :: info
!!$      logical, intent(in), optional        :: append
!!$      integer(psb_ipk_), intent(in), optional        :: iren(:)
!!$      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
!!$      logical, intent(in), optional        :: rscale,cscale
!!$    end subroutine psb_z_csrli_csgetptn
!!$  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_csgetrow
  interface
    subroutine psb_z_csrli_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale,chksz)
      import
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale,chksz
    end subroutine psb_z_csrli_csgetrow
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_cssv
  interface
    subroutine psb_z_csrli_cssv(alpha,a,x,beta,y,info,trans)
      import
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      complex(psb_dpk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_z_csrli_cssv
  end interface
  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_cssm
  interface
    subroutine psb_z_csrli_cssm(alpha,a,x,beta,y,info,trans)
      import
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      complex(psb_dpk_), intent(inout)       :: y(:,:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_z_csrli_cssm
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_csmv
  interface
    subroutine psb_z_csrli_csmv(alpha,a,x,beta,y,info,trans)
      import
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      complex(psb_dpk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_z_csrli_csmv
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_csmm
  interface
    subroutine psb_z_csrli_csmm(alpha,a,x,beta,y,info,trans)
      import
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      complex(psb_dpk_), intent(inout)       :: y(:,:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_z_csrli_csmm
  end interface


  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_maxval
  interface
    function psb_z_csrli_maxval(a) result(res)
      import
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_z_csrli_maxval
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_csnmi
  interface
    function psb_z_csrli_csnmi(a) result(res)
      import
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_z_csrli_csnmi
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_rowsum
  interface
    subroutine psb_z_csrli_rowsum(d,a)
      import
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_csrli_rowsum
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_arwsum
  interface
    subroutine psb_z_csrli_arwsum(d,a)
      import
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_csrli_arwsum
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_colsum
  interface
    subroutine psb_z_csrli_colsum(d,a)
      import
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_csrli_colsum
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_aclsum
  interface
    subroutine psb_z_csrli_aclsum(d,a)
      import
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_csrli_aclsum
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_get_diag
  interface
    subroutine psb_z_csrli_get_diag(a,d,info)
      import
      class(psb_z_csrli_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)     :: d(:)
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_csrli_get_diag
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_scal
  interface
    subroutine psb_z_csrli_scal(d,a,info,side)
      import
      class(psb_z_csrli_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: d(:)
      integer(psb_ipk_), intent(out)            :: info
      character, intent(in), optional :: side
    end subroutine psb_z_csrli_scal
  end interface

  !> \memberof psb_z_csrli_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_scals
  interface
    subroutine psb_z_csrli_scals(d,a,info)
      import
      class(psb_z_csrli_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: d
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_csrli_scals
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


  function z_csrli_get_fmt() result(res)
    implicit none
    character(len=5) :: res
    res = 'CSRLI'
  end function z_csrli_get_fmt


  function z_csrli_get_lambda(a) result(res)
    implicit none
    class(psb_z_csrli_sparse_mat), intent(in) :: a
    complex(psb_dpk_) :: res
    res = a%lambda
  end function z_csrli_get_lambda

  subroutine z_csrli_set_lambda(a,val) 
    implicit none
    class(psb_z_csrli_sparse_mat), intent(inout) :: a
    complex(psb_dpk_), intent(in) :: val
    a%lambda = val
  end subroutine z_csrli_set_lambda

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

  subroutine  z_csrli_free(a)
    implicit none

    class(psb_z_csrli_sparse_mat), intent(inout) :: a

    if (allocated(a%irp)) deallocate(a%irp)
    if (allocated(a%ja)) deallocate(a%ja)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(0_psb_ipk_)
    call a%set_ncols(0_psb_ipk_)
    a%lambda = zzero
    
    return

  end subroutine z_csrli_free


end module psb_z_csrli_mat_mod
