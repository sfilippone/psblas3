!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
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
! package: psb_base_mat_mod
!
! This module contains the definition of the very basic object
! psb_base_sparse_mat holding some information common to all matrix
! type variants, such as number of rows and columns, whether the
! matrix is supposed to be triangular (upper or lower) and with a unit
! (i.e. assumed) diagonal, together with some state variables. This
! base class is in common among all variants of real/complex,
! short/long precision; as such, it only contains information that is
! inherently integer in nature.
!
! The methods associated to this class can be grouped into three sets:
! 1) Fully implemented methods: some methods such as get_nrows or
!    set_nrows can be fully implemented at this level.
! 2) Partially implemented methods: Some methods have an
!    implementation that is split between this level and the leaf
!    level. For example, the matrix transposition can be partially
!    done at this level (swapping of the rows and columns dimensions)
!    but it has to be completed by a method defined at the leaf level
!    (for actually transposing the row and column indices).
! 3) Other methods: There are a number of methods that are defined
!    (i.e their interface is defined) but not implemented at this
!    level. This methods will be overwritten at the leaf level with
!    an actual implementation. If it is not the case, the method
!    defined at this level will raise an error. These methods are
!    defined in the serial/f03/psb_base_mat_impl.f03 file
!
!

module psb_base_mat_mod
  
  use psb_const_mod 
  use psi_serial_mod
  
  !
  !> \namespace  psb_base_mod  \class  psb_base_sparse_mat
  !!  The basic data about your matrix.
  !!   This class is extended twice, to provide the various
  !!   data variations S/D/C/Z and to implement the actual
  !!   storage formats. The grandchild classes are then
  !!   encapsulated to implement the STATE design pattern.
  !!   We have an ambiguity in that the inner class has a
  !!   "state" variable; we hope the context will make it clear. 
  !!
  !!
  !! The methods associated to this class can be grouped into three sets:
  !! -  Fully implemented methods: some methods such as get_nrows or
  !!    set_nrows can be fully implemented at this level.
  !! -  Partially implemented methods: Some methods have an
  !!    implementation that is split between this level and the leaf
  !!    level. For example, the matrix transposition can be partially
  !!    done at this level (swapping of the rows and columns dimensions)
  !!    but it has to be completed by a method defined at the leaf level
  !!    (for actually transposing the row and column indices).
  !! -  Other methods: There are a number of methods that are defined
  !!    (i.e their interface is defined) but not implemented at this
  !!    level. This methods will be overwritten at the leaf level with
  !!    an actual implementation. If it is not the case, the method
  !!    defined at this level will raise an error. These methods are
  !!    defined in the serial/impl/psb_base_mat_impl.f90 file
  !!
  !

  type  :: psb_base_sparse_mat
    !> Row size
    integer(psb_ipk_), private     :: m
    !> Col size
    integer(psb_ipk_), private     :: n 
    !> Matrix state:
    !!    null:   pristine;
    !!    build:  it's being filled with entries;
    !!    assembled: ready to use in computations;
    !!    update: accepts coefficients but only
    !!            in already existing entries.
    !!    The transitions among the states are detailed in
    !!            psb_T_mat_mod. 
    integer(psb_ipk_), private     :: state
    !> How to treat duplicate elements when
    !!            transitioning from the BUILD to the ASSEMBLED state. 
    !!            While many formats would allow for duplicate
    !!            entries, it is much better to constrain the matrices
    !!            NOT to have duplicate entries, except while in the
    !!            BUILD state; in our overall design, only COO matrices
    !!            can ever be in the BUILD state, hence all other formats
    !!            cannot have duplicate entries.
    integer(psb_ipk_), private     :: duplicate 
    !> Is the matrix triangular? (must also be square)
    logical, private     :: triangle
    !> Is the matrix upper or lower? (only if  triangular)
    logical, private     :: upper
    !> Is the matrix diagonal stored or assumed unitary? (only if  triangular)
    logical, private     :: unitd
    !> Are the coefficients sorted ?
    logical, private     :: sorted
    logical, private     :: repeatable_updates=.false.

  contains 

    ! == = =================================
    !
    ! Getters 
    !
    !
    ! == = =================================
    procedure, pass(a) :: get_nrows   => psb_base_get_nrows
    procedure, pass(a) :: get_ncols   => psb_base_get_ncols
    procedure, pass(a) :: get_nzeros  => psb_base_get_nzeros
    procedure, pass(a) :: get_nz_row  => psb_base_get_nz_row
    procedure, pass(a) :: get_size    => psb_base_get_size
    procedure, pass(a) :: get_state   => psb_base_get_state
    procedure, pass(a) :: get_dupl    => psb_base_get_dupl
    procedure, nopass  :: get_fmt     => psb_base_get_fmt
    procedure, nopass  :: has_update  => psb_base_has_update
    procedure, pass(a) :: is_null     => psb_base_is_null
    procedure, pass(a) :: is_bld      => psb_base_is_bld
    procedure, pass(a) :: is_upd      => psb_base_is_upd
    procedure, pass(a) :: is_asb      => psb_base_is_asb
    procedure, pass(a) :: is_sorted   => psb_base_is_sorted
    procedure, pass(a) :: is_upper    => psb_base_is_upper
    procedure, pass(a) :: is_lower    => psb_base_is_lower
    procedure, pass(a) :: is_triangle => psb_base_is_triangle
    procedure, pass(a) :: is_unit     => psb_base_is_unit
    procedure, pass(a) :: is_by_rows  => psb_base_is_by_rows
    procedure, pass(a) :: is_by_cols  => psb_base_is_by_cols
    procedure, pass(a) :: is_repeatable_updates     => psb_base_is_repeatable_updates
    
    ! == = =================================
    !
    ! Setters 
    !
    ! == = =================================
    procedure, pass(a) :: set_nrows    => psb_base_set_nrows
    procedure, pass(a) :: set_ncols    => psb_base_set_ncols
    procedure, pass(a) :: set_dupl     => psb_base_set_dupl
    procedure, pass(a) :: set_state    => psb_base_set_state
    procedure, pass(a) :: set_null     => psb_base_set_null
    procedure, pass(a) :: set_bld      => psb_base_set_bld
    procedure, pass(a) :: set_upd      => psb_base_set_upd
    procedure, pass(a) :: set_asb      => psb_base_set_asb
    procedure, pass(a) :: set_sorted   => psb_base_set_sorted
    procedure, pass(a) :: set_upper    => psb_base_set_upper
    procedure, pass(a) :: set_lower    => psb_base_set_lower
    procedure, pass(a) :: set_triangle => psb_base_set_triangle
    procedure, pass(a) :: set_unit     => psb_base_set_unit

    procedure, pass(a) :: set_repeatable_updates     => psb_base_set_repeatable_updates


    ! == = =================================
    !
    ! Data management
    !
    ! == = =================================  
    procedure, pass(a) :: get_neigh  => psb_base_get_neigh
    procedure, pass(a) :: free       => psb_base_free
    procedure, pass(a) :: asb        => psb_base_mat_asb
    procedure, pass(a) :: trim       => psb_base_trim
    procedure, pass(a) :: reinit     => psb_base_reinit
    procedure, pass(a) :: allocate_mnnz => psb_base_allocate_mnnz
    procedure, pass(a) :: reallocate_nz => psb_base_reallocate_nz
    generic,   public  :: allocate   => allocate_mnnz
    generic,   public  :: reallocate => reallocate_nz


    procedure, pass(a) :: csgetptn => psb_base_csgetptn
    generic, public    :: csget => csgetptn
    procedure, pass(a) :: print => psb_base_sparse_print
    procedure, pass(a) :: sizeof => psb_base_sizeof
    procedure, pass(a) :: transp_1mat => psb_base_transp_1mat
    procedure, pass(a) :: transp_2mat => psb_base_transp_2mat
    generic, public    :: transp => transp_1mat, transp_2mat
    procedure, pass(a) :: transc_1mat => psb_base_transc_1mat
    procedure, pass(a) :: transc_2mat => psb_base_transc_2mat
    generic, public    :: transc => transc_1mat, transc_2mat

    !
    ! Sync: centerpiece of handling of external storage.
    ! Any derived class having extra storage upon sync
    ! will guarantee that both fortran/host side and
    ! external side contain the same data. The base
    ! version is only a placeholder. 
    !
    procedure, pass(a) :: sync          => psb_base_mat_sync
    procedure, pass(a) :: is_host       => psb_base_mat_is_host
    procedure, pass(a) :: is_dev        => psb_base_mat_is_dev
    procedure, pass(a) :: is_sync       => psb_base_mat_is_sync
    procedure, pass(a) :: set_host      => psb_base_mat_set_host
    procedure, pass(a) :: set_dev       => psb_base_mat_set_dev
    procedure, pass(a) :: set_sync      => psb_base_mat_set_sync
 
  end type psb_base_sparse_mat

  !>  Function: psb_base_get_nz_row
  !! \memberof  psb_base_sparse_mat
  !!  Interface for  the get_nz_row method. Equivalent to:
  !!    count(A(idx,:)/=0)
  !!    \param idx   The line we are interested in.
  !
  interface 
    function psb_base_get_nz_row(idx,a) result(res)
      import :: psb_ipk_, psb_long_int_k_, psb_base_sparse_mat
      integer(psb_ipk_), intent(in)                    :: idx
      class(psb_base_sparse_mat), intent(in) :: a
      integer(psb_ipk_) :: res
    end function psb_base_get_nz_row
  end interface
  
  !
  !>  Function: psb_base_get_nzeros
  !! \memberof  psb_base_sparse_mat
  !!  Interface for  the get_nzeros method. Equivalent to: 
  !!    count(A(:,:)/=0) 
  !
  interface 
    function psb_base_get_nzeros(a) result(res)
      import :: psb_ipk_, psb_long_int_k_, psb_base_sparse_mat
      class(psb_base_sparse_mat), intent(in) :: a
      integer(psb_ipk_) :: res
    end function psb_base_get_nzeros
  end interface

  !> Function get_size
  !! \memberof  psb_base_sparse_mat
  !!  how many items can A hold with
  !!           its current space allocation?
  !!           (as opposed to how many are
  !!            currently occupied) 
  !    
  interface 
    function psb_base_get_size(a) result(res)
      import :: psb_ipk_, psb_long_int_k_, psb_base_sparse_mat
      class(psb_base_sparse_mat), intent(in) :: a
      integer(psb_ipk_) :: res
    end function psb_base_get_size
  end interface

  !
  !> Function reinit: transition state from ASB to UPDATE
  !! \memberof  psb_base_sparse_mat
  !!  \param clear [true] explicitly zero out coefficients.
  !    
  interface 
    subroutine psb_base_reinit(a,clear)
      import :: psb_ipk_, psb_long_int_k_, psb_base_sparse_mat
      class(psb_base_sparse_mat), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_base_reinit
  end interface
  

  !
  !> Function
  !! \memberof  psb_base_sparse_mat
  !!  print on file in Matrix Market format. 
  !!  \param iout the output unit
  !!  \param iv(:) [none] renumber both row and column indices
  !!  \param head [none] a descriptive header for the matrix data
  !!  \param ivr(:) [none] renumbering for the rows
  !!  \param ivc(:) [none] renumbering for the cols
  !    
  interface 
    subroutine psb_base_sparse_print(iout,a,iv,head,ivr,ivc)
      import :: psb_ipk_, psb_long_int_k_, psb_base_sparse_mat
      integer(psb_ipk_), intent(in)               :: iout
      class(psb_base_sparse_mat), intent(in) :: a   
      integer(psb_ipk_), intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer(psb_ipk_), intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_base_sparse_print
  end interface


  !
  !> Function  getptn:
  !! \memberof  psb_base_sparse_mat
  !! \brief Get the pattern.
  !!
  !!
  !!         Return a list of NZ pairs
  !!           (IA(i),JA(i))
  !!         each identifying the position of a nonzero in A
  !!         between row indices IMIN:IMAX; 
  !!         IA,JA are reallocated as necessary.
  !!  \param imin  the minimum row index we are interested in 
  !!  \param imax  the minimum row index we are interested in 
  !!  \param nz the number of output coefficients
  !!  \param ia(:)  the output row indices
  !!  \param ja(:)  the output col indices
  !!  \param info  return code
  !!  \param jmin [1] minimum col index 
  !!  \param jmax [a\%get_ncols()] maximum col index 
  !!  \param iren(:) [none] an array to return renumbered indices   (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja 
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !           

  interface 
    subroutine psb_base_csgetptn(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import :: psb_ipk_, psb_long_int_k_, psb_base_sparse_mat
      class(psb_base_sparse_mat), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_base_csgetptn
  end interface
  
  !
  !> Function  get_neigh:
  !! \memberof  psb_base_sparse_mat
  !! \brief Get the neighbours.
  !!
  !!
  !!         Return a list of N indices of neighbours of index idx,
  !!         i.e. the indices of the nonzeros in row idx of matrix A
  !!         \param idx the index we are interested in
  !!         \param neigh(:) the list of indices, reallocated as necessary
  !!         \param n  the number of indices returned
  !!         \param info return code
  !!         \param lev [1] find neighbours recursively for LEV levels,
  !!               i.e. when lev=2 find neighours of neighbours, etc. 
  !           
  interface 
    subroutine psb_base_get_neigh(a,idx,neigh,n,info,lev)
      import :: psb_ipk_, psb_long_int_k_, psb_base_sparse_mat
      class(psb_base_sparse_mat), intent(in) :: a   
      integer(psb_ipk_), intent(in)                :: idx 
      integer(psb_ipk_), intent(out)               :: n   
      integer(psb_ipk_), allocatable, intent(out)  :: neigh(:)
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_ipk_), optional, intent(in)      :: lev 
    end subroutine psb_base_get_neigh
  end interface
  
  !         
  !
  !> Function  allocate_mnnz
  !! \memberof  psb_base_sparse_mat
  !! \brief Three-parameters version of allocate
  !!
  !!  \param m  number of rows
  !!  \param n  number of cols
  !!  \param nz [estimated internally] number of nonzeros to allocate for
  !
  interface 
    subroutine  psb_base_allocate_mnnz(m,n,a,nz) 
      import :: psb_ipk_, psb_long_int_k_, psb_base_sparse_mat
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_base_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in), optional  :: nz
    end subroutine psb_base_allocate_mnnz
  end interface

  
  !         
  !
  !> Function  reallocate_nz
  !! \memberof  psb_base_sparse_mat
  !! \brief One--parameter version of (re)allocate
  !!
  !!  \param nz  number of nonzeros to allocate for
  !
  interface 
    subroutine psb_base_reallocate_nz(nz,a) 
      import :: psb_ipk_, psb_long_int_k_, psb_base_sparse_mat
      integer(psb_ipk_), intent(in) :: nz
      class(psb_base_sparse_mat), intent(inout) :: a
    end subroutine psb_base_reallocate_nz
  end interface

  !         
  !> Function  free
  !! \memberof  psb_base_sparse_mat
  !! \brief destructor
  !
  interface 
    subroutine psb_base_free(a) 
      import :: psb_ipk_, psb_long_int_k_, psb_base_sparse_mat
      class(psb_base_sparse_mat), intent(inout) :: a
    end subroutine psb_base_free
  end interface
  
  !         
  !> Function  trim 
  !! \memberof  psb_base_sparse_mat
  !! \brief Memory trim
  !! Make sure the memory allocation of the sparse matrix is as tight as
  !! possible given the actual number of nonzeros it contains. 
  !
  interface 
    subroutine psb_base_trim(a) 
      import :: psb_ipk_, psb_long_int_k_, psb_base_sparse_mat
      class(psb_base_sparse_mat), intent(inout) :: a
    end subroutine psb_base_trim
  end interface
  
 
contains

  
  !         
  !> Function   sizeof
  !! \memberof  psb_base_sparse_mat
  !! \brief Memory occupation in byes
  !
  function psb_base_sizeof(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    res = 8
  end function psb_base_sizeof
 
  !         
  !> Function  get_fmt
  !! \memberof  psb_base_sparse_mat
  !! \brief return a short descriptive name (e.g. COO CSR etc.)
  !
  function psb_base_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'NULL'
  end function psb_base_get_fmt
  !         
  !> Function  has_update
  !! \memberof  psb_base_sparse_mat
  !! \brief Does the forma have the UPDATE functionality? 
  !
  function psb_base_has_update() result(res)
    implicit none 
    logical  :: res
    res = .true.
  end function psb_base_has_update
  
  !
  ! Standard getter functions: self-explaining. 
  !
  function psb_base_get_dupl(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res
    res = a%duplicate
  end function psb_base_get_dupl
 
 
  function psb_base_get_state(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res
    res = a%state
  end function psb_base_get_state
 
  function psb_base_get_nrows(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res
    res = a%m
  end function psb_base_get_nrows

  function psb_base_get_ncols(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res
    res = a%n
  end function psb_base_get_ncols

  subroutine  psb_base_set_nrows(m,a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    integer(psb_ipk_), intent(in) :: m
    a%m = m
  end subroutine psb_base_set_nrows

  subroutine  psb_base_set_ncols(n,a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    integer(psb_ipk_), intent(in) :: n
    a%n = n
  end subroutine psb_base_set_ncols
  

  subroutine  psb_base_set_state(n,a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    integer(psb_ipk_), intent(in) :: n
    a%state = n
  end subroutine psb_base_set_state


  subroutine  psb_base_set_dupl(n,a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    integer(psb_ipk_), intent(in) :: n
    a%duplicate = n
  end subroutine psb_base_set_dupl

  subroutine  psb_base_set_null(a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a

    a%state = psb_spmat_null_
  end subroutine psb_base_set_null

  subroutine  psb_base_set_bld(a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a

    a%state = psb_spmat_bld_
  end subroutine psb_base_set_bld

  subroutine  psb_base_set_upd(a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a

    a%state = psb_spmat_upd_
  end subroutine psb_base_set_upd

  subroutine  psb_base_set_asb(a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a

    a%state = psb_spmat_asb_
  end subroutine psb_base_set_asb

  subroutine psb_base_set_sorted(a,val) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%sorted = val
    else
      a%sorted = .true.
    end if
  end subroutine psb_base_set_sorted

  subroutine psb_base_set_triangle(a,val) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%triangle = val
    else
      a%triangle = .true.
    end if
  end subroutine psb_base_set_triangle

  subroutine psb_base_set_unit(a,val) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%unitd = val
    else
      a%unitd = .true.
    end if
  end subroutine psb_base_set_unit

  subroutine psb_base_set_lower(a,val) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%upper = .not.val
    else
      a%upper = .false.
    end if
  end subroutine psb_base_set_lower

  subroutine psb_base_set_upper(a,val) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%upper = val
    else
      a%upper = .true.
    end if
  end subroutine psb_base_set_upper

  subroutine psb_base_set_repeatable_updates(a,val) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%repeatable_updates = val
    else
      a%repeatable_updates = .true.
    end if
  end subroutine psb_base_set_repeatable_updates

  function psb_base_is_triangle(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%triangle
  end function psb_base_is_triangle

  function psb_base_is_unit(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%unitd
  end function psb_base_is_unit

  function psb_base_is_upper(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%upper
  end function psb_base_is_upper

  function psb_base_is_lower(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = .not.a%upper
  end function psb_base_is_lower

  function psb_base_is_null(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psb_spmat_null_)
  end function psb_base_is_null

  function psb_base_is_bld(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psb_spmat_bld_)
  end function psb_base_is_bld

  function psb_base_is_upd(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psb_spmat_upd_)
  end function psb_base_is_upd

  function psb_base_is_asb(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psb_spmat_asb_)
  end function psb_base_is_asb

  function psb_base_is_sorted(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%sorted
  end function psb_base_is_sorted


  function psb_base_is_by_rows(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = .false.
  end function psb_base_is_by_rows

  function psb_base_is_by_cols(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = .false.
  end function psb_base_is_by_cols

  function psb_base_is_repeatable_updates(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%repeatable_updates
  end function psb_base_is_repeatable_updates

  
  !
  !  TRANSP: note sorted=.false.
  !    better invoke a fix() too many than
  !    regret it later...
  !
  subroutine psb_base_transp_2mat(a,b)
    implicit none 
    
    class(psb_base_sparse_mat), intent(in)  :: a
    class(psb_base_sparse_mat), intent(out) :: b
    
    b%m         = a%n
    b%n         = a%m
    b%state     = a%state
    b%duplicate = a%duplicate
    b%triangle  = a%triangle
    b%unitd     = a%unitd
    b%upper     = .not.a%upper
    b%sorted    = .false.
    b%repeatable_updates = .false.
    
  end subroutine psb_base_transp_2mat

  subroutine psb_base_transc_2mat(a,b)
    implicit none 
    
    class(psb_base_sparse_mat), intent(in)  :: a
    class(psb_base_sparse_mat), intent(out) :: b

    
    b%m         = a%n
    b%n         = a%m
    b%state     = a%state
    b%duplicate = a%duplicate
    b%triangle  = a%triangle
    b%unitd     = a%unitd
    b%upper     = .not.a%upper
    b%sorted    = .false.
    b%repeatable_updates = .false.

  end subroutine psb_base_transc_2mat

  subroutine psb_base_transp_1mat(a)
    implicit none 
    
    class(psb_base_sparse_mat), intent(inout) :: a
    integer(psb_ipk_) :: itmp

    itmp        = a%m
    a%m         = a%n
    a%n         = itmp
    a%state     = a%state
    a%duplicate = a%duplicate
    a%triangle  = a%triangle
    a%unitd     = a%unitd
    a%upper     = .not.a%upper
    a%sorted    = .false.
    a%repeatable_updates = .false.
    
  end subroutine psb_base_transp_1mat

  subroutine psb_base_transc_1mat(a)
    implicit none 
    
    class(psb_base_sparse_mat), intent(inout) :: a
    
    call a%transp() 
  end subroutine psb_base_transc_1mat



  !
  !> Function  base_asb:
  !! \memberof  psb_base_sparse_mat
  !! \brief Sync: base version calls sync and the set_asb.
  !!           
  !
  subroutine psb_base_mat_asb(a)
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    
    call a%sync()
    call a%set_asb()
  end subroutine psb_base_mat_asb
  !
  ! The base version of SYNC & friends does nothing, it's just
  ! a placeholder.
  ! 
  !
  !> Function  base_sync:
  !! \memberof  psb_base_sparse_mat
  !! \brief Sync: base version is a no-op.
  !!           
  !
  subroutine psb_base_mat_sync(a)
    implicit none 
    class(psb_base_sparse_mat), target, intent(in) :: a
    
  end subroutine psb_base_mat_sync

  !
  !> Function  base_set_host:
  !! \memberof  psb_base_sparse_mat
  !! \brief Set_host: base version is a no-op.
  !!           
  !
  subroutine psb_base_mat_set_host(a)
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    
  end subroutine psb_base_mat_set_host

  !
  !> Function  base_set_dev:
  !! \memberof  psb_base_sparse_mat
  !! \brief Set_dev: base version is a no-op.
  !!           
  !
  subroutine psb_base_mat_set_dev(a)
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    
  end subroutine psb_base_mat_set_dev

  !
  !> Function  base_set_sync:
  !! \memberof  psb_base_sparse_mat
  !! \brief Set_sync: base version is a no-op.
  !!           
  !
  subroutine psb_base_mat_set_sync(a)
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    
  end subroutine psb_base_mat_set_sync

  !
  !> Function  base_is_dev:
  !! \memberof  psb_base_sparse_mat
  !! \brief Is matrix on eaternal device    .
  !!           
  !
  function psb_base_mat_is_dev(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical  :: res
  
    res = .false.
  end function psb_base_mat_is_dev
  
  !
  !> Function  base_is_host
  !! \memberof  psb_base_sparse_mat
  !! \brief Is matrix on standard memory    .
  !!           
  !
  function psb_base_mat_is_host(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical  :: res

    res = .true.
  end function psb_base_mat_is_host

  !
  !> Function  base_is_sync
  !! \memberof  psb_base_sparse_mat
  !! \brief Is matrix on sync               .
  !!           
  !
  function psb_base_mat_is_sync(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical  :: res

    res = .true.
  end function psb_base_mat_is_sync

end module psb_base_mat_mod

