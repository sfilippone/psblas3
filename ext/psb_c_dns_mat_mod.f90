!                Parallel Sparse BLAS   GPU plugin 
!      (C) Copyright 2013
!  
!                         Salvatore Filippone
!                         Alessandro Fanfarillo
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
module psb_c_dns_mat_mod

  use psb_c_base_mat_mod

  type, extends(psb_c_base_sparse_mat) :: psb_c_dns_sparse_mat
    !
    ! DNS format: a very simple dense matrix storage
    !  psb_spk_ : kind for double precision reals
    !  psb_ipk_: kind for normal integers.
    !  psb_sizeof_dp: variable holding size in bytes of
    !                 a double
    !  psb_sizeof_ip: size in bytes of an integer
    !
    !  psb_realloc(n,v,info)    Reallocate: does what it says
    !  psb_realloc(m,n,a,info)  on rank 1 and 2 arrays, may start
    !                           from unallocated
    ! 
    ! 
    integer(psb_ipk_) :: nnz
    complex(psb_spk_), allocatable :: val(:,:)

  contains
    procedure, pass(a) :: get_size     => c_dns_get_size
    procedure, pass(a) :: get_nzeros   => c_dns_get_nzeros
    procedure, nopass  :: get_fmt      => c_dns_get_fmt
    procedure, pass(a) :: sizeof       => c_dns_sizeof
    procedure, pass(a) :: csmv         => psb_c_dns_csmv
    procedure, pass(a) :: csmm         => psb_c_dns_csmm
    procedure, pass(a) :: csnmi        => psb_c_dns_csnmi
    procedure, pass(a) :: reallocate_nz => psb_c_dns_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_c_dns_allocate_mnnz
    procedure, pass(a) :: cp_to_coo    => psb_c_cp_dns_to_coo
    procedure, pass(a) :: cp_from_coo  => psb_c_cp_dns_from_coo
    procedure, pass(a) :: mv_to_coo    => psb_c_mv_dns_to_coo
    procedure, pass(a) :: mv_from_coo  => psb_c_mv_dns_from_coo
    procedure, pass(a) :: get_diag     => psb_c_dns_get_diag
    procedure, pass(a) :: csgetrow     => psb_c_dns_csgetrow
    procedure, pass(a) :: get_nz_row   => c_dns_get_nz_row
    procedure, pass(a) :: trim         => psb_c_dns_trim
    procedure, pass(a) :: free         => c_dns_free
    procedure, pass(a) :: mold         => psb_c_dns_mold

  end type psb_c_dns_sparse_mat

  private :: c_dns_get_nzeros, c_dns_free,  c_dns_get_fmt, &
       & c_dns_get_size, c_dns_sizeof, c_dns_get_nz_row

  !         
  !
  !> Function  reallocate_nz
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief One--parameters version of (re)allocate
  !!
  !!  \param nz  number of nonzeros to allocate for
  !!             i.e. makes sure that the internal storage
  !!             allows for NZ coefficients and their indices. 
  !  
  interface
    subroutine  psb_c_dns_reallocate_nz(nz,a) 
      import :: psb_c_dns_sparse_mat, psb_ipk_
      integer(psb_ipk_), intent(in) :: nz
      class(psb_c_dns_sparse_mat), intent(inout) :: a
    end subroutine psb_c_dns_reallocate_nz
  end interface
  
  !> Function  trim 
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief Memory trim
  !! Make sure the memory allocation of the sparse matrix is as tight as
  !! possible given the actual number of nonzeros it contains. 
  !
  interface
    subroutine  psb_c_dns_trim(a)
      import :: psb_c_dns_sparse_mat
      class(psb_c_dns_sparse_mat), intent(inout) :: a
    end subroutine psb_c_dns_trim
  end interface
  
  !
  !> Function  mold:
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief Allocate a class(psb_c_dns_sparse_mat) with the
  !!     same dynamic type as the input.
  !!     This is equivalent to allocate(  mold=  ) and is provided
  !!     for those compilers not yet supporting mold.
  !!   \param b The output variable
  !!   \param info return code
  ! 
  interface 
    subroutine psb_c_dns_mold(a,b,info) 
      import :: psb_c_dns_sparse_mat, psb_c_base_sparse_mat, psb_epk_, psb_ipk_
      class(psb_c_dns_sparse_mat), intent(in)                  :: a
      class(psb_c_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                             :: info
    end subroutine psb_c_dns_mold
  end interface

  !         
  !
  !> Function  allocate_mnnz
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief Three-parameters version of allocate
  !!
  !!  \param m  number of rows
  !!  \param n  number of cols
  !!  \param nz [estimated internally] number of nonzeros to allocate for
  !
  interface
    subroutine  psb_c_dns_allocate_mnnz(m,n,a,nz) 
      import :: psb_c_dns_sparse_mat, psb_ipk_
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_c_dns_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in), optional :: nz
    end subroutine psb_c_dns_allocate_mnnz
  end interface
  
  !
  !> Function  cp_to_coo:
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief Copy and convert to psb_c_coo_sparse_mat
  !!        Invoked from the source object.
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_c_cp_dns_to_coo(a,b,info) 
      import :: psb_c_coo_sparse_mat, psb_c_dns_sparse_mat, psb_ipk_
      class(psb_c_dns_sparse_mat), intent(in)    :: a
      class(psb_c_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_c_cp_dns_to_coo
  end interface
  
  !
  !> Function  cp_from_coo:
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief Copy and convert from psb_c_coo_sparse_mat
  !!        Invoked from the target object.
  !!   \param b The input variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_c_cp_dns_from_coo(a,b,info) 
      import :: psb_c_dns_sparse_mat, psb_c_coo_sparse_mat, psb_ipk_
      class(psb_c_dns_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_c_cp_dns_from_coo
  end interface
  
  !
  !> Function  mv_to_coo:
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief Convert to psb_c_coo_sparse_mat, freeing the source.
  !!        Invoked from the source object.
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_c_mv_dns_to_coo(a,b,info) 
      import :: psb_c_dns_sparse_mat, psb_c_coo_sparse_mat, psb_ipk_
      class(psb_c_dns_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_c_mv_dns_to_coo
  end interface
  
  !
  !> Function  mv_from_coo:
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief Convert from psb_c_coo_sparse_mat, freeing the source.
  !!        Invoked from the target object.
  !!   \param b The input variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_c_mv_dns_from_coo(a,b,info) 
      import :: psb_c_dns_sparse_mat, psb_c_coo_sparse_mat, psb_ipk_
      class(psb_c_dns_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_c_mv_dns_from_coo
  end interface
  
  !
  !
  !> Function  csgetrow:
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief Get a (subset of) row(s)
  !!        
  !!        getrow is the basic method by which the other (getblk, clip) can
  !!        be implemented.
  !!        
  !!        Returns the set
  !!           NZ, IA(1:nz), JA(1:nz), VAL(1:NZ)
  !!         each identifying the position of a nonzero in A
  !!         i.e.
  !!           VAL(1:NZ) = A(IA(1:NZ),JA(1:NZ))
  !!         with IMIN<=IA(:)<=IMAX
  !!         with JMIN<=JA(:)<=JMAX
  !!         IA,JA are reallocated as necessary.
  !!         
  !!  \param imin  the minimum row index we are interested in 
  !!  \param imax  the minimum row index we are interested in 
  !!  \param nz the number of output coefficients
  !!  \param ia(:)  the output row indices
  !!  \param ja(:)  the output col indices
  !!  \param val(:)  the output coefficients
  !!  \param info  return code
  !!  \param jmin [1] minimum col index 
  !!  \param jmax [a\%get_ncols()] maximum col index 
  !!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja 
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !!           
  !
  interface 
    subroutine psb_c_dns_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale,chksz)
      import :: psb_c_dns_sparse_mat, psb_spk_, psb_ipk_
      class(psb_c_dns_sparse_mat), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional           :: append
      integer(psb_ipk_), intent(in), optional :: iren(:)
      integer(psb_ipk_), intent(in), optional :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale,chksz
    end subroutine psb_c_dns_csgetrow
  end interface

  
  
  !> Function  csmv:
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief Product by a dense rank 1 array.
  !!
  !!        Compute
  !!           Y = alpha*op(A)*X + beta*Y
  !!
  !! \param alpha  Scaling factor for Ax
  !! \param A      the input sparse matrix
  !! \param x(:)   the input dense X
  !! \param beta   Scaling factor for y
  !! \param y(:)   the input/output dense Y
  !! \param info   return code
  !! \param trans  [N] Whether to use A (N), its transpose (T)
  !!               or its conjugate transpose (C)
  !!
  !
  interface 
    subroutine psb_c_dns_csmv(alpha,a,x,beta,y,info,trans) 
      import :: psb_c_dns_sparse_mat, psb_spk_, psb_ipk_
      class(psb_c_dns_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)       :: alpha, beta, x(:)
      complex(psb_spk_), intent(inout)    :: y(:)
      integer(psb_ipk_), intent(out)    :: info
      character, optional, intent(in)   :: trans
    end subroutine psb_c_dns_csmv
  end interface
  
  !> Function  csmm:
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief Product by a dense rank 2 array.
  !!
  !!        Compute
  !!           Y = alpha*op(A)*X + beta*Y
  !!
  !! \param alpha  Scaling factor for Ax
  !! \param A      the input sparse matrix
  !! \param x(:,:)   the input dense X
  !! \param beta   Scaling factor for y
  !! \param y(:,:)   the input/output dense Y
  !! \param info   return code
  !! \param trans  [N] Whether to use A (N), its transpose (T)
  !!               or its conjugate transpose (C)
  !!
  !
  interface 
    subroutine psb_c_dns_csmm(alpha,a,x,beta,y,info,trans) 
      import :: psb_c_dns_sparse_mat, psb_spk_, psb_ipk_
      class(psb_c_dns_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
      complex(psb_spk_), intent(inout)       :: y(:,:)
      integer(psb_ipk_), intent(out)      :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_c_dns_csmm
  end interface
  
  !
  !
  !> Function  csnmi:
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief Operator infinity norm
  !!     CSNMI = MAXVAL(SUM(ABS(A(:,:)),dim=2))
  !! 
  !
  interface 
    function psb_c_dns_csnmi(a) result(res)
      import :: psb_c_dns_sparse_mat, psb_spk_
      class(psb_c_dns_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_c_dns_csnmi
  end interface
  
  !
  !> Function  get_diag:
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief Extract the diagonal of A. 
  !!        
  !!   D(i) = A(i:i), i=1:min(nrows,ncols)
  !!
  !! \param d(:)  The output diagonal
  !! \param info  return code. 
  ! 
  interface 
    subroutine psb_c_dns_get_diag(a,d,info) 
      import :: psb_c_dns_sparse_mat, psb_spk_, psb_ipk_
      class(psb_c_dns_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(out)     :: d(:)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_c_dns_get_diag
  end interface
  

contains 

  !         
  !> Function  sizeof
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief Memory occupation in bytes
  !
  function c_dns_sizeof(a) result(res)
    implicit none 
    class(psb_c_dns_sparse_mat), intent(in) :: a
    integer(psb_epk_) :: res

    res = psb_sizeof_dp  * size(a%val)
    res = res + psb_sizeof_ip 
      
  end function c_dns_sizeof

  !         
  !> Function  get_fmt
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief return a short descriptive name (e.g. COO CSR etc.)
  !
  function c_dns_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'DNS'
  end function c_dns_get_fmt
  
  !         
  !> Function  get_nzeros
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief Current number of nonzero entries
  !
  function c_dns_get_nzeros(a) result(res)
    implicit none 
    class(psb_c_dns_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res
    res = a%nnz
  end function c_dns_get_nzeros

  !         
  !> Function  get_size
  !! \memberof  psb_c_dns_sparse_mat
  !! \brief Maximum number of nonzeros the current structure can hold
  ! this is fixed once you initialize the matrix, with dense storage
  ! you can hold up to MxN entries
  function c_dns_get_size(a) result(res)
    implicit none 
    class(psb_c_dns_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res

    res = size(a%val)

  end function c_dns_get_size


  !
  !> Function get_nz_row.
  !! \memberof  psb_c_coo_sparse_mat
  !! \brief How many nonzeros in a row?
  !!
  !! \param idx  The row to search.
  !!
  !
  function  c_dns_get_nz_row(idx,a) result(res)

    implicit none
    
    class(psb_c_dns_sparse_mat), intent(in) :: a
    integer(psb_ipk_), intent(in)                  :: idx
    integer(psb_ipk_)                              :: res
    
    res = 0 
 
    if ((1<=idx).and.(idx<=a%get_nrows())) then 
      res = count(a%val(idx,:) /= dzero)
    end if
    
  end function c_dns_get_nz_row

  !         
  !> Function  free
  !! \memberof  psb_c_dns_sparse_mat
  !!            Name says all 

  subroutine  c_dns_free(a) 
    implicit none 

    class(psb_c_dns_sparse_mat), intent(inout) :: a

    if (allocated(a%val)) deallocate(a%val)
    a%nnz = 0


    !
    ! Mark the object as empty just in case
    ! 
    call a%set_null()
    call a%set_nrows(izero)
    call a%set_ncols(izero)
    
    return

  end subroutine c_dns_free


end module psb_c_dns_mat_mod
