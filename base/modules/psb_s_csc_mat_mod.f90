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
! package: psb_s_csc_mat_mod
!
! This module contains the definition of the psb_s_csc_sparse_mat type
! which implements an actual storage format (the CSC in this case) for
! a sparse matrix as well as the related methods (those who are
! specific to the type and could not be defined higher in the
! hierarchy). We are at the bottom level of the inheritance chain.
!
! Please refere to psb_s_base_mat_mod for a detailed description
! of the various methods, and to psb_s_csc_impl for implementation details.
! 
module psb_s_csc_mat_mod

  use psb_s_base_mat_mod

  type, extends(psb_s_base_sparse_mat) :: psb_s_csc_sparse_mat

    integer(psb_ipk_), allocatable :: icp(:), ia(:)
    real(psb_spk_), allocatable :: val(:)

  contains
    procedure, pass(a) :: get_size     => s_csc_get_size
    procedure, pass(a) :: get_nzeros   => s_csc_get_nzeros
    procedure, nopass  :: get_fmt      => s_csc_get_fmt
    procedure, pass(a) :: sizeof       => s_csc_sizeof
    procedure, pass(a) :: s_csmm       => psb_s_csc_csmm
    procedure, pass(a) :: s_csmv       => psb_s_csc_csmv
    procedure, pass(a) :: s_inner_cssm => psb_s_csc_cssm
    procedure, pass(a) :: s_inner_cssv => psb_s_csc_cssv
    procedure, pass(a) :: s_scals      => psb_s_csc_scals
    procedure, pass(a) :: s_scal       => psb_s_csc_scal
    procedure, pass(a) :: maxval       => psb_s_csc_maxval
    procedure, pass(a) :: csnmi        => psb_s_csc_csnmi
    procedure, pass(a) :: csnm1        => psb_s_csc_csnm1
    procedure, pass(a) :: rowsum       => psb_s_csc_rowsum
    procedure, pass(a) :: arwsum       => psb_s_csc_arwsum
    procedure, pass(a) :: colsum       => psb_s_csc_colsum
    procedure, pass(a) :: aclsum       => psb_s_csc_aclsum
    procedure, pass(a) :: reallocate_nz => psb_s_csc_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_s_csc_allocate_mnnz
    procedure, pass(a) :: cp_to_coo    => psb_s_cp_csc_to_coo
    procedure, pass(a) :: cp_from_coo  => psb_s_cp_csc_from_coo
    procedure, pass(a) :: cp_to_fmt    => psb_s_cp_csc_to_fmt
    procedure, pass(a) :: cp_from_fmt  => psb_s_cp_csc_from_fmt
    procedure, pass(a) :: mv_to_coo    => psb_s_mv_csc_to_coo
    procedure, pass(a) :: mv_from_coo  => psb_s_mv_csc_from_coo
    procedure, pass(a) :: mv_to_fmt    => psb_s_mv_csc_to_fmt
    procedure, pass(a) :: mv_from_fmt  => psb_s_mv_csc_from_fmt
    procedure, pass(a) :: csput        => psb_s_csc_csput
    procedure, pass(a) :: get_diag     => psb_s_csc_get_diag
    procedure, pass(a) :: csgetptn     => psb_s_csc_csgetptn
    procedure, pass(a) :: s_csgetrow   => psb_s_csc_csgetrow
    procedure, pass(a) :: get_nz_col   => s_csc_get_nz_col
    procedure, pass(a) :: reinit       => psb_s_csc_reinit
    procedure, pass(a) :: trim         => psb_s_csc_trim
    procedure, pass(a) :: print        => psb_s_csc_print
    procedure, pass(a) :: free         => s_csc_free
    procedure, pass(a) :: mold         => psb_s_csc_mold
    procedure, pass(a) :: psb_s_csc_cp_from
    generic, public    :: cp_from => psb_s_csc_cp_from
    procedure, pass(a) :: psb_s_csc_mv_from
    generic, public    :: mv_from => psb_s_csc_mv_from

  end type psb_s_csc_sparse_mat

 private :: s_csc_get_nzeros, s_csc_free,  s_csc_get_fmt, &
       & s_csc_get_size, s_csc_sizeof, s_csc_get_nz_col

  interface
    subroutine  psb_s_csc_reallocate_nz(nz,a) 
      import :: psb_ipk_, psb_s_csc_sparse_mat
      integer(psb_ipk_), intent(in) :: nz
      class(psb_s_csc_sparse_mat), intent(inout) :: a
    end subroutine psb_s_csc_reallocate_nz
  end interface
  
  interface 
    subroutine psb_s_csc_reinit(a,clear)
      import :: psb_ipk_, psb_s_csc_sparse_mat
      class(psb_s_csc_sparse_mat), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_s_csc_reinit
  end interface
  
  interface
    subroutine  psb_s_csc_trim(a)
      import :: psb_ipk_, psb_s_csc_sparse_mat
      class(psb_s_csc_sparse_mat), intent(inout) :: a
    end subroutine psb_s_csc_trim
  end interface
  
  interface
    subroutine  psb_s_csc_allocate_mnnz(m,n,a,nz) 
      import :: psb_ipk_, psb_s_csc_sparse_mat
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_s_csc_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in), optional :: nz
    end subroutine psb_s_csc_allocate_mnnz
  end interface

  interface 
    subroutine psb_s_csc_mold(a,b,info) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_s_base_sparse_mat, psb_long_int_k_
      class(psb_s_csc_sparse_mat), intent(in)               :: a
      class(psb_s_base_sparse_mat), intent(out), allocatable :: b
      integer(psb_ipk_), intent(out)                                 :: info
    end subroutine psb_s_csc_mold
  end interface

  interface
    subroutine psb_s_csc_print(iout,a,iv,head,ivr,ivc)
      import :: psb_ipk_, psb_s_csc_sparse_mat
      integer(psb_ipk_), intent(in)               :: iout
      class(psb_s_csc_sparse_mat), intent(in) :: a   
      integer(psb_ipk_), intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer(psb_ipk_), intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_s_csc_print
  end interface
  
  interface 
    subroutine psb_s_cp_csc_to_coo(a,b,info) 
      import :: psb_ipk_, psb_s_coo_sparse_mat, psb_s_csc_sparse_mat
      class(psb_s_csc_sparse_mat), intent(in) :: a
      class(psb_s_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_s_cp_csc_to_coo
  end interface
  
  interface 
    subroutine psb_s_cp_csc_from_coo(a,b,info) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_s_coo_sparse_mat
      class(psb_s_csc_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_s_cp_csc_from_coo
  end interface
  
  interface 
    subroutine psb_s_cp_csc_to_fmt(a,b,info) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_s_base_sparse_mat
      class(psb_s_csc_sparse_mat), intent(in)   :: a
      class(psb_s_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                       :: info
    end subroutine psb_s_cp_csc_to_fmt
  end interface
  
  interface 
    subroutine psb_s_cp_csc_from_fmt(a,b,info) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_s_base_sparse_mat
      class(psb_s_csc_sparse_mat), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(in)   :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_s_cp_csc_from_fmt
  end interface
  
  interface 
    subroutine psb_s_mv_csc_to_coo(a,b,info) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_s_coo_sparse_mat
      class(psb_s_csc_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(inout)   :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_s_mv_csc_to_coo
  end interface
  
  interface 
    subroutine psb_s_mv_csc_from_coo(a,b,info) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_s_coo_sparse_mat
      class(psb_s_csc_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_s_mv_csc_from_coo
  end interface
  
  interface 
    subroutine psb_s_mv_csc_to_fmt(a,b,info) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_s_base_sparse_mat
      class(psb_s_csc_sparse_mat), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(inout)  :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_s_mv_csc_to_fmt
  end interface
  
  interface 
    subroutine psb_s_mv_csc_from_fmt(a,b,info) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_s_base_sparse_mat
      class(psb_s_csc_sparse_mat), intent(inout)  :: a
      class(psb_s_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                         :: info
    end subroutine psb_s_mv_csc_from_fmt
  end interface
  
  interface 
    subroutine psb_s_csc_cp_from(a,b)
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(inout) :: a
      type(psb_s_csc_sparse_mat), intent(in)   :: b
    end subroutine psb_s_csc_cp_from
  end interface
  
  interface 
    subroutine psb_s_csc_mv_from(a,b)
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(inout)  :: a
      type(psb_s_csc_sparse_mat), intent(inout) :: b
    end subroutine psb_s_csc_mv_from
  end interface
  
  
  interface 
    subroutine psb_s_csc_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in)      :: val(:)
      integer(psb_ipk_), intent(in)             :: nz,ia(:), ja(:),&
           &  imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: gtl(:)
    end subroutine psb_s_csc_csput
  end interface
  
  interface 
    subroutine psb_s_csc_csgetptn(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_s_csc_csgetptn
  end interface
  
  interface 
    subroutine psb_s_csc_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      real(psb_spk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_s_csc_csgetrow
  end interface

  interface 
    subroutine psb_s_csc_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_, psb_s_coo_sparse_mat
      class(psb_s_csc_sparse_mat), intent(in) :: a
      class(psb_s_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_s_csc_csgetblk
  end interface
    
  interface 
    subroutine psb_s_csc_cssv(alpha,a,x,beta,y,info,trans) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)          :: alpha, beta, x(:)
      real(psb_spk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_s_csc_cssv
    subroutine psb_s_csc_cssm(alpha,a,x,beta,y,info,trans) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_spk_), intent(inout)       :: y(:,:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_s_csc_cssm
  end interface
  
  interface 
    subroutine psb_s_csc_csmv(alpha,a,x,beta,y,info,trans) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)          :: alpha, beta, x(:)
      real(psb_spk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_s_csc_csmv
    subroutine psb_s_csc_csmm(alpha,a,x,beta,y,info,trans) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_spk_), intent(inout)       :: y(:,:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_s_csc_csmm
  end interface
  
  
  interface 
    function psb_s_csc_maxval(a) result(res)
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_s_csc_maxval
  end interface
  
  interface 
    function psb_s_csc_csnmi(a) result(res)
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_s_csc_csnmi
  end interface
  
  interface 
    function psb_s_csc_csnm1(a) result(res)
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_s_csc_csnm1
  end interface

  interface 
    subroutine psb_s_csc_rowsum(d,a) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_s_csc_rowsum
  end interface

  interface 
    subroutine psb_s_csc_arwsum(d,a) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_s_csc_arwsum
  end interface
  
  interface 
    subroutine psb_s_csc_colsum(d,a) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_s_csc_colsum
  end interface

  interface 
    subroutine psb_s_csc_aclsum(d,a) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_s_csc_aclsum
  end interface
    
  interface 
    subroutine psb_s_csc_get_diag(a,d,info) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(out)     :: d(:)
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_s_csc_get_diag
  end interface
  
  interface 
    subroutine psb_s_csc_scal(d,a,info,side) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in)      :: d(:)
      integer(psb_ipk_), intent(out)            :: info
      character, intent(in), optional :: side
    end subroutine psb_s_csc_scal
  end interface
  
  interface
    subroutine psb_s_csc_scals(d,a,info) 
      import :: psb_ipk_, psb_s_csc_sparse_mat, psb_spk_
      class(psb_s_csc_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in)      :: d
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_s_csc_scals
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

  
  function s_csc_sizeof(a) result(res)
    implicit none 
    class(psb_s_csc_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    res = 8 
    res = res + psb_sizeof_sp  * size(a%val)
    res = res + psb_sizeof_int * size(a%icp)
    res = res + psb_sizeof_int * size(a%ia)
      
  end function s_csc_sizeof

  function s_csc_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'CSC'
  end function s_csc_get_fmt
  
  function s_csc_get_nzeros(a) result(res)
    implicit none 
    class(psb_s_csc_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res
    res = a%icp(a%get_ncols()+1)-1
  end function s_csc_get_nzeros

  function s_csc_get_size(a) result(res)
    implicit none 
    class(psb_s_csc_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res

    res = 0
    
    if (allocated(a%ia)) then 
      if (res >= 0) then 
        res = min(res,size(a%ia))
      else 
        res = size(a%ia)
      end if
    end if
    if (allocated(a%val)) then 
      if (res >= 0) then 
        res = min(res,size(a%val))
      else 
        res = size(a%val)
      end if
    end if

  end function s_csc_get_size



  function  s_csc_get_nz_col(idx,a) result(res)
    use psb_const_mod
    implicit none
    
    class(psb_s_csc_sparse_mat), intent(in) :: a
    integer(psb_ipk_), intent(in)                  :: idx
    integer(psb_ipk_) :: res
    
    res = 0 
 
    if ((1<=idx).and.(idx<=a%get_ncols())) then 
      res = a%icp(idx+1)-a%icp(idx)
    end if
    
  end function s_csc_get_nz_col



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


  subroutine  s_csc_free(a) 
    implicit none 

    class(psb_s_csc_sparse_mat), intent(inout) :: a

    if (allocated(a%icp)) deallocate(a%icp)
    if (allocated(a%ia)) deallocate(a%ia)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(izero)
    call a%set_ncols(izero)
    
    return

  end subroutine s_csc_free

end module psb_s_csc_mat_mod
