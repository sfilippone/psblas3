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
  

module psb_c_ell_mat_mod

  use psb_c_base_mat_mod

  type, extends(psb_c_base_sparse_mat) :: psb_c_ell_sparse_mat
    !
    ! ITPACK/ELL format, extended.
    ! Based on M. Heroux "A proposal for a sparse BLAS toolkit".
    !  IRN is our addition, should help in transferring to/from
    !  other formats (should come in handy for GPUs). 
    ! Notes:
    !  1. JA holds the column indices, padded with the row index.
    !  2. VAL holds the coefficients, padded with zeros
    !  3. IDIAG hold the position of the diagonal element
    !     or 0 if it is not there, but is only relevant for
    !     triangular matrices. In particular, a unit triangular matrix
    !     will have IDIAG==0.
    !  4. IRN holds the actual number of nonzeros stored in each row
    !  5. Within a row, the indices are sorted for use of SV. 
    !     
    
    integer(psb_ipk_) :: nzt
    integer(psb_ipk_), allocatable :: irn(:), ja(:,:), idiag(:)
    complex(psb_spk_), allocatable :: val(:,:)

  contains
    procedure, pass(a) :: is_by_rows   => c_ell_is_by_rows
    procedure, pass(a) :: get_size     => c_ell_get_size
    procedure, pass(a) :: get_nzeros   => c_ell_get_nzeros
    procedure, nopass  :: get_fmt      => c_ell_get_fmt
    procedure, pass(a) :: sizeof       => c_ell_sizeof
    procedure, pass(a) :: csmm         => psb_c_ell_csmm
    procedure, pass(a) :: csmv         => psb_c_ell_csmv
    procedure, pass(a) :: inner_cssm   => psb_c_ell_cssm
    procedure, pass(a) :: inner_cssv   => psb_c_ell_cssv
    procedure, pass(a) :: scals        => psb_c_ell_scals
    procedure, pass(a) :: scalv        => psb_c_ell_scal
    procedure, pass(a) :: maxval       => psb_c_ell_maxval
    procedure, pass(a) :: csnmi        => psb_c_ell_csnmi
    procedure, pass(a) :: csnm1        => psb_c_ell_csnm1
    procedure, pass(a) :: rowsum       => psb_c_ell_rowsum
    procedure, pass(a) :: arwsum       => psb_c_ell_arwsum
    procedure, pass(a) :: colsum       => psb_c_ell_colsum
    procedure, pass(a) :: aclsum       => psb_c_ell_aclsum
    procedure, pass(a) :: reallocate_nz => psb_c_ell_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_c_ell_allocate_mnnz
    procedure, pass(a) :: cp_to_coo    => psb_c_cp_ell_to_coo
    procedure, pass(a) :: cp_from_coo  => psb_c_cp_ell_from_coo
    procedure, pass(a) :: cp_to_fmt    => psb_c_cp_ell_to_fmt
    procedure, pass(a) :: cp_from_fmt  => psb_c_cp_ell_from_fmt
    procedure, pass(a) :: mv_to_coo    => psb_c_mv_ell_to_coo
    procedure, pass(a) :: mv_from_coo  => psb_c_mv_ell_from_coo
    procedure, pass(a) :: mv_to_fmt    => psb_c_mv_ell_to_fmt
    procedure, pass(a) :: mv_from_fmt  => psb_c_mv_ell_from_fmt
    procedure, pass(a) :: csput_a      => psb_c_ell_csput_a
    procedure, pass(a) :: get_diag     => psb_c_ell_get_diag
    procedure, pass(a) :: csgetptn     => psb_c_ell_csgetptn
    procedure, pass(a) :: csgetrow     => psb_c_ell_csgetrow
    procedure, pass(a) :: get_nz_row   => c_ell_get_nz_row
    procedure, pass(a) :: reinit       => psb_c_ell_reinit
    procedure, pass(a) :: trim         => psb_c_ell_trim
    procedure, pass(a) :: print        => psb_c_ell_print
    procedure, pass(a) :: free         => c_ell_free
    procedure, pass(a) :: mold         => psb_c_ell_mold
    procedure, pass(a) :: get_nrm      => c_ell_get_nrm

  end type psb_c_ell_sparse_mat

  private :: c_ell_get_nzeros, c_ell_free,  c_ell_get_fmt, &
       & c_ell_get_size, c_ell_sizeof, c_ell_get_nz_row, &
       & c_ell_is_by_rows

  interface
    subroutine  psb_c_ell_reallocate_nz(nz,a) 
      import :: psb_c_ell_sparse_mat, psb_ipk_
      integer(psb_ipk_), intent(in) :: nz
      class(psb_c_ell_sparse_mat), intent(inout) :: a
    end subroutine psb_c_ell_reallocate_nz
  end interface
  
  interface 
    subroutine psb_c_ell_reinit(a,clear)
      import :: psb_c_ell_sparse_mat
      class(psb_c_ell_sparse_mat), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_c_ell_reinit
  end interface
  
  interface
    subroutine  psb_c_ell_trim(a)
      import :: psb_c_ell_sparse_mat
      class(psb_c_ell_sparse_mat), intent(inout) :: a
    end subroutine psb_c_ell_trim
  end interface
  
  interface 
    subroutine psb_c_ell_mold(a,b,info) 
      import :: psb_c_ell_sparse_mat, psb_c_base_sparse_mat, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(in)                  :: a
      class(psb_c_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_c_ell_mold
  end interface

  interface
    subroutine  psb_c_ell_allocate_mnnz(m,n,a,nz) 
      import :: psb_c_ell_sparse_mat, psb_ipk_
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_c_ell_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in), optional :: nz
    end subroutine psb_c_ell_allocate_mnnz
  end interface
  
  interface
    subroutine psb_c_ell_print(iout,a,iv,head,ivr,ivc)
      import :: psb_c_ell_sparse_mat, psb_ipk_, psb_lpk_
      integer(psb_ipk_), intent(in)           :: iout
      class(psb_c_ell_sparse_mat), intent(in) :: a   
      integer(psb_lpk_), intent(in), optional :: iv(:)
      character(len=*), optional              :: head
      integer(psb_lpk_), intent(in), optional :: ivr(:), ivc(:)
    end subroutine psb_c_ell_print
  end interface
  
  interface 
    subroutine psb_c_cp_ell_to_coo(a,b,info) 
      import :: psb_c_coo_sparse_mat, psb_c_ell_sparse_mat, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(in)    :: a
      class(psb_c_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_c_cp_ell_to_coo
  end interface
  
  interface 
    subroutine psb_c_cp_ell_from_coo(a,b,info) 
      import :: psb_c_ell_sparse_mat, psb_c_coo_sparse_mat, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_c_cp_ell_from_coo
  end interface
  
  interface 
    subroutine psb_c_cp_ell_to_fmt(a,b,info) 
      import :: psb_c_ell_sparse_mat, psb_c_base_sparse_mat, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(in)     :: a
      class(psb_c_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_c_cp_ell_to_fmt
  end interface
  
  interface 
    subroutine psb_c_cp_ell_from_fmt(a,b,info) 
      import :: psb_c_ell_sparse_mat, psb_c_base_sparse_mat, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(inout) :: a
      class(psb_c_base_sparse_mat), intent(in)   :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_c_cp_ell_from_fmt
  end interface
  
  interface 
    subroutine psb_c_mv_ell_to_coo(a,b,info) 
      import :: psb_c_ell_sparse_mat, psb_c_coo_sparse_mat, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_c_mv_ell_to_coo
  end interface
  
  interface 
    subroutine psb_c_mv_ell_from_coo(a,b,info) 
      import :: psb_c_ell_sparse_mat, psb_c_coo_sparse_mat, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_c_mv_ell_from_coo
  end interface
  
  interface 
    subroutine psb_c_mv_ell_to_fmt(a,b,info) 
      import :: psb_c_ell_sparse_mat, psb_c_base_sparse_mat, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(inout)  :: a
      class(psb_c_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_c_mv_ell_to_fmt
  end interface
  
  interface 
    subroutine psb_c_mv_ell_from_fmt(a,b,info) 
      import :: psb_c_ell_sparse_mat, psb_c_base_sparse_mat, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(inout)  :: a
      class(psb_c_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_c_mv_ell_from_fmt
  end interface
  
  interface 
    subroutine psb_c_ell_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info) 
      import :: psb_c_ell_sparse_mat, psb_spk_, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)      :: val(:)
      integer(psb_ipk_), intent(in)             :: nz,ia(:), ja(:),&
           &  imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_c_ell_csput_a
  end interface
  
  interface 
    subroutine psb_c_ell_csgetptn(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import :: psb_c_ell_sparse_mat, psb_spk_, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(in)        :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_c_ell_csgetptn
  end interface
  
  interface 
    subroutine psb_c_ell_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale,chksz)
      import :: psb_c_ell_sparse_mat, psb_spk_, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(in)        :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale,chksz
    end subroutine psb_c_ell_csgetrow
  end interface

  interface 
    subroutine psb_c_ell_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
      import :: psb_c_ell_sparse_mat, psb_spk_, psb_c_coo_sparse_mat, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(in)    :: a
      class(psb_c_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(in)              :: imin,imax
      integer(psb_ipk_),intent(out)              :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional    :: iren(:)
      integer(psb_ipk_), intent(in), optional    :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_c_ell_csgetblk
  end interface
    
  interface 
    subroutine psb_c_ell_cssv(alpha,a,x,beta,y,info,trans) 
      import :: psb_c_ell_sparse_mat, psb_spk_, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:)
      complex(psb_spk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)           :: info
      character, optional, intent(in)          :: trans
    end subroutine psb_c_ell_cssv
    subroutine psb_c_ell_cssm(alpha,a,x,beta,y,info,trans) 
      import :: psb_c_ell_sparse_mat, psb_spk_, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
      complex(psb_spk_), intent(inout)       :: y(:,:)
      integer(psb_ipk_), intent(out)       :: info
      character, optional, intent(in)      :: trans
    end subroutine psb_c_ell_cssm
  end interface
  
  interface 
    subroutine psb_c_ell_csmv(alpha,a,x,beta,y,info,trans,ivshft) 
      import :: psb_c_ell_sparse_mat, psb_spk_, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:)
      complex(psb_spk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)       :: info
      character, optional, intent(in)      :: trans
      integer(psb_ipk_), optional, intent(in) :: ivshft
    end subroutine psb_c_ell_csmv
    subroutine psb_c_ell_csmm(alpha,a,x,beta,y,info,trans,ivshft) 
      import :: psb_c_ell_sparse_mat, psb_spk_, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
      complex(psb_spk_), intent(inout)       :: y(:,:)
      integer(psb_ipk_), intent(out)       :: info
      character, optional, intent(in)      :: trans
      integer(psb_ipk_), optional, intent(in) :: ivshft
    end subroutine psb_c_ell_csmm
  end interface
  
  
  interface 
    function psb_c_ell_maxval(a) result(res)
      import :: psb_c_ell_sparse_mat, psb_spk_
      class(psb_c_ell_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_c_ell_maxval
  end interface
  
  interface 
    function psb_c_ell_csnmi(a) result(res)
      import :: psb_c_ell_sparse_mat, psb_spk_
      class(psb_c_ell_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_c_ell_csnmi
  end interface
  
  interface 
    function psb_c_ell_csnm1(a) result(res)
      import :: psb_c_ell_sparse_mat, psb_spk_
      class(psb_c_ell_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_c_ell_csnm1
  end interface

  interface 
    subroutine psb_c_ell_rowsum(d,a) 
      import :: psb_c_ell_sparse_mat, psb_spk_
      class(psb_c_ell_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_c_ell_rowsum
  end interface

  interface 
    subroutine psb_c_ell_arwsum(d,a) 
      import :: psb_c_ell_sparse_mat, psb_spk_
      class(psb_c_ell_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_c_ell_arwsum
  end interface
  
  interface 
    subroutine psb_c_ell_colsum(d,a) 
      import :: psb_c_ell_sparse_mat, psb_spk_
      class(psb_c_ell_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_c_ell_colsum
  end interface

  interface 
    subroutine psb_c_ell_aclsum(d,a) 
      import :: psb_c_ell_sparse_mat, psb_spk_
      class(psb_c_ell_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(out)              :: d(:)
    end subroutine psb_c_ell_aclsum
  end interface
    
  interface 
    subroutine psb_c_ell_get_diag(a,d,info) 
      import :: psb_c_ell_sparse_mat, psb_spk_, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(out)     :: d(:)
      integer(psb_ipk_), intent(out)   :: info
    end subroutine psb_c_ell_get_diag
  end interface
  
  interface 
    subroutine psb_c_ell_scal(d,a,info,side) 
      import :: psb_c_ell_sparse_mat, psb_spk_, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)      :: d(:)
      integer(psb_ipk_), intent(out)   :: info
      character, intent(in), optional  :: side
    end subroutine psb_c_ell_scal
  end interface
  
  interface
    subroutine psb_c_ell_scals(d,a,info) 
      import :: psb_c_ell_sparse_mat, psb_spk_, psb_ipk_
      class(psb_c_ell_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)      :: d
      integer(psb_ipk_), intent(out)   :: info
    end subroutine psb_c_ell_scals
  end interface
  
  interface
    subroutine psi_c_convert_ell_from_coo(a,tmp,info,hacksize) 
      import :: psb_c_ell_sparse_mat, psb_c_coo_sparse_mat, psb_ipk_
      implicit none 
      class(psb_c_ell_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(in)    :: tmp
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), intent(in), optional    :: hacksize
    end subroutine psi_c_convert_ell_from_coo
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


  function c_ell_is_by_rows(a) result(res)
    implicit none 
    class(psb_c_ell_sparse_mat), intent(in) :: a
    logical :: res
    res = .true.
  end function c_ell_is_by_rows
  
  function c_ell_sizeof(a) result(res)
    implicit none 
    class(psb_c_ell_sparse_mat), intent(in) :: a
    integer(psb_epk_) :: res
    if (a%is_dev()) call a%sync()
    res = 8 
    res = res + (2*psb_sizeof_sp)  * size(a%val)
    res = res + psb_sizeof_ip * size(a%irn)
    res = res + psb_sizeof_ip * size(a%idiag)
    res = res + psb_sizeof_ip * size(a%ja)
      
  end function c_ell_sizeof

  function c_ell_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'ELL'
  end function c_ell_get_fmt
  
  function c_ell_get_nrm(a) result(res)
    implicit none 
    class(psb_c_ell_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res
    res = size(a%val,2)
  end function c_ell_get_nrm
  
  function c_ell_get_nzeros(a) result(res)
    implicit none 
    class(psb_c_ell_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res
    res = a%nzt
  end function c_ell_get_nzeros

  function c_ell_get_size(a) result(res)
    implicit none 
    class(psb_c_ell_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res

    res = -1
    if (a%is_dev()) call a%sync()
    
    if (allocated(a%ja)) then 
      if (res >= 0) then 
        res = min(res,size(a%ja))
      else 
        res = size(a%ja)
      end if
    end if
    if (allocated(a%val)) then 
      if (res >= 0) then 
        res = min(res,size(a%val))
      else 
        res = size(a%val)
      end if
    end if

  end function c_ell_get_size


  function  c_ell_get_nz_row(idx,a) result(res)

    implicit none
    
    class(psb_c_ell_sparse_mat), intent(in) :: a
    integer(psb_ipk_), intent(in)           :: idx
    integer(psb_ipk_)                       :: res
    
    res = 0 
    if (a%is_dev()) call a%sync()
 
    if ((1<=idx).and.(idx<=a%get_nrows())) then 
      res = a%irn(idx)
    end if
    
  end function c_ell_get_nz_row



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

  subroutine  c_ell_free(a) 
    implicit none 

    class(psb_c_ell_sparse_mat), intent(inout) :: a

    if (allocated(a%idiag)) deallocate(a%idiag)
    if (allocated(a%irn)) deallocate(a%irn)
    if (allocated(a%ja))  deallocate(a%ja)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(izero)
    call a%set_ncols(izero)
    
    return

  end subroutine c_ell_free


end module psb_c_ell_mat_mod
