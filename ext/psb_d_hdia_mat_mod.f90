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
  
module psb_d_hdia_mat_mod

  use psb_d_base_mat_mod


  type, extends(psb_d_base_sparse_mat) :: psb_d_hdia_sparse_mat
    !
    ! HDIA format
    !
    integer(psb_ipk_), allocatable :: hackOffsets(:), diaOffsets(:)
    real(psb_dpk_), allocatable   :: val(:)
    
    
    integer(psb_ipk_) :: nhacks, nzeros
    integer(psb_ipk_) :: hacksize = 32
    integer(psb_epk_) :: dim=0

  contains
    ! procedure, pass(a) :: get_size     => d_hdia_get_size
    procedure, pass(a) :: get_nzeros   => d_hdia_get_nzeros
    procedure, pass(a) :: set_nzeros   => d_hdia_set_nzeros
    procedure, nopass  :: get_fmt      => d_hdia_get_fmt
    procedure, pass(a) :: sizeof       => d_hdia_sizeof
    ! procedure, pass(a) :: csmm         => psb_d_hdia_csmm
    procedure, pass(a) :: csmv         => psb_d_hdia_csmv
    ! procedure, pass(a) :: inner_cssm   => psb_d_hdia_cssm
    ! procedure, pass(a) :: inner_cssv   => psb_d_hdia_cssv
    ! procedure, pass(a) :: scals        => psb_d_hdia_scals
    ! procedure, pass(a) :: scalv        => psb_d_hdia_scal
    ! procedure, pass(a) :: maxval       => psb_d_hdia_maxval
    ! procedure, pass(a) :: csnmi        => psb_d_hdia_csnmi
    ! procedure, pass(a) :: csnm1        => psb_d_hdia_csnm1
    ! procedure, pass(a) :: rowsum       => psb_d_hdia_rowsum
    ! procedure, pass(a) :: arwsum       => psb_d_hdia_arwsum
    ! procedure, pass(a) :: colsum       => psb_d_hdia_colsum
    ! procedure, pass(a) :: aclsum       => psb_d_hdia_aclsum
    ! procedure, pass(a) :: reallocate_nz => psb_d_hdia_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_d_hdia_allocate_mnnz
    procedure, pass(a) :: cp_to_coo    => psb_d_cp_hdia_to_coo
    procedure, pass(a) :: cp_from_coo  => psb_d_cp_hdia_from_coo
    ! procedure, pass(a) :: cp_to_fmt    => psb_d_cp_hdia_to_fmt
    !   procedure, pass(a) :: cp_from_fmt  => psb_d_cp_hdia_from_fmt
    procedure, pass(a) :: mv_to_coo    => psb_d_mv_hdia_to_coo
    procedure, pass(a) :: mv_from_coo  => psb_d_mv_hdia_from_coo
    ! procedure, pass(a) :: mv_to_fmt    => psb_d_mv_hdia_to_fmt
    !    procedure, pass(a) :: mv_from_fmt  => psb_d_mv_hdia_from_fmt
    ! procedure, pass(a) :: csput_a        => psb_d_hdia_csput_a
    ! procedure, pass(a) :: get_diag     => psb_d_hdia_get_diag
    ! procedure, pass(a) :: csgetptn     => psb_d_hdia_csgetptn
    ! procedure, pass(a) :: csgetrow     => psb_d_hdia_csgetrow
    ! procedure, pass(a) :: get_nz_row   => d_hdia_get_nz_row
    ! procedure, pass(a) :: reinit       => psb_d_hdia_reinit
    ! procedure, pass(a) :: trim         => psb_d_hdia_trim
    procedure, pass(a) :: print        => psb_d_hdia_print
    procedure, pass(a) :: free         => d_hdia_free
    procedure, pass(a) :: mold         => psb_d_hdia_mold

  end type psb_d_hdia_sparse_mat

  private :: d_hdia_get_nzeros, d_hdia_set_nzeros, d_hdia_free, &
       & d_hdia_get_fmt, d_hdia_sizeof 
!!$  &
!!$       & d_hdia_get_nz_row d_hdia_get_size, 

!!$  interface
!!$    subroutine  psb_d_hdia_reallocate_nz(nz,a) 
!!$      import :: psb_d_hdia_sparse_mat, psb_ipk_
!!$      integer(psb_ipk_), intent(in) :: nz
!!$      class(psb_d_hdia_sparse_mat), intent(inout) :: a
!!$    end subroutine psb_d_hdia_reallocate_nz
!!$  end interface
!!$  
!!$  interface 
!!$    subroutine psb_d_hdia_reinit(a,clear)
!!$      import :: psb_d_hdia_sparse_mat
!!$      class(psb_d_hdia_sparse_mat), intent(inout) :: a   
!!$      logical, intent(in), optional :: clear
!!$    end subroutine psb_d_hdia_reinit
!!$  end interface
!!$  
!!$  interface
!!$    subroutine  psb_d_hdia_trim(a)
!!$      import :: psb_d_hdia_sparse_mat
!!$      class(psb_d_hdia_sparse_mat), intent(inout) :: a
!!$    end subroutine psb_d_hdia_trim
!!$  end interface
  
  interface 
    subroutine psb_d_hdia_mold(a,b,info) 
      import :: psb_d_hdia_sparse_mat, psb_d_base_sparse_mat, psb_ipk_
      class(psb_d_hdia_sparse_mat), intent(in)                  :: a
      class(psb_d_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_d_hdia_mold
  end interface

  interface
    subroutine  psb_d_hdia_allocate_mnnz(m,n,a,nz) 
      import :: psb_d_hdia_sparse_mat, psb_ipk_
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_d_hdia_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in), optional :: nz
    end subroutine psb_d_hdia_allocate_mnnz
  end interface
  
  interface
    subroutine psb_d_hdia_print(iout,a,iv,head,ivr,ivc)
      import :: psb_d_hdia_sparse_mat, psb_ipk_, psb_lpk_
      integer(psb_ipk_), intent(in)           :: iout
      class(psb_d_hdia_sparse_mat), intent(in) :: a   
      integer(psb_lpk_), intent(in), optional :: iv(:)
      character(len=*), optional              :: head
      integer(psb_lpk_), intent(in), optional :: ivr(:), ivc(:)
    end subroutine psb_d_hdia_print
  end interface
  
  interface 
    subroutine psb_d_cp_hdia_to_coo(a,b,info) 
      import :: psb_d_coo_sparse_mat, psb_d_hdia_sparse_mat, psb_ipk_
      class(psb_d_hdia_sparse_mat), intent(in)    :: a
      class(psb_d_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_d_cp_hdia_to_coo
  end interface
  
  interface 
    subroutine psb_d_cp_hdia_from_coo(a,b,info) 
      import :: psb_d_hdia_sparse_mat, psb_d_coo_sparse_mat, psb_ipk_
      class(psb_d_hdia_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_d_cp_hdia_from_coo
  end interface
  
!!$  interface 
!!$    subroutine psb_d_cp_hdia_to_fmt(a,b,info) 
!!$      import :: psb_d_hdia_sparse_mat, psb_d_base_sparse_mat, psb_ipk_
!!$      class(psb_d_hdia_sparse_mat), intent(in)     :: a
!!$      class(psb_d_base_sparse_mat), intent(inout) :: b
!!$      integer(psb_ipk_), intent(out)              :: info
!!$    end subroutine psb_d_cp_hdia_to_fmt
!!$  end interface
!!$  
!!$  interface 
!!$    subroutine psb_d_cp_hdia_from_fmt(a,b,info) 
!!$      import :: psb_d_hdia_sparse_mat, psb_d_base_sparse_mat, psb_ipk_
!!$      class(psb_d_hdia_sparse_mat), intent(inout) :: a
!!$      class(psb_d_base_sparse_mat), intent(in)   :: b
!!$      integer(psb_ipk_), intent(out)             :: info
!!$    end subroutine psb_d_cp_hdia_from_fmt
!!$  end interface
  
  interface 
    subroutine psb_d_mv_hdia_to_coo(a,b,info) 
      import :: psb_d_hdia_sparse_mat, psb_d_coo_sparse_mat, psb_ipk_
      class(psb_d_hdia_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_d_mv_hdia_to_coo
  end interface
  
  interface 
    subroutine psb_d_mv_hdia_from_coo(a,b,info) 
      import :: psb_d_hdia_sparse_mat, psb_d_coo_sparse_mat, psb_ipk_
      class(psb_d_hdia_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_d_mv_hdia_from_coo
  end interface
  
!!$  interface 
!!$    subroutine psb_d_mv_hdia_to_fmt(a,b,info) 
!!$      import :: psb_d_hdia_sparse_mat, psb_d_base_sparse_mat, psb_ipk_
!!$      class(psb_d_hdia_sparse_mat), intent(inout)  :: a
!!$      class(psb_d_base_sparse_mat), intent(inout) :: b
!!$      integer(psb_ipk_), intent(out)              :: info
!!$    end subroutine psb_d_mv_hdia_to_fmt
!!$  end interface
!!$  
!!$  interface 
!!$    subroutine psb_d_mv_hdia_from_fmt(a,b,info) 
!!$      import :: psb_d_hdia_sparse_mat, psb_d_base_sparse_mat, psb_ipk_
!!$      class(psb_d_hdia_sparse_mat), intent(inout)  :: a
!!$      class(psb_d_base_sparse_mat), intent(inout) :: b
!!$      integer(psb_ipk_), intent(out)              :: info
!!$    end subroutine psb_d_mv_hdia_from_fmt
!!$  end interface
!!$  
!!$  interface 
!!$    subroutine psb_d_hdia_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info) 
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_, psb_ipk_
!!$      class(psb_d_hdia_sparse_mat), intent(inout) :: a
!!$      real(psb_dpk_), intent(in)      :: val(:)
!!$      integer(psb_ipk_), intent(in)             :: nz,ia(:), ja(:),&
!!$           &  imin,imax,jmin,jmax
!!$      integer(psb_ipk_), intent(out)            :: info
!!$    end subroutine psb_d_hdia_csput_a
!!$  end interface
!!$  
!!$  interface 
!!$    subroutine psb_d_hdia_csgetptn(imin,imax,a,nz,ia,ja,info,&
!!$         & jmin,jmax,iren,append,nzin,rscale,cscale)
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_, psb_ipk_
!!$      class(psb_d_hdia_sparse_mat), intent(in)        :: a
!!$      integer(psb_ipk_), intent(in)                  :: imin,imax
!!$      integer(psb_ipk_), intent(out)                 :: nz
!!$      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
!!$      integer(psb_ipk_),intent(out)                  :: info
!!$      logical, intent(in), optional        :: append
!!$      integer(psb_ipk_), intent(in), optional        :: iren(:)
!!$      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
!!$      logical, intent(in), optional        :: rscale,cscale
!!$    end subroutine psb_d_hdia_csgetptn
!!$  end interface
!!$  
!!$  interface 
!!$    subroutine psb_d_hdia_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
!!$         & jmin,jmax,iren,append,nzin,rscale,cscale)
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_, psb_ipk_
!!$      class(psb_d_hdia_sparse_mat), intent(in)        :: a
!!$      integer(psb_ipk_), intent(in)                  :: imin,imax
!!$      integer(psb_ipk_), intent(out)                 :: nz
!!$      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
!!$      real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
!!$      integer(psb_ipk_),intent(out)                  :: info
!!$      logical, intent(in), optional        :: append
!!$      integer(psb_ipk_), intent(in), optional        :: iren(:)
!!$      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
!!$      logical, intent(in), optional        :: rscale,cscale
!!$    end subroutine psb_d_hdia_csgetrow
!!$  end interface
!!$
!!$  interface 
!!$    subroutine psb_d_hdia_csgetblk(imin,imax,a,b,info,&
!!$       & jmin,jmax,iren,append,rscale,cscale)
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_, psb_d_coo_sparse_mat, psb_ipk_
!!$      class(psb_d_hdia_sparse_mat), intent(in)    :: a
!!$      class(psb_d_coo_sparse_mat), intent(inout) :: b
!!$      integer(psb_ipk_), intent(in)              :: imin,imax
!!$      integer(psb_ipk_),intent(out)              :: info
!!$      logical, intent(in), optional        :: append
!!$      integer(psb_ipk_), intent(in), optional    :: iren(:)
!!$      integer(psb_ipk_), intent(in), optional    :: jmin,jmax
!!$      logical, intent(in), optional        :: rscale,cscale
!!$    end subroutine psb_d_hdia_csgetblk
!!$  end interface
!!$    
!!$  interface 
!!$    subroutine psb_d_hdia_cssv(alpha,a,x,beta,y,info,trans) 
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_, psb_ipk_
!!$      class(psb_d_hdia_sparse_mat), intent(in) :: a
!!$      real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
!!$      real(psb_dpk_), intent(inout)       :: y(:)
!!$      integer(psb_ipk_), intent(out)           :: info
!!$      character, optional, intent(in)          :: trans
!!$    end subroutine psb_d_hdia_cssv
!!$    subroutine psb_d_hdia_cssm(alpha,a,x,beta,y,info,trans) 
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_, psb_ipk_
!!$      class(psb_d_hdia_sparse_mat), intent(in) :: a
!!$      real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
!!$      real(psb_dpk_), intent(inout)       :: y(:,:)
!!$      integer(psb_ipk_), intent(out)       :: info
!!$      character, optional, intent(in)      :: trans
!!$    end subroutine psb_d_hdia_cssm
!!$  end interface
  
  interface 
    subroutine psb_d_hdia_csmv(alpha,a,x,beta,y,info,trans,ivshft) 
      import :: psb_d_hdia_sparse_mat, psb_dpk_, psb_ipk_
      class(psb_d_hdia_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      real(psb_dpk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)       :: info
      character, optional, intent(in)      :: trans
      integer(psb_ipk_), optional, intent(in) :: ivshft
    end subroutine psb_d_hdia_csmv
!!$    subroutine psb_d_hdia_csmm(alpha,a,x,beta,y,info,trans) 
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_, psb_ipk_
!!$      class(psb_d_hdia_sparse_mat), intent(in) :: a
!!$      real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
!!$      real(psb_dpk_), intent(inout)       :: y(:,:)
!!$      integer(psb_ipk_), intent(out)       :: info
!!$      character, optional, intent(in)      :: trans
!!$    end subroutine psb_d_hdia_csmm
  end interface
  
  
!!$  interface 
!!$    function psb_d_hdia_maxval(a) result(res)
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_
!!$      class(psb_d_hdia_sparse_mat), intent(in) :: a
!!$      real(psb_dpk_)         :: res
!!$    end function psb_d_hdia_maxval
!!$  end interface
!!$  
!!$  interface 
!!$    function psb_d_hdia_csnmi(a) result(res)
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_
!!$      class(psb_d_hdia_sparse_mat), intent(in) :: a
!!$      real(psb_dpk_)         :: res
!!$    end function psb_d_hdia_csnmi
!!$  end interface
!!$  
!!$  interface 
!!$    function psb_d_hdia_csnm1(a) result(res)
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_
!!$      class(psb_d_hdia_sparse_mat), intent(in) :: a
!!$      real(psb_dpk_)         :: res
!!$    end function psb_d_hdia_csnm1
!!$  end interface
!!$
!!$  interface 
!!$    subroutine psb_d_hdia_rowsum(d,a) 
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_
!!$      class(psb_d_hdia_sparse_mat), intent(in) :: a
!!$      real(psb_dpk_), intent(out)              :: d(:)
!!$    end subroutine psb_d_hdia_rowsum
!!$  end interface
!!$
!!$  interface 
!!$    subroutine psb_d_hdia_arwsum(d,a) 
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_
!!$      class(psb_d_hdia_sparse_mat), intent(in) :: a
!!$      real(psb_dpk_), intent(out)              :: d(:)
!!$    end subroutine psb_d_hdia_arwsum
!!$  end interface
!!$  
!!$  interface 
!!$    subroutine psb_d_hdia_colsum(d,a) 
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_
!!$      class(psb_d_hdia_sparse_mat), intent(in) :: a
!!$      real(psb_dpk_), intent(out)              :: d(:)
!!$    end subroutine psb_d_hdia_colsum
!!$  end interface
!!$
!!$  interface 
!!$    subroutine psb_d_hdia_aclsum(d,a) 
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_
!!$      class(psb_d_hdia_sparse_mat), intent(in) :: a
!!$      real(psb_dpk_), intent(out)              :: d(:)
!!$    end subroutine psb_d_hdia_aclsum
!!$  end interface
!!$    
!!$  interface 
!!$    subroutine psb_d_hdia_get_diag(a,d,info) 
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_, psb_ipk_
!!$      class(psb_d_hdia_sparse_mat), intent(in) :: a
!!$      real(psb_dpk_), intent(out)     :: d(:)
!!$      integer(psb_ipk_), intent(out)   :: info
!!$    end subroutine psb_d_hdia_get_diag
!!$  end interface
!!$  
!!$  interface 
!!$    subroutine psb_d_hdia_scal(d,a,info,side) 
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_, psb_ipk_
!!$      class(psb_d_hdia_sparse_mat), intent(inout) :: a
!!$      real(psb_dpk_), intent(in)      :: d(:)
!!$      integer(psb_ipk_), intent(out)   :: info
!!$      character, intent(in), optional  :: side
!!$    end subroutine psb_d_hdia_scal
!!$  end interface
  
!!$  interface
!!$    subroutine psb_d_hdia_scals(d,a,info) 
!!$      import :: psb_d_hdia_sparse_mat, psb_dpk_, psb_ipk_
!!$      class(psb_d_hdia_sparse_mat), intent(inout) :: a
!!$      real(psb_dpk_), intent(in)      :: d
!!$      integer(psb_ipk_), intent(out)   :: info
!!$    end subroutine psb_d_hdia_scals
!!$  end interface
!!$  


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

  
   function d_hdia_sizeof(a) result(res)
     use psb_realloc_mod, only : psb_size
     implicit none 
     class(psb_d_hdia_sparse_mat), intent(in) :: a
     integer(psb_epk_) :: res
     integer(psb_ipk_) :: i

     if (a%is_dev()) call a%sync()
     res = 0
     
     res = res + psb_size(a%hackOffsets)*psb_sizeof_ip
     res = res + psb_size(a%diaOffsets)*psb_sizeof_ip
     res = res + psb_size(a%val) * psb_sizeof_dp 

   end function d_hdia_sizeof

  function d_hdia_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'HDIA'
  end function d_hdia_get_fmt
  
  function d_hdia_get_nzeros(a) result(res)
    implicit none 
    class(psb_d_hdia_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res
    res = a%nzeros
  end function d_hdia_get_nzeros

  subroutine  d_hdia_set_nzeros(a,nz) 
    implicit none 
    class(psb_d_hdia_sparse_mat), intent(inout) :: a
    integer(psb_ipk_), intent(in) :: nz
    a%nzeros = nz
  end subroutine d_hdia_set_nzeros
  
  ! function d_hdia_get_size(a) result(res)
  !   implicit none 
  !   class(psb_d_hdia_sparse_mat), intent(in) :: a
  !   integer(psb_ipk_) :: res

  !   res = -1
    
  !   if (allocated(a%ja)) then 
  !     if (res >= 0) then 
  !       res = min(res,size(a%ja))
  !     else 
  !       res = size(a%ja)
  !     end if
  !   end if
  !   if (allocated(a%val)) then 
  !     if (res >= 0) then 
  !       res = min(res,size(a%val))
  !     else 
  !       res = size(a%val)
  !     end if
  !   end if

  ! end function d_hdia_get_size


  ! function  d_hdia_get_nz_row(idx,a) result(res)

  !   implicit none
    
  !   class(psb_d_hdia_sparse_mat), intent(in) :: a
  !   integer(psb_ipk_), intent(in)           :: idx
  !   integer(psb_ipk_)                       :: res
    
  !   res = 0 
 
  !   if ((1<=idx).and.(idx<=a%get_nrows())) then 
  !     res = a%irn(idx)
  !   end if
    
  ! end function d_hdia_get_nz_row



  ! ! == ===================================
  ! !
  ! !
  ! !
  ! ! Data management
  ! !
  ! !
  ! !
  ! !
  ! !
  ! ! == ===================================  

  subroutine d_hdia_free(a) 
    implicit none 

    class(psb_d_hdia_sparse_mat), intent(inout) :: a
    integer(psb_ipk_) :: i, info

    
    if (allocated(a%hackOffsets))&
         &  deallocate(a%hackOffsets,stat=info)
    if (allocated(a%diaOffsets))&
         &  deallocate(a%diaOffsets,stat=info)
    if (allocated(a%val))&
         &  deallocate(a%val,stat=info)
    a%nhacks=0

    call a%set_null()
    call a%set_nrows(izero)
    call a%set_ncols(izero)
    
    return

  end subroutine d_hdia_free


end module psb_d_hdia_mat_mod
