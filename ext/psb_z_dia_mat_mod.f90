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
  

module psb_z_dia_mat_mod

  use psb_z_base_mat_mod

  type, extends(psb_z_base_sparse_mat) :: psb_z_dia_sparse_mat
    !
    ! DIA format, extended.
    !     
    
    integer(psb_ipk_), allocatable :: offset(:)
    integer(psb_ipk_) :: nzeros
    complex(psb_dpk_), allocatable :: data(:,:)

  contains
    ! procedure, pass(a) :: get_size     => z_dia_get_size
    procedure, pass(a) :: get_nzeros   => z_dia_get_nzeros
    procedure, nopass  :: get_fmt      => z_dia_get_fmt
    procedure, pass(a) :: sizeof       => z_dia_sizeof
    procedure, pass(a) :: csmm         => psb_z_dia_csmm
    procedure, pass(a) :: csmv         => psb_z_dia_csmv
    ! procedure, pass(a) :: inner_cssm   => psb_z_dia_cssm
    ! procedure, pass(a) :: inner_cssv   => psb_z_dia_cssv
    procedure, pass(a) :: scals        => psb_z_dia_scals
    procedure, pass(a) :: scalv        => psb_z_dia_scal
    procedure, pass(a) :: maxval       => psb_z_dia_maxval
    procedure, pass(a) :: rowsum       => psb_z_dia_rowsum
    procedure, pass(a) :: arwsum       => psb_z_dia_arwsum
    procedure, pass(a) :: colsum       => psb_z_dia_colsum
    procedure, pass(a) :: aclsum       => psb_z_dia_aclsum
    procedure, pass(a) :: reallocate_nz => psb_z_dia_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_z_dia_allocate_mnnz
    procedure, pass(a) :: cp_to_coo    => psb_z_cp_dia_to_coo
    procedure, pass(a) :: cp_from_coo  => psb_z_cp_dia_from_coo
    ! procedure, pass(a) :: mv_to_coo    => psb_z_mv_dia_to_coo
    procedure, pass(a) :: mv_from_coo  => psb_z_mv_dia_from_coo
    ! procedure, pass(a) :: mv_to_fmt    => psb_z_mv_dia_to_fmt
    ! procedure, pass(a) :: mv_from_fmt  => psb_z_mv_dia_from_fmt
    ! procedure, pass(a) :: csput_a      => psb_z_dia_csput_a
    procedure, pass(a) :: get_diag     => psb_z_dia_get_diag
    procedure, pass(a) :: csgetptn     => psb_z_dia_csgetptn
    procedure, pass(a) :: csgetrow     => psb_z_dia_csgetrow
    ! procedure, pass(a) :: get_nz_row   => z_dia_get_nz_row
    procedure, pass(a) :: reinit       => psb_z_dia_reinit
    ! procedure, pass(a) :: trim         => psb_z_dia_trim
    procedure, pass(a) :: print        => psb_z_dia_print
    procedure, pass(a) :: free         => z_dia_free
    procedure, pass(a) :: mold         => psb_z_dia_mold

  end type psb_z_dia_sparse_mat

  private :: z_dia_get_nzeros, z_dia_free,  z_dia_get_fmt, &
       & z_dia_sizeof !, z_dia_get_size, z_dia_get_nz_row

  interface
    subroutine  psb_z_dia_reallocate_nz(nz,a) 
      import :: psb_z_dia_sparse_mat, psb_ipk_
      integer(psb_ipk_), intent(in) :: nz
      class(psb_z_dia_sparse_mat), intent(inout) :: a
    end subroutine psb_z_dia_reallocate_nz
  end interface
  
  interface 
    subroutine psb_z_dia_reinit(a,clear)
      import :: psb_z_dia_sparse_mat
      class(psb_z_dia_sparse_mat), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_z_dia_reinit
  end interface
  
  interface
    subroutine  psb_z_dia_trim(a)
      import :: psb_z_dia_sparse_mat
      class(psb_z_dia_sparse_mat), intent(inout) :: a
    end subroutine psb_z_dia_trim
  end interface
  
  interface 
    subroutine psb_z_dia_mold(a,b,info) 
      import :: psb_z_dia_sparse_mat, psb_z_base_sparse_mat, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(in)                  :: a
      class(psb_z_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_z_dia_mold
  end interface

  interface
    subroutine  psb_z_dia_allocate_mnnz(m,n,a,nz) 
      import :: psb_z_dia_sparse_mat, psb_ipk_
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_z_dia_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in), optional :: nz
    end subroutine psb_z_dia_allocate_mnnz
  end interface
  
  interface
    subroutine psb_z_dia_print(iout,a,iv,head,ivr,ivc)
      import :: psb_z_dia_sparse_mat, psb_ipk_, psb_lpk_
      integer(psb_ipk_), intent(in)           :: iout
      class(psb_z_dia_sparse_mat), intent(in) :: a   
      integer(psb_lpk_), intent(in), optional :: iv(:)
      character(len=*), optional              :: head
      integer(psb_lpk_), intent(in), optional :: ivr(:), ivc(:)
    end subroutine psb_z_dia_print
  end interface
  
  interface 
    subroutine psb_z_cp_dia_to_coo(a,b,info) 
      import :: psb_z_coo_sparse_mat, psb_z_dia_sparse_mat, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(in)    :: a
      class(psb_z_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_z_cp_dia_to_coo
  end interface
  
  interface 
    subroutine psb_z_cp_dia_from_coo(a,b,info) 
      import :: psb_z_dia_sparse_mat, psb_z_coo_sparse_mat, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_z_cp_dia_from_coo
  end interface
  
  interface 
    subroutine psb_z_cp_dia_to_fmt(a,b,info) 
      import :: psb_z_dia_sparse_mat, psb_z_base_sparse_mat, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(in)     :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_z_cp_dia_to_fmt
  end interface
  
  interface 
    subroutine psb_z_cp_dia_from_fmt(a,b,info) 
      import :: psb_z_dia_sparse_mat, psb_z_base_sparse_mat, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(inout) :: a
      class(psb_z_base_sparse_mat), intent(in)   :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_z_cp_dia_from_fmt
  end interface
  
  interface 
    subroutine psb_z_mv_dia_to_coo(a,b,info) 
      import :: psb_z_dia_sparse_mat, psb_z_coo_sparse_mat, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_z_mv_dia_to_coo
  end interface
  
  interface 
    subroutine psb_z_mv_dia_from_coo(a,b,info) 
      import :: psb_z_dia_sparse_mat, psb_z_coo_sparse_mat, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_z_mv_dia_from_coo
  end interface
  
  interface 
    subroutine psb_z_mv_dia_to_fmt(a,b,info) 
      import :: psb_z_dia_sparse_mat, psb_z_base_sparse_mat, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(inout)  :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_z_mv_dia_to_fmt
  end interface
  
  interface 
    subroutine psb_z_mv_dia_from_fmt(a,b,info) 
      import :: psb_z_dia_sparse_mat, psb_z_base_sparse_mat, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(inout)  :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_z_mv_dia_from_fmt
  end interface
  
  interface 
    subroutine psb_z_dia_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info) 
      import :: psb_z_dia_sparse_mat, psb_dpk_, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: val(:)
      integer(psb_ipk_), intent(in)             :: nz,ia(:), ja(:),&
           &  imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_dia_csput_a
  end interface
  
  interface 
    subroutine psb_z_dia_csgetptn(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import :: psb_z_dia_sparse_mat, psb_dpk_, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(in)        :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_z_dia_csgetptn
  end interface
  
  interface 
    subroutine psb_z_dia_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale,chksz)
      import :: psb_z_dia_sparse_mat, psb_dpk_, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(in)        :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale,chksz
    end subroutine psb_z_dia_csgetrow
  end interface

  interface 
    subroutine psb_z_dia_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
      import :: psb_z_dia_sparse_mat, psb_dpk_, psb_z_coo_sparse_mat, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(in)    :: a
      class(psb_z_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(in)              :: imin,imax
      integer(psb_ipk_),intent(out)              :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional    :: iren(:)
      integer(psb_ipk_), intent(in), optional    :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_z_dia_csgetblk
  end interface
    
  interface 
    subroutine psb_z_dia_cssv(alpha,a,x,beta,y,info,trans) 
      import :: psb_z_dia_sparse_mat, psb_dpk_, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      complex(psb_dpk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)           :: info
      character, optional, intent(in)          :: trans
    end subroutine psb_z_dia_cssv
    subroutine psb_z_dia_cssm(alpha,a,x,beta,y,info,trans) 
      import :: psb_z_dia_sparse_mat, psb_dpk_, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      complex(psb_dpk_), intent(inout)       :: y(:,:)
      integer(psb_ipk_), intent(out)       :: info
      character, optional, intent(in)      :: trans
    end subroutine psb_z_dia_cssm
  end interface
  
  interface 
    subroutine psb_z_dia_csmv(alpha,a,x,beta,y,info,trans) 
      import :: psb_z_dia_sparse_mat, psb_dpk_, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      complex(psb_dpk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)       :: info
      character, optional, intent(in)      :: trans
    end subroutine psb_z_dia_csmv
    subroutine psb_z_dia_csmm(alpha,a,x,beta,y,info,trans) 
      import :: psb_z_dia_sparse_mat, psb_dpk_, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      complex(psb_dpk_), intent(inout)       :: y(:,:)
      integer(psb_ipk_), intent(out)       :: info
      character, optional, intent(in)      :: trans
    end subroutine psb_z_dia_csmm
  end interface
  
  
  interface 
    function psb_z_dia_maxval(a) result(res)
      import :: psb_z_dia_sparse_mat, psb_dpk_
      class(psb_z_dia_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_z_dia_maxval
  end interface
  
  interface 
    function psb_z_dia_csnmi(a) result(res)
      import :: psb_z_dia_sparse_mat, psb_dpk_
      class(psb_z_dia_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_z_dia_csnmi
  end interface
  
  interface 
    function psb_z_dia_csnm1(a) result(res)
      import :: psb_z_dia_sparse_mat, psb_dpk_
      class(psb_z_dia_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_z_dia_csnm1
  end interface

  interface 
    subroutine psb_z_dia_rowsum(d,a) 
      import :: psb_z_dia_sparse_mat, psb_dpk_
      class(psb_z_dia_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_dia_rowsum
  end interface

  interface 
    subroutine psb_z_dia_arwsum(d,a) 
      import :: psb_z_dia_sparse_mat, psb_dpk_
      class(psb_z_dia_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_dia_arwsum
  end interface
  
  interface 
    subroutine psb_z_dia_colsum(d,a) 
      import :: psb_z_dia_sparse_mat, psb_dpk_
      class(psb_z_dia_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_dia_colsum
  end interface

  interface 
    subroutine psb_z_dia_aclsum(d,a) 
      import :: psb_z_dia_sparse_mat, psb_dpk_
      class(psb_z_dia_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_dia_aclsum
  end interface
    
  interface 
    subroutine psb_z_dia_get_diag(a,d,info) 
      import :: psb_z_dia_sparse_mat, psb_dpk_, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)     :: d(:)
      integer(psb_ipk_), intent(out)   :: info
    end subroutine psb_z_dia_get_diag
  end interface
  
  interface 
    subroutine psb_z_dia_scal(d,a,info,side) 
      import :: psb_z_dia_sparse_mat, psb_dpk_, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: d(:)
      integer(psb_ipk_), intent(out)   :: info
      character, intent(in), optional  :: side
    end subroutine psb_z_dia_scal
  end interface
  
  interface
    subroutine psb_z_dia_scals(d,a,info) 
      import :: psb_z_dia_sparse_mat, psb_dpk_, psb_ipk_
      class(psb_z_dia_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: d
      integer(psb_ipk_), intent(out)   :: info
    end subroutine psb_z_dia_scals
  end interface
  
  interface  psi_convert_dia_from_coo
    subroutine psi_z_convert_dia_from_coo(a,tmp,info)
      import :: psb_z_dia_sparse_mat, psb_ipk_, psb_z_coo_sparse_mat
      implicit none 
      class(psb_z_dia_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(in)    :: tmp
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psi_z_convert_dia_from_coo
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

  
   function z_dia_sizeof(a) result(res)
     implicit none 
     class(psb_z_dia_sparse_mat), intent(in) :: a
     integer(psb_epk_) :: res
     if (a%is_dev()) call a%sync()
     res = 8 
     res = res + (2*psb_sizeof_dp)  * size(a%data)
     res = res + psb_sizeof_ip * size(a%offset)
     
   end function z_dia_sizeof

  function z_dia_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'DIA'
  end function z_dia_get_fmt
  
  function z_dia_get_nzeros(a) result(res)
    implicit none 
    class(psb_z_dia_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res
    res = a%nzeros
  end function z_dia_get_nzeros

  ! function z_dia_get_size(a) result(res)
  !   implicit none 
  !   class(psb_z_dia_sparse_mat), intent(in) :: a
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

  ! end function z_dia_get_size


  ! function  z_dia_get_nz_row(idx,a) result(res)

  !   implicit none
    
  !   class(psb_z_dia_sparse_mat), intent(in) :: a
  !   integer(psb_ipk_), intent(in)           :: idx
  !   integer(psb_ipk_)                       :: res
    
  !   res = 0 
 
  !   if ((1<=idx).and.(idx<=a%get_nrows())) then 
  !     res = a%irn(idx)
  !   end if
    
  ! end function z_dia_get_nz_row



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

  subroutine  z_dia_free(a) 
    implicit none 

    class(psb_z_dia_sparse_mat), intent(inout) :: a

    if (allocated(a%data)) deallocate(a%data)
    if (allocated(a%offset)) deallocate(a%offset)
    call a%set_null()
    call a%set_nrows(izero)
    call a%set_ncols(izero)
    
    return

  end subroutine z_dia_free


end module psb_z_dia_mat_mod
