module psi_d_serial_mod
  use psb_const_mod, only : psb_ipk_, psb_dpk_

  interface psb_gelp
    ! 2-D version
    subroutine psb_dgelp(trans,iperm,x,info)
      import :: psb_ipk_, psb_dpk_
      implicit none 
      real(psb_dpk_), intent(inout)     ::  x(:,:)
      integer(psb_ipk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_dgelp
    subroutine psb_dgelpv(trans,iperm,x,info)
      import :: psb_ipk_, psb_dpk_
      implicit none 
      real(psb_dpk_), intent(inout)     ::  x(:)
      integer(psb_ipk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_dgelpv
  end interface psb_gelp
  
  interface psb_geaxpby 
    subroutine psi_daxpby(m,n,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_dpk_
      implicit none 
      integer(psb_ipk_), intent(in)      :: m, n
      real(psb_dpk_), intent (in)       ::  x(:,:)
      real(psb_dpk_), intent (inout)    ::  y(:,:)
      real(psb_dpk_), intent (in)       ::  alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_daxpby
    subroutine psi_daxpbyv(m,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_dpk_
      implicit none 
      integer(psb_ipk_), intent(in)      :: m
      real(psb_dpk_), intent (in)       ::  x(:)
      real(psb_dpk_), intent (inout)    ::  y(:)
      real(psb_dpk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_daxpbyv
  end interface psb_geaxpby

  interface psi_gth
    subroutine psi_dgthmv(n,k,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      real(psb_dpk_) :: x(:,:), y(:),alpha,beta
    end subroutine psi_dgthmv
    subroutine psi_dgthv(n,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, idx(:)
      real(psb_dpk_) :: x(:), y(:),alpha,beta
    end subroutine psi_dgthv
    subroutine psi_dgthzmv(n,k,idx,x,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      real(psb_dpk_) :: x(:,:), y(:)
      
    end subroutine psi_dgthzmv
    subroutine psi_dgthzv(n,idx,x,y)
      import :: psb_ipk_, psb_dpk_
      implicit none 
      integer(psb_ipk_) :: n, idx(:)
      real(psb_dpk_) :: x(:), y(:)
    end subroutine psi_dgthzv
  end interface psi_gth

  interface psi_sct
    subroutine psi_dsctmv(n,k,idx,x,beta,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      real(psb_dpk_) :: beta, x(:), y(:,:)
    end subroutine psi_dsctmv
    subroutine psi_dsctv(n,idx,x,beta,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      
      integer(psb_ipk_) :: n, idx(:)
      real(psb_dpk_) :: beta, x(:), y(:)
    end subroutine psi_dsctv
  end interface psi_sct

end module psi_d_serial_mod
