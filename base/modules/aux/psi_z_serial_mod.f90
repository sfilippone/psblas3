module psi_z_serial_mod
  use psb_const_mod, only : psb_ipk_, psb_dpk_

  interface psb_gelp
    ! 2-D version
    subroutine psb_zgelp(trans,iperm,x,info)
      import :: psb_ipk_, psb_dpk_
      implicit none 
      complex(psb_dpk_), intent(inout)     ::  x(:,:)
      integer(psb_ipk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_zgelp
    subroutine psb_zgelpv(trans,iperm,x,info)
      import :: psb_ipk_, psb_dpk_
      implicit none 
      complex(psb_dpk_), intent(inout)     ::  x(:)
      integer(psb_ipk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_zgelpv
  end interface psb_gelp
  
  interface psb_geaxpby 
    subroutine psi_zaxpby(m,n,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_dpk_
      implicit none 
      integer(psb_ipk_), intent(in)      :: m, n
      complex(psb_dpk_), intent (in)       ::  x(:,:)
      complex(psb_dpk_), intent (inout)    ::  y(:,:)
      complex(psb_dpk_), intent (in)       ::  alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_zaxpby
    subroutine psi_zaxpbyv(m,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_dpk_
      implicit none 
      integer(psb_ipk_), intent(in)      :: m
      complex(psb_dpk_), intent (in)       ::  x(:)
      complex(psb_dpk_), intent (inout)    ::  y(:)
      complex(psb_dpk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_zaxpbyv
  end interface psb_geaxpby

  interface psi_gth
    subroutine psi_zgthmv(n,k,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      complex(psb_dpk_) :: x(:,:), y(:),alpha,beta
    end subroutine psi_zgthmv
    subroutine psi_zgthv(n,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, idx(:)
      complex(psb_dpk_) :: x(:), y(:),alpha,beta
    end subroutine psi_zgthv
    subroutine psi_zgthzmv(n,k,idx,x,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      complex(psb_dpk_) :: x(:,:), y(:)
      
    end subroutine psi_zgthzmv
    subroutine psi_zgthzmm(n,k,idx,x,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      complex(psb_dpk_) :: x(:,:), y(:,:)
      
    end subroutine psi_zgthzmm
    subroutine psi_zgthzv(n,idx,x,y)
      import :: psb_ipk_, psb_dpk_
      implicit none 
      integer(psb_ipk_) :: n, idx(:)
      complex(psb_dpk_) :: x(:), y(:)
    end subroutine psi_zgthzv
  end interface psi_gth

  interface psi_sct
    subroutine psi_zsctmm(n,k,idx,x,beta,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      complex(psb_dpk_) :: beta, x(:,:), y(:,:)
    end subroutine psi_zsctmm
    subroutine psi_zsctmv(n,k,idx,x,beta,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      complex(psb_dpk_) :: beta, x(:), y(:,:)
    end subroutine psi_zsctmv
    subroutine psi_zsctv(n,idx,x,beta,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      
      integer(psb_ipk_) :: n, idx(:)
      complex(psb_dpk_) :: beta, x(:), y(:)
    end subroutine psi_zsctv
  end interface psi_sct

end module psi_z_serial_mod
