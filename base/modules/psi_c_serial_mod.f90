module psi_c_serial_mod
  use psb_const_mod, only : psb_ipk_, psb_spk_

  interface psb_gelp
    ! 2-D version
    subroutine psb_cgelp(trans,iperm,x,info)
      import :: psb_ipk_, psb_spk_
      implicit none 
      complex(psb_spk_), intent(inout)     ::  x(:,:)
      integer(psb_ipk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_cgelp
    subroutine psb_cgelpv(trans,iperm,x,info)
      import :: psb_ipk_, psb_spk_
      implicit none 
      complex(psb_spk_), intent(inout)     ::  x(:)
      integer(psb_ipk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_cgelpv
  end interface psb_gelp
  
  interface psb_geaxpby 
    subroutine psi_caxpby(m,n,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_spk_
      implicit none 
      integer(psb_ipk_), intent(in)      :: m, n
      complex(psb_spk_), intent (in)       ::  x(:,:)
      complex(psb_spk_), intent (inout)    ::  y(:,:)
      complex(psb_spk_), intent (in)       ::  alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_caxpby
    subroutine psi_caxpbyv(m,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_spk_
      implicit none 
      integer(psb_ipk_), intent(in)      :: m
      complex(psb_spk_), intent (in)       ::  x(:)
      complex(psb_spk_), intent (inout)    ::  y(:)
      complex(psb_spk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_caxpbyv
  end interface psb_geaxpby

  interface psi_gth
    subroutine psi_cgthmv(n,k,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      complex(psb_spk_) :: x(:,:), y(:),alpha,beta
    end subroutine psi_cgthmv
    subroutine psi_cgthv(n,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, idx(:)
      complex(psb_spk_) :: x(:), y(:),alpha,beta
    end subroutine psi_cgthv
    subroutine psi_cgthzmv(n,k,idx,x,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      complex(psb_spk_) :: x(:,:), y(:)
      
    end subroutine psi_cgthzmv
    subroutine psi_cgthzv(n,idx,x,y)
      import :: psb_ipk_, psb_spk_
      implicit none 
      integer(psb_ipk_) :: n, idx(:)
      complex(psb_spk_) :: x(:), y(:)
    end subroutine psi_cgthzv
  end interface psi_gth

  interface psi_sct
    subroutine psi_csctmv(n,k,idx,x,beta,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      complex(psb_spk_) :: beta, x(:), y(:,:)
    end subroutine psi_csctmv
    subroutine psi_csctv(n,idx,x,beta,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      
      integer(psb_ipk_) :: n, idx(:)
      complex(psb_spk_) :: beta, x(:), y(:)
    end subroutine psi_csctv
  end interface psi_sct

end module psi_c_serial_mod
