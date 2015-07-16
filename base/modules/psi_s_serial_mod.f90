module psi_s_serial_mod
  use psb_const_mod, only : psb_ipk_, psb_spk_

  interface psb_gelp
    ! 2-D version
    subroutine psb_sgelp(trans,iperm,x,info)
      import :: psb_ipk_, psb_spk_
      implicit none 
      real(psb_spk_), intent(inout)     ::  x(:,:)
      integer(psb_ipk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_sgelp
    subroutine psb_sgelpv(trans,iperm,x,info)
      import :: psb_ipk_, psb_spk_
      implicit none 
      real(psb_spk_), intent(inout)     ::  x(:)
      integer(psb_ipk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_sgelpv
  end interface psb_gelp
  
  interface psb_geaxpby 
    subroutine psi_saxpby(m,n,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_spk_
      implicit none 
      integer(psb_ipk_), intent(in)      :: m, n
      real(psb_spk_), intent (in)       ::  x(:,:)
      real(psb_spk_), intent (inout)    ::  y(:,:)
      real(psb_spk_), intent (in)       ::  alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_saxpby
    subroutine psi_saxpbyv(m,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_spk_
      implicit none 
      integer(psb_ipk_), intent(in)      :: m
      real(psb_spk_), intent (in)       ::  x(:)
      real(psb_spk_), intent (inout)    ::  y(:)
      real(psb_spk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_saxpbyv
  end interface psb_geaxpby

  interface psi_gth
    subroutine psi_sgthmv(n,k,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      real(psb_spk_) :: x(:,:), y(:),alpha,beta
    end subroutine psi_sgthmv
    subroutine psi_sgthv(n,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, idx(:)
      real(psb_spk_) :: x(:), y(:),alpha,beta
    end subroutine psi_sgthv
    subroutine psi_sgthzmv(n,k,idx,x,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      real(psb_spk_) :: x(:,:), y(:)
      
    end subroutine psi_sgthzmv
    subroutine psi_sgthzmm(n,k,idx,x,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      real(psb_spk_) :: x(:,:), y(:,:)
      
    end subroutine psi_sgthzmm
    subroutine psi_sgthzv(n,idx,x,y)
      import :: psb_ipk_, psb_spk_
      implicit none 
      integer(psb_ipk_) :: n, idx(:)
      real(psb_spk_) :: x(:), y(:)
    end subroutine psi_sgthzv
  end interface psi_gth

  interface psi_sct
    subroutine psi_ssctmm(n,k,idx,x,beta,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      real(psb_spk_) :: beta, x(:,:), y(:,:)
    end subroutine psi_ssctmm
    subroutine psi_ssctmv(n,k,idx,x,beta,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      real(psb_spk_) :: beta, x(:), y(:,:)
    end subroutine psi_ssctmv
    subroutine psi_ssctv(n,idx,x,beta,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      
      integer(psb_ipk_) :: n, idx(:)
      real(psb_spk_) :: beta, x(:), y(:)
    end subroutine psi_ssctv
  end interface psi_sct

end module psi_s_serial_mod
