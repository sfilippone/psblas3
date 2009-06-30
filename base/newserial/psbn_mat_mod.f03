
module psbn_d_mat_mod

  use psbn_d_base_mat_mod
  
  type :: psbn_d_sparse_mat

    class(psbn_d_base_sparse_mat), allocatable  :: a 
    
  contains
    
    procedure, pass(a) :: d_csmv
    procedure, pass(a) :: d_csmm
    generic, public    :: psbn_csmm => d_csmm, d_csmv

    procedure, pass(a) :: d_cssv
    procedure, pass(a) :: d_cssm
    generic, public    :: psbn_cssm => d_cssm, d_cssv
    
  end type psbn_d_sparse_mat

contains 

  subroutine d_csmm(alpha,a,x,beta,y,info,trans) 
    class(psbn_d_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:,:)
    real(kind(1.d0)), intent(inout) :: y(:,:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans

    call a%a%psbn_csmm(alpha,x,beta,y,info,trans) 

  end subroutine d_csmm

  subroutine d_csmv(alpha,a,x,beta,y,info,trans) 
    class(psbn_d_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:)
    real(kind(1.d0)), intent(inout) :: y(:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans

    call a%a%psbn_csmm(alpha,x,beta,y,info,trans) 

  end subroutine d_csmv

  subroutine d_cssm(alpha,a,x,beta,y,info,trans) 
    class(psbn_d_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:,:)
    real(kind(1.d0)), intent(inout) :: y(:,:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    
    call a%a%psbn_cssm(alpha,x,beta,y,info,trans) 

  end subroutine d_cssm

  subroutine d_cssv(alpha,a,x,beta,y,info,trans) 
    class(psbn_d_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:)
    real(kind(1.d0)), intent(inout) :: y(:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans

    call a%a%psbn_cssm(alpha,x,beta,y,info,trans) 

  end subroutine d_cssv

end module psbn_d_mat_mod

