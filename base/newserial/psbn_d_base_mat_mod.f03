
module psbn_d_base_mat_mod

  use psbn_base_mat_mod

  type, extends(psbn_base_sparse_mat) :: psbn_d_base_sparse_mat
  contains
    procedure, pass(a) :: d_base_csmv
    procedure, pass(a) :: d_base_csmm
    generic, public    :: psbn_csmm => d_base_csmm, d_base_csmv
    procedure, pass(a) :: d_base_cssv
    procedure, pass(a) :: d_base_cssm
    generic, public    :: psbn_cssm => d_base_cssm, d_base_cssv
    procedure, pass(a) :: csins
    
  end type psbn_d_base_sparse_mat

contains 

  subroutine csins(nz,val,ia,ja,a,imin,imax,jmin,jmax,info,gtl) 
    use psb_error_mod
    use psb_realloc_mod
    class(psbn_d_base_sparse_mat), intent(inout) :: a
    real(psb_dpk_), intent(in)      :: val(:)
    integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
    integer, intent(out)            :: info
    integer, intent(in), optional   :: gtl(:)

    Integer :: err_act
    character(len=20)  :: name='csins'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine csins

  subroutine d_base_csmm(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    class(psbn_d_base_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
    real(psb_dpk_), intent(inout) :: y(:,:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans

    Integer :: err_act
    character(len=20)  :: name='d_base_csmm'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
         
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine d_base_csmm

  subroutine d_base_csmv(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    class(psbn_d_base_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:)
    real(kind(1.d0)), intent(inout) :: y(:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans

    Integer :: err_act
    character(len=20)  :: name='d_base_csmv'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return


  end subroutine d_base_csmv

  subroutine d_base_cssm(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    class(psbn_d_base_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:,:)
    real(kind(1.d0)), intent(inout) :: y(:,:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    
    Integer :: err_act
    character(len=20)  :: name='d_base_cssm'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine d_base_cssm

  subroutine d_base_cssv(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    class(psbn_d_base_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:)
    real(kind(1.d0)), intent(inout) :: y(:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans

    Integer :: err_act
    character(len=20)  :: name='d_base_cssv'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return


  end subroutine d_base_cssv

end module psbn_d_base_mat_mod
