module psb_d_diagprec

  use psb_d_base_prec_mod
  
  type, extends(psb_d_base_prec_type) :: psb_d_diag_prec_type
    real(psb_dpk_), allocatable     :: d(:)
    type(psb_d_vect_type), allocatable :: dv
  contains
    procedure, pass(prec) :: d_apply_v => psb_d_diag_apply_vect
    procedure, pass(prec) :: d_apply   => psb_d_diag_apply
    procedure, pass(prec) :: precbld   => psb_d_diag_precbld
    procedure, pass(prec) :: precinit  => psb_d_diag_precinit  
    procedure, pass(prec) :: precseti  => psb_d_diag_precseti
    procedure, pass(prec) :: precsetr  => psb_d_diag_precsetr
    procedure, pass(prec) :: precsetc  => psb_d_diag_precsetc
    procedure, pass(prec) :: precfree  => psb_d_diag_precfree
    procedure, pass(prec) :: precdescr => psb_d_diag_precdescr
    procedure, pass(prec) :: sizeof    => psb_d_diag_sizeof
    procedure, pass(prec) :: get_nzeros => psb_d_diag_get_nzeros
  end type psb_d_diag_prec_type

  private :: psb_d_diag_precseti,&
       & psb_d_diag_precsetr, psb_d_diag_precsetc, psb_d_diag_sizeof,&
       & psb_d_diag_precinit, psb_d_diag_precfree, psb_d_diag_precdescr,&
       & psb_d_diag_get_nzeros
  
  
  
  interface  psb_d_diag_apply_vect
    subroutine psb_d_diag_apply_vect(alpha,prec,x,beta,y,desc_data,info,trans,work)
      import :: psb_desc_type, psb_d_diag_prec_type, psb_d_vect_type, psb_dpk_
      type(psb_desc_type),intent(in)    :: desc_data
      class(psb_d_diag_prec_type), intent(inout)  :: prec
      type(psb_d_vect_type),intent(inout)   :: x
      real(psb_dpk_),intent(in)         :: alpha, beta
      type(psb_d_vect_type),intent(inout)   :: y
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
      real(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine psb_d_diag_apply_vect
  end interface psb_d_diag_apply_vect
  
  interface  psb_d_diag_apply
    subroutine psb_d_diag_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
      import :: psb_desc_type, psb_d_diag_prec_type, psb_d_vect_type, psb_dpk_
      type(psb_desc_type),intent(in)    :: desc_data
      class(psb_d_diag_prec_type), intent(in)  :: prec
      real(psb_dpk_),intent(inout)      :: x(:)
      real(psb_dpk_),intent(in)         :: alpha, beta
      real(psb_dpk_),intent(inout)      :: y(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
      real(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine psb_d_diag_apply
  end interface psb_d_diag_apply
  
  interface  psb_d_diag_precbld
    subroutine psb_d_diag_precbld(a,desc_a,prec,info,upd,amold,afmt,vmold)
      import :: psb_desc_type, psb_d_diag_prec_type, psb_d_vect_type, psb_dpk_, &
           & psb_dspmat_type, psb_d_base_sparse_mat, psb_d_base_vect_type
      type(psb_dspmat_type), intent(in), target :: a
      type(psb_desc_type), intent(in), target   :: desc_a
      class(psb_d_diag_prec_type),intent(inout) :: prec
      integer, intent(out)                      :: info
      character, intent(in), optional           :: upd
      character(len=*), intent(in), optional    :: afmt
      class(psb_d_base_sparse_mat), intent(in), optional :: amold
      class(psb_d_base_vect_type), intent(in), optional  :: vmold
    end subroutine psb_d_diag_precbld
  end interface psb_d_diag_precbld
  

contains
  

  subroutine psb_d_diag_precinit(prec,info)
    Implicit None
    
    class(psb_d_diag_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='d_diag_precinit'

    call psb_erractionsave(err_act)

    info = psb_success_

    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_d_diag_precinit


  subroutine psb_d_diag_precseti(prec,what,val,info)
    Implicit None
    
    class(psb_d_diag_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    integer, intent(in)                      :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='d_diag_precset'

    call psb_erractionsave(err_act)

    info = psb_success_
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_d_diag_precseti

  subroutine psb_d_diag_precsetr(prec,what,val,info)
    
    Implicit None
    
    class(psb_d_diag_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    real(psb_dpk_), intent(in)               :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='d_diag_precset'

    call psb_erractionsave(err_act)

    info = psb_success_
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_d_diag_precsetr

  subroutine psb_d_diag_precsetc(prec,what,val,info)
    Implicit None
    
    class(psb_d_diag_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    character(len=*), intent(in)             :: val
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='d_diag_precset'

    call psb_erractionsave(err_act)

    info = psb_success_
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_d_diag_precsetc

  subroutine psb_d_diag_precfree(prec,info)
    
    Implicit None

    class(psb_d_diag_prec_type), intent(inout) :: prec
    integer, intent(out)                :: info
    
    Integer :: err_act, nrow
    character(len=20)  :: name='d_diag_precset'
    
    call psb_erractionsave(err_act)
    
    info = psb_success_

    if (allocated(prec%dv)) call prec%dv%free(info)
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine psb_d_diag_precfree
  

  subroutine psb_d_diag_precdescr(prec,iout)
    Implicit None

    class(psb_d_diag_prec_type), intent(in) :: prec
    integer, intent(in), optional    :: iout

    Integer :: err_act, nrow, info
    character(len=20)  :: name='d_diag_precdescr'

    integer :: iout_

    call psb_erractionsave(err_act)

    info = psb_success_
   
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if

    write(iout_,*) 'Diagonal scaling'

    call psb_erractionsave(err_act)

    info = psb_success_
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine psb_d_diag_precdescr

  function psb_d_diag_sizeof(prec) result(val)
    class(psb_d_diag_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    
    val = psb_sizeof_dp * prec%get_nzeros()
    return
  end function psb_d_diag_sizeof

  function psb_d_diag_get_nzeros(prec) result(val)
    class(psb_d_diag_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val

    val = 0
    if (allocated(prec%dv)) val = val + prec%dv%get_nrows()
    return
  end function psb_d_diag_get_nzeros


end module psb_d_diagprec
