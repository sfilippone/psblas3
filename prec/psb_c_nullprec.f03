module psb_c_nullprec
  use psb_prec_type

  
  type, extends(psb_c_base_prec_type) :: psb_c_null_prec_type
  contains
    procedure, pass(prec) :: apply     => c_null_apply
    procedure, pass(prec) :: precbld   => c_null_precbld
    procedure, pass(prec) :: precinit  => c_null_precinit
    procedure, pass(prec) :: precseti  => c_null_precseti
    procedure, pass(prec) :: precsetr  => c_null_precsetr
    procedure, pass(prec) :: precsetc  => c_null_precsetc
    procedure, pass(prec) :: precfree  => c_null_precfree
    procedure, pass(prec) :: precdescr => c_null_precdescr
    procedure, pass(prec) :: sizeof    => c_null_sizeof
  end type psb_c_null_prec_type


contains
  

  subroutine c_null_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
    use psb_sparse_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_c_null_prec_type), intent(in)  :: prec
    complex(psb_spk_),intent(in)         :: x(:)
    complex(psb_spk_),intent(in)         :: alpha, beta
    complex(psb_spk_),intent(inout)      :: y(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans
    complex(psb_spk_),intent(inout), optional, target :: work(:)
    Integer :: err_act, nrow
    character(len=20)  :: name='c_null_prec_apply'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 0 
    
    nrow = psb_cd_get_local_rows(desc_data)
    if (size(x) < nrow) then 
      info = 36
      call psb_errpush(info,name,i_err=(/2,nrow,0,0,0/))
      goto 9999
    end if
    if (size(y) < nrow) then 
      info = 36
      call psb_errpush(info,name,i_err=(/3,nrow,0,0,0/))
      goto 9999
    end if

    call psb_geaxpby(alpha,x,beta,y,desc_data,info)
    if (info /= 0 ) then 
      info = 4010
      call psb_errpush(infoi,name,a_err="psb_geaxpby")
      goto 9999
    end if
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine c_null_apply


  subroutine c_null_precinit(prec,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_c_null_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='c_null_precinit'

    call psb_erractionsave(err_act)

    info = 0

    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine c_null_precinit

  subroutine c_null_precbld(a,desc_a,prec,info,upd)
    
    use psb_sparse_mod
    Implicit None
    
    type(psb_c_sparse_mat), intent(in), target :: a
    type(psb_desc_type), intent(in), target  :: desc_a
    class(psb_c_null_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    character, intent(in), optional          :: upd
    Integer :: err_act, nrow
    character(len=20)  :: name='c_null_precbld'

    call psb_erractionsave(err_act)

    info = 0

    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine c_null_precbld

  subroutine c_null_precseti(prec,what,val,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_c_null_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    integer, intent(in)                      :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='c_null_precset'

    call psb_erractionsave(err_act)

    info = 0
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine c_null_precseti

  subroutine c_null_precsetr(prec,what,val,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_c_null_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    real(psb_spk_), intent(in)               :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='c_null_precset'

    call psb_erractionsave(err_act)

    info = 0
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine c_null_precsetr

  subroutine c_null_precsetc(prec,what,val,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_c_null_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    character(len=*), intent(in)             :: val
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='c_null_precset'

    call psb_erractionsave(err_act)

    info = 0
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine c_null_precsetc

  subroutine c_null_precfree(prec,info)
    
    use psb_sparse_mod
    Implicit None

    class(psb_c_null_prec_type), intent(inout) :: prec
    integer, intent(out)                :: info
    
    Integer :: err_act, nrow
    character(len=20)  :: name='c_null_precset'
    
    call psb_erractionsave(err_act)
    
    info = 0
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine c_null_precfree
  

  subroutine c_null_precdescr(prec,iout)
    
    use psb_sparse_mod
    Implicit None

    class(psb_c_null_prec_type), intent(in) :: prec
    integer, intent(in), optional    :: iout

    Integer :: err_act, nrow, info
    character(len=20)  :: name='c_null_precset'
    integer :: iout_

    call psb_erractionsave(err_act)

    info = 0
   
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if

    write(iout_,*) 'No preconditioning'

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine c_null_precdescr

  function c_null_sizeof(prec) result(val)
    use psb_sparse_mod
    class(psb_c_null_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    
    val = 0

    return
  end function c_null_sizeof

end module psb_c_nullprec
