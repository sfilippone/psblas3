module psb_z_diagprec
  use psb_prec_type

  
  type, extends(psb_z_base_prec_type) :: psb_z_diag_prec_type
    complex(psb_dpk_), allocatable :: d(:)
  contains
    procedure, pass(prec) :: apply     => z_diag_apply
    procedure, pass(prec) :: precbld   => z_diag_precbld
    procedure, pass(prec) :: precinit  => z_diag_precinit  
    procedure, pass(prec) :: precseti  => z_diag_precseti
    procedure, pass(prec) :: precsetr  => z_diag_precsetr
    procedure, pass(prec) :: precsetc  => z_diag_precsetc
    procedure, pass(prec) :: precfree  => z_diag_precfree
    procedure, pass(prec) :: precdescr => z_diag_precdescr
    procedure, pass(prec) :: sizeof    => z_diag_sizeof
  end type psb_z_diag_prec_type


contains
  

  subroutine z_diag_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
    use psb_sparse_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_z_diag_prec_type), intent(in)  :: prec
    complex(psb_dpk_),intent(in)         :: x(:)
    complex(psb_dpk_),intent(in)         :: alpha, beta
    complex(psb_dpk_),intent(inout)      :: y(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans
    complex(psb_dpk_),intent(inout), optional, target :: work(:)
    Integer :: err_act, nrow
    character :: trans_
    character(len=20)  :: name='z_diag_prec_apply'
    complex(psb_dpk_), pointer :: ww(:)

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the DIAG preonditioner???
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
    if (.not.allocated(prec%d)) then
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner: D")
      goto 9999
    end if
    if (size(prec%d) < nrow) then
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner: D")
      goto 9999
    end if
    if (present(trans)) then 
      trans_ = psb_toupper(trans)
    else
      trans_='N'
    end if

    select case(trans_)
    case('N')
    case('T','C')
    case default
      info=40
      call psb_errpush(info,name,&
           & i_err=(/6,0,0,0,0/),a_err=trans_)
      goto 9999
    end select
    
    if (size(work) >= size(x)) then 
      ww => work
    else
      allocate(ww(size(x)),stat=info)
      if (info /= 0) then 
        call psb_errpush(4025,name,&
             & i_err=(/size(x),0,0,0,0/),a_err='complex(psb_dpk_)')
        goto 9999      
      end if
    end if


    if (trans_=='C') then 
      ww(1:nrow) = x(1:nrow)*conjg(prec%d(1:nrow))
    else
      ww(1:nrow) = x(1:nrow)*prec%d(1:nrow)
    endif
    call psb_geaxpby(alpha,ww,beta,y,desc_data,info)

    if (size(work) < size(x)) then 
      deallocate(ww,stat=info)
      if (info /= 0) then 
        call psb_errpush(4010,name,a_err='Deallocate')
        goto 9999      
      end if
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

  end subroutine z_diag_apply

  subroutine z_diag_precinit(prec,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_z_diag_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='z_diag_precinit'

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
  end subroutine z_diag_precinit


  subroutine z_diag_precbld(a,desc_a,prec,info,upd)
    
    use psb_sparse_mod
    Implicit None
    
    type(psb_z_sparse_mat), intent(in), target :: a
    type(psb_desc_type), intent(in), target  :: desc_a
    class(psb_z_diag_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    character, intent(in), optional          :: upd
    Integer :: err_act, nrow,i
    character(len=20)  :: name='z_diag_precbld'

    call psb_erractionsave(err_act)

    info = 0
    nrow = psb_cd_get_local_cols(desc_a)
    if (allocated(prec%d)) then 
      if (size(prec%d) < nrow) then 
        deallocate(prec%d,stat=info)
      end if
    end if
    if ((info == 0).and.(.not.allocated(prec%d))) then 
      allocate(prec%d(nrow), stat=info)
    end if
    if (info /= 0) then 
      info = 4000
      call psb_errpush(info,name)
      goto 9999
    end if

    call a%get_diag(prec%d,info) 
    if (info /= 0) then 
      info = 4010
      call psb_errpush(info,name, a_err='get_diag')
      goto 9999
    end if
    
    do i=1,nrow
      if (prec%d(i) == dzero) then
        prec%d(i) = done
      else
        prec%d(i) = done/prec%d(i)
      endif
    end do

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine z_diag_precbld

  subroutine z_diag_precseti(prec,what,val,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_z_diag_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    integer, intent(in)                      :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='z_diag_precset'

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
  end subroutine z_diag_precseti

  subroutine z_diag_precsetr(prec,what,val,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_z_diag_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    real(psb_dpk_), intent(in)               :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='z_diag_precset'

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
  end subroutine z_diag_precsetr

  subroutine z_diag_precsetc(prec,what,val,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_z_diag_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    character(len=*), intent(in)             :: val
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='z_diag_precset'

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
  end subroutine z_diag_precsetc

  subroutine z_diag_precfree(prec,info)
    
    use psb_sparse_mod
    Implicit None

    class(psb_z_diag_prec_type), intent(inout) :: prec
    integer, intent(out)                :: info
    
    Integer :: err_act, nrow
    character(len=20)  :: name='z_diag_precset'
    
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
    
  end subroutine z_diag_precfree
  

  subroutine z_diag_precdescr(prec,iout)
    
    use psb_sparse_mod
    Implicit None

    class(psb_z_diag_prec_type), intent(in) :: prec
    integer, intent(in), optional    :: iout

    Integer :: err_act, nrow, info
    character(len=20)  :: name='z_diag_precdescr'

    integer :: iout_

    call psb_erractionsave(err_act)

    info = 0
   
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if

    write(iout_,*) 'Diagonal scaling'

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
    
  end subroutine z_diag_precdescr

  function z_diag_sizeof(prec) result(val)
    use psb_sparse_mod
    class(psb_z_diag_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    
    val = 0
    val = val + 2*psb_sizeof_dp * size(prec%d)
    return
  end function z_diag_sizeof

end module psb_z_diagprec
