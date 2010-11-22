module psb_z_bjacprec

  use psb_z_base_prec_mod
  
  type, extends(psb_z_base_prec_type) :: psb_z_bjac_prec_type
    integer, allocatable                :: iprcparm(:)
    type(psb_zspmat_type), allocatable  :: av(:)
    complex(psb_dpk_), allocatable      :: d(:)
  contains
    procedure, pass(prec) :: apply     => psb_z_bjac_apply
    procedure, pass(prec) :: precbld   => psb_z_bjac_precbld
    procedure, pass(prec) :: precinit  => psb_z_bjac_precinit
    procedure, pass(prec) :: precseti  => psb_z_bjac_precseti
    procedure, pass(prec) :: precsetr  => psb_z_bjac_precsetr
    procedure, pass(prec) :: precsetc  => psb_z_bjac_precsetc
    procedure, pass(prec) :: precfree  => psb_z_bjac_precfree
    procedure, pass(prec) :: precdescr => psb_z_bjac_precdescr
    procedure, pass(prec) :: sizeof    => psb_z_bjac_sizeof
  end type psb_z_bjac_prec_type

  private :: psb_z_bjac_apply, psb_z_bjac_precbld, psb_z_bjac_precseti,&
       & psb_z_bjac_precsetr, psb_z_bjac_precsetc, psb_z_bjac_sizeof,&
       & psb_z_bjac_precinit, psb_z_bjac_precfree, psb_z_bjac_precdescr
  

  character(len=15), parameter, private :: &
       &  fact_names(0:2)=(/'None          ','ILU(n)        ',&
       &  'ILU(eps)      '/)

contains
  

  subroutine psb_z_bjac_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
    use psb_sparse_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_z_bjac_prec_type), intent(in)  :: prec
    complex(psb_dpk_),intent(in)         :: alpha,beta
    complex(psb_dpk_),intent(in)         :: x(:)
    complex(psb_dpk_),intent(inout)      :: y(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans
    complex(psb_dpk_),intent(inout), optional, target :: work(:)

    ! Local variables
    integer :: n_row,n_col
    complex(psb_dpk_), pointer :: ww(:), aux(:)
    integer :: ictxt,np,me, err_act, int_err(5)
    integer            :: debug_level, debug_unit
    character          :: trans_
    character(len=20)  :: name='z_bjac_prec_apply'
    character(len=20)  :: ch_err

    info = psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    ictxt       = psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)
    
    
    trans_ = psb_toupper(trans)
    select case(trans_)
    case('N','T','C')
      ! Ok
    case default
      call psb_errpush(psb_err_iarg_invalid_i_,name)
      goto 9999
    end select
    
    
    n_row = psb_cd_get_local_rows(desc_data)
    n_col = psb_cd_get_local_cols(desc_data)

    if (size(x) < n_row) then 
      info = 36
      call psb_errpush(info,name,i_err=(/2,n_row,0,0,0/))
      goto 9999
    end if
    if (size(y) < n_row) then 
      info = 36
      call psb_errpush(info,name,i_err=(/3,n_row,0,0,0/))
      goto 9999
    end if
    if (.not.allocated(prec%d)) then
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner: D")
      goto 9999
    end if
    if (size(prec%d) < n_row) then
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner: D")
      goto 9999
    end if

    
    if (n_col <= size(work)) then 
      ww => work(1:n_col)
      if ((4*n_col+n_col) <= size(work)) then 
        aux => work(n_col+1:)
      else
        allocate(aux(4*n_col),stat=info)
        if (info /= psb_success_) then 
          call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
          goto 9999      
        end if
        
      endif
    else
      allocate(ww(n_col),aux(4*n_col),stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
        goto 9999      
      end if
    endif
    
    
    select case(prec%iprcparm(psb_f_type_))
    case(psb_f_ilu_n_) 
      
      select case(trans_)
      case('N')
        call psb_spsm(zone,prec%av(psb_l_pr_),x,zzero,ww,desc_data,info,&
             & trans=trans_,scale='L',diag=prec%d,choice=psb_none_,work=aux)
        if(info == psb_success_) call psb_spsm(alpha,prec%av(psb_u_pr_),ww,beta,y,desc_data,info,&
             & trans=trans_,scale='U',choice=psb_none_, work=aux)
        
      case('T')
        call psb_spsm(zone,prec%av(psb_u_pr_),x,zzero,ww,desc_data,info,&
             & trans=trans_,scale='L',diag=prec%d,choice=psb_none_, work=aux)
        if(info == psb_success_)  call psb_spsm(alpha,prec%av(psb_l_pr_),ww,beta,y,desc_data,info,&
           & trans=trans_,scale='U',choice=psb_none_,work=aux)

      case('C')
        call psb_spsm(zone,prec%av(psb_u_pr_),x,zzero,ww,desc_data,info,&
             & trans=trans_,scale='L',diag=conjg(prec%d),choice=psb_none_, work=aux)
        if(info == psb_success_)  call psb_spsm(alpha,prec%av(psb_l_pr_),ww,beta,y,desc_data,info,&
           & trans=trans_,scale='U',choice=psb_none_,work=aux)
        
      end select
      if (info /= psb_success_) then 
        ch_err="psb_spsm"
        goto 9999
      end if
      
      
    case default
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='Invalid factorization')
      goto 9999
    end select
    
    call psb_halo(y,desc_data,info,data=psb_comm_mov_)
    
    if (n_col <= size(work)) then 
      if ((4*n_col+n_col) <= size(work)) then 
      else
        deallocate(aux)
      endif
    else
      deallocate(ww,aux)
    endif
    
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_errpush(info,name,i_err=int_err,a_err=ch_err)
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return


  end subroutine psb_z_bjac_apply

  subroutine psb_z_bjac_precinit(prec,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_z_bjac_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='z_bjac_precinit'

    call psb_erractionsave(err_act)

    info = psb_success_
    call psb_realloc(psb_ifpsz,prec%iprcparm,info)
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_
      call psb_Errpush(info,name)
      goto 9999
    end if
    
    prec%iprcparm(:)                = 0
    prec%iprcparm(psb_p_type_)      = psb_bjac_
    prec%iprcparm(psb_f_type_)      = psb_f_ilu_n_
    prec%iprcparm(psb_ilu_fill_in_) = 0


    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_z_bjac_precinit


  subroutine psb_z_bjac_precbld(a,desc_a,prec,info,upd,mold,afmt)

    use psb_sparse_mod
    use psb_prec_mod
    Implicit None

    type(psb_zspmat_type), intent(in), target :: a
    type(psb_desc_type), intent(in), target   :: desc_a
    class(psb_z_bjac_prec_type),intent(inout) :: prec
    integer, intent(out)                      :: info
    character, intent(in), optional           :: upd
    character(len=*), intent(in), optional    :: afmt
    class(psb_z_base_sparse_mat), intent(in), optional :: mold

    !     .. Local Scalars ..                                                       
    integer  ::    i, m
    integer  ::    int_err(5)
    character ::        trans, unitd
    type(psb_z_csr_sparse_mat), allocatable  :: lf, uf
    integer   nztota,  err_act, n_row, nrow_a,n_col, nhalo
    integer :: ictxt,np,me
    character(len=20)  :: name='z_bjac_precbld'
    character(len=20)  :: ch_err


    if(psb_get_errstatus() /= 0) return 
    info = psb_success_

    call psb_erractionsave(err_act)

    ictxt=psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)

    m = a%get_nrows()
    if (m < 0) then
      info = psb_err_iarg_neg_
      int_err(1) = 1
      int_err(2) = m
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    endif
    trans = 'N'
    unitd = 'U'

    select case(prec%iprcparm(psb_f_type_))

    case(psb_f_ilu_n_) 

      if (allocated(prec%av)) then 
        if (size(prec%av) < psb_bp_ilu_avsz) then 
          do i=1,size(prec%av) 
            call prec%av(i)%free()
          enddo
          deallocate(prec%av,stat=info)
        endif
      end if
      if (.not.allocated(prec%av)) then 
        allocate(prec%av(psb_max_avsz),stat=info)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_alloc_dealloc_,name)
          goto 9999
        end if
      endif

      nrow_a = psb_cd_get_local_rows(desc_a)
      nztota = a%get_nzeros()

      n_col  = psb_cd_get_local_cols(desc_a)
      nhalo  = n_col-nrow_a
      n_row  = nrow_a

      allocate(lf,uf,stat=info)
      if (info == psb_success_) call lf%allocate(n_row,n_row,nztota)
      if (info == psb_success_) call uf%allocate(n_row,n_row,nztota)

      if(info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='psb_sp_all'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

      if (allocated(prec%d)) then 
        if (size(prec%d) < n_row) then 
          deallocate(prec%d)
        endif
      endif
      if (.not.allocated(prec%d)) then 
        allocate(prec%d(n_row),stat=info)
        if (info /= psb_success_) then 
          call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
          goto 9999      
        end if

      endif
      ! This is where we have no renumbering, thus no need 
      call psb_ilu_fct(a,lf,uf,prec%d,info)

      if(info == psb_success_) then
        call prec%av(psb_l_pr_)%mv_from(lf)
        call prec%av(psb_u_pr_)%mv_from(uf)
        call prec%av(psb_l_pr_)%set_asb()
        call prec%av(psb_u_pr_)%set_asb()
        call prec%av(psb_l_pr_)%trim()
        call prec%av(psb_u_pr_)%trim()
      else
        info=psb_err_from_subroutine_
        ch_err='psb_ilu_fct'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
      
    case(psb_f_none_) 
      info=psb_err_from_subroutine_
      ch_err='Inconsistent prec  psb_f_none_'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999

    case default
      info=psb_err_from_subroutine_
      ch_err='Unknown psb_f_type_'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end select

    if (present(mold)) then 
      call prec%av(psb_l_pr_)%cscnv(info,mold=mold)
      call prec%av(psb_u_pr_)%cscnv(info,mold=mold)
    else if (present(afmt)) then 
      call prec%av(psb_l_pr_)%cscnv(info,type=afmt)
      call prec%av(psb_u_pr_)%cscnv(info,type=afmt)
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

  end subroutine psb_z_bjac_precbld

  subroutine psb_z_bjac_precseti(prec,what,val,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_z_bjac_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    integer, intent(in)                      :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='z_bjac_precset'

    call psb_erractionsave(err_act)

    info = psb_success_
    if (.not.allocated(prec%iprcparm)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
      goto 9999
    end if
    
    select case(what)
    case (psb_f_type_) 
      if (prec%iprcparm(psb_p_type_) /= psb_bjac_) then 
        write(psb_err_unit,*) 'WHAT is invalid for current preconditioner ',prec%iprcparm(psb_p_type_),&
             & 'ignoring user specification'
        return
      endif
      prec%iprcparm(psb_f_type_)     = val
      
    case (psb_ilu_fill_in_) 
      if ((prec%iprcparm(psb_p_type_) /= psb_bjac_).or.(prec%iprcparm(psb_f_type_) /= psb_f_ilu_n_)) then 
        write(psb_err_unit,*) 'WHAT is invalid for current preconditioner ',prec%iprcparm(psb_p_type_),&
             & 'ignoring user specification'
        return
      endif
      prec%iprcparm(psb_ilu_fill_in_) = val
      
    case default
      write(psb_err_unit,*) 'WHAT is invalid, ignoring user specification'
      
    end select
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_z_bjac_precseti

  subroutine psb_z_bjac_precsetr(prec,what,val,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_z_bjac_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    real(psb_dpk_), intent(in)               :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='z_bjac_precset'

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
  end subroutine psb_z_bjac_precsetr

  subroutine psb_z_bjac_precsetc(prec,what,val,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_z_bjac_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    character(len=*), intent(in)             :: val
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='c_bjac_precset'

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
  end subroutine psb_z_bjac_precsetc

  subroutine psb_z_bjac_precfree(prec,info)
    
    use psb_sparse_mod
    Implicit None

    class(psb_z_bjac_prec_type), intent(inout) :: prec
    integer, intent(out)                :: info
    
    Integer :: err_act, i
    character(len=20)  :: name='z_bjac_precfree'
    
    call psb_erractionsave(err_act)
    
    info = psb_success_
    if (allocated(prec%av)) then 
      do i=1,size(prec%av) 
        call prec%av(i)%free()
      enddo
      deallocate(prec%av,stat=info)
    end if
    if (allocated(prec%d)) then 
      deallocate(prec%d,stat=info)
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
    
  end subroutine psb_z_bjac_precfree
  

  subroutine psb_z_bjac_precdescr(prec,iout)
    
    use psb_sparse_mod
    Implicit None

    class(psb_z_bjac_prec_type), intent(in) :: prec
    integer, intent(in), optional    :: iout

    Integer :: err_act, nrow, info
    character(len=20)  :: name='z_bjac_precdescr'
    integer :: iout_

    call psb_erractionsave(err_act)

    info = psb_success_
   
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if

    if (.not.allocated(prec%iprcparm)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
      goto 9999
    end if
    
    write(iout_,*) 'Block Jacobi with: ',&
         &  fact_names(prec%iprcparm(psb_f_type_))

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
    
  end subroutine psb_z_bjac_precdescr

  function psb_z_bjac_sizeof(prec) result(val)
    use psb_sparse_mod
    class(psb_z_bjac_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    
    val = 0
    if (allocated(prec%d)) then 
      val = val + 2*psb_sizeof_dp * size(prec%d)
    endif
    if (allocated(prec%av)) then 
      val = val + psb_sizeof(prec%av(psb_l_pr_))
      val = val + psb_sizeof(prec%av(psb_u_pr_))
    endif
    return
  end function psb_z_bjac_sizeof

end module psb_z_bjacprec
