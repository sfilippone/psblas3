
subroutine psb_dprecset(p,ptype,iv,rs,rv,info)

  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  implicit none

  type(psb_dprec_type), intent(inout)    :: p
  character(len=*), intent(in)           :: ptype
  integer, optional, intent(in)          :: iv(:)
  real(kind(1.d0)), optional, intent(in) :: rs
  real(kind(1.d0)), optional, intent(in) :: rv(:)
  integer, optional, intent(out)         :: info

  type(psb_dbase_prec), pointer          :: bpv(:)=>null()
  character(len=len(ptype))              :: typeup
  integer                                :: isz, err

  if (present(info)) info = 0

  if (.not.associated(p%baseprecv)) then 
    allocate(p%baseprecv(1),stat=err)
    call psb_nullify_baseprec(p%baseprecv(1))
  endif

  if (.not.associated(p%baseprecv(1)%iprcparm)) then 
    allocate(p%baseprecv(1)%iprcparm(ifpsz),stat=err)
    if (err/=0) then 
      write(0,*)'Precset Memory Failure',err
    endif
  end if

  select case(toupper(ptype(1:len_trim(ptype))))
  case ('NONE','NOPREC') 
    p%baseprecv(1)%iprcparm(p_type_)     = noprec_
    p%baseprecv(1)%iprcparm(f_type_)     = f_none_
    p%baseprecv(1)%iprcparm(restr_)      = psb_none_
    p%baseprecv(1)%iprcparm(prol_)       = psb_none_
    p%baseprecv(1)%iprcparm(iren_)       = 0
    p%baseprecv(1)%iprcparm(n_ovr_)      = 0
    p%baseprecv(1)%iprcparm(jac_sweeps_) = 1

  case ('DIAG','DIAGSC')
    p%baseprecv(1)%iprcparm(p_type_)     = diagsc_
    p%baseprecv(1)%iprcparm(f_type_)     = f_none_
    p%baseprecv(1)%iprcparm(restr_)      = psb_none_
    p%baseprecv(1)%iprcparm(prol_)       = psb_none_
    p%baseprecv(1)%iprcparm(iren_)       = 0 
    p%baseprecv(1)%iprcparm(n_ovr_)      = 0
    p%baseprecv(1)%iprcparm(jac_sweeps_) = 1

  case ('BJA','ILU') 
    p%baseprecv(1)%iprcparm(p_type_)      = bja_
    p%baseprecv(1)%iprcparm(f_type_)      = f_ilu_n_
    p%baseprecv(1)%iprcparm(restr_)       = psb_none_
    p%baseprecv(1)%iprcparm(prol_)        = psb_none_
    p%baseprecv(1)%iprcparm(iren_)        = 0
    p%baseprecv(1)%iprcparm(n_ovr_)       = 0
    p%baseprecv(1)%iprcparm(ilu_fill_in_) = 0
    p%baseprecv(1)%iprcparm(jac_sweeps_)  = 1

  case ('ASM','AS')
    ! Defaults first 
    p%baseprecv(1)%iprcparm(p_type_)      = asm_
    p%baseprecv(1)%iprcparm(f_type_)      = f_ilu_n_
    p%baseprecv(1)%iprcparm(restr_)       = psb_halo_
    p%baseprecv(1)%iprcparm(prol_)        = psb_none_
    p%baseprecv(1)%iprcparm(iren_)        = 0
    p%baseprecv(1)%iprcparm(n_ovr_)       = 1
    p%baseprecv(1)%iprcparm(ilu_fill_in_) = 0
    p%baseprecv(1)%iprcparm(jac_sweeps_)  = 1
    if (present(iv)) then 
      isz = size(iv) 
      if (isz >= 1) p%baseprecv(1)%iprcparm(n_ovr_)  = iv(1)
      if (isz >= 2) p%baseprecv(1)%iprcparm(restr_)  = iv(2)
      if (isz >= 3) p%baseprecv(1)%iprcparm(prol_)   = iv(3)
      if (isz >= 4) p%baseprecv(1)%iprcparm(f_type_) = iv(4) 
      ! Do not consider renum for the time being. 
!!$      if (isz >= 5) p%baseprecv(1)%iprcparm(iren_) = iv(5)
    end if


  case ('ML', '2LEV')

    select case (size(p%baseprecv)) 
    case(1)
      ! Reallocate
      allocate(bpv(2),stat=err)
      if (err/=0) then 
        write(0,*)'Precset Memory Failure 2l:1',err
      endif
      bpv(1) = p%baseprecv(1)
      call psb_nullify_baseprec(bpv(2))
      deallocate(p%baseprecv)
      p%baseprecv => bpv
      nullify(bpv)

    case(2)
      ! Do nothing

    case default
      ! Error

    end select

    allocate(p%baseprecv(2)%iprcparm(ifpsz),stat=err)
    if (err/=0) then 
      write(0,*)'Precset Memory Failure 2l:2',err
    endif
    allocate(p%baseprecv(2)%dprcparm(dfpsz),stat=err)
    if (err/=0) then 
      write(0,*)'Precset Memory Failure 2l:3',err
    endif



    p%baseprecv(2)%iprcparm(p_type_)       = bja_
    p%baseprecv(2)%iprcparm(ml_type_)      = mult_ml_prec_
    p%baseprecv(2)%iprcparm(aggr_alg_)     = loc_aggr_
    p%baseprecv(2)%iprcparm(smth_kind_)    = smth_omg_
    p%baseprecv(2)%iprcparm(coarse_mat_)   = mat_distr_
    p%baseprecv(2)%iprcparm(smth_pos_)     = post_smooth_
    p%baseprecv(2)%iprcparm(glb_smth_)     = 1
    p%baseprecv(2)%iprcparm(om_choice_)    = lib_choice_
    p%baseprecv(2)%iprcparm(f_type_)       = f_ilu_n_
    p%baseprecv(2)%iprcparm(ilu_fill_in_)  = 0
    p%baseprecv(2)%dprcparm(smooth_omega_) = 4.d0/3.d0         
    p%baseprecv(2)%iprcparm(jac_sweeps_)   = 1


    if (present(iv)) then 
      isz = size(iv)
      if (isz >= 1) p%baseprecv(2)%iprcparm(ml_type_)      = iv(1)
      if (isz >= 2) p%baseprecv(2)%iprcparm(aggr_alg_)     = iv(2) 
      if (isz >= 3) p%baseprecv(2)%iprcparm(smth_kind_)    = iv(3) 
      if (isz >= 4) p%baseprecv(2)%iprcparm(coarse_mat_)   = iv(4) 
      if (isz >= 5) p%baseprecv(2)%iprcparm(smth_pos_)     = iv(5)
      if (isz >= 6) p%baseprecv(2)%iprcparm(glb_smth_)     = iv(6)
      if (isz >= 7) p%baseprecv(2)%iprcparm(f_type_)       = iv(7)
      if (isz >= 8) p%baseprecv(2)%iprcparm(jac_sweeps_)   = iv(8)

    end if

    if (present(rs)) then 
      p%baseprecv(2)%iprcparm(om_choice_)    = user_choice_
      p%baseprecv(2)%dprcparm(smooth_omega_) = rs      
    end if


  case default
    write(0,*) 'Unknown preconditioner type request "',ptype,'"'
    err = 2

  end select

  if (present(info)) info = err

end subroutine psb_dprecset
