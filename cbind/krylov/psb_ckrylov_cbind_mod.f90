module psb_ckrylov_cbind_mod

  use psb_base_krylov_cbind_mod

contains 

  function  psb_c_ckrylov(methd,&
       & ah,ph,bh,xh,cdh,options) bind(c) result(res)
    use psb_base_mod
    use psb_prec_mod
    use psb_krylov_mod
    use psb_objhandle_mod
    use psb_prec_cbind_mod
    use psb_base_string_cbind_mod
    implicit none 
    integer(psb_c_ipk)      :: res
    type(psb_c_cspmat)    :: ah
    type(psb_c_descriptor)  :: cdh
    type(psb_c_cprec)       :: ph
    type(psb_c_cvector)     :: bh,xh
    character(c_char)       :: methd(*)
    type(solveroptions)     :: options

    res= psb_c_ckrylov_opt(methd, ah, ph, bh, xh, options%eps,cdh,  &
         & itmax=options%itmax, iter=options%iter,&
         & itrace=options%itrace, istop=options%istop,&
         & irst=options%irst, err=options%err)
    
  end function psb_c_ckrylov


  function  psb_c_ckrylov_opt(methd,&
       & ah,ph,bh,xh,eps,cdh,itmax,iter,err,itrace,irst,istop) bind(c) result(res)
    use psb_base_mod
    use psb_prec_mod
    use psb_krylov_mod
    use psb_objhandle_mod
    use psb_prec_cbind_mod
    use psb_base_string_cbind_mod
    implicit none 
    integer(psb_c_ipk)      :: res
    type(psb_c_cspmat)    :: ah
    type(psb_c_descriptor)  :: cdh
    type(psb_c_cprec)       :: ph
    type(psb_c_cvector)     :: bh,xh
    integer(psb_c_ipk), value :: itmax,itrace,irst,istop
    real(c_double), value :: eps
    integer(psb_c_ipk)    :: iter
    real(c_double)        :: err
    character(c_char)       :: methd(*)
    type(solveroptions)     :: options
    type(psb_desc_type), pointer   :: descp
    type(psb_cspmat_type), pointer :: ap
    type(psb_cprec_type), pointer  :: precp
    type(psb_c_vect_type), pointer :: xp, bp

    integer               :: info,fitmax,fitrace,first,fistop,fiter
    character(len=20)     :: fmethd
    real(psb_spk_)      :: feps,ferr

    res = -1
    if (c_associated(cdh%item)) then 
      call c_f_pointer(cdh%item,descp)
    else
      return 
    end if
    if (c_associated(xh%item)) then 
      call c_f_pointer(xh%item,xp)
    else
      return 
    end if
    if (c_associated(bh%item)) then 
      call c_f_pointer(bh%item,bp)
    else
      return 
    end if
    if (c_associated(ah%item)) then 
      call c_f_pointer(ah%item,ap)
    else
      return 
    end if
    if (c_associated(ph%item)) then 
      call c_f_pointer(ph%item,precp)
    else
      return 
    end if

    
    call stringc2f(methd,fmethd)
    feps    = eps
    fitmax  = itmax
    fitrace = itrace
    first   = irst
    fistop  = istop

    call psb_krylov(fmethd, ap, precp, bp, xp, feps, &
         & descp, info,&
         & itmax=fitmax,iter=fiter,itrace=fitrace,istop=fistop,&
         & irst=first, err=ferr)
    iter = fiter
    err  = ferr
    res = info
    
  end function psb_c_ckrylov_opt

end module psb_ckrylov_cbind_mod
