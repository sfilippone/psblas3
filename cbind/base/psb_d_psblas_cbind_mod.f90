module psb_d_psblas_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_string_cbind_mod
  
contains
  
  function psb_c_dgeaxpby(alpha,xh,beta,yh,cdh) bind(c) result(res)
    implicit none 
    integer(psb_c_int) :: res

    type(psb_c_dvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    real(c_double), value   :: alpha,beta
    
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    integer                 :: info
    

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
    if (c_associated(yh%item)) then 
      call c_f_pointer(yh%item,yp)
    else
      return 
    end if

    call psb_geaxpby(alpha,xp,beta,yp,descp,info)

    res = info

  end function psb_c_dgeaxpby

  function psb_c_dgenrm2(xh,cdh) bind(c) result(res)
    implicit none 
    real(c_double) :: res

    type(psb_c_dvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp
    integer                :: info

    res = -1.0

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
    
    res = psb_genrm2(xp,descp,info)

  end function psb_c_dgenrm2
  
  function psb_c_dgeamax(xh,cdh) bind(c) result(res)
    implicit none 
    real(c_double) :: res

    type(psb_c_dvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp
    integer                 :: info

    res = -1.0
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
    
    res = psb_geamax(xp,descp,info)

  end function psb_c_dgeamax
  
  function psb_c_dgeasum(xh,cdh) bind(c) result(res)
    implicit none 
    real(c_double) :: res

    type(psb_c_dvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer   :: descp
    type(psb_d_vect_type), pointer :: xp
    integer                 :: info

    res = -1.0

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
    
    res = psb_geasum(xp,descp,info)

  end function psb_c_dgeasum

  
  function psb_c_dspnrmi(ah,cdh) bind(c) result(res)
    implicit none 
    real(c_double) :: res

    type(psb_c_dspmat)   :: ah
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_dspmat_type), pointer :: ap
    integer                 ::  info

    res = -1.0
    if (c_associated(cdh%item)) then 
      call c_f_pointer(cdh%item,descp)
    else
      return 
    end if
    if (c_associated(ah%item)) then 
      call c_f_pointer(ah%item,ap)
    else
      return 
    end if

    res = psb_spnrmi(ap,descp,info)

  end function psb_c_dspnrmi

  function psb_c_dgedot(xh,yh,cdh) bind(c) result(res)
    implicit none 
    real(c_double) :: res

    type(psb_c_dvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    integer               :: info

    res = -1.0
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
    if (c_associated(yh%item)) then 
      call c_f_pointer(yh%item,yp)
    else
      return 
    end if
    res = psb_gedot(xp,yp,descp,info)

  end function psb_c_dgedot


  function psb_c_dspmm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
    implicit none 
    integer(psb_c_int) :: res

    type(psb_c_dspmat) :: ah
    type(psb_c_dvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    real(c_double), value :: alpha, beta
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    type(psb_dspmat_type), pointer :: ap
    integer               :: info

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
    if (c_associated(yh%item)) then 
      call c_f_pointer(yh%item,yp)
    else
      return 
    end if
    if (c_associated(ah%item)) then 
      call c_f_pointer(ah%item,ap)
    else
      return 
    end if

    call psb_spmm(alpha,ap,xp,beta,yp,descp,info)

    res = info

  end function psb_c_dspmm


  function psb_c_dspmm_opt(alpha,ah,xh,beta,yh,cdh,trans,doswap) bind(c) result(res)
    implicit none 
    integer(psb_c_int) :: res

    type(psb_c_dspmat) :: ah
    type(psb_c_dvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    real(c_double), value :: alpha, beta
    character(c_char)      :: trans
    logical(c_bool), value :: doswap

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    type(psb_dspmat_type), pointer :: ap
    character :: ftrans
    logical   :: fdoswap
    integer               :: info

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
    if (c_associated(yh%item)) then 
      call c_f_pointer(yh%item,yp)
    else
      return 
    end if
    if (c_associated(ah%item)) then 
      call c_f_pointer(ah%item,ap)
    else
      return 
    end if

    fdoswap = doswap
    ftrans  = trans
    call psb_spmm(alpha,ap,xp,beta,yp,descp,info,trans=ftrans,doswap=fdoswap)

    res = info

  end function psb_c_dspmm_opt
  

  function psb_c_dspsm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
    implicit none 
    integer(psb_c_int) :: res

    type(psb_c_dspmat) :: ah
    type(psb_c_dvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    real(c_double), value :: alpha, beta
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    type(psb_dspmat_type), pointer :: ap
    integer               :: info

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
    if (c_associated(yh%item)) then 
      call c_f_pointer(yh%item,yp)
    else
      return 
    end if
    if (c_associated(ah%item)) then 
      call c_f_pointer(ah%item,ap)
    else
      return 
    end if

    call psb_spsm(alpha,ap,xp,beta,yp,descp,info)

    res = info

  end function psb_c_dspsm
  

end module psb_d_psblas_cbind_mod
