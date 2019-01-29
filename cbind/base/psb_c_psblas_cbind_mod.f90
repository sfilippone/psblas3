module psb_c_psblas_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_string_cbind_mod
  
contains
  
  function psb_c_cgeaxpby(alpha,xh,beta,yh,cdh) bind(c) result(res)
    implicit none 
    integer(psb_c_int) :: res

    type(psb_c_cvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    complex(c_float_complex), value   :: alpha,beta
    
    type(psb_desc_type), pointer :: descp
    type(psb_c_vect_type), pointer :: xp,yp
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

  end function psb_c_cgeaxpby

  function psb_c_cgenrm2(xh,cdh) bind(c) result(res)
    implicit none 
    real(c_float_complex) :: res

    type(psb_c_cvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_c_vect_type), pointer :: xp
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

  end function psb_c_cgenrm2
  
  function psb_c_cgeamax(xh,cdh) bind(c) result(res)
    implicit none 
    real(c_float_complex) :: res

    type(psb_c_cvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_c_vect_type), pointer :: xp
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

  end function psb_c_cgeamax
  
  function psb_c_cgeasum(xh,cdh) bind(c) result(res)
    implicit none 
    real(c_float_complex) :: res

    type(psb_c_cvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer   :: descp
    type(psb_c_vect_type), pointer :: xp
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

  end function psb_c_cgeasum

  
  function psb_c_cspnrmi(ah,cdh) bind(c) result(res)
    implicit none 
    real(c_float_complex) :: res

    type(psb_c_cspmat)   :: ah
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_cspmat_type), pointer :: ap
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

  end function psb_c_cspnrmi

  function psb_c_cgedot(xh,yh,cdh) bind(c) result(res)
    implicit none 
    complex(c_float_complex) :: res

    type(psb_c_cvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_c_vect_type), pointer :: xp,yp
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

  end function psb_c_cgedot


  function psb_c_cspmm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
    implicit none 
    integer(psb_c_int) :: res

    type(psb_c_cspmat) :: ah
    type(psb_c_cvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    complex(c_float_complex), value :: alpha, beta
    type(psb_desc_type), pointer :: descp
    type(psb_c_vect_type), pointer :: xp,yp
    type(psb_cspmat_type), pointer :: ap
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

  end function psb_c_cspmm


  function psb_c_cspmm_opt(alpha,ah,xh,beta,yh,cdh,trans,doswap) bind(c) result(res)
    implicit none 
    integer(psb_c_int) :: res

    type(psb_c_cspmat) :: ah
    type(psb_c_cvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    complex(c_float_complex), value :: alpha, beta
    character(c_char)      :: trans
    logical(c_bool), value :: doswap

    type(psb_desc_type), pointer :: descp
    type(psb_c_vect_type), pointer :: xp,yp
    type(psb_cspmat_type), pointer :: ap
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

  end function psb_c_cspmm_opt
  

  function psb_c_cspsm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
    implicit none 
    integer(psb_c_int) :: res

    type(psb_c_cspmat) :: ah
    type(psb_c_cvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    complex(c_float_complex), value :: alpha, beta
    type(psb_desc_type), pointer :: descp
    type(psb_c_vect_type), pointer :: xp,yp
    type(psb_cspmat_type), pointer :: ap
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

  end function psb_c_cspsm
  

end module psb_c_psblas_cbind_mod
