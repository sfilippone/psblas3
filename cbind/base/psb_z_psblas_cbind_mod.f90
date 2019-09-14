module psb_z_psblas_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_string_cbind_mod
  
contains
  
  function psb_c_zgeaxpby(alpha,xh,beta,yh,cdh) bind(c) result(res)
    implicit none 
    integer(psb_c_ipk_) :: res

    type(psb_c_zvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    complex(c_double_complex), value   :: alpha,beta
    
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
    integer(psb_c_ipk_)          :: info
    

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

  end function psb_c_zgeaxpby

  function psb_c_zgenrm2(xh,cdh) bind(c) result(res)
    implicit none 
    real(c_double_complex) :: res

    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
    integer(psb_c_ipk_)               :: info

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

  end function psb_c_zgenrm2
  
  function psb_c_zgeamax(xh,cdh) bind(c) result(res)
    implicit none 
    real(c_double_complex) :: res

    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
    integer(psb_c_ipk_)          :: info

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

  end function psb_c_zgeamax
  
  function psb_c_zgeasum(xh,cdh) bind(c) result(res)
    implicit none 
    real(c_double_complex) :: res

    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer   :: descp
    type(psb_z_vect_type), pointer :: xp
    integer(psb_c_ipk_)          :: info

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

  end function psb_c_zgeasum

  
  function psb_c_zspnrmi(ah,cdh) bind(c) result(res)
    implicit none 
    real(c_double_complex) :: res

    type(psb_c_zspmat)   :: ah
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_zspmat_type), pointer :: ap
    integer(psb_c_ipk_)              ::  info

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

  end function psb_c_zspnrmi

  function psb_c_zgedot(xh,yh,cdh) bind(c) result(res)
    implicit none 
    complex(c_double_complex) :: res

    type(psb_c_zvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
    integer(psb_c_ipk_)               :: info

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

  end function psb_c_zgedot


  function psb_c_zspmm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
    implicit none 
    integer(psb_c_ipk_) :: res

    type(psb_c_zspmat) :: ah
    type(psb_c_zvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    complex(c_double_complex), value :: alpha, beta
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
    type(psb_zspmat_type), pointer :: ap
    integer(psb_c_ipk_)               :: info

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

  end function psb_c_zspmm


  function psb_c_zspmm_opt(alpha,ah,xh,beta,yh,cdh,trans,doswap) bind(c) result(res)
    implicit none 
    integer(psb_c_ipk_) :: res

    type(psb_c_zspmat) :: ah
    type(psb_c_zvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    complex(c_double_complex), value :: alpha, beta
    character(c_char)      :: trans
    logical(c_bool), value :: doswap

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
    type(psb_zspmat_type), pointer :: ap
    character :: ftrans
    logical   :: fdoswap
    integer(psb_c_ipk_)   :: info

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

  end function psb_c_zspmm_opt
  

  function psb_c_zspsm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
    implicit none 
    integer(psb_c_ipk_) :: res

    type(psb_c_zspmat) :: ah
    type(psb_c_zvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    complex(c_double_complex), value :: alpha, beta
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
    type(psb_zspmat_type), pointer :: ap
    integer(psb_c_ipk_)               :: info

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

  end function psb_c_zspsm
  

end module psb_z_psblas_cbind_mod
