module psb_s_serial_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_string_cbind_mod
  use psb_base_tools_cbind_mod

contains

  
  function psb_c_svect_get_nrows(xh) bind(c) result(res)
    implicit none 

    integer(psb_c_int) :: res   
    type(psb_c_svector) :: xh

    type(psb_s_vect_type), pointer :: vp
    integer               :: info

    res = -1

    if (c_associated(xh%item)) then 
      call c_f_pointer(xh%item,vp)
      res = vp%get_nrows()
    end if

  end function psb_c_svect_get_nrows
  
  function psb_c_svect_f_get_cpy(v,xh) bind(c) result(res)
    implicit none 

    integer(psb_c_int)    :: res   
    real(c_float)    :: v(*)
    type(psb_c_svector) :: xh
    
    type(psb_s_vect_type), pointer :: vp
    real(psb_spk_), allocatable :: fv(:)
    integer               :: info, sz

    res = -1

    if (c_associated(xh%item)) then 
      call c_f_pointer(xh%item,vp)
      fv = vp%get_vect()
      sz = size(fv)
      v(1:sz) = fv(1:sz)
    end if

  end function psb_c_svect_f_get_cpy

  
  function psb_c_svect_zero(xh) bind(c) result(res)
    implicit none 

    integer(psb_c_int)    :: res   
    type(psb_c_svector) :: xh
    
    type(psb_s_vect_type), pointer :: vp
    integer               :: info

    res = -1

    if (c_associated(xh%item)) then 
      call c_f_pointer(xh%item,vp)
      call vp%zero()
    end if

  end function psb_c_svect_zero

  
  function psb_c_smat_get_nrows(mh) bind(c) result(res)
    use psb_base_mod
    use psb_objhandle_mod
    use psb_base_string_cbind_mod
    implicit none 
    integer(psb_c_int) :: res

    type(psb_c_sspmat) :: mh
    type(psb_sspmat_type), pointer :: ap
    integer                 ::  info

    res = 0
    if (c_associated(mh%item)) then 
      call c_f_pointer(mh%item,ap)
    else
      return 
    end if
    
    res = ap%get_nrows()

  end function psb_c_smat_get_nrows

  
  function psb_c_smat_get_ncols(mh) bind(c) result(res)
    use psb_base_mod
    use psb_objhandle_mod
    use psb_base_string_cbind_mod
    implicit none 
    integer(psb_c_int) :: res

    type(psb_c_sspmat) :: mh
    type(psb_sspmat_type), pointer :: ap
    integer                 ::  info

    res = 0
    if (c_associated(mh%item)) then 
      call c_f_pointer(mh%item,ap)
    else
      return 
    end if
    
    res = ap%get_ncols()

  end function psb_c_smat_get_ncols

  
  function psb_c_smat_name_print(mh,name) bind(c) result(res)
    use psb_base_mod
    use psb_objhandle_mod
    use psb_base_string_cbind_mod
    implicit none 
    integer(psb_c_int) :: res

    character(c_char)        :: name(*)
    type(psb_c_sspmat) :: mh
    type(psb_sspmat_type), pointer :: ap
    integer                 ::  info
    character(1024)         :: fname

    res = 0
    if (c_associated(mh%item)) then 
      call c_f_pointer(mh%item,ap)
    else
      return 
    end if
    call stringc2f(name,fname)
    
    call ap%print(fname,head='PSBLAS Cbinding Interface')

  end function psb_c_smat_name_print
  

end module psb_s_serial_cbind_mod

