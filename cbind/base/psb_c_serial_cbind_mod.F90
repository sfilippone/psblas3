module psb_c_serial_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_string_cbind_mod
  use psb_base_tools_cbind_mod

contains


  function psb_c_cvect_get_nrows(xh) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_cvector) :: xh

    type(psb_c_vect_type), pointer :: vp
    integer(psb_c_ipk_)               :: info

    res = -1

    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
      res = vp%get_nrows()
    end if

  end function psb_c_cvect_get_nrows

  function psb_c_cvect_f_get_cpy(v,xh) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_)    :: res
    complex(c_float_complex)    :: v(*)
    type(psb_c_cvector) :: xh

    type(psb_c_vect_type), pointer :: vp
    complex(psb_spk_), allocatable :: fv(:)
    integer(psb_c_ipk_)           :: info, sz

    res = -1

    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
      fv = vp%get_vect()
      sz = size(fv)
      v(1:sz) = fv(1:sz)
    end if

  end function psb_c_cvect_f_get_cpy


  function psb_c_cvect_zero(xh) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_)    :: res
    type(psb_c_cvector) :: xh

    type(psb_c_vect_type), pointer :: vp
    integer(psb_c_ipk_)               :: info

    res = -1

    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
      call vp%zero()
    end if

  end function psb_c_cvect_zero

  function psb_c_cvect_f_get_pnt(xh) bind(c) result(res)
    implicit none

    type(c_ptr)        :: res
    type(psb_c_cvector) :: xh

    type(psb_c_vect_type), pointer :: vp

    res = c_null_ptr

    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
      if(vp%is_dev()) call vp%sync()
      res = c_loc(vp%v%v)
    end if

  end function psb_c_cvect_f_get_pnt


  function psb_c_cmat_get_nrows(mh) bind(c) result(res)
    use psb_base_mod
    use psb_objhandle_mod
    use psb_base_string_cbind_mod
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_cspmat) :: mh
    type(psb_cspmat_type), pointer :: ap
    integer(psb_c_ipk_)               ::  info

    res = 0
    if (c_associated(mh%item)) then
      call c_f_pointer(mh%item,ap)
    else
      return
    end if

    res = ap%get_nrows()

  end function psb_c_cmat_get_nrows


  function psb_c_cmat_get_ncols(mh) bind(c) result(res)
    use psb_base_mod
    use psb_objhandle_mod
    use psb_base_string_cbind_mod
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_cspmat) :: mh
    type(psb_cspmat_type), pointer :: ap
    integer(psb_c_ipk_)               ::  info

    res = 0
    if (c_associated(mh%item)) then
      call c_f_pointer(mh%item,ap)
    else
      return
    end if

    res = ap%get_ncols()

  end function psb_c_cmat_get_ncols


  function psb_c_cmat_name_print(mh,name) bind(c) result(res)
    use psb_base_mod
    use psb_objhandle_mod
    use psb_base_string_cbind_mod
    implicit none
    integer(psb_c_ipk_) :: res

    character(c_char)        :: name(*)
    type(psb_c_cspmat) :: mh
    type(psb_cspmat_type), pointer :: ap
    integer(psb_c_ipk_)      ::  info
    character(1024)         :: fname

    res = 0
    if (c_associated(mh%item)) then
      call c_f_pointer(mh%item,ap)
    else
      return
    end if
    call stringc2f(name,fname)

    call ap%print(fname,head='PSBLAS Cbinding Interface')

  end function psb_c_cmat_name_print


end module psb_c_serial_cbind_mod
