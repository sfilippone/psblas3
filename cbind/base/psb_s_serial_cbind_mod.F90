module psb_s_serial_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_string_cbind_mod
  use psb_base_tools_cbind_mod

contains


  function psb_c_svect_get_nrows(xh) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_svector) :: xh

    type(psb_s_vect_type), pointer :: vp
    integer(psb_c_ipk_)               :: info

    res = -1

    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
      res = vp%get_nrows()
    end if

  end function psb_c_svect_get_nrows

  function psb_c_svect_f_get_cpy(v,xh) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_)    :: res
    real(c_float)    :: v(*)
    type(psb_c_svector) :: xh

    type(psb_s_vect_type), pointer :: vp
    real(psb_spk_), allocatable :: fv(:)
    integer(psb_c_ipk_)           :: info, sz

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

    integer(psb_c_ipk_)    :: res
    type(psb_c_svector) :: xh

    type(psb_s_vect_type), pointer :: vp
    integer(psb_c_ipk_)               :: info

    res = -1

    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
      call vp%zero()
    end if

  end function psb_c_svect_zero

  function psb_c_svect_f_get_pnt(xh) bind(c) result(res)
    implicit none

    type(c_ptr)        :: res
    type(psb_c_svector) :: xh

    type(psb_s_vect_type), pointer :: vp

    res = c_null_ptr

    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
      if(vp%is_dev()) call vp%sync()
      res = c_loc(vp%v%v)
    end if

  end function psb_c_svect_f_get_pnt


  function psb_c_smat_get_nrows(mh) bind(c) result(res)
    use psb_base_mod
    use psb_objhandle_mod
    use psb_base_string_cbind_mod
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_sspmat) :: mh
    type(psb_sspmat_type), pointer :: ap
    integer(psb_c_ipk_)               ::  info

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
    integer(psb_c_ipk_) :: res

    type(psb_c_sspmat) :: mh
    type(psb_sspmat_type), pointer :: ap
    integer(psb_c_ipk_)               ::  info

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
    integer(psb_c_ipk_) :: res

    character(c_char)        :: name(*)
    type(psb_c_sspmat) :: mh
    type(psb_sspmat_type), pointer :: ap
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

  end function psb_c_smat_name_print

  function psb_c_svect_set_scal(x,val) result(info)
    use psb_base_mod
    implicit none

    type(psb_c_svector) :: x
    type(psb_s_vect_type), pointer :: xp
    integer(psb_c_ipk_) :: info
    real(c_float), value :: val

    info = -1;

    if (c_associated(x%item)) then
      call c_f_pointer(x%item,xp)
    else
      return
    end if

    call xp%set(val)

    info = 0

  end function psb_c_svect_set_scal

  function psb_c_svect_set_vect(x,val,n) result(info)
    use psb_base_mod
    implicit none

    type(psb_c_svector) :: x
    type(psb_s_vect_type), pointer :: xp
    integer(psb_c_ipk_) :: info
    integer(psb_c_ipk_), value :: n
    real(c_float)    :: val(*)

    info = -1;

    if (c_associated(x%item)) then
      call c_f_pointer(x%item,xp)
    else
      return
    end if

    call xp%set(val(1:n))

    info = 0

  end function psb_c_svect_set_vect


end module psb_s_serial_cbind_mod
