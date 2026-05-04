submodule (psb_d_serial_cbind_mod) psb_d_serial_cbind_impl
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_tools_cbind_mod

contains


  module function psb_c_dvect_get_nrows(xh) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_dvector) :: xh

    type(psb_d_vect_type), pointer :: vp
    integer(psb_c_ipk_)               :: info

    res = -1

    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
      res = vp%get_nrows()
    end if

  end function psb_c_dvect_get_nrows

  module function psb_c_dvect_f_get_cpy(v,xh) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_)    :: res
    real(c_double)    :: v(*)
    type(psb_c_dvector) :: xh

    type(psb_d_vect_type), pointer :: vp
    real(psb_dpk_), allocatable :: fv(:)
    integer(psb_c_ipk_)           :: info, sz

    res = -1

    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
      fv = vp%get_vect()
      sz = size(fv)
      v(1:sz) = fv(1:sz)
    end if

  end function psb_c_dvect_f_get_cpy


  module function psb_c_dvect_zero(xh) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_)    :: res
    type(psb_c_dvector) :: xh

    type(psb_d_vect_type), pointer :: vp
    integer(psb_c_ipk_)               :: info

    res = -1

    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
      call vp%zero()
    end if

  end function psb_c_dvect_zero

  module function psb_c_dvect_f_get_pnt(xh) bind(c) result(res)
    implicit none

    type(c_ptr)        :: res
    type(psb_c_dvector) :: xh

    type(psb_d_vect_type), pointer :: vp

    res = c_null_ptr

    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
      if(vp%is_dev()) call vp%sync()
      res = c_loc(vp%v%v)
    end if

  end function psb_c_dvect_f_get_pnt


  module function psb_c_dmat_get_nrows(mh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_dspmat) :: mh
    type(psb_dspmat_type), pointer :: ap
    integer(psb_c_ipk_)               ::  info

    res = 0
    if (c_associated(mh%item)) then
      call c_f_pointer(mh%item,ap)
    else
      return
    end if

    res = ap%get_nrows()

  end function psb_c_dmat_get_nrows


  module function psb_c_dmat_get_ncols(mh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_dspmat) :: mh
    type(psb_dspmat_type), pointer :: ap
    integer(psb_c_ipk_)               ::  info

    res = 0
    if (c_associated(mh%item)) then
      call c_f_pointer(mh%item,ap)
    else
      return
    end if

    res = ap%get_ncols()

  end function psb_c_dmat_get_ncols

  module function psb_c_dmat_name_print(mh,name) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res
    character(c_char)        :: name(*)

    type(psb_c_dspmat) :: mh
    type(psb_dspmat_type), pointer :: ap
    integer(psb_c_ipk_)      ::  info
    character(1024)         :: fname

    res = 0
    if (c_associated(mh%item)) then
      call c_f_pointer(mh%item,ap)
    else
      return
    end if
    call psb_stringc2f(name,fname)

    call ap%print(fname,head='PSBLAS Cbinding Interface')

  end function psb_c_dmat_name_print

  module function psb_c_dvect_set_scal(x,val) bind(c) result(info)
    implicit none

    type(psb_c_dvector) :: x
    type(psb_d_vect_type), pointer :: xp
    integer(psb_c_ipk_) :: info
    real(c_double), value :: val

    info = -1;

    if (c_associated(x%item)) then
      call c_f_pointer(x%item,xp)
    else
      return
    end if

    call xp%set(val)

    info = 0

  end function psb_c_dvect_set_scal

  module function psb_c_dvect_set_scal_bound(x,val,ifirst,ilast) bind(c) result(info)
    implicit none

    type(psb_c_dvector) :: x
    type(psb_d_vect_type), pointer :: xp
    integer(psb_c_ipk_) :: info
    integer(psb_c_ipk_), value :: ifirst, ilast 
    real(c_double)    :: val

    info = -1;

    if (c_associated(x%item)) then
      call c_f_pointer(x%item,xp)
    else
      return
    end if

    call xp%set(val,first=ifirst,last=ilast)

    info = 0

  end function psb_c_dvect_set_scal_bound

  module function psb_c_dvect_set_vect(x,val,n) bind(c) result(info)
    implicit none
    type(psb_c_dvector) :: x
    type(psb_d_vect_type), pointer :: xp
    integer(psb_c_ipk_) :: info
    integer(psb_c_ipk_), value :: n
    real(c_double)    :: val(*)

    info = -1;

    if (c_associated(x%item)) then
      call c_f_pointer(x%item,xp)
    else
      return
    end if

    call xp%set(val(1:n))

    info = 0

  end function psb_c_dvect_set_vect

  module function psb_c_dvect_set_entry(x,index,val) bind(c) result(info)
    implicit none

    type(psb_c_dvector) :: x
    type(psb_d_vect_type), pointer :: xp
    integer(psb_c_ipk_) :: info
    integer(psb_c_ipk_), value :: index
    real(c_double), value :: val
    integer(psb_c_ipk_) :: ixb

    info = -1;

    if (c_associated(x%item)) then
      call c_f_pointer(x%item,xp)
    else
      return
    end if

    ixb = psb_c_get_index_base()
    call xp%set_entry((index+(1-ixb)),val)
    info = 0

  end function psb_c_dvect_set_entry
 
  module function psb_c_dvect_get_entry(x,index) bind(c) result(res)
    implicit none

    type(psb_c_dvector) :: x
    type(psb_d_vect_type), pointer :: xp
    integer(psb_c_ipk_), value :: index
    real(c_double) :: res
    integer(psb_c_ipk_) :: ixb
    
    if (c_associated(x%item)) then
      call c_f_pointer(x%item,xp)
    else
      return
    end if

    ixb = psb_c_get_index_base()
    res = xp%get_entry((index+(1-ixb)))
  end function psb_c_dvect_get_entry

  module function psb_c_dvect_clone(xh,yh) bind(c) result(info)
    implicit none

    integer(psb_c_ipk_) :: info
    type(psb_c_dvector) :: xh,yh

    type(psb_d_vect_type), pointer :: xp,yp

    info = -1

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
    call xp%clone(yp,info)
    
  end function psb_c_dvect_clone

end submodule psb_d_serial_cbind_impl
