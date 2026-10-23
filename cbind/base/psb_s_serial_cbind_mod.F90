module psb_s_serial_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_tools_cbind_mod


  interface
    module function psb_c_svect_get_nrows(xh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_svector) :: xh
      
      type(psb_s_vect_type), pointer :: vp
      integer(psb_c_ipk_)               :: info
    end function psb_c_svect_get_nrows
  end interface
  
  interface
    module function psb_c_svect_f_get_cpy(v,xh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      real(c_float)    :: v(*)
      type(psb_c_svector) :: xh
    end function psb_c_svect_f_get_cpy
  end interface
  

  interface
    module function psb_c_svect_zero(xh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_svector) :: xh
    end function psb_c_svect_zero
  end interface
  
  interface
    module function psb_c_svect_f_get_pnt(xh) bind(c) result(res)
      type(c_ptr)        :: res
      type(psb_c_svector) :: xh
    end function psb_c_svect_f_get_pnt
  end interface

  interface
    module function psb_c_smat_get_nrows(mh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_sspmat) :: mh
    end function psb_c_smat_get_nrows
  end interface

  interface
    module function psb_c_smat_get_ncols(mh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_sspmat) :: mh
    end function psb_c_smat_get_ncols
  end interface

  interface
    module function psb_c_smat_name_print(mh,name) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_sspmat) :: mh
      character(c_char)        :: name(*)
    end function psb_c_smat_name_print
  end interface

  interface
    module function psb_c_svect_set_scal(x,val) bind(c) result(info)
      type(psb_c_svector) :: x
      integer(psb_c_ipk_) :: info
      real(c_float), value :: val
    end function psb_c_svect_set_scal
  end interface

  interface
    module function psb_c_svect_set_scal_bound(x,val,ifirst,ilast) bind(c) result(info)
      type(psb_c_svector) :: x
      integer(psb_c_ipk_) :: info
      integer(psb_c_ipk_), value :: ifirst, ilast 
      real(c_float)    :: val
    end function psb_c_svect_set_scal_bound
  end interface

  interface
    module function psb_c_svect_set_vect(x,val,n) bind(c) result(info)
      type(psb_c_svector) :: x
      integer(psb_c_ipk_) :: info
      integer(psb_c_ipk_), value :: n
      real(c_float)    :: val(*)
    end function psb_c_svect_set_vect
  end interface

  interface
    module function psb_c_svect_set_entry(x,index,val) bind(c) result(info)
      type(psb_c_svector) :: x
      integer(psb_c_ipk_) :: info
      integer(psb_c_ipk_), value :: index
      real(c_float), value :: val
    end function psb_c_svect_set_entry
  end interface
 
  interface
    module function psb_c_svect_get_entry(x,index) bind(c) result(res)
      type(psb_c_svector) :: x
      integer(psb_c_ipk_), value :: index
      real(c_float) :: res
    end function psb_c_svect_get_entry
  end interface

  interface
    module function psb_c_svect_clone(xh,yh) bind(c) result(info)
      integer(psb_c_ipk_) :: info
      type(psb_c_svector) :: xh,yh
    end function psb_c_svect_clone
  end interface

    function psb_c_snnz(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_sspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_sspmat_type), pointer  :: ap
    integer(psb_c_ipk_)               :: info

    res = 0

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

    res = psb_nnz(ap,descp,info)

  end function psb_c_snnz

  function psb_c_sis_matupd(ah,cdh) bind(c) result(res)
    implicit none
    logical(c_bool) :: res

    type(psb_c_sspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_sspmat_type), pointer  :: ap
    integer(psb_c_ipk_)               :: info

    res = .false.

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

    res = ap%is_upd()
  end function

  function psb_c_sis_matasb(ah,cdh) bind(c) result(res)
    implicit none
    logical(c_bool) :: res

    type(psb_c_sspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_sspmat_type), pointer  :: ap
    integer(psb_c_ipk_)               :: info

    res = .false.

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

    res = ap%is_asb()
  end function

  function psb_c_sis_matbld(ah,cdh) bind(c) result(res)
    implicit none
    logical(c_bool) :: res

    type(psb_c_sspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_sspmat_type), pointer  :: ap
    integer(psb_c_ipk_)               :: info

    res = .false.

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

    res = ap%is_bld()
  end function

  function psb_c_sset_matupd(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_sspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_sspmat_type), pointer  :: ap

    res = -1

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

    call ap%set_upd()

    res = psb_success_
  end function

  function psb_c_sset_matasb(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_sspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_sspmat_type), pointer  :: ap

    res = -1;

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

    call ap%set_asb()

    res = psb_success_

  end function

  function psb_c_sset_matbld(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_sspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_sspmat_type), pointer  :: ap

    res = -1

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

    call ap%set_bld()

    res = psb_success_
  end function

  function psb_c_scopy_mat(ah,bh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_sspmat)   :: ah,bh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_sspmat_type), pointer  :: ap,bp
    integer(psb_c_ipk_)               :: info

    res = -1

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
    if (c_associated(bh%item)) then
      call c_f_pointer(bh%item,bp)
    else
      return
    end if

    call ap%clone(bp,info)

    res = info
  end function
  
end module psb_s_serial_cbind_mod
