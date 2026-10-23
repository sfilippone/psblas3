module psb_c_serial_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_tools_cbind_mod


  interface
    module function psb_c_cvect_get_nrows(xh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cvector) :: xh
      
      type(psb_c_vect_type), pointer :: vp
      integer(psb_c_ipk_)               :: info
    end function psb_c_cvect_get_nrows
  end interface
  
  interface
    module function psb_c_cvect_f_get_cpy(v,xh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      complex(c_float_complex)    :: v(*)
      type(psb_c_cvector) :: xh
    end function psb_c_cvect_f_get_cpy
  end interface
  

  interface
    module function psb_c_cvect_zero(xh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cvector) :: xh
    end function psb_c_cvect_zero
  end interface
  
  interface
    module function psb_c_cvect_f_get_pnt(xh) bind(c) result(res)
      type(c_ptr)        :: res
      type(psb_c_cvector) :: xh
    end function psb_c_cvect_f_get_pnt
  end interface

  interface
    module function psb_c_cmat_get_nrows(mh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cspmat) :: mh
    end function psb_c_cmat_get_nrows
  end interface

  interface
    module function psb_c_cmat_get_ncols(mh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cspmat) :: mh
    end function psb_c_cmat_get_ncols
  end interface

  interface
    module function psb_c_cmat_name_print(mh,name) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cspmat) :: mh
      character(c_char)        :: name(*)
    end function psb_c_cmat_name_print
  end interface

  interface
    module function psb_c_cvect_set_scal(x,val) bind(c) result(info)
      type(psb_c_cvector) :: x
      integer(psb_c_ipk_) :: info
      complex(c_float_complex), value :: val
    end function psb_c_cvect_set_scal
  end interface

  interface
    module function psb_c_cvect_set_scal_bound(x,val,ifirst,ilast) bind(c) result(info)
      type(psb_c_cvector) :: x
      integer(psb_c_ipk_) :: info
      integer(psb_c_ipk_), value :: ifirst, ilast 
      complex(c_float_complex)    :: val
    end function psb_c_cvect_set_scal_bound
  end interface

  interface
    module function psb_c_cvect_set_vect(x,val,n) bind(c) result(info)
      type(psb_c_cvector) :: x
      integer(psb_c_ipk_) :: info
      integer(psb_c_ipk_), value :: n
      complex(c_float_complex)    :: val(*)
    end function psb_c_cvect_set_vect
  end interface

  interface
    module function psb_c_cvect_set_entry(x,index,val) bind(c) result(info)
      type(psb_c_cvector) :: x
      integer(psb_c_ipk_) :: info
      integer(psb_c_ipk_), value :: index
      complex(c_float_complex), value :: val
    end function psb_c_cvect_set_entry
  end interface
 
  interface
    module function psb_c_cvect_get_entry(x,index) bind(c) result(res)
      type(psb_c_cvector) :: x
      integer(psb_c_ipk_), value :: index
      complex(c_float_complex) :: res
    end function psb_c_cvect_get_entry
  end interface

  interface
    module function psb_c_cvect_clone(xh,yh) bind(c) result(info)
      integer(psb_c_ipk_) :: info
      type(psb_c_cvector) :: xh,yh
    end function psb_c_cvect_clone
  end interface

    function psb_c_cnnz(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_cspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_cspmat_type), pointer  :: ap
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

  end function psb_c_cnnz

  function psb_c_cis_matupd(ah,cdh) bind(c) result(res)
    implicit none
    logical(c_bool) :: res

    type(psb_c_cspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_cspmat_type), pointer  :: ap
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

  function psb_c_cis_matasb(ah,cdh) bind(c) result(res)
    implicit none
    logical(c_bool) :: res

    type(psb_c_cspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_cspmat_type), pointer  :: ap
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

  function psb_c_cis_matbld(ah,cdh) bind(c) result(res)
    implicit none
    logical(c_bool) :: res

    type(psb_c_cspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_cspmat_type), pointer  :: ap
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

  function psb_c_cset_matupd(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_cspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_cspmat_type), pointer  :: ap

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

  function psb_c_cset_matasb(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_cspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_cspmat_type), pointer  :: ap

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

  function psb_c_cset_matbld(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_cspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_cspmat_type), pointer  :: ap

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

  function psb_c_ccopy_mat(ah,bh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_cspmat)   :: ah,bh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_cspmat_type), pointer  :: ap,bp
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
  
end module psb_c_serial_cbind_mod
