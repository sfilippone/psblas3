module psb_z_serial_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_tools_cbind_mod

contains


  function psb_c_zvect_get_nrows(xh) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_zvector) :: xh

    type(psb_z_vect_type), pointer :: vp
    integer(psb_c_ipk_)               :: info

    res = -1

    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
      res = vp%get_nrows()
    end if

  end function psb_c_zvect_get_nrows

  function psb_c_zvect_f_get_cpy(v,xh) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_)    :: res
    complex(c_double_complex)    :: v(*)
    type(psb_c_zvector) :: xh

    type(psb_z_vect_type), pointer :: vp
    complex(psb_dpk_), allocatable :: fv(:)
    integer(psb_c_ipk_)           :: info, sz

    res = -1

    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
      fv = vp%get_vect()
      sz = size(fv)
      v(1:sz) = fv(1:sz)
    end if

  end function psb_c_zvect_f_get_cpy


  function psb_c_zvect_zero(xh) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_)    :: res
    type(psb_c_zvector) :: xh

    type(psb_z_vect_type), pointer :: vp
    integer(psb_c_ipk_)               :: info

    res = -1

    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
      call vp%zero()
    end if

  end function psb_c_zvect_zero

  function psb_c_zvect_f_get_pnt(xh) bind(c) result(res)
    implicit none

    type(c_ptr)        :: res
    type(psb_c_zvector) :: xh

    type(psb_z_vect_type), pointer :: vp

    res = c_null_ptr

    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
      if(vp%is_dev()) call vp%sync()
      res = c_loc(vp%v%v)
    end if

  end function psb_c_zvect_f_get_pnt


  function psb_c_zmat_get_nrows(mh) bind(c) result(res)
    use psb_base_mod
    use psb_objhandle_mod
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_zspmat) :: mh
    type(psb_zspmat_type), pointer :: ap
    integer(psb_c_ipk_)               ::  info

    res = 0
    if (c_associated(mh%item)) then
      call c_f_pointer(mh%item,ap)
    else
      return
    end if

    res = ap%get_nrows()

  end function psb_c_zmat_get_nrows


  function psb_c_zmat_get_ncols(mh) bind(c) result(res)
    use psb_base_mod
    use psb_objhandle_mod
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_zspmat) :: mh
    type(psb_zspmat_type), pointer :: ap
    integer(psb_c_ipk_)               ::  info

    res = 0
    if (c_associated(mh%item)) then
      call c_f_pointer(mh%item,ap)
    else
      return
    end if

    res = ap%get_ncols()

  end function psb_c_zmat_get_ncols


  function psb_c_zmat_name_print(mh,name) bind(c) result(res)
    use psb_base_mod
    use psb_objhandle_mod
    implicit none
    integer(psb_c_ipk_) :: res

    character(c_char)        :: name(*)
    type(psb_c_zspmat) :: mh
    type(psb_zspmat_type), pointer :: ap
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

  end function psb_c_zmat_name_print

  function psb_c_zvect_set_scal(x,val) bind(c) result(info)
    use psb_base_mod
    implicit none

    type(psb_c_zvector) :: x
    type(psb_z_vect_type), pointer :: xp
    integer(psb_c_ipk_) :: info
    complex(c_double_complex), value :: val

    info = -1;

    if (c_associated(x%item)) then
      call c_f_pointer(x%item,xp)
    else
      return
    end if

    call xp%set(val)

    info = 0

  end function psb_c_zvect_set_scal

  function psb_c_zvect_set_scal_bound(x,val,ifirst,ilast) bind(c) result(info)
    use psb_base_mod
    implicit none

    type(psb_c_zvector) :: x
    type(psb_z_vect_type), pointer :: xp
    integer(psb_c_ipk_) :: info
    integer(psb_c_ipk_), value :: ifirst, ilast 
    complex(c_double_complex)    :: val

    info = -1;

    if (c_associated(x%item)) then
      call c_f_pointer(x%item,xp)
    else
      return
    end if

    call xp%set(val,first=ifirst,last=ilast)

    info = 0

  end function psb_c_zvect_set_scal_bound

  function psb_c_zvect_set_vect(x,val,n) bind(c) result(info)
    use psb_base_mod
    implicit none

    type(psb_c_zvector) :: x
    type(psb_z_vect_type), pointer :: xp
    integer(psb_c_ipk_) :: info
    integer(psb_c_ipk_), value :: n
    complex(c_double_complex)    :: val(*)

    info = -1;

    if (c_associated(x%item)) then
      call c_f_pointer(x%item,xp)
    else
      return
    end if

    call xp%set(val(1:n))

    info = 0

  end function psb_c_zvect_set_vect

  function psb_c_zvect_set_entry(x,index,val) bind(c) result(info)
    use psb_base_mod
    implicit none
    type(psb_c_zvector) :: x
    type(psb_z_vect_type), pointer :: xp
    integer(psb_c_ipk_) :: info
    integer(psb_c_ipk_), value :: index
    complex(c_double_complex), value :: val
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
    
  end function psb_c_zvect_set_entry
  
  
  function psb_c_zvect_get_entry(x,index) bind(c) result(res)
    use psb_base_mod
    implicit none
    
    type(psb_c_zvector) :: x
    type(psb_z_vect_type), pointer :: xp
    integer(psb_c_ipk_), value :: index
    complex(c_double_complex) :: res
    integer(psb_c_ipk_) :: ixb

    if (c_associated(x%item)) then
      call c_f_pointer(x%item,xp)
    else
      return
    end if

    ixb = psb_c_get_index_base()
    res = xp%get_entry((index+(1-ixb)))
  end function psb_c_zvect_get_entry

  function psb_c_zvect_clone(xh,yh) bind(c) result(info)
    implicit none

    integer(psb_c_ipk_) :: info
    type(psb_c_zvector) :: xh,yh

    type(psb_z_vect_type), pointer :: xp,yp

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
    
  end function psb_c_zvect_clone

  function psb_c_znnz(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res
    
    type(psb_c_zspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_zspmat_type), pointer  :: ap
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

  end function psb_c_znnz

  function psb_c_zis_matupd(ah,cdh) bind(c) result(res)
    implicit none
    logical(c_bool) :: res

    type(psb_c_zspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_zspmat_type), pointer  :: ap
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

  function psb_c_zis_matasb(ah,cdh) bind(c) result(res)
    implicit none
    logical(c_bool) :: res

    type(psb_c_zspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_zspmat_type), pointer  :: ap
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

  function psb_c_zis_matbld(ah,cdh) bind(c) result(res)
    implicit none
    logical(c_bool) :: res

    type(psb_c_zspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_zspmat_type), pointer  :: ap
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

  function psb_c_zset_matupd(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_zspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_zspmat_type), pointer  :: ap

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

  function psb_c_zset_matasb(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_zspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_zspmat_type), pointer  :: ap

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

  function psb_c_zset_matbld(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_zspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_zspmat_type), pointer  :: ap

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

  function psb_c_zcopy_mat(ah,bh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_zspmat)   :: ah,bh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_zspmat_type), pointer  :: ap,bp
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
  
end module psb_z_serial_cbind_mod
