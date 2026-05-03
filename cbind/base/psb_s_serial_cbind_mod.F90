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

end module psb_s_serial_cbind_mod
