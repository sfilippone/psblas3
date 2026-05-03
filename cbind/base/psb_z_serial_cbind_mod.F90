module psb_z_serial_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_tools_cbind_mod


  interface
    module function psb_c_zvect_get_nrows(xh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zvector) :: xh
      
      type(psb_z_vect_type), pointer :: vp
      integer(psb_c_ipk_)               :: info
    end function psb_c_zvect_get_nrows
  end interface
  
  interface
    module function psb_c_zvect_f_get_cpy(v,xh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      complex(c_double_complex)    :: v(*)
      type(psb_c_zvector) :: xh
    end function psb_c_zvect_f_get_cpy
  end interface
  

  interface
    module function psb_c_zvect_zero(xh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zvector) :: xh
    end function psb_c_zvect_zero
  end interface
  
  interface
    module function psb_c_zvect_f_get_pnt(xh) bind(c) result(res)
      type(c_ptr)        :: res
      type(psb_c_zvector) :: xh
    end function psb_c_zvect_f_get_pnt
  end interface

  interface
    module function psb_c_zmat_get_nrows(mh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zspmat) :: mh
    end function psb_c_zmat_get_nrows
  end interface

  interface
    module function psb_c_zmat_get_ncols(mh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zspmat) :: mh
    end function psb_c_zmat_get_ncols
  end interface

  interface
    module function psb_c_zmat_name_print(mh,name) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zspmat) :: mh
      character(c_char)        :: name(*)
    end function psb_c_zmat_name_print
  end interface

  interface
    module function psb_c_zvect_set_scal(x,val) bind(c) result(info)
      type(psb_c_zvector) :: x
      integer(psb_c_ipk_) :: info
      complex(c_double_complex), value :: val
    end function psb_c_zvect_set_scal
  end interface

  interface
    module function psb_c_zvect_set_scal_bound(x,val,ifirst,ilast) bind(c) result(info)
      type(psb_c_zvector) :: x
      integer(psb_c_ipk_) :: info
      integer(psb_c_ipk_), value :: ifirst, ilast 
      complex(c_double_complex)    :: val
    end function psb_c_zvect_set_scal_bound
  end interface

  interface
    module function psb_c_zvect_set_vect(x,val,n) bind(c) result(info)
      type(psb_c_zvector) :: x
      integer(psb_c_ipk_) :: info
      integer(psb_c_ipk_), value :: n
      complex(c_double_complex)    :: val(*)
    end function psb_c_zvect_set_vect
  end interface

  interface
    module function psb_c_zvect_set_entry(x,index,val) bind(c) result(info)
      type(psb_c_zvector) :: x
      integer(psb_c_ipk_) :: info
      integer(psb_c_ipk_), value :: index
      complex(c_double_complex), value :: val
    end function psb_c_zvect_set_entry
  end interface
 
  interface
    module function psb_c_zvect_get_entry(x,index) bind(c) result(res)
      type(psb_c_zvector) :: x
      integer(psb_c_ipk_), value :: index
      complex(c_double_complex) :: res
    end function psb_c_zvect_get_entry
  end interface

  interface
    module function psb_c_zvect_clone(xh,yh) bind(c) result(info)
      integer(psb_c_ipk_) :: info
      type(psb_c_zvector) :: xh,yh
    end function psb_c_zvect_clone
  end interface

end module psb_z_serial_cbind_mod
