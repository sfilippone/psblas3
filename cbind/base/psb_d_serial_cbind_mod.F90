module psb_d_serial_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_tools_cbind_mod


  interface
    module function psb_c_dvect_get_nrows(xh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dvector) :: xh
      
      type(psb_d_vect_type), pointer :: vp
      integer(psb_c_ipk_)               :: info
    end function psb_c_dvect_get_nrows
  end interface
  
  interface
    module function psb_c_dvect_f_get_cpy(v,xh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      real(c_double)    :: v(*)
      type(psb_c_dvector) :: xh
    end function psb_c_dvect_f_get_cpy
  end interface
  

  interface
    module function psb_c_dvect_zero(xh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dvector) :: xh
    end function psb_c_dvect_zero
  end interface
  
  interface
    module function psb_c_dvect_f_get_pnt(xh) bind(c) result(res)
      type(c_ptr)        :: res
      type(psb_c_dvector) :: xh
    end function psb_c_dvect_f_get_pnt
  end interface

  interface
    module function psb_c_dmat_get_nrows(mh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dspmat) :: mh
    end function psb_c_dmat_get_nrows
  end interface

  interface
    module function psb_c_dmat_get_ncols(mh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dspmat) :: mh
    end function psb_c_dmat_get_ncols
  end interface

  interface
    module function psb_c_dmat_name_print(mh,name) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dspmat) :: mh
      character(c_char)        :: name(*)
    end function psb_c_dmat_name_print
  end interface

  interface
    module function psb_c_dvect_set_scal(x,val) bind(c) result(info)
      type(psb_c_dvector) :: x
      integer(psb_c_ipk_) :: info
      real(c_double), value :: val
    end function psb_c_dvect_set_scal
  end interface

  interface
    module function psb_c_dvect_set_scal_bound(x,val,ifirst,ilast) bind(c) result(info)
      type(psb_c_dvector) :: x
      integer(psb_c_ipk_) :: info
      integer(psb_c_ipk_), value :: ifirst, ilast 
      real(c_double)    :: val
    end function psb_c_dvect_set_scal_bound
  end interface

  interface
    module function psb_c_dvect_set_vect(x,val,n) bind(c) result(info)
      type(psb_c_dvector) :: x
      integer(psb_c_ipk_) :: info
      integer(psb_c_ipk_), value :: n
      real(c_double)    :: val(*)
    end function psb_c_dvect_set_vect
  end interface

  interface
    module function psb_c_dvect_set_entry(x,index,val) bind(c) result(info)
      type(psb_c_dvector) :: x
      integer(psb_c_ipk_) :: info
      integer(psb_c_ipk_), value :: index
      real(c_double), value :: val
    end function psb_c_dvect_set_entry
  end interface
 
  interface
    module function psb_c_dvect_get_entry(x,index) bind(c) result(res)
      type(psb_c_dvector) :: x
      integer(psb_c_ipk_), value :: index
      real(c_double) :: res
    end function psb_c_dvect_get_entry
  end interface

  interface
    module function psb_c_dvect_clone(xh,yh) bind(c) result(info)
      integer(psb_c_ipk_) :: info
      type(psb_c_dvector) :: xh,yh
    end function psb_c_dvect_clone
  end interface

end module psb_d_serial_cbind_mod
