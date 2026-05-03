module psb_c_comm_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod

  interface
    module function psb_c_covrl(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_covrl
  end interface

  interface
    module function psb_c_covrl_opt(xh,cdh,update,mode) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: update, mode
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_covrl_opt
  end interface


  interface
    module function psb_c_chalo(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_chalo
  end interface

  interface
    module function psb_c_chalo_opt(xh,cdh,tran,data,mode) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: data, mode
      character(c_char)      :: tran
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_chalo_opt
  end interface

  interface
    module function psb_c_cvscatter(ng,gx,xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      integer(psb_c_lpk_), value :: ng
      complex(c_float_complex), target :: gx(*)
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cvscatter
  end interface

  interface
    module function psb_c_cvgather_f(v,xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      complex(c_float_complex), target :: v(*)
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cvgather_f
  end interface

  interface
    module function psb_c_cspgather_f(gah,ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cspmat)   :: ah, gah
      type(psb_c_descriptor) :: cdh
    end function psb_c_cspgather_f
  end interface
  
end module psb_c_comm_cbind_mod
