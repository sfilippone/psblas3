module psb_s_comm_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod

  interface
    module function psb_c_sovrl(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sovrl
  end interface

  interface
    module function psb_c_sovrl_opt(xh,cdh,update,mode) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: update, mode
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sovrl_opt
  end interface


  interface
    module function psb_c_shalo(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_shalo
  end interface

  interface
    module function psb_c_shalo_opt(xh,cdh,tran,data,mode) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: data, mode
      character(c_char)      :: tran
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_shalo_opt
  end interface

  interface
    module function psb_c_svscatter(ng,gx,xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      integer(psb_c_lpk_), value :: ng
      real(c_float), target :: gx(*)
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_svscatter
  end interface

  interface
    module function psb_c_svgather_f(v,xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      real(c_float), target :: v(*)
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_svgather_f
  end interface

  interface
    module function psb_c_sspgather_f(gah,ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_sspmat)   :: ah, gah
      type(psb_c_descriptor) :: cdh
    end function psb_c_sspgather_f
  end interface
  
end module psb_s_comm_cbind_mod
