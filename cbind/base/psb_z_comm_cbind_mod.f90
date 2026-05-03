module psb_z_comm_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod

  interface
    module function psb_c_zovrl(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zovrl
  end interface

  interface
    module function psb_c_zovrl_opt(xh,cdh,update,mode) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: update, mode
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zovrl_opt
  end interface


  interface
    module function psb_c_zhalo(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zhalo
  end interface

  interface
    module function psb_c_zhalo_opt(xh,cdh,tran,data,mode) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: data, mode
      character(c_char)      :: tran
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zhalo_opt
  end interface

  interface
    module function psb_c_zvscatter(ng,gx,xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      integer(psb_c_lpk_), value :: ng
      complex(c_double_complex), target :: gx(*)
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zvscatter
  end interface

  interface
    module function psb_c_zvgather_f(v,xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      complex(c_double_complex), target :: v(*)
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zvgather_f
  end interface

  interface
    module function psb_c_zspgather_f(gah,ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zspmat)   :: ah, gah
      type(psb_c_descriptor) :: cdh
    end function psb_c_zspgather_f
  end interface
  
end module psb_z_comm_cbind_mod
