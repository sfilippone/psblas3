module psb_d_comm_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod

  interface
    module function psb_c_dovrl(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dovrl
  end interface

  interface
    module function psb_c_dovrl_opt(xh,cdh,update,mode) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: update, mode
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dovrl_opt
  end interface


  interface
    module function psb_c_dhalo(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dhalo
  end interface

  interface
    module function psb_c_dhalo_opt(xh,cdh,tran,data,mode) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: data, mode
      character(c_char)      :: tran
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dhalo_opt
  end interface

  interface
    module function psb_c_dvscatter(ng,gx,xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      integer(psb_c_lpk_), value :: ng
      real(c_double), target :: gx(*)
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dvscatter
  end interface

  interface
    module function psb_c_dvgather_f(v,xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      real(c_double), target :: v(*)
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dvgather_f
  end interface

  interface
    module function psb_c_dspgather_f(gah,ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dspmat)   :: ah, gah
      type(psb_c_descriptor) :: cdh
    end function psb_c_dspgather_f
  end interface
  
end module psb_d_comm_cbind_mod
