module psb_z_tools_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_cpenv_mod
  use psb_objhandle_mod
  use psb_base_tools_cbind_mod
#ifdef PSB_HAVE_CUDA
  use psb_cuda_mod
#endif
  
  ! Should define   geall_opt  with DUPL argument
  interface
    module function psb_c_zgeall(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgeall
  end interface

  interface
    module function psb_c_zgeall_remote(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
      
    end function psb_c_zgeall_remote
  end interface

  interface
    module function psb_c_zgeall_remote_options(xh,cdh,bldmode,dupl) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
      integer(psb_c_ipk_), value :: dupl
      integer(psb_c_ipk_), value :: bldmode
    end function psb_c_zgeall_remote_options
  end interface

  interface
    module function psb_c_zgeasb(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgeasb
  end interface
  
  interface
    module function psb_c_zgeasb_options(xh,cdh,dupl) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
      integer(psb_c_ipk_), value :: dupl
    end function psb_c_zgeasb_options
  end interface
  
  interface
    module function psb_c_zgeasb_options_format(xh,cdh,dupl,format) bind(c) result(res)
      ! Takes into account format argument as a c string, and uses it to call the appropriate psb_geasb
      ! with mold argument
      integer(psb_c_ipk_) :: res
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
      character(kind=c_char), dimension(*) :: format
      integer(psb_c_ipk_), value :: dupl
    end function psb_c_zgeasb_options_format
  end interface
   
  interface
    module function psb_c_zgefree(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgefree
  end interface

  interface
    module function psb_c_zgeins(nz,irw,val,xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: nz
      integer(psb_c_lpk_)        :: irw(*)
      complex(c_double_complex)        :: val(*)
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgeins
  end interface

  interface
    module function psb_c_zspall(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zspall
  end interface
  
  interface
    module function psb_c_zspall_remote(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zspall_remote
  end interface

  interface
    module function psb_c_zspasb(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zspasb
  end interface

  interface
    module function psb_c_zspfree(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zspfree
  end interface

  interface
    module function psb_c_zspasb_opt(mh,cdh,afmt,upd,dupl) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zspmat) :: mh
      type(psb_c_descriptor) :: cdh
      integer(psb_c_ipk_), value :: upd,dupl
      character(c_char)     :: afmt(*)
    end function psb_c_zspasb_opt
  end interface

  interface
    module function psb_c_zspins(nz,irw,icl,val,mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: nz
      integer(psb_c_lpk_)      :: irw(*), icl(*)
      complex(c_double_complex)        :: val(*)
      type(psb_c_zspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zspins
  end interface

  interface
    module function psb_c_zsprn(mh,cdh,clear) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      logical(c_bool), value :: clear
      type(psb_c_zspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zsprn
  end interface
!!$
!!$  module function psb_c_zspprint(mh) bind(c) result(res)
!!$
!!$    implicit none
!!$    integer(psb_c_ipk_) :: res
!!$    integer(psb_c_ipk_),  value :: mh
!!$    integer(psb_c_ipk_)         :: info
!!$
!!$
!!$    res = -1
!!$    call psb_check_double_spmat_handle(mh,info)
!!$    if (info < 0) return
!!$
!!$    call psb_csprt(0,double_spmat_pool(mh)%item,head='Debug mat')
!!$
!!$    res = 0
!!$
!!$    return
!!$  end function psb_c_zspprint

  interface 
    module function psb_c_zgetelem(xh,index,cdh) bind(c) result(res)
      type(psb_c_zvector)      :: xh
      integer(psb_c_lpk_), value :: index
      type(psb_c_descriptor)     :: cdh
      complex(c_double_complex)           :: res
    end function psb_c_zgetelem
  end interface
  
  interface 
    module function psb_c_zsetelem(index,val,xh,cdh) bind(c) result(res)
      type(psb_c_zvector)      :: xh
      integer(psb_c_lpk_), value :: index
      type(psb_c_descriptor)     :: cdh
      complex(c_double_complex), value    :: val
      integer(psb_c_ipk_) :: res
    end function psb_c_zsetelem
  end interface

  interface 
    module function psb_c_zmatgetelem(ah,rowindex,colindex,cdh) bind(c) result(res)
      type(psb_c_zspmat)      :: ah
      integer(psb_c_lpk_), value :: rowindex, colindex
      type(psb_c_descriptor)     :: cdh
      complex(c_double_complex)           :: res
    end function psb_c_zmatgetelem
  end interface

end module psb_z_tools_cbind_mod
