module psb_c_tools_cbind_mod
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
    module function psb_c_cgeall(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgeall
  end interface

  interface
    module function psb_c_cgeall_remote(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
      
    end function psb_c_cgeall_remote
  end interface

  interface
    module function psb_c_cgeall_remote_options(xh,cdh,bldmode,dupl) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
      integer(psb_c_ipk_), value :: dupl
      integer(psb_c_ipk_), value :: bldmode
    end function psb_c_cgeall_remote_options
  end interface

  interface
    module function psb_c_cgeasb(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgeasb
  end interface
  
  interface
    module function psb_c_cgeasb_options(xh,cdh,dupl) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
      integer(psb_c_ipk_), value :: dupl
    end function psb_c_cgeasb_options
  end interface
  
  interface
    module function psb_c_cgeasb_options_format(xh,cdh,dupl,format) bind(c) result(res)
      ! Takes into account format argument as a c string, and uses it to call the appropriate psb_geasb
      ! with mold argument
      integer(psb_c_ipk_) :: res
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
      character(kind=c_char), dimension(*) :: format
      integer(psb_c_ipk_), value :: dupl
    end function psb_c_cgeasb_options_format
  end interface
   
  interface
    module function psb_c_cgefree(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgefree
  end interface

  interface
    module function psb_c_cgeins(nz,irw,val,xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: nz
      integer(psb_c_lpk_)        :: irw(*)
      complex(c_float_complex)        :: val(*)
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgeins
  end interface

  interface
    module function psb_c_cspall(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cspall
  end interface
  
  interface
    module function psb_c_cspall_remote(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cspall_remote
  end interface

  interface
    module function psb_c_cspasb(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cspasb
  end interface

  interface
    module function psb_c_cspfree(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cspfree
  end interface

  interface
    module function psb_c_cspasb_opt(mh,cdh,afmt,upd,dupl) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cspmat) :: mh
      type(psb_c_descriptor) :: cdh
      integer(psb_c_ipk_), value :: upd,dupl
      character(c_char)     :: afmt(*)
    end function psb_c_cspasb_opt
  end interface

  interface
    module function psb_c_cspins(nz,irw,icl,val,mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: nz
      integer(psb_c_lpk_)      :: irw(*), icl(*)
      complex(c_float_complex)        :: val(*)
      type(psb_c_cspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cspins
  end interface

  interface
    module function psb_c_csprn(mh,cdh,clear) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      logical(c_bool), value :: clear
      type(psb_c_cspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_csprn
  end interface
!!$
!!$  module function psb_c_cspprint(mh) bind(c) result(res)
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
!!$  end function psb_c_cspprint

  interface 
    module function psb_c_cgetelem(xh,index,cdh) bind(c) result(res)
      type(psb_c_cvector)      :: xh
      integer(psb_c_lpk_), value :: index
      type(psb_c_descriptor)     :: cdh
      complex(c_float_complex)           :: res
    end function psb_c_cgetelem
  end interface
  
  interface 
    module function psb_c_csetelem(index,val,xh,cdh) bind(c) result(res)
      type(psb_c_cvector)      :: xh
      integer(psb_c_lpk_), value :: index
      type(psb_c_descriptor)     :: cdh
      complex(c_float_complex), value    :: val
      integer(psb_c_ipk_) :: res
    end function psb_c_csetelem
  end interface

  interface 
    module function psb_c_cmatgetelem(ah,rowindex,colindex,cdh) bind(c) result(res)
      type(psb_c_cspmat)      :: ah
      integer(psb_c_lpk_), value :: rowindex, colindex
      type(psb_c_descriptor)     :: cdh
      complex(c_float_complex)           :: res
    end function psb_c_cmatgetelem
  end interface

end module psb_c_tools_cbind_mod
