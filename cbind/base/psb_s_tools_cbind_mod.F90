module psb_s_tools_cbind_mod
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
    module function psb_c_sgeall(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgeall
  end interface

  interface
    module function psb_c_sgeall_remote(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
      
    end function psb_c_sgeall_remote
  end interface

  interface
    module function psb_c_sgeall_remote_options(xh,cdh,bldmode,dupl) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
      integer(psb_c_ipk_), value :: dupl
      integer(psb_c_ipk_), value :: bldmode
    end function psb_c_sgeall_remote_options
  end interface

  interface
    module function psb_c_sgeasb(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgeasb
  end interface
  
  interface
    module function psb_c_sgeasb_options(xh,cdh,dupl) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
      integer(psb_c_ipk_), value :: dupl
    end function psb_c_sgeasb_options
  end interface
  
  interface
    module function psb_c_sgeasb_options_format(xh,cdh,dupl,format) bind(c) result(res)
      ! Takes into account format argument as a c string, and uses it to call the appropriate psb_geasb
      ! with mold argument
      integer(psb_c_ipk_) :: res
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
      character(kind=c_char), dimension(*) :: format
      integer(psb_c_ipk_), value :: dupl
    end function psb_c_sgeasb_options_format
  end interface
   
  interface
    module function psb_c_sgefree(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgefree
  end interface

  interface
    module function psb_c_sgeins(nz,irw,val,xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: nz
      integer(psb_c_lpk_)        :: irw(*)
      real(c_float)        :: val(*)
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgeins
  end interface

  interface
    module function psb_c_sspall(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_sspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sspall
  end interface
  
  interface
    module function psb_c_sspall_remote(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_sspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sspall_remote
  end interface

  interface
    module function psb_c_sspasb(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_sspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sspasb
  end interface

  interface
    module function psb_c_sspfree(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_sspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sspfree
  end interface

  interface
    module function psb_c_sspasb_opt(mh,cdh,afmt,upd,dupl) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_sspmat) :: mh
      type(psb_c_descriptor) :: cdh
      integer(psb_c_ipk_), value :: upd,dupl
      character(c_char)     :: afmt(*)
    end function psb_c_sspasb_opt
  end interface

  interface
    module function psb_c_sspins(nz,irw,icl,val,mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: nz
      integer(psb_c_lpk_)      :: irw(*), icl(*)
      real(c_float)        :: val(*)
      type(psb_c_sspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sspins
  end interface

  interface
    module function psb_c_ssprn(mh,cdh,clear) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      logical(c_bool), value :: clear
      type(psb_c_sspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_ssprn
  end interface
!!$
!!$  module function psb_c_sspprint(mh) bind(c) result(res)
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
!!$  end function psb_c_sspprint

  interface 
    module function psb_c_sgetelem(xh,index,cdh) bind(c) result(res)
      type(psb_c_svector)      :: xh
      integer(psb_c_lpk_), value :: index
      type(psb_c_descriptor)     :: cdh
      real(c_float)           :: res
    end function psb_c_sgetelem
  end interface
  
  interface 
    module function psb_c_ssetelem(index,val,xh,cdh) bind(c) result(res)
      type(psb_c_svector)      :: xh
      integer(psb_c_lpk_), value :: index
      type(psb_c_descriptor)     :: cdh
      real(c_float), value    :: val
      integer(psb_c_ipk_) :: res
    end function psb_c_ssetelem
  end interface

  interface 
    module function psb_c_smatgetelem(ah,rowindex,colindex,cdh) bind(c) result(res)
      type(psb_c_sspmat)      :: ah
      integer(psb_c_lpk_), value :: rowindex, colindex
      type(psb_c_descriptor)     :: cdh
      real(c_float)           :: res
    end function psb_c_smatgetelem
  end interface

end module psb_s_tools_cbind_mod
