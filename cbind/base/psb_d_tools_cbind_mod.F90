module psb_d_tools_cbind_mod
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
    module function psb_c_dgeall(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgeall
  end interface

  interface
    module function psb_c_dgeall_remote(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
      
    end function psb_c_dgeall_remote
  end interface

  interface
    module function psb_c_dgeall_remote_options(xh,cdh,bldmode,dupl) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
      integer(psb_c_ipk_), value :: dupl
      integer(psb_c_ipk_), value :: bldmode
    end function psb_c_dgeall_remote_options
  end interface

  interface
    module function psb_c_dgeasb(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgeasb
  end interface
  
  interface
    module function psb_c_dgeasb_options(xh,cdh,dupl) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
      integer(psb_c_ipk_), value :: dupl
    end function psb_c_dgeasb_options
  end interface
  
  interface
    module function psb_c_dgeasb_options_format(xh,cdh,dupl,format) bind(c) result(res)
      ! Takes into account format argument as a c string, and uses it to call the appropriate psb_geasb
      ! with mold argument
      integer(psb_c_ipk_) :: res
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
      character(kind=c_char), dimension(*) :: format
      integer(psb_c_ipk_), value :: dupl
    end function psb_c_dgeasb_options_format
  end interface
   
  interface
    module function psb_c_dgefree(xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgefree
  end interface

  interface
    module function psb_c_dgeins(nz,irw,val,xh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: nz
      integer(psb_c_lpk_)        :: irw(*)
      real(c_double)        :: val(*)
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgeins
  end interface

  interface
    module function psb_c_dspall(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dspall
  end interface
  
  interface
    module function psb_c_dspall_remote(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dspall_remote
  end interface

  interface
    module function psb_c_dspasb(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dspasb
  end interface

  interface
    module function psb_c_dspfree(mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dspfree
  end interface

  interface
    module function psb_c_dspasb_opt(mh,cdh,afmt,upd,dupl) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dspmat) :: mh
      type(psb_c_descriptor) :: cdh
      integer(psb_c_ipk_), value :: upd,dupl
      character(c_char)     :: afmt(*)
    end function psb_c_dspasb_opt
  end interface

  interface
    module function psb_c_dspins(nz,irw,icl,val,mh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value :: nz
      integer(psb_c_lpk_)      :: irw(*), icl(*)
      real(c_double)        :: val(*)
      type(psb_c_dspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dspins
  end interface

  interface
    module function psb_c_dsprn(mh,cdh,clear) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      logical(c_bool), value :: clear
      type(psb_c_dspmat) :: mh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dsprn
  end interface
!!$
!!$  module function psb_c_dspprint(mh) bind(c) result(res)
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
!!$  end function psb_c_dspprint

  interface 
    module function psb_c_dgetelem(xh,index,cdh) bind(c) result(res)
      type(psb_c_dvector)      :: xh
      integer(psb_c_lpk_), value :: index
      type(psb_c_descriptor)     :: cdh
      real(c_double)           :: res
    end function psb_c_dgetelem
  end interface
  
  interface 
    module function psb_c_dsetelem(index,val,xh,cdh) bind(c) result(res)
      type(psb_c_dvector)      :: xh
      integer(psb_c_lpk_), value :: index
      type(psb_c_descriptor)     :: cdh
      real(c_double), value    :: val
      integer(psb_c_ipk_) :: res
    end function psb_c_dsetelem
  end interface

  interface 
    module function psb_c_dmatgetelem(ah,rowindex,colindex,cdh) bind(c) result(res)
      type(psb_c_dspmat)      :: ah
      integer(psb_c_lpk_), value :: rowindex, colindex
      type(psb_c_descriptor)     :: cdh
      real(c_double)           :: res
    end function psb_c_dmatgetelem
  end interface

end module psb_d_tools_cbind_mod
