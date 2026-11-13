module psb_base_tools_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_cpenv_mod
#ifdef PSB_HAVE_CUDA
  use psb_cuda_mod
#endif

  interface 
    module function psb_c_error() bind(c) result(res)
      integer(psb_c_ipk_) :: res
    end function psb_c_error
  end interface

  interface 
    module function psb_c_clean_errstack() bind(c) result(res)
      integer(psb_c_ipk_) :: res
    end function psb_c_clean_errstack
  end interface

  interface 
    module function psb_c_cdall_vg(ng,vg,cctxt,cdh) bind(c,name='psb_c_cdall_vg') result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_lpk_), value :: ng
      type(psb_c_object_type), value :: cctxt
      integer(psb_c_ipk_)        :: vg(*)
      type(psb_c_object_type) :: cdh
    end function psb_c_cdall_vg
  end interface

  interface 
    module function psb_c_cdall_vl(nl,vl,cctxt,cdh) bind(c,name='psb_c_cdall_vl') result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_object_type), value :: cctxt
      integer(psb_c_ipk_), value :: nl
      integer(psb_c_lpk_)        :: vl(*)
      type(psb_c_object_type) :: cdh
    end function psb_c_cdall_vl
  end interface
  
  interface 
    module function psb_c_cdall_vl_opt(nl,vl,cctxt,cdh) bind(c,name='psb_c_cdall_vl_opt') result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_object_type), value :: cctxt
      integer(psb_c_ipk_), value :: nl
      integer(psb_c_lpk_)        :: vl(*)
      type(psb_c_object_type) :: cdh
    end function psb_c_cdall_vl_opt
  end interface

  interface 
    module function psb_c_cdall_vl_lidx(nl,vl,lidx,cctxt,cdh) bind(c,name='psb_c_cdall_vl_lidx') result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_object_type), value :: cctxt
      integer(psb_c_ipk_), value :: nl
      integer(psb_c_lpk_)        :: vl(*)
      integer(psb_c_ipk_)        :: lidx(*)
      type(psb_c_object_type) :: cdh
    end function psb_c_cdall_vl_lidx
  end interface

  interface 
    module function psb_c_cdall_nl(nl,cctxt,cdh) bind(c,name='psb_c_cdall_nl') result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_object_type), value :: cctxt
      integer(psb_c_ipk_), value :: nl
      type(psb_c_object_type) :: cdh
    end function psb_c_cdall_nl
  end interface
  
  interface 
    module function psb_c_cdall_repl(n,cctxt,cdh) bind(c,name='psb_c_cdall_repl') result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_lpk_), value :: n
      type(psb_c_object_type), value :: cctxt
      type(psb_c_object_type) :: cdh
    end function psb_c_cdall_repl
  end interface

  interface 
    module function psb_c_cdasb(cdh) bind(c,name='psb_c_cdasb') result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_object_type) :: cdh
    end function psb_c_cdasb
  end interface

  interface 
    module function psb_c_cdasb_format(cdh,format) bind(c,name='psb_c_cdasb_format') result(res)
      ! Takes as input the desired format bewten CPU or GPU, and assembles accordingly
      ! via the mold parameter of psb_cdasb
      integer(psb_c_ipk_) :: res
      type(psb_c_object_type) :: cdh
      character(c_char), dimension(*) :: format
    end function psb_c_cdasb_format
  end interface

  interface 
    module function psb_c_cdfree(cdh) bind(c,name='psb_c_cdfree') result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_object_type) :: cdh
    end function psb_c_cdfree
  end interface

  interface
    module function psb_c_cdins(nz,ia,ja,cdh) bind(c,name='psb_c_cdins') result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value   :: nz
      type(psb_c_object_type) :: cdh
      integer(psb_c_lpk_)          :: ia(*),ja(*)
    end function psb_c_cdins
  end interface
  
  interface 
    module function psb_c_cdins_lidx(nz,ja,lidx,cdh) bind(c,name='psb_c_cdins_lidx') result(res)
      integer(psb_c_ipk_) :: res
      integer(psb_c_ipk_), value   :: nz
      type(psb_c_object_type) :: cdh
      integer(psb_c_lpk_)          :: ja(*)
      integer(psb_c_ipk_) :: lidx(*)
    end function psb_c_cdins_lidx
  end interface

  interface 
    module function psb_c_cd_is_asb(cdh) bind(c,name='psb_c_cd_is_asb') result(res)
      logical(c_bool) :: res
      type(psb_c_object_type) :: cdh
    end function psb_c_cd_is_asb
  end interface
  
  interface 
    module function psb_c_cd_check_addr(cdh) &
         & bind(c,name='psb_c_cd_check_addr') result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_object_type) :: cdh
    end function psb_c_cd_check_addr
  end interface
  
  interface
    module function psb_c_cd_get_local_rows(cdh) bind(c,name='psb_c_cd_get_local_rows') result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_object_type) :: cdh
    end function psb_c_cd_get_local_rows
  end interface
  
  interface 
    module function psb_c_cd_get_local_cols(cdh) bind(c,name='psb_c_cd_get_local_cols') result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_object_type) :: cdh
    end function psb_c_cd_get_local_cols
  end interface

  interface 
    module function psb_c_cd_get_global_rows(cdh) bind(c,name='psb_c_cd_get_global_rows') result(res)
      integer(psb_c_lpk_) :: res
      type(psb_c_object_type) :: cdh
    end function psb_c_cd_get_global_rows
  end interface

  interface 
    module function psb_c_cd_get_global_cols(cdh) bind(c,name='psb_c_cd_get_global_cols') result(res)
      integer(psb_c_lpk_) :: res
      type(psb_c_object_type) :: cdh
    end function psb_c_cd_get_global_cols
  end interface
  
  interface 
    module function psb_c_cd_get_global_indices(idx,nidx,owned,cdh) &
         & bind(c,name='psb_c_cd_get_global_indices') result(res)
      integer(psb_c_ipk_)            :: res
      type(psb_c_object_type)        :: cdh
      integer(psb_c_ipk_), value     :: nidx
      integer(psb_c_lpk_)            :: idx(nidx)
      logical(c_bool), value         :: owned
    end function psb_c_cd_get_global_indices
  end interface
  
  interface 
    module function psb_c_is_owned(x,cdh) bind(c,name='psb_c_is_owned') result(res)
      type(psb_c_object_type)        :: cdh
      integer(psb_c_lpk_), value     :: x
      logical(c_bool) :: res
    end function psb_c_is_owned
  end interface

end module psb_base_tools_cbind_mod
