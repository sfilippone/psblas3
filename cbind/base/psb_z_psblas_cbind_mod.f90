module psb_z_psblas_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod


  interface
    module function psb_c_zgeaxpby(alpha,xh,beta,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zvector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      complex(c_double_complex), value   :: alpha,beta
    end function psb_c_zgeaxpby
  end interface

  interface
    module function psb_c_zgeaxpbyz(alpha,xh,beta,yh,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zvector) :: xh,yh,zh
      type(psb_c_descriptor) :: cdh
      complex(c_double_complex), value   :: alpha,beta
    end function psb_c_zgeaxpbyz
  end interface

  interface
    module function psb_c_zgemlt(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgemlt
  end interface

  interface
    module function psb_c_zgemlt2(alpha,xh,yh,beta,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zvector)  :: xh,yh, zh
      type(psb_c_descriptor) :: cdh
      complex(psb_dpk_), intent(in), value  :: alpha,beta
    end function psb_c_zgemlt2
  end interface

  interface
    module function psb_c_zgediv(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgediv
  end interface
  
  interface
    module function psb_c_zgediv2(xh,yh,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zvector)  :: xh,yh,zh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgediv2
  end interface

  interface
    module function psb_c_zgediv_check(xh,yh,cdh, flag) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
      logical(c_bool), value :: flag
    end function psb_c_zgediv_check
  end interface
  
  interface
    module function psb_c_zgediv2_check(xh,yh,zh,cdh, flag) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zvector)  :: xh,yh,zh
      type(psb_c_descriptor) :: cdh
      logical(c_bool), value :: flag
    end function psb_c_zgediv2_check
  end interface
  
  interface
    module function psb_c_zgeinv(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgeinv
  end interface

  interface
    module function psb_c_zgeinv_check(xh,yh,cdh, flag) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
      logical(c_bool), value :: flag
    end function psb_c_zgeinv_check
  end interface

  interface
    module function psb_c_zgeabs(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgeabs
  end interface

  interface
    module function psb_c_zgecmp(xh,ch,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zvector)  :: xh,zh
      type(psb_c_descriptor) :: cdh
      real(c_double_complex), value :: ch
    end function psb_c_zgecmp
  end interface

  interface
    module function psb_c_zgecmpmat(ah,bh,tol,cdh) bind(c) result(res)
      logical(c_bool)       :: res
      type(psb_c_zspmat)  :: ah,bh
      type(psb_c_descriptor) :: cdh
      real(c_double_complex), value :: tol
    end function psb_c_zgecmpmat
  end interface

  interface
    module function psb_c_zgecmpmat_val(ah,val,tol,cdh) bind(c) result(res)
      logical(c_bool)   :: res
      type(psb_c_zspmat)  :: ah
      type(psb_c_descriptor) :: cdh
      complex(c_double_complex), value :: val
      real(c_double_complex), value :: tol
    end function psb_c_zgecmpmat_val
  end interface

  interface
    module function psb_c_zgeaddconst(xh,bh,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zvector)  :: xh,zh
      type(psb_c_descriptor) :: cdh
      real(c_double_complex), value :: bh
    end function psb_c_zgeaddconst
  end interface


  interface
    module function psb_c_zgenrm2(xh,cdh) bind(c) result(res)
      real(c_double_complex) :: res
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgenrm2
  end interface

  interface
    module function psb_c_zgenrmi(xh,cdh) bind(c) result(res)
      real(c_double_complex) :: res
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgenrmi
  end interface

  interface
    module function psb_c_zgenrm2_weight(xh,wh,cdh) bind(c) result(res)
      real(c_double_complex) :: res
      type(psb_c_zvector) :: xh, wh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgenrm2_weight
  end interface
  
  interface
    module function psb_c_zgenrm2_weightmask(xh,wh,idvh,cdh) bind(c) result(res)
      real(c_double_complex) :: res
      type(psb_c_zvector) :: xh, wh, idvh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgenrm2_weightmask
  end interface

  interface
    module function psb_c_zgeamax(xh,cdh) bind(c) result(res)
      real(c_double_complex) :: res
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgeamax
  end interface


  interface
    module function psb_c_zgeasum(xh,cdh) bind(c) result(res)
      real(c_double_complex) :: res
      type(psb_c_zvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgeasum
  end interface

  interface
    module function psb_c_zspnrmi(ah,cdh) bind(c) result(res)
      real(c_double_complex) :: res
      type(psb_c_zspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_zspnrmi
  end interface

  interface
    module function psb_c_zgedot(xh,yh,cdh) bind(c) result(res)
      complex(c_double_complex) :: res
      type(psb_c_zvector) :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zgedot
  end interface

  interface
    module function psb_c_zspmm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zspmat) :: ah
      type(psb_c_zvector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      complex(c_double_complex), value :: alpha, beta
    end function psb_c_zspmm
  end interface

  interface
    module function psb_c_zspmm_opt(alpha,ah,xh,beta,yh,&
      & cdh,trans,doswap) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zspmat) :: ah
      type(psb_c_zvector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      complex(c_double_complex), value :: alpha, beta
      character(c_char)      :: trans
      logical(c_bool), value :: doswap
    end function psb_c_zspmm_opt
  end interface

  interface
    module function psb_c_zspsm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zspmat) :: ah
      type(psb_c_zvector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      complex(c_double_complex), value :: alpha, beta
      type(psb_desc_type), pointer :: descp
    end function psb_c_zspsm
  end interface
  
  interface
    module function psb_c_znnz(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_znnz
  end interface

  interface
    module function psb_c_zis_matupd(ah,cdh) bind(c) result(res)
      logical(c_bool) :: res
      type(psb_c_zspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_zis_matasb(ah,cdh) bind(c) result(res)
      logical(c_bool) :: res
      type(psb_c_zspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_zis_matbld(ah,cdh) bind(c) result(res)
      logical(c_bool) :: res
      type(psb_c_zspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_zset_matupd(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_zspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_zset_matasb(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface
  
  interface
    module function psb_c_zset_matbld(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_zcopy_mat(ah,bh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_zspmat)   :: ah,bh
      type(psb_c_descriptor) :: cdh
    end function
  end interface
  
  interface
    module function psb_c_zspscal(alpha,ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)              ::  res
      complex(c_double_complex), value :: alpha
      type(psb_c_zspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_zspscal
  end interface

  interface
    module function psb_c_zspscalpid(alpha,ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)              ::  res
      complex(c_double_complex), value :: alpha
      type(psb_c_zspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_zspscalpid
  end interface

  interface
    module function psb_c_zspaxpby(alpha,ah,beta,bh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)              ::  res
      complex(c_double_complex), value :: alpha
      type(psb_c_zspmat)   :: ah
      complex(c_double_complex), value :: beta
      type(psb_c_zspmat)   :: bh
      type(psb_c_descriptor) :: cdh
    end function psb_c_zspaxpby
  end interface

end module psb_z_psblas_cbind_mod
