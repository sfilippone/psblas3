module psb_c_psblas_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod


  interface
    module function psb_c_cgeaxpby(alpha,xh,beta,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cvector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      complex(c_float_complex), value   :: alpha,beta
    end function psb_c_cgeaxpby
  end interface

  interface
    module function psb_c_cgeaxpbyz(alpha,xh,beta,yh,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cvector) :: xh,yh,zh
      type(psb_c_descriptor) :: cdh
      complex(c_float_complex), value   :: alpha,beta
    end function psb_c_cgeaxpbyz
  end interface

  interface
    module function psb_c_cgemlt(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgemlt
  end interface

  interface
    module function psb_c_cgemlt2(alpha,xh,yh,beta,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cvector)  :: xh,yh, zh
      type(psb_c_descriptor) :: cdh
      complex(psb_spk_), intent(in), value  :: alpha,beta
    end function psb_c_cgemlt2
  end interface

  interface
    module function psb_c_cgediv(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgediv
  end interface
  
  interface
    module function psb_c_cgediv2(xh,yh,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cvector)  :: xh,yh,zh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgediv2
  end interface

  interface
    module function psb_c_cgediv_check(xh,yh,cdh, flag) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
      logical(c_bool), value :: flag
    end function psb_c_cgediv_check
  end interface
  
  interface
    module function psb_c_cgediv2_check(xh,yh,zh,cdh, flag) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cvector)  :: xh,yh,zh
      type(psb_c_descriptor) :: cdh
      logical(c_bool), value :: flag
    end function psb_c_cgediv2_check
  end interface
  
  interface
    module function psb_c_cgeinv(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgeinv
  end interface

  interface
    module function psb_c_cgeinv_check(xh,yh,cdh, flag) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
      logical(c_bool), value :: flag
    end function psb_c_cgeinv_check
  end interface

  interface
    module function psb_c_cgeabs(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgeabs
  end interface

  interface
    module function psb_c_cgecmp(xh,ch,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cvector)  :: xh,zh
      type(psb_c_descriptor) :: cdh
      real(c_float_complex), value :: ch
    end function psb_c_cgecmp
  end interface

  interface
    module function psb_c_cgecmpmat(ah,bh,tol,cdh) bind(c) result(res)
      logical(c_bool)       :: res
      type(psb_c_cspmat)  :: ah,bh
      type(psb_c_descriptor) :: cdh
      real(c_float_complex), value :: tol
    end function psb_c_cgecmpmat
  end interface

  interface
    module function psb_c_cgecmpmat_val(ah,val,tol,cdh) bind(c) result(res)
      logical(c_bool)   :: res
      type(psb_c_cspmat)  :: ah
      type(psb_c_descriptor) :: cdh
      complex(c_float_complex), value :: val
      real(c_float_complex), value :: tol
    end function psb_c_cgecmpmat_val
  end interface

  interface
    module function psb_c_cgeaddconst(xh,bh,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cvector)  :: xh,zh
      type(psb_c_descriptor) :: cdh
      real(c_float_complex), value :: bh
    end function psb_c_cgeaddconst
  end interface


  interface
    module function psb_c_cgenrm2(xh,cdh) bind(c) result(res)
      real(c_float_complex) :: res
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgenrm2
  end interface

  interface
    module function psb_c_cgenrmi(xh,cdh) bind(c) result(res)
      real(c_float_complex) :: res
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgenrmi
  end interface

  interface
    module function psb_c_cgenrm2_weight(xh,wh,cdh) bind(c) result(res)
      real(c_float_complex) :: res
      type(psb_c_cvector) :: xh, wh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgenrm2_weight
  end interface
  
  interface
    module function psb_c_cgenrm2_weightmask(xh,wh,idvh,cdh) bind(c) result(res)
      real(c_float_complex) :: res
      type(psb_c_cvector) :: xh, wh, idvh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgenrm2_weightmask
  end interface

  interface
    module function psb_c_cgeamax(xh,cdh) bind(c) result(res)
      real(c_float_complex) :: res
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgeamax
  end interface


  interface
    module function psb_c_cgeasum(xh,cdh) bind(c) result(res)
      real(c_float_complex) :: res
      type(psb_c_cvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgeasum
  end interface

  interface
    module function psb_c_cspnrmi(ah,cdh) bind(c) result(res)
      real(c_float_complex) :: res
      type(psb_c_cspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_cspnrmi
  end interface

  interface
    module function psb_c_cgedot(xh,yh,cdh) bind(c) result(res)
      complex(c_float_complex) :: res
      type(psb_c_cvector) :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cgedot
  end interface

  interface
    module function psb_c_cspmm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cspmat) :: ah
      type(psb_c_cvector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      complex(c_float_complex), value :: alpha, beta
    end function psb_c_cspmm
  end interface

  interface
    module function psb_c_cspmm_opt(alpha,ah,xh,beta,yh,&
      & cdh,trans,doswap) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cspmat) :: ah
      type(psb_c_cvector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      complex(c_float_complex), value :: alpha, beta
      character(c_char)      :: trans
      logical(c_bool), value :: doswap
    end function psb_c_cspmm_opt
  end interface

  interface
    module function psb_c_cspsm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cspmat) :: ah
      type(psb_c_cvector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      complex(c_float_complex), value :: alpha, beta
      type(psb_desc_type), pointer :: descp
    end function psb_c_cspsm
  end interface
  
  interface
    module function psb_c_cnnz(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_cnnz
  end interface

  interface
    module function psb_c_cis_matupd(ah,cdh) bind(c) result(res)
      logical(c_bool) :: res
      type(psb_c_cspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_cis_matasb(ah,cdh) bind(c) result(res)
      logical(c_bool) :: res
      type(psb_c_cspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_cis_matbld(ah,cdh) bind(c) result(res)
      logical(c_bool) :: res
      type(psb_c_cspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_cset_matupd(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_cspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_cset_matasb(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface
  
  interface
    module function psb_c_cset_matbld(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_ccopy_mat(ah,bh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_cspmat)   :: ah,bh
      type(psb_c_descriptor) :: cdh
    end function
  end interface
  
  interface
    module function psb_c_cspscal(alpha,ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)              ::  res
      complex(c_float_complex), value :: alpha
      type(psb_c_cspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_cspscal
  end interface

  interface
    module function psb_c_cspscalpid(alpha,ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)              ::  res
      complex(c_float_complex), value :: alpha
      type(psb_c_cspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_cspscalpid
  end interface

  interface
    module function psb_c_cspaxpby(alpha,ah,beta,bh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)              ::  res
      complex(c_float_complex), value :: alpha
      type(psb_c_cspmat)   :: ah
      complex(c_float_complex), value :: beta
      type(psb_c_cspmat)   :: bh
      type(psb_c_descriptor) :: cdh
    end function psb_c_cspaxpby
  end interface

end module psb_c_psblas_cbind_mod
