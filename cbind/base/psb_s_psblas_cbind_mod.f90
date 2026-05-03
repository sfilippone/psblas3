module psb_s_psblas_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod


  interface
    module function psb_c_sgeaxpby(alpha,xh,beta,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_svector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      real(c_float), value   :: alpha,beta
    end function psb_c_sgeaxpby
  end interface

  interface
    module function psb_c_sgeaxpbyz(alpha,xh,beta,yh,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_svector) :: xh,yh,zh
      type(psb_c_descriptor) :: cdh
      real(c_float), value   :: alpha,beta
    end function psb_c_sgeaxpbyz
  end interface

  interface
    module function psb_c_sgemlt(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_svector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgemlt
  end interface

  interface
    module function psb_c_sgemlt2(alpha,xh,yh,beta,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_svector)  :: xh,yh, zh
      type(psb_c_descriptor) :: cdh
      real(psb_spk_), intent(in), value  :: alpha,beta
    end function psb_c_sgemlt2
  end interface

  interface
    module function psb_c_sgediv(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_svector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgediv
  end interface
  
  interface
    module function psb_c_sgediv2(xh,yh,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_svector)  :: xh,yh,zh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgediv2
  end interface

  interface
    module function psb_c_sgediv_check(xh,yh,cdh, flag) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_svector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
      logical(c_bool), value :: flag
    end function psb_c_sgediv_check
  end interface
  
  interface
    module function psb_c_sgediv2_check(xh,yh,zh,cdh, flag) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_svector)  :: xh,yh,zh
      type(psb_c_descriptor) :: cdh
      logical(c_bool), value :: flag
    end function psb_c_sgediv2_check
  end interface
  
  interface
    module function psb_c_sgeinv(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_svector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgeinv
  end interface

  interface
    module function psb_c_sgeinv_check(xh,yh,cdh, flag) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_svector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
      logical(c_bool), value :: flag
    end function psb_c_sgeinv_check
  end interface

  interface
    module function psb_c_sgeabs(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_svector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgeabs
  end interface

  interface
    module function psb_c_sgecmp(xh,ch,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_svector)  :: xh,zh
      type(psb_c_descriptor) :: cdh
      real(c_float), value :: ch
    end function psb_c_sgecmp
  end interface

  interface
    module function psb_c_sgecmpmat(ah,bh,tol,cdh) bind(c) result(res)
      logical(c_bool)       :: res
      type(psb_c_sspmat)  :: ah,bh
      type(psb_c_descriptor) :: cdh
      real(c_float), value :: tol
    end function psb_c_sgecmpmat
  end interface

  interface
    module function psb_c_sgecmpmat_val(ah,val,tol,cdh) bind(c) result(res)
      logical(c_bool)   :: res
      type(psb_c_sspmat)  :: ah
      type(psb_c_descriptor) :: cdh
      real(c_float), value :: val
      real(c_float), value :: tol
    end function psb_c_sgecmpmat_val
  end interface

  interface
    module function psb_c_sgeaddconst(xh,bh,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_svector)  :: xh,zh
      type(psb_c_descriptor) :: cdh
      real(c_float), value :: bh
    end function psb_c_sgeaddconst
  end interface

  interface
    module function psb_c_smask(ch,xh,mh,t,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_svector)  :: ch,xh,mh
      type(psb_c_descriptor) :: cdh
      logical(c_bool)        :: t
    end function psb_c_smask
  end interface
  
  interface
    module function psb_c_sminquotient(xh,yh,cdh) bind(c) result(res)
      real(psb_spk_)       :: res
      type(psb_c_svector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sminquotient
  end interface

  interface
    module function psb_c_sgenrm2(xh,cdh) bind(c) result(res)
      real(c_float) :: res
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgenrm2
  end interface

  interface
    module function psb_c_sgenrmi(xh,cdh) bind(c) result(res)
      real(c_float) :: res
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgenrmi
  end interface

  interface
    module function psb_c_sgenrm2_weight(xh,wh,cdh) bind(c) result(res)
      real(c_float) :: res
      type(psb_c_svector) :: xh, wh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgenrm2_weight
  end interface
  
  interface
    module function psb_c_sgenrm2_weightmask(xh,wh,idvh,cdh) bind(c) result(res)
      real(c_float) :: res
      type(psb_c_svector) :: xh, wh, idvh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgenrm2_weightmask
  end interface

  interface
    module function psb_c_sgeamax(xh,cdh) bind(c) result(res)
      real(c_float) :: res
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgeamax
  end interface

  interface
    module function psb_c_sgemin(xh,cdh) bind(c) result(res)
      real(c_float) :: res
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgemin
  end interface

  interface
    module function psb_c_sgeasum(xh,cdh) bind(c) result(res)
      real(c_float) :: res
      type(psb_c_svector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgeasum
  end interface

  interface
    module function psb_c_sspnrmi(ah,cdh) bind(c) result(res)
      real(c_float) :: res
      type(psb_c_sspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_sspnrmi
  end interface

  interface
    module function psb_c_sgedot(xh,yh,cdh) bind(c) result(res)
      real(c_float) :: res
      type(psb_c_svector) :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sgedot
  end interface

  interface
    module function psb_c_sspmm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_sspmat) :: ah
      type(psb_c_svector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      real(c_float), value :: alpha, beta
    end function psb_c_sspmm
  end interface

  interface
    module function psb_c_sspmm_opt(alpha,ah,xh,beta,yh,&
      & cdh,trans,doswap) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_sspmat) :: ah
      type(psb_c_svector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      real(c_float), value :: alpha, beta
      character(c_char)      :: trans
      logical(c_bool), value :: doswap
    end function psb_c_sspmm_opt
  end interface

  interface
    module function psb_c_sspsm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_sspmat) :: ah
      type(psb_c_svector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      real(c_float), value :: alpha, beta
      type(psb_desc_type), pointer :: descp
    end function psb_c_sspsm
  end interface
  
  interface
    module function psb_c_snnz(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_sspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_snnz
  end interface

  interface
    module function psb_c_sis_matupd(ah,cdh) bind(c) result(res)
      logical(c_bool) :: res
      type(psb_c_sspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_sis_matasb(ah,cdh) bind(c) result(res)
      logical(c_bool) :: res
      type(psb_c_sspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_sis_matbld(ah,cdh) bind(c) result(res)
      logical(c_bool) :: res
      type(psb_c_sspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_sset_matupd(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_sspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_sset_matasb(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_sspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface
  
  interface
    module function psb_c_sset_matbld(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_sspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_scopy_mat(ah,bh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_sspmat)   :: ah,bh
      type(psb_c_descriptor) :: cdh
    end function
  end interface
  
  interface
    module function psb_c_sspscal(alpha,ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)              ::  res
      real(c_float), value :: alpha
      type(psb_c_sspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_sspscal
  end interface

  interface
    module function psb_c_sspscalpid(alpha,ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)              ::  res
      real(c_float), value :: alpha
      type(psb_c_sspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_sspscalpid
  end interface

  interface
    module function psb_c_sspaxpby(alpha,ah,beta,bh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)              ::  res
      real(c_float), value :: alpha
      type(psb_c_sspmat)   :: ah
      real(c_float), value :: beta
      type(psb_c_sspmat)   :: bh
      type(psb_c_descriptor) :: cdh
    end function psb_c_sspaxpby
  end interface

end module psb_s_psblas_cbind_mod
