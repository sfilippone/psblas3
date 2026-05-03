module psb_d_psblas_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod


  interface
    module function psb_c_dgeaxpby(alpha,xh,beta,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dvector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      real(c_double), value   :: alpha,beta
    end function psb_c_dgeaxpby
  end interface

  interface
    module function psb_c_dgeaxpbyz(alpha,xh,beta,yh,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dvector) :: xh,yh,zh
      type(psb_c_descriptor) :: cdh
      real(c_double), value   :: alpha,beta
    end function psb_c_dgeaxpbyz
  end interface

  interface
    module function psb_c_dgemlt(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgemlt
  end interface

  interface
    module function psb_c_dgemlt2(alpha,xh,yh,beta,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dvector)  :: xh,yh, zh
      type(psb_c_descriptor) :: cdh
      real(psb_dpk_), intent(in), value  :: alpha,beta
    end function psb_c_dgemlt2
  end interface

  interface
    module function psb_c_dgediv(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgediv
  end interface
  
  interface
    module function psb_c_dgediv2(xh,yh,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dvector)  :: xh,yh,zh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgediv2
  end interface

  interface
    module function psb_c_dgediv_check(xh,yh,cdh, flag) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
      logical(c_bool), value :: flag
    end function psb_c_dgediv_check
  end interface
  
  interface
    module function psb_c_dgediv2_check(xh,yh,zh,cdh, flag) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dvector)  :: xh,yh,zh
      type(psb_c_descriptor) :: cdh
      logical(c_bool), value :: flag
    end function psb_c_dgediv2_check
  end interface
  
  interface
    module function psb_c_dgeinv(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgeinv
  end interface

  interface
    module function psb_c_dgeinv_check(xh,yh,cdh, flag) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
      logical(c_bool), value :: flag
    end function psb_c_dgeinv_check
  end interface

  interface
    module function psb_c_dgeabs(xh,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgeabs
  end interface

  interface
    module function psb_c_dgecmp(xh,ch,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dvector)  :: xh,zh
      type(psb_c_descriptor) :: cdh
      real(c_double), value :: ch
    end function psb_c_dgecmp
  end interface

  interface
    module function psb_c_dgecmpmat(ah,bh,tol,cdh) bind(c) result(res)
      logical(c_bool)       :: res
      type(psb_c_dspmat)  :: ah,bh
      type(psb_c_descriptor) :: cdh
      real(c_double), value :: tol
    end function psb_c_dgecmpmat
  end interface

  interface
    module function psb_c_dgecmpmat_val(ah,val,tol,cdh) bind(c) result(res)
      logical(c_bool)   :: res
      type(psb_c_dspmat)  :: ah
      type(psb_c_descriptor) :: cdh
      real(c_double), value :: val
      real(c_double), value :: tol
    end function psb_c_dgecmpmat_val
  end interface

  interface
    module function psb_c_dgeaddconst(xh,bh,zh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dvector)  :: xh,zh
      type(psb_c_descriptor) :: cdh
      real(c_double), value :: bh
    end function psb_c_dgeaddconst
  end interface

  interface
    module function psb_c_dmask(ch,xh,mh,t,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dvector)  :: ch,xh,mh
      type(psb_c_descriptor) :: cdh
      logical(c_bool)        :: t
    end function psb_c_dmask
  end interface
  
  interface
    module function psb_c_dminquotient(xh,yh,cdh) bind(c) result(res)
      real(psb_dpk_)       :: res
      type(psb_c_dvector)  :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dminquotient
  end interface

  interface
    module function psb_c_dgenrm2(xh,cdh) bind(c) result(res)
      real(c_double) :: res
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgenrm2
  end interface

  interface
    module function psb_c_dgenrmi(xh,cdh) bind(c) result(res)
      real(c_double) :: res
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgenrmi
  end interface

  interface
    module function psb_c_dgenrm2_weight(xh,wh,cdh) bind(c) result(res)
      real(c_double) :: res
      type(psb_c_dvector) :: xh, wh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgenrm2_weight
  end interface
  
  interface
    module function psb_c_dgenrm2_weightmask(xh,wh,idvh,cdh) bind(c) result(res)
      real(c_double) :: res
      type(psb_c_dvector) :: xh, wh, idvh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgenrm2_weightmask
  end interface

  interface
    module function psb_c_dgeamax(xh,cdh) bind(c) result(res)
      real(c_double) :: res
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgeamax
  end interface

  interface
    module function psb_c_dgemin(xh,cdh) bind(c) result(res)
      real(c_double) :: res
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgemin
  end interface

  interface
    module function psb_c_dgeasum(xh,cdh) bind(c) result(res)
      real(c_double) :: res
      type(psb_c_dvector) :: xh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgeasum
  end interface

  interface
    module function psb_c_dspnrmi(ah,cdh) bind(c) result(res)
      real(c_double) :: res
      type(psb_c_dspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_dspnrmi
  end interface

  interface
    module function psb_c_dgedot(xh,yh,cdh) bind(c) result(res)
      real(c_double) :: res
      type(psb_c_dvector) :: xh,yh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dgedot
  end interface

  interface
    module function psb_c_dspmm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dspmat) :: ah
      type(psb_c_dvector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      real(c_double), value :: alpha, beta
    end function psb_c_dspmm
  end interface

  interface
    module function psb_c_dspmm_opt(alpha,ah,xh,beta,yh,&
      & cdh,trans,doswap) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dspmat) :: ah
      type(psb_c_dvector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      real(c_double), value :: alpha, beta
      character(c_char)      :: trans
      logical(c_bool), value :: doswap
    end function psb_c_dspmm_opt
  end interface

  interface
    module function psb_c_dspsm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dspmat) :: ah
      type(psb_c_dvector) :: xh,yh
      type(psb_c_descriptor) :: cdh
      real(c_double), value :: alpha, beta
      type(psb_desc_type), pointer :: descp
    end function psb_c_dspsm
  end interface
  
  interface
    module function psb_c_dnnz(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_dnnz
  end interface

  interface
    module function psb_c_dis_matupd(ah,cdh) bind(c) result(res)
      logical(c_bool) :: res
      type(psb_c_dspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_dis_matasb(ah,cdh) bind(c) result(res)
      logical(c_bool) :: res
      type(psb_c_dspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_dis_matbld(ah,cdh) bind(c) result(res)
      logical(c_bool) :: res
      type(psb_c_dspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_dset_matupd(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_) :: res
      type(psb_c_dspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_dset_matasb(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface
  
  interface
    module function psb_c_dset_matbld(ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function
  end interface

  interface
    module function psb_c_dcopy_mat(ah,bh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)    :: res
      type(psb_c_dspmat)   :: ah,bh
      type(psb_c_descriptor) :: cdh
    end function
  end interface
  
  interface
    module function psb_c_dspscal(alpha,ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)              ::  res
      real(c_double), value :: alpha
      type(psb_c_dspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_dspscal
  end interface

  interface
    module function psb_c_dspscalpid(alpha,ah,cdh) bind(c) result(res)
      integer(psb_c_ipk_)              ::  res
      real(c_double), value :: alpha
      type(psb_c_dspmat)   :: ah
      type(psb_c_descriptor) :: cdh
    end function psb_c_dspscalpid
  end interface

  interface
    module function psb_c_dspaxpby(alpha,ah,beta,bh,cdh) bind(c) result(res)
      integer(psb_c_ipk_)              ::  res
      real(c_double), value :: alpha
      type(psb_c_dspmat)   :: ah
      real(c_double), value :: beta
      type(psb_c_dspmat)   :: bh
      type(psb_c_descriptor) :: cdh
    end function psb_c_dspaxpby
  end interface

end module psb_d_psblas_cbind_mod
