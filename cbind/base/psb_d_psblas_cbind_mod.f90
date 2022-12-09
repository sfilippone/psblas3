module psb_d_psblas_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_string_cbind_mod

contains

  function psb_c_dgeaxpby(alpha,xh,beta,yh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_dvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    real(c_double), value   :: alpha,beta

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    integer(psb_c_ipk_)          :: info


    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if

    call psb_geaxpby(alpha,xp,beta,yp,descp,info)

    res = info

  end function psb_c_dgeaxpby

  function psb_c_dgeaxpbyz(alpha,xh,beta,yh,zh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_dvector) :: xh,yh,zh
    type(psb_c_descriptor) :: cdh
    real(c_double), value   :: alpha,beta

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp,zp
    integer(psb_c_ipk_)          :: info


    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if

    if (c_associated(zh%item)) then
      call c_f_pointer(zh%item,zp)
    else
      return
    end if

    call psb_geaxpby(alpha,xp,beta,yp,zp,descp,info)

    res = info

  end function psb_c_dgeaxpbyz

  function psb_c_dgemlt(xh,yh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_dvector)  :: xh,yh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    integer(psb_c_ipk_)          :: info

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if

    call psb_gemlt(xp,yp,descp,info)

    res = info

  end function psb_c_dgemlt

  function psb_c_dgemlt2(alpha,xh,yh,beta,zh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_dvector)  :: xh,yh, zh
    type(psb_c_descriptor) :: cdh


    type(psb_desc_type), pointer        :: descp
    type(psb_d_vect_type), pointer    :: xp,yp,zp
    integer(psb_c_ipk_)                 :: info
    real(psb_dpk_), intent(in), value  :: alpha,beta

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if
    if (c_associated(zh%item)) then
      call c_f_pointer(zh%item,zp)
    else
      return
    end if

    call psb_gemlt(alpha,xp,yp,beta,zp,descp,info)

    res = info

  end function psb_c_dgemlt2

  function psb_c_dgediv(xh,yh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_dvector)  :: xh,yh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    integer(psb_c_ipk_)          :: info

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if

    call psb_gediv(xp,yp,descp,info)

    res = info

  end function psb_c_dgediv

  function psb_c_dgediv2(xh,yh,zh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_dvector)  :: xh,yh,zh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp,zp
    integer(psb_c_ipk_)          :: info

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if
    if (c_associated(zh%item)) then
      call c_f_pointer(zh%item,zp)
    else
      return
    end if

    call psb_gediv(xp,yp,zp,descp,info)

    res = info

  end function psb_c_dgediv2

  function psb_c_dgediv_check(xh,yh,cdh, flag) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_dvector)  :: xh,yh
    type(psb_c_descriptor) :: cdh
    logical(c_bool), value :: flag

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    integer(psb_c_ipk_)          :: info
    logical   :: fflag

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if

    fflag = flag
    call psb_gediv(xp,yp,descp,info,fflag)

    res = info

  end function psb_c_dgediv_check

  function psb_c_dgediv2_check(xh,yh,zh,cdh, flag) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_dvector)  :: xh,yh,zh
    type(psb_c_descriptor) :: cdh
    logical(c_bool), value :: flag

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp,zp
    integer(psb_c_ipk_)          :: info
    logical   :: fflag

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if
    if (c_associated(zh%item)) then
      call c_f_pointer(zh%item,zp)
    else
      return
    end if

    fflag = flag
    call psb_gediv(xp,yp,zp,descp,info,fflag)

    res = info

  end function psb_c_dgediv2_check

  function psb_c_dgeinv(xh,yh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_dvector)  :: xh,yh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    integer(psb_c_ipk_)          :: info

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if

    call psb_geinv(xp,yp,descp,info)

    res = info

  end function psb_c_dgeinv

  function psb_c_dgeinv_check(xh,yh,cdh, flag) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_dvector)  :: xh,yh
    type(psb_c_descriptor) :: cdh
    logical(c_bool), value :: flag

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    integer(psb_c_ipk_)          :: info
    logical   :: fflag

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if

    fflag = flag
    call psb_geinv(xp,yp,descp,info,fflag)

    res = info

  end function psb_c_dgeinv_check

  function psb_c_dgeabs(xh,yh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_dvector)  :: xh,yh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    integer(psb_c_ipk_)          :: info

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if

    call psb_geabs(xp,yp,descp,info)

    res = info

  end function psb_c_dgeabs

  function psb_c_dgecmp(xh,ch,zh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_dvector)  :: xh,zh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,zp
    integer(psb_c_ipk_)          :: info
    real(c_double), value :: ch

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(zh%item)) then
      call c_f_pointer(zh%item,zp)
    else
      return
    end if

    call psb_gecmp(xp,ch,zp,descp,info)

    res = info

  end function psb_c_dgecmp

  function psb_c_dgecmpmat(ah,bh,tol,cdh) bind(c) result(res)
    implicit none
    logical(c_bool)       :: res

    type(psb_c_dspmat)  :: ah,bh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_dspmat_type), pointer :: ap,bp
    integer(psb_c_ipk_)          :: info
    real(c_double), value :: tol
    logical                  :: isequal

    res = .false.

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if
    if (c_associated(bh%item)) then
      call c_f_pointer(bh%item,bp)
    else
      return
    end if

    call psb_gecmp(ap,bp,tol,descp,isequal,info)

    res = isequal

  end function psb_c_dgecmpmat

  function psb_c_dgecmpmat_val(ah,val,tol,cdh) bind(c) result(res)
    implicit none
    logical(c_bool)   :: res

    type(psb_c_dspmat)  :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_dspmat_type), pointer :: ap
    integer(psb_c_ipk_)          :: info
    real(c_double), value :: val
    real(c_double), value :: tol
    logical                  :: isequal

    res = .false.

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if

    call psb_gecmp(ap,val,tol,descp,isequal,info)

    res = isequal

  end function psb_c_dgecmpmat_val

  function psb_c_dgeaddconst(xh,bh,zh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_dvector)  :: xh,zh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,zp
    integer(psb_c_ipk_)          :: info
    real(c_double), value :: bh

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(zh%item)) then
      call c_f_pointer(zh%item,zp)
    else
      return
    end if

    call psb_geaddconst(xp,bh,zp,descp,info)

    res = info

  end function psb_c_dgeaddconst

  function psb_c_dgescal(xh,ch,zh,cdh) bind(c) result(res)
  implicit none
  integer(psb_c_ipk_)    :: res

  type(psb_c_dvector)  :: xh,zh
  type(psb_c_descriptor) :: cdh

  type(psb_desc_type), pointer :: descp
  type(psb_d_vect_type), pointer :: xp,zp
  integer(psb_c_ipk_)          :: info
  real(c_double) :: ch

  res = -1

  if (c_associated(cdh%item)) then
    call c_f_pointer(cdh%item,descp)
  else
    return
  end if
  if (c_associated(xh%item)) then
    call c_f_pointer(xh%item,xp)
  else
    return
  end if
  if (c_associated(zh%item)) then
    call c_f_pointer(zh%item,zp)
  else
    return
  end if

  call psb_gescal(xp,ch,zp,descp,info)

  res = info

end function psb_c_dgescal

  function psb_c_dmask(ch,xh,mh,t,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_dvector)  :: ch,xh,mh
    type(psb_c_descriptor) :: cdh
    logical(c_bool)        :: t

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: cp,xp,mp
    integer(psb_c_ipk_)          :: info
    logical                      :: fp

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ch%item)) then
      call c_f_pointer(ch%item,cp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(mh%item)) then
      call c_f_pointer(mh%item,mp)
    else
      return
    end if

    call psb_mask(cp,xp,mp,fp,descp,info)

    t = fp
    res = info

  end function psb_c_dmask

  function psb_c_dminquotient(xh,yh,cdh) bind(c) result(res)
    implicit none
    real(psb_dpk_)       :: res

    type(psb_c_dvector)  :: xh,yh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    integer(psb_c_ipk_)          :: info

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if

    res = psb_minquotient(xp,yp,descp,info)


  end function psb_c_dminquotient

  function psb_c_dgenrm2(xh,cdh) bind(c) result(res)
    implicit none
    real(c_double) :: res

    type(psb_c_dvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp
    integer(psb_c_ipk_)               :: info

    res = -1.0

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if

    res = psb_genrm2(xp,descp,info)

  end function psb_c_dgenrm2

  function psb_c_dgenrmi(xh,cdh) bind(c) result(res)
    implicit none
    real(c_double) :: res

    type(psb_c_dvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp
    type(psb_d_vect_type) ::  yp
    integer(psb_c_ipk_)               :: info

    res = -1.0

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if

    call psb_geall(yp,descp,info)
    call psb_geabs(xp,yp,descp,info)
    res = psb_geasum(yp,descp,info)
    call psb_gefree(yp,descp,info)

  end function psb_c_dgenrmi

  function psb_c_dgenrm2_weight(xh,wh,cdh) bind(c) result(res)
    implicit none
    real(c_double) :: res

    type(psb_c_dvector) :: xh, wh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp, wp
    integer(psb_c_ipk_)               :: info

    res = -1.0

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(wh%item)) then
      call c_f_pointer(wh%item,wp)
    else
      return
    end if

    res = psb_genrm2(xp,wp,descp,info)

  end function psb_c_dgenrm2_weight

  function psb_c_dgenrm2_weightmask(xh,wh,idvh,cdh) bind(c) result(res)
    implicit none
    real(c_double) :: res

    type(psb_c_dvector) :: xh, wh, idvh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp, wp, idvp
    integer(psb_c_ipk_)               :: info

    res = -1.0

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(wh%item)) then
      call c_f_pointer(wh%item,wp)
    else
      return
    end if
    if (c_associated(idvh%item)) then
      call c_f_pointer(idvh%item,idvp)
    else
      return
    end if

    res = psb_genrm2(xp,wp,idvp,descp,info)

  end function psb_c_dgenrm2_weightmask

  function psb_c_dgeamax(xh,cdh) bind(c) result(res)
    implicit none
    real(c_double) :: res

    type(psb_c_dvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp
    integer(psb_c_ipk_)          :: info

    res = -1.0
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if

    res = psb_geamax(xp,descp,info)

  end function psb_c_dgeamax

  function psb_c_dgemin(xh,cdh) bind(c) result(res)
    implicit none
    real(c_double) :: res

    type(psb_c_dvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp
    integer(psb_c_ipk_)          :: info

    res = -1.0
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if

    res = psb_gemin(xp,descp,info)

  end function psb_c_dgemin

  function psb_c_dgeasum(xh,cdh) bind(c) result(res)
    implicit none
    real(c_double) :: res

    type(psb_c_dvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer   :: descp
    type(psb_d_vect_type), pointer :: xp
    integer(psb_c_ipk_)          :: info

    res = -1.0

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if

    res = psb_geasum(xp,descp,info)

  end function psb_c_dgeasum


  function psb_c_dspnrmi(ah,cdh) bind(c) result(res)
    implicit none
    real(c_double) :: res

    type(psb_c_dspmat)   :: ah
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_dspmat_type), pointer :: ap
    integer(psb_c_ipk_)              ::  info

    res = -1.0
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if

    res = psb_spnrmi(ap,descp,info)

  end function psb_c_dspnrmi

  function psb_c_dgedot(xh,yh,cdh) bind(c) result(res)
    implicit none
    real(c_double) :: res

    type(psb_c_dvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    integer(psb_c_ipk_)               :: info

    res = -1.0
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if
    res = psb_gedot(xp,yp,descp,info)

  end function psb_c_dgedot


  function psb_c_dspmm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_dspmat) :: ah
    type(psb_c_dvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    real(c_double), value :: alpha, beta
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    type(psb_dspmat_type), pointer :: ap
    integer(psb_c_ipk_)               :: info

    res = -1
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if

    call psb_spmm(alpha,ap,xp,beta,yp,descp,info)

    res = info

  end function psb_c_dspmm


  function psb_c_dspmm_opt(alpha,ah,xh,beta,yh,cdh,trans,doswap) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_dspmat) :: ah
    type(psb_c_dvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    real(c_double), value :: alpha, beta
    character(c_char)      :: trans
    logical(c_bool), value :: doswap

    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    type(psb_dspmat_type), pointer :: ap
    character :: ftrans
    logical   :: fdoswap
    integer(psb_c_ipk_)   :: info

    res = -1
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if

    fdoswap = doswap
    ftrans  = trans
    call psb_spmm(alpha,ap,xp,beta,yp,descp,info,trans=ftrans,doswap=fdoswap)

    res = info

  end function psb_c_dspmm_opt


  function psb_c_dspsm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_dspmat) :: ah
    type(psb_c_dvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    real(c_double), value :: alpha, beta
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp,yp
    type(psb_dspmat_type), pointer :: ap
    integer(psb_c_ipk_)               :: info

    res = -1
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,xp)
    else
      return
    end if
    if (c_associated(yh%item)) then
      call c_f_pointer(yh%item,yp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if

    call psb_spsm(alpha,ap,xp,beta,yp,descp,info)

    res = info

  end function psb_c_dspsm

  function psb_c_dnnz(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_dspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_dspmat_type), pointer  :: ap
    integer(psb_c_ipk_)               :: info

    res = 0

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if

    res = psb_nnz(ap,descp,info)

  end function psb_c_dnnz

  function psb_c_dis_matupd(ah,cdh) bind(c) result(res)
    implicit none
    logical(c_bool) :: res

    type(psb_c_dspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_dspmat_type), pointer  :: ap
    integer(psb_c_ipk_)               :: info

    res = .false.

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if

    res = ap%is_upd()
  end function

  function psb_c_dis_matasb(ah,cdh) bind(c) result(res)
    implicit none
    logical(c_bool) :: res

    type(psb_c_dspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_dspmat_type), pointer  :: ap
    integer(psb_c_ipk_)               :: info

    res = .false.

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if

    res = ap%is_asb()
  end function

  function psb_c_dis_matbld(ah,cdh) bind(c) result(res)
    implicit none
    logical(c_bool) :: res

    type(psb_c_dspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_dspmat_type), pointer  :: ap
    integer(psb_c_ipk_)               :: info

    res = .false.

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if

    res = ap%is_bld()
  end function

  function psb_c_dset_matupd(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_dspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_dspmat_type), pointer  :: ap

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if

    call ap%set_upd()

    res = psb_success_
  end function

  function psb_c_dset_matasb(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_dspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_dspmat_type), pointer  :: ap

    res = -1;

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if

    call ap%set_asb()

    res = psb_success_

  end function

  function psb_c_dset_matbld(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_dspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_dspmat_type), pointer  :: ap

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if

    call ap%set_bld()

    res = psb_success_
  end function

  function psb_c_dcopy_mat(ah,bh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_dspmat)   :: ah,bh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_dspmat_type), pointer  :: ap,bp
    integer(psb_c_ipk_)               :: info

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if
    if (c_associated(bh%item)) then
      call c_f_pointer(bh%item,bp)
    else
      return
    end if

    call ap%clone(bp,info)

    res = info
  end function

  function psb_c_dspscal(alpha,ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)              ::  res

    real(c_double), value :: alpha
    type(psb_c_dspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_dspmat_type), pointer :: ap
    integer(psb_c_ipk_)              ::  info

    res = -1
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if

    call ap%scal(alpha,info)

    res = info

  end function psb_c_dspscal

  function psb_c_dspscalpid(alpha,ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)              ::  res

    real(c_double), value :: alpha
    type(psb_c_dspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_dspmat_type), pointer :: ap
    integer(psb_c_ipk_)              ::  info

    res = -1
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if

    call ap%scalpid(alpha,info)

    res = info

  end function psb_c_dspscalpid

  function psb_c_dspaxpby(alpha,ah,beta,bh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)              ::  res

    real(c_double), value :: alpha
    type(psb_c_dspmat)   :: ah
    real(c_double), value :: beta
    type(psb_c_dspmat)   :: bh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_dspmat_type), pointer :: ap,bp
    integer(psb_c_ipk_)              ::  info

    res = -1
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(ah%item)) then
      call c_f_pointer(ah%item,ap)
    else
      return
    end if
    if (c_associated(bh%item)) then
      call c_f_pointer(bh%item,bp)
    else
      return
    end if

    call ap%spaxpby(alpha,beta,bp,info)

    res = info
  end function psb_c_dspaxpby

end module psb_d_psblas_cbind_mod
