module psb_z_psblas_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_string_cbind_mod

contains

  function psb_c_zgeaxpby(alpha,xh,beta,yh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_zvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    complex(c_double_complex), value   :: alpha,beta

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
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

  end function psb_c_zgeaxpby

  function psb_c_zgeaxpbyz(alpha,xh,beta,yh,zh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_zvector) :: xh,yh,zh
    type(psb_c_descriptor) :: cdh
    complex(c_double_complex), value   :: alpha,beta

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp,zp
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

  end function psb_c_zgeaxpbyz

  function psb_c_zgemlt(xh,yh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_zvector)  :: xh,yh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
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

  end function psb_c_zgemlt

  function psb_c_zgemlt2(alpha,xh,yh,beta,zh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_zvector)  :: xh,yh, zh
    type(psb_c_descriptor) :: cdh


    type(psb_desc_type), pointer        :: descp
    type(psb_z_vect_type), pointer    :: xp,yp,zp
    integer(psb_c_ipk_)                 :: info
    complex(psb_dpk_), intent(in), value  :: alpha,beta

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

  end function psb_c_zgemlt2

  function psb_c_zgediv(xh,yh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_zvector)  :: xh,yh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
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

  end function psb_c_zgediv

  function psb_c_zgediv2(xh,yh,zh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_zvector)  :: xh,yh,zh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp,zp
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

  end function psb_c_zgediv2

  function psb_c_zgediv_check(xh,yh,cdh, flag) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_zvector)  :: xh,yh
    type(psb_c_descriptor) :: cdh
    logical(c_bool), value :: flag

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
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

  end function psb_c_zgediv_check

  function psb_c_zgediv2_check(xh,yh,zh,cdh, flag) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_zvector)  :: xh,yh,zh
    type(psb_c_descriptor) :: cdh
    logical(c_bool), value :: flag

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp,zp
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

  end function psb_c_zgediv2_check

  function psb_c_zgeinv(xh,yh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_zvector)  :: xh,yh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
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

  end function psb_c_zgeinv

  function psb_c_zgeinv_check(xh,yh,cdh, flag) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_zvector)  :: xh,yh
    type(psb_c_descriptor) :: cdh
    logical(c_bool), value :: flag

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
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

  end function psb_c_zgeinv_check

  function psb_c_zgeabs(xh,yh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_zvector)  :: xh,yh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
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

  end function psb_c_zgeabs

  function psb_c_zgecmp(xh,ch,zh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_zvector)  :: xh,zh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,zp
    integer(psb_c_ipk_)          :: info
    real(c_double_complex), value :: ch

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

  end function psb_c_zgecmp

  function psb_c_zgeaddconst(xh,bh,zh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_)    :: res

    type(psb_c_zvector)  :: xh,zh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,zp
    integer(psb_c_ipk_)          :: info
    real(c_double_complex), value :: bh

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

  end function psb_c_zgeaddconst


  function psb_c_zgenrm2(xh,cdh) bind(c) result(res)
    implicit none
    real(c_double_complex) :: res

    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
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

  end function psb_c_zgenrm2

  function psb_c_zgenrmi(xh,cdh) bind(c) result(res)
    implicit none
    real(c_double_complex) :: res

    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
    type(psb_z_vect_type) ::  yp
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

  end function psb_c_zgenrmi

  function psb_c_zgenrm2_weight(xh,wh,cdh) bind(c) result(res)
    implicit none
    real(c_double_complex) :: res

    type(psb_c_zvector) :: xh, wh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp, wp
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

  end function psb_c_zgenrm2_weight

  function psb_c_zgenrm2_weightmask(xh,wh,idvh,cdh) bind(c) result(res)
    implicit none
    real(c_double_complex) :: res

    type(psb_c_zvector) :: xh, wh, idvh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp, wp, idvp
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

  end function psb_c_zgenrm2_weightmask

  function psb_c_zgeamax(xh,cdh) bind(c) result(res)
    implicit none
    real(c_double_complex) :: res

    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
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

  end function psb_c_zgeamax


  function psb_c_zgeasum(xh,cdh) bind(c) result(res)
    implicit none
    real(c_double_complex) :: res

    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer   :: descp
    type(psb_z_vect_type), pointer :: xp
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

  end function psb_c_zgeasum


  function psb_c_zspnrmi(ah,cdh) bind(c) result(res)
    implicit none
    real(c_double_complex) :: res

    type(psb_c_zspmat)   :: ah
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_zspmat_type), pointer :: ap
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

  end function psb_c_zspnrmi

  function psb_c_zgedot(xh,yh,cdh) bind(c) result(res)
    implicit none
    complex(c_double_complex) :: res

    type(psb_c_zvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
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

  end function psb_c_zgedot


  function psb_c_zspmm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_zspmat) :: ah
    type(psb_c_zvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    complex(c_double_complex), value :: alpha, beta
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
    type(psb_zspmat_type), pointer :: ap
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

  end function psb_c_zspmm


  function psb_c_zspmm_opt(alpha,ah,xh,beta,yh,cdh,trans,doswap) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_zspmat) :: ah
    type(psb_c_zvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    complex(c_double_complex), value :: alpha, beta
    character(c_char)      :: trans
    logical(c_bool), value :: doswap

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
    type(psb_zspmat_type), pointer :: ap
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

  end function psb_c_zspmm_opt


  function psb_c_zspsm(alpha,ah,xh,beta,yh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_zspmat) :: ah
    type(psb_c_zvector) :: xh,yh
    type(psb_c_descriptor) :: cdh
    complex(c_double_complex), value :: alpha, beta
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp,yp
    type(psb_zspmat_type), pointer :: ap
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

  end function psb_c_zspsm

  function psb_c_znnz(ah,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_zspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_zspmat_type), pointer  :: ap
    integer(psb_c_ipk_)               :: info

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

  end function psb_c_znnz

  function psb_c_zis_matupd(ah,cdh) bind(c) result(res)
    implicit none
    logical :: res

    type(psb_c_zspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_zspmat_type), pointer  :: ap
    integer(psb_c_ipk_)               :: info

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

    res = psb_is_matupd(ap,descp,info)
  end function

  function psb_c_zis_matasb(ah,cdh) bind(c) result(res)
    implicit none
    logical :: res

    type(psb_c_zspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_zspmat_type), pointer  :: ap
    integer(psb_c_ipk_)               :: info

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

    res = psb_is_matasb(ap,descp,info)
  end function

  function psb_c_zis_matbld(ah,cdh) bind(c) result(res)
    implicit none
    logical :: res

    type(psb_c_zspmat)   :: ah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer      :: descp
    type(psb_zspmat_type), pointer  :: ap
    integer(psb_c_ipk_)               :: info

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

    res = psb_is_matbld(ap,descp,info)
  end function

end module psb_z_psblas_cbind_mod
