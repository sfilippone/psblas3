module psb_c_comm_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_string_cbind_mod

contains

  function psb_c_covrl(xh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_cvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_c_vect_type), pointer :: xp
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

    call psb_ovrl(xp,descp,info)

    res = info

  end function psb_c_covrl

  function psb_c_covrl_opt(xh,cdh,update,mode) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res
    integer(psb_c_ipk_), value :: update, mode

    type(psb_c_cvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_c_vect_type), pointer :: xp
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

    call psb_ovrl(xp,descp,info,update=update,mode=mode)

    res = info

  end function psb_c_covrl_opt


  function psb_c_chalo(xh,cdh) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res

    type(psb_c_cvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_c_vect_type), pointer :: xp
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

    call psb_halo(xp,descp,info)

    res = info

  end function psb_c_chalo

  function psb_c_chalo_opt(xh,cdh,tran,data,mode) bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res
    integer(psb_c_ipk_), value :: data, mode
    character(c_char)      :: tran


    type(psb_c_cvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_c_vect_type), pointer :: xp
    character :: ftran
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

    ftran  = tran
    call psb_halo(xp,descp,info,data=data,mode=mode,tran=ftran)

    res = info

  end function psb_c_chalo_opt


  function psb_c_cvscatter(ng,gx,xh,cdh) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_)    :: res
    integer(psb_c_lpk_), value :: ng
    complex(c_float_complex), target :: gx(*)
    type(psb_c_cvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_c_vect_type), pointer :: vp
    complex(psb_spk_), pointer :: pgx(:)
    integer(psb_c_ipk_)       :: info, sz

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
    else
      return
    end if

    pgx => gx(1:ng)

    call psb_scatter(pgx,vp,descp,info)
    res = info

  end function psb_c_cvscatter

  function psb_c_cvgather_f(v,xh,cdh) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_)    :: res
    complex(c_float_complex), target :: v(*)
    type(psb_c_cvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_c_vect_type), pointer :: vp
    complex(psb_spk_), allocatable :: fv(:)
    integer(psb_c_ipk_)           :: info, sz

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(xh%item)) then
      call c_f_pointer(xh%item,vp)
    else
      return
    end if

    call psb_gather(fv,vp,descp,info)
    res = info
    if (res /=0) return
    sz = size(fv)
    v(1:sz) = fv(1:sz)
  end function psb_c_cvgather_f

  function psb_c_cspgather_f(gah,ah,cdh) bind(c) result(res)
    implicit none

    integer(psb_c_ipk_)    :: res
    type(psb_c_cspmat)   :: ah, gah
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_cspmat_type), pointer :: ap, gap
    integer(psb_c_ipk_)               :: info, sz

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
    if (c_associated(gah%item)) then
      call c_f_pointer(gah%item,gap)
    else
      return
    end if
    call psb_gather(gap,ap,descp,info)
    res = info
  end function psb_c_cspgather_f

end module psb_c_comm_cbind_mod
