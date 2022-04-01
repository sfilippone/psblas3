module psb_c_tools_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_cpenv_mod
  use psb_objhandle_mod
  use psb_base_string_cbind_mod
  use psb_base_tools_cbind_mod
  
contains

  ! Should define   geall_opt  with DUPL argument
  function psb_c_cgeall(xh,cdh) bind(c) result(res)

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
      return
    end if
    allocate(xp)
    call psb_geall(xp,descp,info)
    xh%item = c_loc(xp)
    res = min(0,info)

    return
  end function psb_c_cgeall

  function psb_c_cgeasb(xh,cdh) bind(c) result(res)

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

    call psb_geasb(xp,descp,info)
    res = min(0,info)

    return
  end function psb_c_cgeasb

  function psb_c_cgefree(xh,cdh) bind(c) result(res)

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

    call psb_gefree(xp,descp,info)
    res = min(0,info)
    deallocate(xp,stat=info)
    res = min(0,info)
    xh%item = c_null_ptr

    return
  end function psb_c_cgefree


  function psb_c_cgeins(nz,irw,val,xh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    integer(psb_c_ipk_), value :: nz
    integer(psb_c_lpk_)        :: irw(*)
    complex(c_float_complex)        :: val(*)
    type(psb_c_cvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_c_vect_type), pointer :: xp
    integer(psb_c_ipk_)               :: ixb, info

    res = -1
    info = 0 
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

    ixb = psb_c_get_index_base()
    if (ixb == 1) then
      call psb_geins(nz,irw(1:nz),val(1:nz),&
           & xp,descp,info)
    else
      call psb_geins(nz,(irw(1:nz)+(1-ixb)),val(1:nz),&
           & xp,descp,info)
    end if

    res = min(0,info)

    return
  end function psb_c_cgeins

 function psb_c_cspall(mh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_cspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_cspmat_type), pointer :: ap
    integer(psb_c_ipk_)               :: info,n

    res = -1
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(mh%item)) then
      return
    end if
    allocate(ap)
    call psb_spall(ap,descp,info)
    mh%item = c_loc(ap)
    res = min(0,info)

    return
  end function psb_c_cspall


 function psb_c_cspall_remote(mh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_cspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_cspmat_type), pointer :: ap
    integer(psb_c_ipk_)               :: info,n

    res = -1
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(mh%item)) then
      return
    end if
    allocate(ap)
    call psb_spall(ap,descp,info,bldmode=psb_matbld_remote_)
    mh%item = c_loc(ap)
    res = min(0,info)

    return
  end function psb_c_cspall_remote

 function psb_c_cspasb(mh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_cspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_cspmat_type), pointer :: ap
    integer(psb_c_ipk_)               :: info,n

    res = -1
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(mh%item)) then
      call c_f_pointer(mh%item,ap)
    else
      return
    end if

    call psb_spasb(ap,descp,info)
    res = min(0,info)
    return
  end function psb_c_cspasb

  function psb_c_cspfree(mh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_cspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_cspmat_type), pointer :: ap
    integer(psb_c_ipk_)               :: info,n

    res = -1
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(mh%item)) then
      call c_f_pointer(mh%item,ap)
    else
      return
    end if

    call psb_spfree(ap,descp,info)
    res = min(0,info)
    deallocate(ap,stat=info)
    mh%item=c_null_ptr
    return
  end function psb_c_cspfree


#if 0

  function psb_c_cspasb_opt(mh,cdh,afmt,upd) bind(c) result(res)

#ifdef HAVE_LIBRSB
    use psb_c_rsb_mat_mod
#endif
    implicit none
    integer(psb_c_ipk_) :: res
    integer(psb_c_ipk_), value :: cdh, mh,upd,dupl
    character(c_char)     :: afmt(*)
    integer(psb_c_ipk_)    :: info,n
    character(len=5)      :: fafmt
#ifdef HAVE_LIBRSB
    type(psb_c_rsb_sparse_mat) :: arsb
#endif

    res = -1
    call psb_check_descriptor_handle(cdh,info)
    if (info < 0) return
    call psb_check_double_spmat_handle(mh,info)
    if (info < 0) return

    call stringc2f(afmt,fafmt)
    select case(fafmt)
#ifdef HAVE_LIBRSB
    case('RSB')
      call psb_spasb(double_spmat_pool(mh)%item,descriptor_pool(cdh)%item,info,&
           & upd=upd,mold=arsb)
#endif
    case default
      call psb_spasb(double_spmat_pool(mh)%item,descriptor_pool(cdh)%item,info,&
           & afmt=fafmt,upd=upd)
    end select

    res = min(0,info)

    return
  end function psb_c_cspasb_opt
#endif

 function psb_c_cspins(nz,irw,icl,val,mh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    integer(psb_c_ipk_), value :: nz
    integer(psb_c_lpk_)      :: irw(*), icl(*)
    complex(c_float_complex)        :: val(*)
    type(psb_c_cspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_cspmat_type), pointer :: ap
    integer(psb_c_ipk_)               :: ixb,info,n

    res = -1
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(mh%item)) then
      call c_f_pointer(mh%item,ap)
    else
      return
    end if

    ixb = psb_c_get_index_base()
    if (ixb == 1) then
      call psb_spins(nz,irw(1:nz),icl(1:nz),val(1:nz),ap,descp,info)
    else
      call psb_spins(nz,(irw(1:nz)+(1-ixb)),(icl(1:nz)+(1-ixb)),val(1:nz),ap,descp,info)
    end if
    res = min(0,info)
    return
  end function psb_c_cspins


  function psb_c_csprn(mh,cdh,clear) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    logical(c_bool), value :: clear
    type(psb_c_cspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_cspmat_type), pointer :: ap
    integer(psb_c_ipk_)     :: info
    logical                :: fclear

    res = -1
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
    else
      return
    end if
    if (c_associated(mh%item)) then
      call c_f_pointer(mh%item,ap)
    else
      return
    end if

    fclear = clear
    call psb_sprn(ap,descp,info,clear=fclear)
    res = min(0,info)

    return
  end function psb_c_csprn
!!$
!!$  function psb_c_cspprint(mh) bind(c) result(res)
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
!!$  end function psb_c_cspprint

  function psb_c_cgetelem(xh,index,cdh) bind(c) result(res)
    implicit none

    type(psb_c_cvector)      :: xh
    integer(psb_c_lpk_), value :: index
    type(psb_c_descriptor)     :: cdh
    complex(c_float_complex)           :: res

    type(psb_c_vect_type), pointer :: xp
    type(psb_desc_type), pointer     :: descp
    integer(psb_c_ipk_)              :: info, ixb

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

    ixb = psb_c_get_index_base()
    if (ixb == 1) then
      res = psb_getelem(xp,index,descp,info)
    else
      res = psb_getelem(xp,index+(1-ixb),descp,info)
    end if

    return

  end function psb_c_cgetelem

end module psb_c_tools_cbind_mod
