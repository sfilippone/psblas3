module psb_z_tools_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_cpenv_mod
  use psb_objhandle_mod
  use psb_base_string_cbind_mod
  use psb_base_tools_cbind_mod

contains

 function psb_c_zgeall(xh,cdh) bind(c) result(res)

    implicit none 
    integer(psb_c_ipk) :: res   
    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
    integer(psb_c_ipk)               :: info

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
  end function psb_c_zgeall

  function psb_c_zgeasb(xh,cdh) bind(c) result(res)

    implicit none 
    integer(psb_c_ipk) :: res   
    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
    integer(psb_c_ipk)               :: info

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
  end function psb_c_zgeasb
  
  function psb_c_zgefree(xh,cdh) bind(c) result(res)

    implicit none 
    integer(psb_c_ipk) :: res   
    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
    integer(psb_c_ipk)               :: info

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
    xh%item = c_null_ptr
    
    return
  end function psb_c_zgefree
  

 function psb_c_zgeins(nz,irw,val,xh,cdh) bind(c) result(res)

    implicit none 
    integer(psb_c_ipk) :: res   
    integer(psb_c_ipk), value :: nz
    integer(psb_c_lpk)        :: irw(*)
    complex(c_double_complex)        :: val(*)
    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
    integer(psb_c_ipk)               :: ixb, info

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
      call psb_geins(nz,irw(1:nz),val(1:nz),&
           & xp,descp,info, dupl=psb_dupl_ovwrt_)
    else
      call psb_geins(nz,(irw(1:nz)+(1-ixb)),val(1:nz),&
           & xp,descp,info, dupl=psb_dupl_ovwrt_)
    end if

    res = min(0,info)

    return
  end function psb_c_zgeins


 function psb_c_zgeins_add(nz,irw,val,xh,cdh) bind(c) result(res)

    implicit none 
    integer(psb_c_ipk) :: res   
    integer(psb_c_ipk), value :: nz
    integer(psb_c_lpk)       :: irw(*)
    complex(c_double_complex)        :: val(*)
    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
    integer(psb_c_ipk)               :: ixb, info

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
      call psb_geins(nz,irw(1:nz),val(1:nz),&
           & xp,descp,info, dupl=psb_dupl_add_)
    else
      call psb_geins(nz,(irw(1:nz)+(1-ixb)),val(1:nz),&
           & xp,descp,info, dupl=psb_dupl_add_)
    end if
    res = min(0,info)

    return
  end function psb_c_zgeins_add


 function psb_c_zspall(mh,cdh) bind(c) result(res)

    implicit none 
    integer(psb_c_ipk) :: res   
    type(psb_c_zspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_zspmat_type), pointer :: ap
    integer(psb_c_ipk)               :: info,n

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
  end function psb_c_zspall



 function psb_c_zspasb(mh,cdh) bind(c) result(res)
   
    implicit none 
    integer(psb_c_ipk) :: res   
    type(psb_c_zspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_zspmat_type), pointer :: ap
    integer(psb_c_ipk)               :: info,n

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
  end function psb_c_zspasb


  function psb_c_zspfree(mh,cdh) bind(c) result(res)
   
    implicit none 
    integer(psb_c_ipk) :: res   
    type(psb_c_zspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_zspmat_type), pointer :: ap
    integer(psb_c_ipk)               :: info,n

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
    mh%item=c_null_ptr
    return
  end function psb_c_zspfree


#if 0

  function psb_c_zspasb_opt(mh,cdh,afmt,upd,dupl) bind(c) result(res)

#ifdef HAVE_LIBRSB
    use psb_z_rsb_mat_mod
#endif
    implicit none 
    integer(psb_c_ipk) :: res   
    integer(psb_c_ipk), value :: cdh, mh,upd,dupl
    character(c_char)     :: afmt(*)
    integer(psb_c_ipk)    :: info,n, fdupl
    character(len=5)      :: fafmt
#ifdef HAVE_LIBRSB
    type(psb_z_rsb_sparse_mat) :: arsb
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
           & upd=upd,dupl=dupl,mold=arsb)
#endif
    case default
      call psb_spasb(double_spmat_pool(mh)%item,descriptor_pool(cdh)%item,info,&
           & afmt=fafmt,upd=upd,dupl=dupl)
    end select
    
    res = min(0,info)

    return
  end function psb_c_zspasb_opt
#endif

 function psb_c_zspins(nz,irw,icl,val,mh,cdh) bind(c) result(res)

    implicit none 
    integer(psb_c_ipk) :: res   
    integer(psb_c_ipk), value :: nz
    integer(psb_c_lpk)      :: irw(*), icl(*) 
    complex(c_double_complex)        :: val(*)
    type(psb_c_zspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_zspmat_type), pointer :: ap
    integer(psb_c_ipk)               :: ixb,info,n

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
  end function psb_c_zspins


  function psb_c_zsprn(mh,cdh,clear) bind(c) result(res)

    implicit none 
    integer(psb_c_ipk) :: res   
    logical(c_bool), value :: clear
    type(psb_c_zspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_zspmat_type), pointer :: ap
    integer(psb_c_ipk)     :: info
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
  end function psb_c_zsprn
!!$
!!$  function psb_c_zspprint(mh) bind(c) result(res)
!!$
!!$    implicit none 
!!$    integer(psb_c_ipk) :: res   
!!$    integer(psb_c_ipk),  value :: mh
!!$    integer(psb_c_ipk)         :: info
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
!!$  end function psb_c_zspprint


end module psb_z_tools_cbind_mod

