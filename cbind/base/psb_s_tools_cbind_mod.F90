module psb_s_tools_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_string_cbind_mod
  use psb_base_tools_cbind_mod

contains

 function psb_c_sgeall(xh,cdh) bind(c) result(res)

    implicit none 
    integer(psb_c_int) :: res   
    type(psb_c_svector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_s_vect_type), pointer :: xp
    integer               :: info

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
  end function psb_c_sgeall

  function psb_c_sgeasb(xh,cdh) bind(c) result(res)

    implicit none 
    integer(psb_c_int) :: res   
    type(psb_c_svector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_s_vect_type), pointer :: xp
    integer               :: info

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
  end function psb_c_sgeasb
  
  function psb_c_sgefree(xh,cdh) bind(c) result(res)

    implicit none 
    integer(psb_c_int) :: res   
    type(psb_c_svector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_s_vect_type), pointer :: xp
    integer               :: info

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
  end function psb_c_sgefree
  

 function psb_c_sgeins(nz,irw,val,xh,cdh) bind(c) result(res)

    implicit none 
    integer(psb_c_int) :: res   
    integer(psb_c_int), value :: nz
    integer(psb_c_int)        :: irw(*)
    real(c_float)        :: val(*)
    type(psb_c_svector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_s_vect_type), pointer :: xp
    integer               :: info

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

    call psb_geins(nz,irw(1:nz),val(1:nz),&
         & xp,descp,info, dupl=psb_dupl_ovwrt_)
    res = min(0,info)

    return
  end function psb_c_sgeins


 function psb_c_sgeins_add(nz,irw,val,xh,cdh) bind(c) result(res)

    implicit none 
    integer(psb_c_int) :: res   
    integer(psb_c_int), value :: nz
    integer(psb_c_int)        :: irw(*)
    real(c_float)        :: val(*)
    type(psb_c_svector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_s_vect_type), pointer :: xp
    integer               :: info

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

    call psb_geins(nz,irw(1:nz),val(1:nz),&
         & xp,descp,info, dupl=psb_dupl_add_)
    res = min(0,info)

    return
  end function psb_c_sgeins_add


 function psb_c_sspall(mh,cdh) bind(c) result(res)

    implicit none 
    integer(psb_c_int) :: res   
    type(psb_c_sspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_sspmat_type), pointer :: ap
    integer               :: info,n

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
  end function psb_c_sspall



 function psb_c_sspasb(mh,cdh) bind(c) result(res)
   
    implicit none 
    integer(psb_c_int) :: res   
    type(psb_c_sspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_sspmat_type), pointer :: ap
    integer               :: info,n

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
  end function psb_c_sspasb


  function psb_c_sspfree(mh,cdh) bind(c) result(res)
   
    implicit none 
    integer(psb_c_int) :: res   
    type(psb_c_sspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_sspmat_type), pointer :: ap
    integer               :: info,n

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
  end function psb_c_sspfree


#if 0

  function psb_c_sspasb_opt(mh,cdh,afmt,upd,dupl) bind(c) result(res)

#ifdef HAVE_LIBRSB
    use psb_s_rsb_mat_mod
#endif
    implicit none 
    integer(psb_c_int) :: res   
    integer(psb_c_int), value :: cdh, mh,upd,dupl
    character(c_char)     :: afmt(*)
    integer               :: info,n, fdupl
    character(len=5)      :: fafmt
#ifdef HAVE_LIBRSB
    type(psb_s_rsb_sparse_mat) :: arsb
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
  end function psb_c_sspasb_opt
#endif

 function psb_c_sspins(nz,irw,icl,val,mh,cdh) bind(c) result(res)

    implicit none 
    integer(psb_c_int) :: res   
    integer(psb_c_int), value :: nz
    integer(psb_c_int)        :: irw(*), icl(*) 
    real(c_float)        :: val(*)
    type(psb_c_sspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_sspmat_type), pointer :: ap
    integer               :: info,n

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


    call psb_spins(nz,irw(1:nz),icl(1:nz),val(1:nz),ap,descp,info)
    res = min(0,info)
    return
  end function psb_c_sspins


  function psb_c_ssprn(mh,cdh,clear) bind(c) result(res)

    implicit none 
    integer(psb_c_int) :: res   
    logical(c_bool), value :: clear
    type(psb_c_sspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_sspmat_type), pointer :: ap
    integer                :: info
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
  end function psb_c_ssprn
!!$
!!$  function psb_c_sspprint(mh) bind(c) result(res)
!!$
!!$    implicit none 
!!$    integer(psb_c_int) :: res   
!!$    integer(psb_c_int),  value :: mh
!!$    integer                :: info
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
!!$  end function psb_c_sspprint


end module psb_s_tools_cbind_mod

