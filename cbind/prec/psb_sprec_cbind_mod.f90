module psb_sprec_cbind_mod

  use iso_c_binding
  use psb_prec_mod, only : psb_sprec_type
  use psb_objhandle_mod
  use psb_base_string_cbind_mod

  type, bind(c) :: psb_c_sprec
    type(c_ptr) :: item = c_null_ptr
  end type psb_c_sprec
  
  
contains 


  function  psb_c_sprecinit(ictxt,ph,ptype) bind(c) result(res)
    use psb_base_mod
    use psb_prec_mod
    use psb_base_string_cbind_mod
    implicit none 
    integer(psb_c_ipk_)          :: res
    integer(psb_c_ipk_), value :: ictxt   

    type(psb_c_sprec) :: ph
    character(c_char)       :: ptype(*)
    type(psb_sprec_type), pointer :: precp
    integer(psb_c_ipk_)              :: info
    character(len=80)       :: fptype

    res = -1
    if (c_associated(ph%item)) then 
      return 
    end if

    allocate(precp,stat=info)
    if (info /= 0) return
    ph%item = c_loc(precp)

    call stringc2f(ptype,fptype)
    
    call psb_precinit(ictxt,precp,fptype,info) 
    
    res = min(0,info)
    return
  end function psb_c_sprecinit
    


  function  psb_c_sprecbld(ah,cdh,ph) bind(c) result(res)
    use psb_base_mod
    use psb_prec_mod
    use psb_objhandle_mod
    use psb_base_string_cbind_mod
    implicit none 
    
    integer(psb_c_ipk_) :: res
    type(psb_c_sspmat) :: ah
    type(psb_c_sprec) :: ph
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_sspmat_type), pointer :: ap
    type(psb_sprec_type), pointer :: precp
    integer(psb_c_ipk_)              :: info

    res = -1
!!$    write(*,*) 'Entry:   ', psb_c_cd_get_local_rows(cdh)
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
    if (c_associated(ph%item)) then 
      call c_f_pointer(ph%item,precp)
    else
      return 
    end if

    call psb_precbld(ap,descp, precp, info)

    res = min(info,0)
    
  end function psb_c_sprecbld


  function  psb_c_sprecfree(ph) bind(c) result(res)
    use psb_base_mod
    use psb_prec_mod
    use psb_objhandle_mod
    use psb_base_string_cbind_mod
    implicit none 
    
    integer(psb_c_ipk_) :: res
    type(psb_c_sprec) :: ph
    type(psb_sprec_type), pointer :: precp
    integer(psb_c_ipk_)              :: info

    res = -1
    if (c_associated(ph%item)) then 
      call c_f_pointer(ph%item,precp)
    else
      return 
    end if
    
    call psb_precfree(precp, info)

    res = min(info,0)
    
  end function psb_c_sprecfree


end module psb_sprec_cbind_mod
