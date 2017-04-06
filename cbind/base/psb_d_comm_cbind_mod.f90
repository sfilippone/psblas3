module psb_d_comm_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_string_cbind_mod
  use psb_base_tools_cbind_mod 

contains

  
  function psb_c_dhalo(xh,cdh) bind(c) result(res)
    implicit none 

    integer(psb_c_int)    :: res   
    type(psb_c_dvector) :: xh
    type(psb_c_descriptor) :: cdh
    
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: vp
    integer               :: info, sz

    res = -1

    if (c_associated(cdh%item)) then 
      call c_f_pointer(cdh%item,descp)
    else
      return 
    end if
    if (c_associated(xh%item)) then 
      call c_f_pointer(xh%item,vp)
      call psb_halo(vp,descp,info)
      res = info 
    end if
    
  end function psb_c_dhalo
  
  function psb_c_dhalo_opt(xh,cdh,trans,mode) bind(c) result(res)
    implicit none 

    integer(psb_c_int)    :: res   
    type(psb_c_dvector) :: xh
    type(psb_c_descriptor) :: cdh
    character(c_char)   :: trans
    integer(psb_c_int), value  :: mode
    
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: vp
    character             :: trans_
    integer(psb_ipk_)     :: mode_
    integer               :: info, sz

    res = -1

    if (c_associated(cdh%item)) then 
      call c_f_pointer(cdh%item,descp)
    else
      return 
    end if
    trans_ = trans
    mode_  = mode
    if (c_associated(xh%item)) then 
      call c_f_pointer(xh%item,vp)
      call psb_halo(vp,descp,info,tran=trans_, mode=mode_)
      res = info 
    end if
    
  end function psb_c_dhalo_opt
  
  function psb_c_dovrl(xh,cdh) bind(c) result(res)
    implicit none 

    integer(psb_c_int)    :: res   
    type(psb_c_dvector) :: xh
    type(psb_c_descriptor) :: cdh
    
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: vp
    integer               :: info, sz

    res = -1

    if (c_associated(cdh%item)) then 
      call c_f_pointer(cdh%item,descp)
    else
      return 
    end if
    if (c_associated(xh%item)) then 
      call c_f_pointer(xh%item,vp)
      call psb_ovrl(vp,descp,info)
      res = info 
    end if
    
  end function psb_c_dovrl
  
  function psb_c_dovrl_opt(xh,cdh,update,mode) bind(c) result(res)
    implicit none 

    integer(psb_c_int)    :: res   
    type(psb_c_dvector) :: xh
    type(psb_c_descriptor) :: cdh
    integer(psb_c_int), value :: update, mode
    
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: vp
    integer(psb_ipk_)     :: mode_, update_
    integer               :: info, sz

    res = -1

    if (c_associated(cdh%item)) then 
      call c_f_pointer(cdh%item,descp)
    else
      return 
    end if
    update_ = update
    mode_   = mode
    if (c_associated(xh%item)) then 
      call c_f_pointer(xh%item,vp)
      call psb_ovrl(vp,descp,info,update=update_,mode=mode_)
      res = info 
    end if
    
  end function psb_c_dovrl_opt
  
  function psb_c_dvscatter(ng,gx,xh,cdh) bind(c) result(res)
    implicit none 

    integer(psb_c_int)    :: res
    integer(psb_c_int), value :: ng
    real(c_double), target :: gx(*)
    type(psb_c_dvector) :: xh
    type(psb_c_descriptor) :: cdh
    
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: vp
    real(psb_dpk_), pointer :: pgx(:)
    integer               :: info, sz

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

  end function psb_c_dvscatter
  
  function psb_c_dvgather_f(v,xh,cdh) bind(c) result(res)
    implicit none 

    integer(psb_c_int)    :: res   
    real(c_double)    :: v(*)
    type(psb_c_dvector) :: xh
    type(psb_c_descriptor) :: cdh
    
    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: vp
    real(psb_dpk_), allocatable :: fv(:)
    integer               :: info, sz

    res = -1

    if (c_associated(cdh%item)) then 
      call c_f_pointer(cdh%item,descp)
    else
      return 
    end if
    if (c_associated(xh%item)) then 
      call c_f_pointer(xh%item,vp)
      call psb_gather(fv,vp,descp,info)
      res = info 
      if (res /=0) return          
      sz = size(fv)
      v(1:sz) = fv(1:sz)
    end if
  end function psb_c_dvgather_f
    
  function psb_c_dspgather_f(gah,ah,cdh) bind(c) result(res)
    implicit none 

    integer(psb_c_int)    :: res   
    type(psb_c_dspmat)   :: ah, gah
    type(psb_c_descriptor) :: cdh
    
    type(psb_desc_type), pointer :: descp
    type(psb_dspmat_type), pointer :: ap, gap
    integer               :: info, sz

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
  end function psb_c_dspgather_f
    
  
end module psb_d_comm_cbind_mod

