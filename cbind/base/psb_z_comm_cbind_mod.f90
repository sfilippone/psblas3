module psb_z_comm_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_base_string_cbind_mod
  
contains

  function psb_c_z_ovrl(xh,cdh) bind(c) result(res)
    implicit none 
    integer(psb_c_int) :: res

    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
    integer                 :: info
    

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

  end function psb_c_z_ovrl
 
  function psb_c_z_ovrl_opt(xh,cdh,update,mode) bind(c) result(res)
    implicit none 
    integer(psb_c_int) :: res
    integer(psb_c_int), value :: update, mode

    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
    integer                 :: info
    

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

  end function psb_c_z_ovrl_opt

 
  function psb_c_z_halo(xh,cdh) bind(c) result(res)
    implicit none 
    integer(psb_c_int) :: res

    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
    integer                 :: info
    

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

  end function psb_c_z_halo
 
  function psb_c_z_halo_opt(xh,cdh,tran,data,mode) bind(c) result(res)
    implicit none 
    integer(psb_c_int) :: res
    integer(psb_c_int), value :: data, mode
    character(c_char)      :: tran
        

    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
    character :: ftran
    integer                 :: info
    

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
    
  end function psb_c_z_halo_opt

  
  function psb_c_z_vscatter(ng,gx,xh,cdh) bind(c) result(res)
    implicit none 

    integer(psb_c_int)    :: res
    integer(psb_c_int), value :: ng
    complex(c_double_complex), target :: gx(*)
    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: vp
    complex(psb_dpk_), pointer :: pgx(:)
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

  end function psb_c_z_vscatter
  
  function psb_c_zvgather(v,xh,cdh) bind(c) result(res)
    implicit none 

    integer(psb_c_int)    :: res   
    complex(c_double_complex), target :: v(*)
    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: vp
    complex(psb_dpk_), allocatable :: fv(:)
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

    call psb_gather(fv,vp,descp,info)
    res = info 
    if (res /=0) return          
    sz = size(fv)
    v(1:sz) = fv(1:sz)
  end function psb_c_zvgather
    
  function psb_c_zspgather(gah,ah,cdh) bind(c) result(res)
    implicit none 

    integer(psb_c_int)    :: res   
    type(psb_c_zspmat)   :: ah, gah
    type(psb_c_descriptor) :: cdh
    
    type(psb_desc_type), pointer :: descp
    type(psb_zspmat_type), pointer :: ap, gap
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
  end function psb_c_zspgather
     
end module psb_z_comm_cbind_mod
