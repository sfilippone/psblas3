module psb_z_tools_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_cpenv_mod
  use psb_objhandle_mod
  use psb_base_tools_cbind_mod
#ifdef PSB_HAVE_CUDA
  use psb_cuda_mod
#endif
  
contains

  ! Should define   geall_opt  with DUPL argument
  function psb_c_zgeall(xh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
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
  end function psb_c_zgeall

  function psb_c_zgeall_remote(xh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
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
    call psb_geall(xp,descp,info,bldmode=psb_matbld_remote_,dupl=psb_dupl_add_)
    xh%item = c_loc(xp)
    res = min(0,info)

    return
  end function psb_c_zgeall_remote

  function psb_c_zgeall_remote_options(xh,cdh,bldmode,dupl) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    integer(psb_c_ipk_), value :: dupl
    integer(psb_c_ipk_), value :: bldmode


    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
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
    call psb_geall(xp,descp,info,bldmode=bldmode,dupl=dupl)
    xh%item = c_loc(xp)
    res = min(0,info)

    return
  end function psb_c_zgeall_remote_options

  function psb_c_zgeasb(xh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
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
  end function psb_c_zgeasb

  function psb_c_zgeasb_options(xh,cdh,dupl) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    integer(psb_c_ipk_), value :: dupl


    type(psb_desc_type), pointer :: descp
    type(psb_d_vect_type), pointer :: xp
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

    call psb_geasb(xp,descp,info,dupl=dupl)
    res = min(0,info)

    return
  end function psb_c_zgeasb_options

  function psb_c_zgeasb_options_format(xh,cdh,dupl,format) bind(c) result(res)
    ! Takes into account format argument as a c string, and uses it to call the appropriate psb_geasb
    ! with mold argument
    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh
    character(kind=c_char), dimension(*) :: format
    integer(psb_c_ipk_), value :: dupl

    ! Local variables
    character(len=6) :: fformat
    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
    integer(psb_c_ipk_)               :: info
    ! mold variables
#ifdef PSB_HAVE_CUDA
    type(psb_z_vect_cuda), target   :: vgpu
#endif
    type(psb_z_base_vect_type), target   :: vect
    class(psb_z_base_vect_type), pointer :: vmold

    ! Select mold based on format
    call psb_stringc2f(format,fformat)

    select case (psb_toupper(fformat))
#ifdef PSB_HAVE_CUDA  
    case('GPU','DEVICE')
      vmold => vgpu
#endif  
    case('CPU','HOST')
      vmold => vect
    case default
      write(psb_out_unit,*) 'psb_c_zgeasb_options_format: Unknown format ',fformat
      vmold => vect
    end select
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

    call psb_geasb(xp,descp,info,dupl=dupl,mold=vmold)
    res = min(0,info)

    return
  end function psb_c_zgeasb_options_format

 
  function psb_c_zgefree(xh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
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
  end function psb_c_zgefree


  function psb_c_zgeins(nz,irw,val,xh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    integer(psb_c_ipk_), value :: nz
    integer(psb_c_lpk_)        :: irw(*)
    complex(c_double_complex)        :: val(*)
    type(psb_c_zvector) :: xh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_z_vect_type), pointer :: xp
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
    select case(ixb)
    case (0)  
      !write(0,*) 'C_GEINS:  IDX_BASE',ixb,' :',irw(1:nz)+(1-ixb),val(1:nz)
      call psb_geins(nz,(irw(1:nz)+(1-ixb)),val(1:nz),&
           & xp,descp,info)
    case(1)
      !write(0,*) 'C_GEINS:  IDX_BASE',ixb,' :',irw(1:nz),val(1:nz)
      call psb_geins(nz,irw(1:nz),val(1:nz),&
           & xp,descp,info)
    case default
      write(0,*) 'C_GEINS: Unkonwn inndex base ',ixb
      info =-2
    end select

    res = min(0,info)

    return
  end function psb_c_zgeins

 function psb_c_zspall(mh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_zspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_zspmat_type), pointer :: ap
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
  end function psb_c_zspall


 function psb_c_zspall_remote(mh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_zspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_zspmat_type), pointer :: ap
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
    call psb_spall(ap,descp,info,bldmode=psb_matbld_remote_,dupl=psb_dupl_add_)
    mh%item = c_loc(ap)
    res = min(0,info)

    return
  end function psb_c_zspall_remote

 function psb_c_zspasb(mh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_zspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_zspmat_type), pointer :: ap
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
  end function psb_c_zspasb

  function psb_c_zspfree(mh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_zspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_zspmat_type), pointer :: ap
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
  end function psb_c_zspfree




  function psb_c_zspasb_opt(mh,cdh,afmt,upd,dupl) bind(c) result(res)

#if 0
#ifdef PSB_HAVE_LIBRSB 
    use psb_z_rsb_mat_mod
#endif
#endif
#if defined(PSB_HAVE_CUDA)
    use psb_cuda_mod
#endif
    use psb_ext_mod
    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_zspmat) :: mh
    type(psb_c_descriptor) :: cdh
    integer(psb_c_ipk_), value :: upd,dupl
    character(c_char)     :: afmt(*)
    integer(psb_c_ipk_)    :: info,n
    character(len=5)      :: fafmt
    integer(psb_ipk_), parameter :: hksz = 32
    ! mold variables
#if 0
#ifdef PSB_HAVE_LIBRSB
    type(psb_z_rsb_sparse_mat) :: arsb
#endif
#endif
    type(psb_z_ell_sparse_mat), target    :: aell
    type(psb_z_csr_sparse_mat), target   :: acsr
    type(psb_z_csc_sparse_mat), target   :: acsc
    type(psb_z_coo_sparse_mat), target   :: acoo
    type(psb_z_hll_sparse_mat), target    :: ahll
    type(psb_z_hdia_sparse_mat), target  :: ahdia
    type(psb_z_dns_sparse_mat), target   :: adns
#if defined(PSB_HAVE_CUDA)
    type(psb_z_cuda_hlg_sparse_mat), target   :: ahlg
    type(psb_z_cuda_csrg_sparse_mat), target   :: acsrg
    type(psb_z_cuda_elg_sparse_mat), target    :: aelg
#endif
    class(psb_z_base_sparse_mat), pointer :: amold
    !Local variables
    type(psb_desc_type), pointer :: descp
    type(psb_zspmat_type), pointer :: ap

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
    call psb_stringc2f(afmt,fafmt)

    ! Set the mold variable based on afmt
    select case (psb_toupper(fafmt))
#if defined(PSB_HAVE_CUDA)
    case('ELG')
      amold => aelg
    case('HLG')
      call psi_set_hksz(hksz)
      amold => ahlg
    case('CSRG')
      amold => acsrg
    case('ELL')
      amold => aell
    case('HLL')
      call psi_set_hksz(hksz)
      amold => ahll
    case('CSR')
      amold => acsr
    case('CSC')
      amold => acsc
    case('DNS')
      amold => adns
    case default
      write(*,*) 'Unknown format defaulting to HLG'
      amold => ahlg
#else
    case('ELL')
      amold => aell
    case('HLL')
      call psi_set_hksz(hksz)
      amold => ahll
      amold => ahdia
    case('CSR')
      amold => acsr
    case('CSC')
      amold => acsc
    case('DNS')
      amold => adns
    case default
      write(*,*) 'Unknown format defaulting to CSR'
      amold => acsr
#endif
  end select

    select case(fafmt)
#if 0
#ifdef PSB_HAVE_LIBRSB
    case('RSB')
      call psb_spasb(double_spmat_pool(mh)%item,descriptor_pool(cdh)%item,info,&
           & upd=upd,mold=arsb)
#endif
#endif
    case('ELL','HLL','CSR','DNS','CSC')
      call psb_spasb(ap,descp,info,upd=upd,mold=amold)
#if defined(PSB_HAVE_CUDA)
    case('ELG','HLG','CSRG')
      call psb_spasb(ap,descp,info,upd=upd,mold=amold)
#endif  
    case default
      write(psb_out_unit,*) 'psb_c_zspasb_opt: Unknown format ',fafmt
      call psb_spasb(ap,descp,info,afmt=fafmt,upd=upd,dupl=dupl)
    end select

    res = min(0,info)

    return
  end function psb_c_zspasb_opt


  function psb_c_zspins(nz,irw,icl,val,mh,cdh) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    integer(psb_c_ipk_), value :: nz
    integer(psb_c_lpk_)      :: irw(*), icl(*)
    complex(c_double_complex)        :: val(*)
    type(psb_c_zspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_zspmat_type), pointer :: ap
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
  end function psb_c_zspins


  function psb_c_zsprn(mh,cdh,clear) bind(c) result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    logical(c_bool), value :: clear
    type(psb_c_zspmat) :: mh
    type(psb_c_descriptor) :: cdh

    type(psb_desc_type), pointer :: descp
    type(psb_zspmat_type), pointer :: ap
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
  end function psb_c_zsprn
!!$
!!$  function psb_c_zspprint(mh) bind(c) result(res)
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
!!$  end function psb_c_zspprint

  function psb_c_zgetelem(xh,index,cdh) bind(c) result(res)
    implicit none

    type(psb_c_zvector)      :: xh
    integer(psb_c_lpk_), value :: index
    type(psb_c_descriptor)     :: cdh
    complex(c_double_complex)           :: res

    type(psb_z_vect_type), pointer :: xp
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

  end function psb_c_zgetelem

  function psb_c_zsetelem(index,val,xh,cdh) bind(c) result(res)
    implicit none

    type(psb_c_zvector)      :: xh
    integer(psb_c_lpk_), value :: index
    type(psb_c_descriptor)     :: cdh
    complex(c_double_complex), value    :: val
    integer(psb_c_ipk_) :: res

    type(psb_z_vect_type), pointer :: xp
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
      call  psb_setelem(index,val,xp,descp,info)
    else
      call psb_setelem(index+(1-ixb),val,xp,descp,info)
    end if
    res=info
    return

  end function psb_c_zsetelem

  function psb_c_zmatgetelem(ah,rowindex,colindex,cdh) bind(c) result(res)
    implicit none

    type(psb_c_zspmat)      :: ah
    integer(psb_c_lpk_), value :: rowindex, colindex
    type(psb_c_descriptor)     :: cdh
    complex(c_double_complex)           :: res

    type(psb_zspmat_type), pointer :: ap
    type(psb_desc_type), pointer     :: descp
    integer(psb_c_ipk_)              :: info, ixb

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

    ixb = psb_c_get_index_base()
    if (ixb == 1) then
      res = psb_getelem(ap,rowindex,colindex,descp,info)
    else
      res = psb_getelem(ap,rowindex+(1-ixb),colindex+(1-ixb),descp,info)
    end if

    return

  end function psb_c_zmatgetelem

end module psb_z_tools_cbind_mod
