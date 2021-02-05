module psb_base_tools_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_cpenv_mod
  use psb_base_string_cbind_mod

contains

  ! Aggiungere funzione per estrarre comunicatore

  function psb_c_error() bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res
    res = 0
    call psb_error()
  end function psb_c_error

  function psb_c_clean_errstack() bind(c) result(res)
    implicit none
    integer(psb_c_ipk_) :: res
    res = 0
    call psb_clean_errstack()
  end function psb_c_clean_errstack

  function psb_c_cdall_vg(ng,vg,cctxt,cdh) bind(c,name='psb_c_cdall_vg') result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    integer(psb_c_lpk_), value :: ng
    type(psb_c_object_type), value :: cctxt
    integer(psb_c_ipk_)        :: vg(*)
    type(psb_c_object_type) :: cdh
    type(psb_desc_type), pointer :: descp
    integer(psb_c_ipk_)           :: info
    type(psb_ctxt_type) :: ctxt
    ctxt = psb_c2f_ctxt(cctxt)

    res = -1
    if (ng <=0) then
      write(0,*) 'Invalid size'
      return
    end if

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
      call descp%free(info)
      if (info == 0) deallocate(descp,stat=info)
      if (info /= 0) return
    end if

    allocate(descp,stat=info)
    if (info < 0) return

    call psb_cdall(ctxt,descp,info,vg=vg(1:ng))
    cdh%item = c_loc(descp)
    res = info

  end function psb_c_cdall_vg


  function psb_c_cdall_vl(nl,vl,cctxt,cdh) bind(c,name='psb_c_cdall_vl') result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type), value :: cctxt
    integer(psb_c_ipk_), value :: nl
    integer(psb_c_lpk_)        :: vl(*)
    type(psb_c_object_type) :: cdh
    type(psb_desc_type), pointer :: descp
    integer(psb_c_ipk_)           :: info, ixb
    type(psb_ctxt_type) :: ctxt
    ctxt = psb_c2f_ctxt(cctxt)

    res = -1
    if (nl <=0) then
      write(0,*) 'Invalid size'
      return
    end if

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
      call descp%free(info)
      if (info == 0) deallocate(descp,stat=info)
      if (info /= 0) return
    end if

    allocate(descp,stat=info)
    if (info < 0) return

    ixb = psb_c_get_index_base()

    if (ixb == 1) then
      call psb_cdall(ctxt,descp,info,vl=vl(1:nl))
    else
      call psb_cdall(ctxt,descp,info,vl=(vl(1:nl)+(1-ixb)))
    end if
    cdh%item = c_loc(descp)
    res = info

  end function psb_c_cdall_vl

  function psb_c_cdall_nl(nl,cctxt,cdh) bind(c,name='psb_c_cdall_nl') result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type), value :: cctxt
    integer(psb_c_ipk_), value :: nl
    type(psb_c_object_type) :: cdh
    type(psb_desc_type), pointer :: descp
    integer(psb_c_ipk_)           :: info
    type(psb_ctxt_type) :: ctxt
    ctxt = psb_c2f_ctxt(cctxt)

    res = -1
    if (nl <=0) then
      write(0,*) 'Invalid size'
      return
    end if

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
      call descp%free(info)
      if (info == 0) deallocate(descp,stat=info)
      if (info /= 0) return
    end if

    allocate(descp,stat=info)
    if (info < 0) return

    call psb_cdall(ctxt,descp,info,nl=nl)
    cdh%item = c_loc(descp)
    res = info

  end function psb_c_cdall_nl

  function psb_c_cdall_repl(n,cctxt,cdh) bind(c,name='psb_c_cdall_repl') result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    integer(psb_c_lpk_), value :: n
    type(psb_c_object_type), value :: cctxt
    type(psb_c_object_type) :: cdh
    type(psb_desc_type), pointer :: descp
    integer(psb_c_ipk_)           :: info
    type(psb_ctxt_type) :: ctxt
    ctxt = psb_c2f_ctxt(cctxt)


    res = -1
    if (n <=0) then
      write(0,*) 'Invalid size'
      return
    end if

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
      call descp%free(info)
      if (info == 0) deallocate(descp,stat=info)
      if (info /= 0) return
    end if

    allocate(descp,stat=info)
    if (info < 0) return

    call psb_cdall(ctxt,descp,info,mg=n,repl=.true.)
    cdh%item = c_loc(descp)
    res = info

  end function psb_c_cdall_repl

  function psb_c_cdasb(cdh) bind(c,name='psb_c_cdasb') result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: cdh
    type(psb_desc_type), pointer :: descp
    integer(psb_c_ipk_)           :: info

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
      call psb_cdasb(descp,info)
      res = info
    end if

  end function psb_c_cdasb


 function psb_c_cdfree(cdh) bind(c,name='psb_c_cdfree') result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: cdh
    type(psb_desc_type), pointer :: descp
    integer(psb_c_ipk_)           :: info

    res = -1
    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
      call descp%free(info)
      if (info == 0) deallocate(descp,stat=info)
      if (info /= 0) return
      cdh%item = c_null_ptr
    end if

    res = info
    return
  end function psb_c_cdfree

  function psb_c_cdins(nz,ia,ja,cdh) bind(c,name='psb_c_cdins') result(res)

    implicit none
    integer(psb_c_ipk_) :: res
    integer(psb_c_ipk_), value   :: nz
    type(psb_c_object_type) :: cdh
    integer(psb_c_lpk_)          :: ia(*),ja(*)

    type(psb_desc_type), pointer :: descp
    integer(psb_c_ipk_)           :: ixb,info

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
      ixb = psb_c_get_index_base()
      if (ixb == 1) then
        call psb_cdins(nz,ia(1:nz),ja(1:nz),descp,info)
      else
        call psb_cdins(nz,(ia(1:nz)+(1-ixb)),(ja(1:nz)+(1-ixb)),descp,info)
      end if
        
      res = info
    end if
    return
  end function psb_c_cdins



  function psb_c_cd_get_local_rows(cdh) bind(c,name='psb_c_cd_get_local_rows') result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: cdh

    type(psb_desc_type), pointer :: descp
    integer               :: info

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
      res = descp%get_local_rows()

    end if

  end function psb_c_cd_get_local_rows



  function psb_c_cd_get_local_cols(cdh) bind(c,name='psb_c_cd_get_local_cols') result(res)
    implicit none

    integer(psb_c_ipk_) :: res
    type(psb_c_object_type) :: cdh

    type(psb_desc_type), pointer :: descp
    integer               :: info

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
      res = descp%get_local_cols()

    end if

  end function psb_c_cd_get_local_cols

  function psb_c_cd_get_global_rows(cdh) bind(c,name='psb_c_cd_get_global_rows') result(res)
    implicit none

    integer(psb_c_lpk_) :: res
    type(psb_c_object_type) :: cdh

    type(psb_desc_type), pointer :: descp
    integer               :: info

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
      res = descp%get_global_rows()

    end if

  end function psb_c_cd_get_global_rows



  function psb_c_cd_get_global_cols(cdh) bind(c,name='psb_c_cd_get_global_cols') result(res)
    implicit none

    integer(psb_c_lpk_) :: res
    type(psb_c_object_type) :: cdh

    type(psb_desc_type), pointer :: descp
    integer               :: info

    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)
      res = descp%get_global_cols()

    end if

  end function psb_c_cd_get_global_cols

  function psb_c_cd_get_global_indices(idx,nidx,owned,cdh) &
       & bind(c,name='psb_c_cd_get_global_indices') result(res)
    implicit none

    integer(psb_c_ipk_)            :: res
    type(psb_c_object_type)        :: cdh

    integer(psb_c_lpk_)            :: idx(nidx)
    integer(psb_c_ipk_), value     :: nidx
    logical(c_bool), value         :: owned


    type(psb_desc_type), pointer   :: descp
    integer(psb_lpk_), allocatable :: myidx(:)
    integer(psb_c_ipk_)            :: ixb
    logical                        :: fowned
    res = -1

    if (c_associated(cdh%item)) then
      call c_f_pointer(cdh%item,descp)

      fowned = owned
      myidx = descp%get_global_indices(owned=fowned)
      ixb = psb_c_get_index_base()
      idx(1:nidx) = myidx(1:nidx) - (1-ixb)
      res = 0

    end if

  end function psb_c_cd_get_global_indices


end module psb_base_tools_cbind_mod
