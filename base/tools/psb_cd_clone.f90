  subroutine psb_cd_clone(desc, desc_out, info)

    use psb_error_mod
    use psb_penv_mod
    use psb_realloc_mod
    use psb_desc_mod, psb_protect_name => psb_cd_clone
    implicit none
    !....parameters...

    class(psb_desc_type), intent(inout), target :: desc
    class(psb_desc_type), intent(inout)         :: desc_out
    integer(psb_ipk_), intent(out)              :: info
    !locals
    integer(psb_ipk_) :: np,me,ictxt, err_act
    integer(psb_ipk_) :: debug_level, debug_unit
    character(len=20)   :: name

    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    if (psb_get_errstatus() /= 0) return 
    info = psb_success_
    call psb_erractionsave(err_act)
    name = 'psb_cdcpy'

    if (desc%is_valid()) then 
      ictxt = desc%get_context()

      ! check on blacs grid 
      call psb_info(ictxt, me, np)
      if (debug_level >= psb_debug_ext_) &
           & write(debug_unit,*) me,' ',trim(name),': Entered'
      if (np == -1) then
        info = psb_err_context_error_
        call psb_errpush(info,name)
        goto 9999
      endif

      desc_out%base_desc  => desc%base_desc
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%halo_index,desc_out%halo_index,info)
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%ext_index,desc_out%ext_index,info)
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%ovrlap_index,&
           & desc_out%ovrlap_index,info)
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%bnd_elem,desc_out%bnd_elem,info)
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%ovrlap_elem,desc_out%ovrlap_elem,info)
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%ovr_mst_idx,desc_out%ovr_mst_idx,info)
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%lprm,desc_out%lprm,info)
      if (info == psb_success_)&
           & call psb_safe_ab_cpy(desc%idx_space,desc_out%idx_space,info)
!!$      if ((info == psb_success_).and.(allocated(desc%indxmap))) &
!!$           & call desc%indxmap%clone(desc_out%indxmap,info)
!!$      associate(indxin => desc%indxmap) 
!!$        if ((info == psb_success_).and.(allocated(desc%indxmap))) &
!!$             & call indxin%clone(desc_out%indxmap,info)
!!$      end associate
      if ((info == psb_success_).and.(allocated(desc%indxmap))) &
           & allocate(desc_out%indxmap,source=desc%indxmap,stat=info)

    else
      call desc_out%free(info)
    end if
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info,name)
      goto 9999
    endif
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),': Done'

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_ret_) then
      return
    else
      call psb_error(ictxt)
    end if
    return

  end subroutine psb_cd_clone
