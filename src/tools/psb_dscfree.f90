! File: psb_dscfree.f90
!
! Subroutine: psb_dscfree
!   Frees a descriptor data structure.
! 
! Parameters: 
!    desc_a   - type(<psb_desc_type>).         The communication descriptor to be freed.
!    info     - integer.                       Eventually returns an error code.
subroutine psb_dscfree(desc_a,info)
  !...free descriptor structure...
  use psb_descriptor_type
  use psb_const_mod
  use psb_error_mod
  implicit none
  !....parameters...
  type(psb_desc_type), intent(inout) :: desc_a
  integer, intent(out)               :: info
  !...locals....
  integer             :: int_err(5)
  integer             :: temp(1)
  real(kind(1.d0))    :: real_err(5)
  integer             :: icontxt,nprow,npcol,me,mypcol, err_act
  character(len=20)   :: name, char_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)
  name = 'psb_dscfree'


  if (.not.associated(desc_a%matrix_data)) then
     info=295
     call psb_errpush(info,name)
     return
  end if

  icontxt=desc_a%matrix_data(psb_ctxt_)
  deallocate(desc_a%matrix_data)
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
  !     ....verify blacs grid correctness..
  if (nprow.eq.-1) then
     info = 2010
     call psb_errpush(info,name)
     goto 9999
  else if (npcol.ne.1) then
     info = 2030
     int_err(1) = npcol
     call psb_errpush(info,name,int_err)
     goto 9999
  endif

  !...deallocate desc_a....
  if(.not.associated(desc_a%loc_to_glob)) then
     info=295
     call psb_errpush(info,name)
     goto 9999
  end if

    !deallocate loc_to_glob  field
    deallocate(desc_a%loc_to_glob,stat=info)
    if (info /= 0) then
       info=2051
       call psb_errpush(info,name)
       goto 9999
    end if

  if (.not.associated(desc_a%glob_to_loc)) then
     info=295
     call psb_errpush(info,name)
     goto 9999
  end if

    !deallocate glob_to_loc field
    deallocate(desc_a%glob_to_loc,stat=info)
    if (info /= 0) then
       info=2052
       call psb_errpush(info,name)
       goto 9999
    end if

  if (.not.associated(desc_a%halo_index)) then
     info=295
     call psb_errpush(info,name)
     goto 9999
  end if

    !deallocate halo_index field
    deallocate(desc_a%halo_index,stat=info)
    if (info /= 0) then
       info=2053
       call psb_errpush(info,name)
       goto 9999
    end if

  if (.not.associated(desc_a%bnd_elem)) then
     info=296
     call psb_errpush(info,name)
     goto 9999
  end if

    !deallocate halo_index field
  deallocate(desc_a%bnd_elem,stat=info)
  if (info /= 0) then
     info=2054
     call psb_errpush(info,name)
     goto 9999
  end if

  if (.not.associated(desc_a%ovrlap_index)) then
     info=295
     call psb_errpush(info,name)
     goto 9999
  end if

  !deallocate ovrlap_index  field
  deallocate(desc_a%ovrlap_index,stat=info)
  if (info /= 0) then
     info=2055
     call psb_errpush(info,name)
     goto 9999
  end if

  !deallocate ovrlap_elem  field
    deallocate(desc_a%ovrlap_elem,stat=info)
  if (info /= 0) then 
    info=2056
    call psb_errpush(info,name)
    goto 9999
  end if

    !deallocate ovrlap_index  field
    deallocate(desc_a%lprm,stat=info)
  if (info /= 0) then 
    info=2057
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_nullify_desc(desc_a)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_ret) then
     return
  else
     call psb_error(icontxt)
  end if
  return

end subroutine psb_dscfree
