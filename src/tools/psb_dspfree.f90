! File: psb_dspfree.f90
!
! Subroutine: psb_dspfree
!    Frees a sparse matrix structure.
! 
! Parameters: 
!    a        - type(<psb_dspmat_type>).          The sparse matrix to be freed.      
!    desc_a   - type(<psb_desc_type>).            The communication descriptor.
!    info     - integer.                          Eventually returns an error code.
!
subroutine psb_dspfree(a, desc_a,info)
  !...free sparse matrix structure...
  use psb_descriptor_type
  use psb_spmat_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  implicit none

  !....parameters...
  type(psb_desc_type), intent(in) :: desc_a
  type(psb_dspmat_type), intent(inout)       ::a
  integer, intent(out)        :: info
  !...locals....
  integer             :: int_err(5)
  integer             :: temp(1)
  real(kind(1.d0))    :: real_err(5)
  integer             :: icontxt,nprow,npcol,me,mypcol,err, err_act
  integer,parameter   :: ione=1
  character(len=20)   :: name, ch_err

  info=0
  name = 'psb_dspfree'
  call psb_erractionsave(err_act)

  if (.not.associated(desc_a%matrix_data)) then 
     info=295
     call psb_errpush(info,name)
     return
  else
     icontxt=desc_a%matrix_data(psb_ctxt_)
  end if

  !...deallocate a....

  if ((info.eq.0).and.(.not.associated(a%pr))) info=2951
  if (info.eq.0) then
     !deallocate pr field
     deallocate(a%pr,stat=info)
     if (info.ne.0) info=2045
  end if
  if ((info.eq.0).and.(.not.associated(a%pl))) info=2952
  !deallocate pl  field
  if (info.eq.0) then 
     deallocate(a%pl,stat=info)
     if (info.ne.0) info=2046
  end if
  if ((info.eq.0).and.(.not.associated(a%ia2))) info=2953
  if (info.eq.0) then
     !deallocate ia2 field
     deallocate(a%ia2,stat=info)
     if (info.ne.0) info=2047
  end if
  if ((info.eq.0).and.(.not.associated(a%ia1))) info=2954
  if (info.eq.0) then
     !deallocate ia1  field
     deallocate(a%ia1,stat=info)
     if (info.ne.0) info=2048
  endif
  if ((info.eq.0).and.(.not.associated(a%aspk))) info=2955
  if (info.eq.0) then
     !deallocate aspk field
     deallocate(a%aspk,stat=info)
     if (info.ne.0) info=2049
  endif
  if (info.eq.0) call psb_nullify_sp(a)

  if(info.ne.0) then
     call psb_errpush(info,name)
     goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_dspfree



subroutine psb_dspfrees(a, info)
  !...free sparse matrix structure...
  use psb_descriptor_type
  use psb_spmat_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  implicit none

  !....parameters...
  type(psb_dspmat_type), intent(inout)       ::a
  integer, intent(out)        :: info
  !...locals....
  integer             :: int_err(5)
  integer             :: temp(1)
  real(kind(1.d0))    :: real_err(5)
  integer             :: icontxt,nprow,npcol,me,mypcol,err, err_act
  integer,parameter   :: ione=1
  character(len=20)   :: name, ch_err

  info=0
  name = 'psb_dspfrees'
  call psb_erractionsave(err_act)

  !...deallocate a....

  if ((info.eq.0).and.(.not.associated(a%pr))) info=2951
  if (info.eq.0) then
     !deallocate pr field
     deallocate(a%pr,stat=info)
     if (info.ne.0) info=2045
  end if
  if ((info.eq.0).and.(.not.associated(a%pl))) info=2952
  !deallocate pl  field
  if (info.eq.0) then 
     deallocate(a%pl,stat=info)
     if (info.ne.0) info=2046
  end if
  if ((info.eq.0).and.(.not.associated(a%ia2))) info=2953
  if (info.eq.0) then
     !deallocate ia2 field
     deallocate(a%ia2,stat=info)
     if (info.ne.0) info=2047
  end if
  if ((info.eq.0).and.(.not.associated(a%ia1))) info=2954
  if (info.eq.0) then
     !deallocate ia1  field
     deallocate(a%ia1,stat=info)
     if (info.ne.0) info=2048
  endif
  if ((info.eq.0).and.(.not.associated(a%aspk))) info=2955
  if (info.eq.0) then
     !deallocate aspk field
     deallocate(a%aspk,stat=info)
     if (info.ne.0) info=2049
  endif
  if (info.eq.0) call psb_nullify_sp(a)

  if(info.ne.0) then
     call psb_errpush(info,name)
     goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error()
     return
  end if
  return

end subroutine psb_dspfrees
