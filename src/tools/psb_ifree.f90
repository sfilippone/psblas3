! File: psb_ifree.f90
!
! Subroutine: psb_ifree
!    frees a dense integer matrix structure
! 
! Parameters: 
!    x        - integer, pointer, dimension(:,:).    The dense matrix to be freed.
!    desc_a   - type(<psb_desc_type>).               The communication descriptor.
!    info     - integer.                             Eventually returns an error code
subroutine psb_ifree(x, desc_a, info)
  !...free dense matrix structure...
  use psb_const_mod
  use psb_descriptor_type
  use psb_error_mod
  implicit none

  !....parameters...
  integer, pointer                :: x(:,:)
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(out)            :: info
  
  !...locals....
  integer             :: int_err(5)
  integer             :: temp(1)
  real(kind(1.d0))    :: real_err(5)
  integer             :: icontxt,nprow,npcol,me,mypcol,err_act
  integer,parameter   :: ione=1
  character(len=20)   :: name, ch_err

  info=0
  call psb_erractionsave(err_act)
  name = 'psb_ifree'
  
  if (.not.associated(desc_a%matrix_data)) then
     info=295
     call psb_errpush(info,name)
     return
  end if

  icontxt=desc_a%matrix_data(psb_ctxt_)
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

  if (.not.associated(x)) then
     info=290
     call psb_errpush(info,name)
     goto 9999
  end if
  
  !deallocate x
  deallocate(x,stat=info)
  if (info.ne.0) then
     info=2045
     call psb_errpush(info,name)
     goto 9999
  else
     nullify(x)
  endif
  
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_ifree



! Subroutine: psb_ifreev
!    frees a dense integer matrix structure
! 
! Parameters: 
!    x        - integer, pointer, dimension(:).      The dense matrix to be freed.
!    desc_a   - type(<psb_desc_type>).               The communication descriptor.
!    info     - integer.                             Eventually returns an error code
subroutine psb_ifreev(x, desc_a,info)
  !...free dense matrix structure...
  use psb_const_mod
  use psb_descriptor_type
  use psb_error_mod
  implicit none
  !....parameters...
  integer, pointer                :: x(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(out)            :: info 
  !...locals....
  integer             :: int_err(5)
  integer             :: temp(1)
  real(kind(1.d0))    :: real_err(5)
  integer             :: icontxt,nprow,npcol,me,mypcol,err_act
  integer,parameter   :: ione=1
  character(len=20)   :: name, ch_err

  info=0
  call psb_erractionsave(err_act)
  name = 'psb_ifreev'

  
  if (.not.associated(desc_a%matrix_data)) then
     info=295
     call psb_errpush(info,name)
     return
  end if

  icontxt=desc_a%matrix_data(psb_ctxt_)
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

  if (.not.associated(x)) then
     info=290
     call psb_errpush(info,name,int_err)
     goto 9999
  end if
  
  !deallocate x
  deallocate(x,stat=info)
  if (info.ne.0) then 
     info=2045
     call psb_errpush(info,name,int_err)
     goto 9999
  else
     nullify(x)
  endif
  
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_ifreev
