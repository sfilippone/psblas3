! File: psb_dfree.f90
!
! Subroutine: psb_dfree
!    frees a dense matrix structure
! 
! Parameters: 
!    x        - real, pointer, dimension(:,:).    The dense matrix to be freed.
!    desc_a   - type(<psb_desc_type>).            The communication descriptor.
!    info     - integer.                          Eventually returns an error code
subroutine psb_dfree(x, desc_a, info)
  !...free dense matrix structure...
  use psb_const_mod
  use psb_descriptor_type
  use psb_error_mod
  implicit none

  !....parameters...
  real(kind(1.d0)),pointer    :: x(:,:)
  type(psb_desc_type), intent(in) :: desc_a
  integer                     :: info

  !...locals....
  integer             :: int_err(5)
  integer             :: icontxt,nprow,npcol,me,mypcol,err, err_act
  integer,parameter   :: ione=1
  character(len=20)   :: name


  info=0
  call psb_erractionsave(err_act)
  name='psb_dfree'

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

  if (.not.associated(desc_a%matrix_data)) then
     info=295
     call psb_errpush(info,name)
     goto 9999
  end if

  if (.not.associated(x)) then
     info=295
     call psb_errpush(info,name)
     goto 9999
  end if

  !deallocate x
  deallocate(x,stat=info)
  if (info.ne.no_err) then
     info=4000
     call psb_errpush(info,name)
     goto 9999
  else
     nullify(x)
  endif
  

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

end subroutine psb_dfree



! Subroutine: psb_dfreev
!    frees a dense matrix structure
! 
! Parameters: 
!    x        - real, pointer, dimension(:).    The dense matrix to be freed.
!    desc_a   - type(<psb_desc_type>).          The communication descriptor.
!    info     - integer.                        Eventually returns an error code
subroutine psb_dfreev(x, desc_a, info)
  !...free dense matrix structure...
  use psb_const_mod
  use psb_descriptor_type
  use psb_error_mod

  implicit none
  !....parameters...
  real(kind(1.d0)),pointer    :: x(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer                     :: info

  !...locals....
  integer             :: int_err(5)
  integer             :: icontxt,nprow,npcol,me,mypcol,err, err_act
  integer,parameter   :: ione=1
  character(len=20)   :: name


  info=0
  call psb_erractionsave(err_act)
  name='psb_dfreev'

  icontxt=desc_a%matrix_data(psb_ctxt_)

  if (.not.associated(desc_a%matrix_data)) then
     info=295
     call psb_errpush(info,name)
     goto 9999
  end if

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
     info=295
     call psb_errpush(info,name)
     goto 9999
  end if

  !deallocate x
  deallocate(x,stat=info)
  if (info.ne.no_err) then
     info=4000
     call psb_errpush(info,name)
  else
     nullify(x)
  endif
  
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

end subroutine psb_dfreev
