subroutine psb_init(ictxt,np)
  use psb_const_mod
  use psb_blacs_mod
  use psb_error_mod
  integer, intent(out) :: ictxt
  integer, intent(in), optional :: np

  integer :: np_, npavail, iam, info
  character(len=20), parameter :: name='psb_init'

  call blacs_pinfo(iam, npavail)
  call blacs_get(izero, izero, ictxt)
  
  if (present(np)) then 
    np_ = max(1,min(np,npavail))
  else
    np_ = npavail
  endif

  call blacs_gridinit(ictxt, 'R', np_, ione)
 
  if (present(np)) then 
    if (np_ < np) then 
      info = 2011
      call psb_errpush(info,name)
      call psb_error(ictxt)
    endif
  endif

  
end subroutine psb_init
