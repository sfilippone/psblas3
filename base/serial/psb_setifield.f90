subroutine psb_setifield(val,field,info,isize,ierr)
  integer :: val,field,isize,ierr
  integer :: info(*) 

  ierr = 0

  if ((field < 1).or.(field > isize)) then
    ierr = -1
    return
  endif

  info(field) = val
  return
end subroutine psb_setifield
