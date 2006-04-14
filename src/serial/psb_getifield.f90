subroutine psb_getifield(val,field,info,isize,ierr)
  integer :: val,field,isize,ierr
  integer :: info(*) 

  ierr = 0
  val = -1 
  if ((field < 1).or.(field > isize)) then
    ierr = -1
    return
  endif

  val = info(field) 
  return
end subroutine psb_getifield
