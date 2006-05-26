subroutine psb_barrier(ictxt)
  integer, intent(in) :: ictxt
  
  call blacs_barrier(ictxt,'All')

end subroutine psb_barrier
