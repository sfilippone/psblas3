subroutine psb_set_coher(icontxt,isvch)
  integer :: icontxt, isvch
  ! Ensure global coherence for convergence checks.
  Call blacs_get(icontxt,16,isvch)
  Call blacs_set(icontxt,16,1)
end subroutine psb_set_coher
subroutine psb_restore_coher(icontxt,isvch)
  integer :: icontxt, isvch
  ! Ensure global coherence for convergence checks.
  Call blacs_set(icontxt,16,isvch)
end subroutine psb_restore_coher
