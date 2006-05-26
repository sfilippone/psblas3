subroutine psb_set_coher(ictxt,isvch)
  integer :: ictxt, isvch
  ! Ensure global coherence for convergence checks.
  ! Do nothing: ESSL does coherence by default,
  ! and does not handle req=16 
!!$  Call blacs_get(ictxt,16,isvch)
!!$  Call blacs_set(ictxt,16,1)
end subroutine psb_set_coher
subroutine psb_restore_coher(ictxt,isvch)
  integer :: ictxt, isvch
  ! Ensure global coherence for convergence checks.
!!$  Call blacs_set(ictxt,16,isvch)
end subroutine psb_restore_coher
