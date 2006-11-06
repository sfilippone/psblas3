subroutine psb_set_coher(ictxt,isvch)
  integer :: ictxt, isvch
  ! Ensure global coherence for convergence checks.
#ifdef NORMAL
  Call blacs_get(ictxt,16,isvch)
  Call blacs_set(ictxt,16,1)
#endif
#ifdef HAVE_ESSL
  ! Do nothing: ESSL does coherence by default,
  ! and does not handle req=16 
#endif
end subroutine psb_set_coher
subroutine psb_restore_coher(ictxt,isvch)
  integer :: ictxt, isvch
  ! Ensure global coherence for convergence checks.
#ifdef NORMAL
  Call blacs_set(ictxt,16,isvch)
#endif
#ifdef HAVE_ESSL
  ! Do nothing: ESSL does coherence by default,
  ! and does not handle req=16 
#endif
end subroutine psb_restore_coher
subroutine psb_get_mpicomm(ictxt,comm)
  integer :: ictxt, comm
  call blacs_get(ictxt,10,comm)
end subroutine psb_get_mpicomm
subroutine psb_get_rank(rank,ictxt,id)
  integer :: rank,ictxt, id
  integer :: blacs_pnum
  rank =  blacs_pnum(ictxt,id,0)
end subroutine psb_get_rank

#ifdef HAVE_ESSL
!
! Need these, as they are not in the ESSL implementation 
! of the BLACS. 
!
integer function krecvid(contxt,proc_to_comm,myrow)
  integer contxt,proc_to_comm,myrow
  
  krecvid=32766
  
  return
end function krecvid
integer function ksendid(contxt,proc_to_comm,myrow)
  integer contxt,proc_to_comm,myrow

  ksendid=32766

  return
end function ksendid
#endif
