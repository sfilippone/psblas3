subroutine psb_exit(ictxt)
  integer, intent(in) :: ictxt
  
  integer  :: nprow, npcol, myprow, mypcol

  call blacs_gridinfo(ictxt, nprow, npcol, myprow, mypcol)
  if ((myprow >=0).and.(mypcol>=0)) then
    call blacs_gridexit(ictxt)
  end if
  call blacs_exit(0)

end subroutine psb_exit
