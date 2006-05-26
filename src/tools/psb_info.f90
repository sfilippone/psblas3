subroutine psb_info(ictxt,iam,np)
  integer, intent(in)  :: ictxt
  integer, intent(out) :: iam, np
  integer              :: nprow, npcol, myprow, mypcol

  call blacs_gridinfo(ictxt, nprow, npcol, myprow, mypcol)

  iam = myprow
  np  = nprow

end subroutine psb_info
