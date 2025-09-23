subroutine psb_saplusat(ain,aout,info)
  use psb_s_mat_mod

  implicit none
  type(psb_sspmat_type), intent(inout) :: ain
  type(psb_sspmat_type), intent(out)   :: aout
  integer(psb_ipk_) :: info

  type(psb_s_coo_sparse_mat) :: acoo1, acoo2
  integer(psb_ipk_) :: nr, nc, nz1, nz2
  integer(psb_ipk_) :: err_act
  character(len=20) :: name, ch_err

  name='psb_saplusat'
  info  = psb_success_
  call psb_erractionsave(err_act)


  
  nr = ain%get_nrows()
  nc = ain%get_ncols()

  if (nr /= nc) then
    info=psb_err_internal_error_
    call psb_errpush(info,name)
    goto 9999
  end if

  call ain%cp_to(acoo1)
  call acoo1%cp_to_coo(acoo2,info)
  nz1 = acoo1%get_nzeros()
  nz2 = acoo2%get_nzeros()
  call acoo1%reallocate(nz1+nz2)
  acoo1%ia(nz1+1:nz1+nz2) = acoo2%ja(1:nz2)
  acoo1%ja(nz1+1:nz1+nz2) = acoo2%ia(1:nz2)
  acoo1%val(nz1+1:nz1+nz2) = acoo2%val(1:nz2)
  call acoo1%set_nrows(nr)
  call acoo1%set_ncols(nr)
  call acoo1%set_nzeros(nz1+nz2)
  call aout%cp_from(acoo1)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return


end subroutine psb_saplusat
