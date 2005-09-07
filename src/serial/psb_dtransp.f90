! File:  psb_dtransp.f90 
! Subroutine: 
! Parameters:

subroutine psb_dtransp(a,b,c,fmt)
  use psb_spmat_type
  use psb_serial_mod, only : ipcoo2csr, ipcsr2coo, fixcoo
  implicit none

  type(psb_dspmat_type)      :: a,b
  integer, optional          :: c
  character(len=*), optional :: fmt

  character(len=5)           :: fmt_
  integer  ::c_, info, nz 
  integer, pointer :: itmp(:)=>null()
  if (present(c)) then 
    c_=c
  else
    c_=1
  endif
  if (present(fmt)) then 
    fmt_ = fmt
  else 
    fmt_='CSR'
  endif
  if (associated(b%aspk)) call psb_spfree(b,info)
  call psb_spclone(a,b,info)
  
  if (b%fida=='CSR') then 
    call psb_ipcsr2coo(b,info)
  else if (b%fida=='COO') then 
    ! do nothing 
  else
    write(0,*) 'Unimplemented case in TRANSP '
  endif
!!$  nz = b%infoa(nnz_)
!!$  write(0,*) 'TRANSP CHECKS:',a%m,a%k,&
!!$       &minval(b%ia1(1:nz)),maxval(b%ia1(1:nz)),&
!!$       &minval(b%ia2(1:nz)),maxval(b%ia2(1:nz))
  itmp  => b%ia1
  b%ia1 => b%ia2
  b%ia2 => itmp

  b%m = a%k 
  b%k = a%m
!!$  write(0,*) 'Calling IPCOO2CSR from transp90 ',b%m,b%k
  if (fmt_=='CSR') then 
    call psb_ipcoo2csr(b,info)
    b%fida='CSR'
  else if (fmt_=='COO') then 
    call psb_fixcoo(b,info)
    b%fida='COO'
  else
    write(0,*) 'Unknown FMT in TRANSP : "',fmt_,'"'
  endif

  return
end subroutine psb_dtransp
