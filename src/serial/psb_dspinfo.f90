! File:  psb_dspinfo.f90 
! Subroutine: 
! Parameters:

!*****************************************************************************
!*                                                                           *
!* Extract info from sparse matrix A. The required info is always a single   *
!* integer.                            Input FIDA might be anything, once    *
!* we get to actually write the code.....                                    *
!*                                                                           *
!*****************************************************************************
subroutine psb_dspinfo(ireq,a,ires,info,iaux)
  use psb_spmat_type
  use psb_const_mod
  use psb_error_mod
  implicit none

  type(psb_dspmat_type), intent(in) :: a
  integer, intent(in)               :: ireq
  integer, intent(out)              :: ires, info
  integer, intent(in), optional     :: iaux

  integer :: i,j,k,ip,jp,nr,irw,nz, err_act
  character(len=20)                 :: name, ch_err

  name='psb_dspinfo'
  info  = 0
  call psb_erractionsave(err_act)


  if (ireq == psb_nztotreq_) then 
     if (a%fida == 'CSR') then 
        nr   = a%m
        ires = a%ia2(nr+1)-1
     else if ((a%fida == 'COO').or.(a%fida == 'COI')) then 
        ires = a%infoa(psb_nnz_)
     else if (a%fida == 'JAD') then 
        ires=-1
        info=135
        ch_err=a%fida(1:3)
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     else
        ires=-1
        info=136
        ch_err=a%fida(1:3)
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

  else if (ireq == psb_nzrowreq_) then 
     if (.not.present(iaux)) then 
        write(0,*) 'Need IAUX when ireq=nzrowreq'
        ires=-1
        return
     endif
     irw = iaux
     if (a%fida == 'CSR') then 
        ires = a%ia2(irw+1)-a%ia2(irw)
     else if ((a%fida == 'COO').or.(a%fida == 'COI')) then 

        if (a%infoa(psb_srtd_) == psb_isrtdcoo_) then 
!!$      write(0,*) 'Gtrow_: srtd coo',irw
           ! In this case we can do a binary search. 
           nz = a%infoa(psb_nnz_)
           call ibsrch(ip,irw,nz,a%ia1)
           jp = ip
           ! expand [ip,jp] to contain all row entries.
           do 
              if (ip < 2) exit
              if (a%ia1(ip-1) == irw) then  
                 ip = ip -1 
              else 
                 exit
              end if
           end do

           do
              if (jp > nz) exit
              if (a%ia1(jp) == irw) then
                 jp =jp + 1
              else
                 exit
              endif
           end do
           ires = jp-ip
        else
           ires = count(a%ia1(1:a%infoa(psb_nnz_))==irw)
        endif
!!$      ires = 0
!!$      do i=1, a%infoa(psb_nnz_) 
!!$        if (a%ia1(i) == irw) ires = ires + 1
!!$      enddo
     else if (a%fida == 'JAD') then 
        ires=-1
        info=135
        ch_err=a%fida(1:3)
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     else
        ires=-1
        info=136
        ch_err=a%fida(1:3)
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

  else  if (ireq == psb_nzsizereq_) then 
     if (a%fida == 'CSR') then 
        ires = size(a%aspk)
     else if ((a%fida == 'COO').or.(a%fida == 'COI')) then 
        ires = size(a%aspk)
     else if (a%fida == 'JAD') then 
        ires = a%infoa(psb_nnz_)
     else
        ires=-1
        info=136
        ch_err=a%fida(1:3)
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

  else 
     write(0,*) 'Unknown request into SPINFO'
     ires=-1
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error()
     return
  end if
  return

end subroutine psb_dspinfo
