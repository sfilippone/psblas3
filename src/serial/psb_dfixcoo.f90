! File:  psb_dfixcoo.f90 
! Subroutine: 
! Parameters:

Subroutine psb_dfixcoo(A,INFO,idir)
  use psb_spmat_type
  use psb_const_mod
  implicit none

  !....Parameters...
  Type(psb_dspmat_type), intent(inout) :: A
  Integer, intent(out)                 :: info
  integer, intent(in), optional :: idir

  integer, allocatable :: iaux(:)
  !locals
  Integer              :: nza, nzl,iret,idir_
  integer              :: i,j, irw, icl
  logical, parameter   :: debug=.false.

  info  = 0
  if(debug) write(0,*)'fixcoo: ',size(a%ia1),size(a%ia2)
  if (a%fida /= 'COO') then 
    write(0,*) 'Fixcoo Invalid input ',a%fida
    info = -1
    return
  end if
  if (present(idir)) then 
    idir_ = idir
  else
    idir_ = 0
  endif

  nza = a%infoa(psb_nnz_)
  if (nza < 2) return

  allocate(iaux(nza+2),stat=info) 
  if (info /= 0) return

  select case(idir_) 

  case(0) !  Row major order

    call mrgsrt(nza,a%ia1,iaux,iret)
    if (iret.eq.0) call reordvn(nza,a%aspk,a%ia1,a%ia2,iaux)
    i    = 1
    j    = i
    do while (i.le.nza)
      do while ((a%ia1(j).eq.a%ia1(i)))
        j = j+1
        if (j > nza) exit
      enddo
      nzl = j - i
      call mrgsrt(nzl,a%ia2(i),iaux,iret)
      if (iret.eq.0) &
           & call reordvn(nzl,a%aspk(i),a%ia1(i),a%ia2(i),iaux)
      i = j
    enddo

    i = 1
    irw = a%ia1(i)
    icl = a%ia2(i)
    j = 1
    do 
      j = j + 1
      if (j > nza) exit
      if ((a%ia1(j) == irw).and.(a%ia2(j) == icl)) then 
        a%aspk(i) = a%aspk(i) + a%aspk(j)
      else
        i = i+1
        a%aspk(i) = a%aspk(j)
        a%ia1(i) = a%ia1(j)
        a%ia2(i) = a%ia2(j)
        irw = a%ia1(i) 
        icl = a%ia2(i) 
      endif
    enddo
    a%infoa(psb_nnz_) = i    
    a%infoa(psb_srtd_) = psb_isrtdcoo_

    if(debug) write(0,*)'FIXCOO: end second loop'

  case(1) !  Col major order

    call mrgsrt(nza,a%ia2,iaux,iret)
    if (iret.eq.0) call reordvn(nza,a%aspk,a%ia1,a%ia2,iaux)
    i    = 1
    j    = i
    do while (i.le.nza)
      do while ((a%ia2(j).eq.a%ia2(i)))
        j = j+1
        if (j > nza) exit
      enddo
      nzl = j - i
      call mrgsrt(nzl,a%ia1(i:i+nzl-1),iaux,iret)
      if (iret.eq.0) &
           & call reordvn(nzl,a%aspk(i:i+nzl-1),a%ia1(i:i+nzl-1),a%ia2(i:i+nzl-1),iaux)
      i = j
    enddo

    i = 1
    irw = a%ia1(i)
    icl = a%ia2(i)
    j = 1
    do 
      j = j + 1
      if (j > nza) exit
      if ((a%ia1(j) == irw).and.(a%ia2(j) == icl)) then 
        a%aspk(i) = a%aspk(i) + a%aspk(j)
      else
        i = i+1
        a%aspk(i) = a%aspk(j)
        a%ia1(i) = a%ia1(j)
        a%ia2(i) = a%ia2(j)
        irw = a%ia1(i) 
        icl = a%ia2(i) 
      endif
    enddo
    a%infoa(psb_nnz_) = i    
    a%infoa(psb_srtd_) = psb_isrtdcoo_

    if(debug) write(0,*)'FIXCOO: end second loop'
  case default
    write(0,*) 'Fixcoo: unknown direction ',idir_
  end select

  deallocate(iaux)
  return
end Subroutine psb_dfixcoo
