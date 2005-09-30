! File:  psb_dspgtrow.f90 
! Subroutine: psb_dspgtrow
!    Gets one or more rows from a sparse matrix. 
! Parameters:

!*****************************************************************************
!*                                                                           *
!* Takes a specified row from matrix A and copies into matrix B (possibly    *
!*  appending to B). Output is always COO. Input might be anything, once     *
!* we get to actually write the code.....                                    *
!*                                                                           *
!*****************************************************************************
subroutine psb_dspgtrow(irw,a,b,info,append,iren,lrw)
  ! Output is always in  COO format  into B, irrespective of 
  ! the input format 
  use psb_spmat_type
  use psb_const_mod
  implicit none

  type(psb_dspmat_type), intent(in)     :: a
  integer, intent(in)                   :: irw
  type(psb_dspmat_type), intent(inout)  :: b
  integer,intent(out)                   :: info
  logical, intent(in), optional         :: append
  integer, intent(in), target, optional :: iren(:)
  integer, intent(in), optional         :: lrw

  logical :: append_ 
  integer, pointer :: iren_(:)
  integer :: i,j,k,ip,jp,nr,idx, nz,iret,nzb, nza, lrw_, irw_, err_act
  character(len=20)                 :: name, ch_err

  name='psb_dspgtrow'
  info  = 0
  call psb_erractionsave(err_act)

  irw_ = irw 
  if (present(lrw)) then
    lrw_ = lrw
  else
    lrw_ = irw
  endif
  if (lrw_ < irw) then
    write(0,*) 'SPGTROW input error: fixing lrw',irw,lrw_
    lrw_ = irw
  end if
  if (present(append)) then
    append_=append
  else
    append_=.false.
  endif
  if (present(iren)) then
    iren_=>iren
  else 
    iren_ => null()
  end if


  if (append_) then 
    nzb = b%infoa(psb_nnz_)
  else
    nzb = 0
    b%m = 0 
    b%k = 0
  endif

  if (a%fida == 'CSR') then 
     call csr_dspgtrow(irw_,a,b,append_,iren_,lrw_)
     
  else if (a%fida == 'COO') then 
     call coo_dspgtrow(irw_,a,b,append_,iren_,lrw_)
     
  else if (a%fida == 'JAD') then 
     info=135
     ch_err=a%fida(1:3)
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  else
     info=136
     ch_err=a%fida(1:3)
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if
  
  
  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error()
     return
  end if
  return
  
contains
  
  subroutine csr_dspgtrow(irw,a,b,append,iren,lrw)

    use psb_spmat_type
    use psb_const_mod
    implicit none
    
    type(psb_dspmat_type), intent(in)     :: a
    integer                               :: irw
    type(psb_dspmat_type), intent(inout)  :: b
    logical, intent(in)                   :: append
    integer, pointer                      :: iren(:)
    integer                               :: lrw

    integer :: idx,i,j ,nr,nz,nzb

    if (a%pl(1) /= 0) then
      write(0,*) 'Fatal error in SPGTROW: do not feed a permuted mat so far!',&
           & a%pl(1)
      idx = -1 
    else
      idx = irw
    endif
!!$    write(0,*) 'csr_gtrow: ',irw,lrw,a%pl(1),idx
    if (idx<0) then 
      write(0,*) ' spgtrow Error : idx no good ',idx
      return
    end if
    nr = lrw - irw + 1 
    nz = a%ia2(idx+nr) - a%ia2(idx)
    if (append) then 
      nzb = b%infoa(psb_nnz_)
    else
      nzb = 0
    endif
    if (min(size(b%ia1),size(b%ia2),size(b%aspk)) < nzb+nz) then 
      call psb_spreall(b,nzb+nz,iret)
    endif
    b%fida='COO'
!!$    write(0,*) 'csr_gtrow: ',out_,b%fida,nzb
    if (associated(iren)) then 
      k=0
      do i=irw,lrw
        do j=a%ia2(i),a%ia2(i+1)-1
          k             = k + 1
          b%aspk(nzb+k) = a%aspk(j)
          b%ia1(nzb+k)  = iren(i)
          b%ia2(nzb+k)  = iren(a%ia1(j))
        end do
      enddo
    else
      k=0
!!$      write(0,*) 'csr_gtrow: ilp',irw,lrw
      do i=irw,lrw
!!$        write(0,*) 'csr_gtrow: jlp',a%ia2(i),a%ia2(i+1)-1
        do j=a%ia2(i),a%ia2(i+1)-1
          k             = k + 1
          b%aspk(nzb+k) = a%aspk(j)
          b%ia1(nzb+k)  = i
          b%ia2(nzb+k)  = a%ia1(j)
!!$          write(0,*) 'csr_gtrow: in:',a%aspk(j),i,a%ia1(j)
        end do
      enddo
    end if
    b%infoa(psb_nnz_) = nzb+nz
    if (a%pr(1) /= 0) then
      write(0,*) 'Feeling lazy today, Right Permutation will have to wait'
    endif
    b%m = b%m+lrw-irw+1
    b%k = max(b%k,a%k)

  end subroutine csr_dspgtrow

  subroutine coo_dspgtrow(irw,a,b,append,iren,lrw)

    use psb_spmat_type
    use psb_const_mod
    implicit none

    type(psb_dspmat_type), intent(in)     :: a
    integer                               :: irw
    type(psb_dspmat_type), intent(inout)  :: b
    logical, intent(in)                   :: append
    integer, pointer                      :: iren(:)
    integer                               :: lrw

    nza = a%infoa(psb_nnz_)
    if (a%pl(1) /= 0) then
      write(0,*) 'Fatal error in SPGTROW: do not feed a permuted mat so far!'
      idx = -1 
    else
      idx = irw
    endif
    if (idx<0) then 
      write(0,*) ' spgtrow Error : idx no good ',idx
      return
    end if

    if (a%infoa(psb_srtd_) == psb_isrtdcoo_) then 
!!$      write(0,*) 'Gtrow_: srtd coo',irw
      ! In this case we can do a binary search. 
      do
        call ibsrch(ip,irw,nza,a%ia1)
        if (ip /= -1) exit
        irw = irw + 1
        if (irw > lrw) then
          write(0,*) 'Warning : did not find any rows. Is this an error?'
          exit
        end if
      end do
      
      if (ip /= -1) then 
        ! expand [ip,jp] to contain all row entries.
        do 
          if (ip < 2) exit
          if (a%ia1(ip-1) == irw) then  
            ip = ip -1 
          else 
            exit
          end if
        end do

      end if
      
      do
        call ibsrch(jp,lrw,nza,a%ia1)
        if (jp /= -1) exit
        lrw = lrw - 1
        if (irw > lrw) then
          write(0,*) 'Warning : did not find any rows. Is this an error?'
          exit
        end if
      end do
      
      if (jp /= -1) then 
        ! expand [ip,jp] to contain all row entries.
        do 
          if (jp == nza) exit
          if (a%ia1(jp+1) == lrw) then  
            jp = jp + 1
          else 
            exit
          end if
        end do
      end if
      if ((ip /= -1) .and.(jp /= -1)) then 
        ! Now do the copy.
        nz = jp - ip +1 
        if (size(b%ia1) < nzb+nz) then 
          call psb_spreall(b,nzb+nz,iret)
        endif
        b%fida='COO'        
        if (associated(iren)) then 
          do i=ip,jp
            nzb = nzb + 1
            b%aspk(nzb) = a%aspk(i)
            b%ia1(nzb)  = iren(a%ia1(i))
            b%ia2(nzb)  = iren(a%ia2(i))
          enddo
        else
          do i=ip,jp
            nzb = nzb + 1
            b%aspk(nzb) = a%aspk(i)
            b%ia1(nzb)  = a%ia1(i)
            b%ia2(nzb)  = a%ia2(i)
          enddo
        end if
      end if

    else

      nz = (nza*(lrw-irw+1))/max(a%m,1)
      
      if (size(b%ia1) < nzb+nz) then 
        call psb_spreall(b,nzb+nz,iret)
      endif
      
      if (associated(iren)) then 
        k = 0 
        do i=1,a%infoa(psb_nnz_)
          if ((a%ia1(i)>=irw).and.(a%ia1(i)<=lrw)) then 
            k = k + 1 
            if (k > nz) then
              nz = k 
              call psb_spreall(b,nzb+nz,iret)
            end if
            b%aspk(nzb+k) = a%aspk(i)
            b%ia1(nzb+k)  = iren(a%ia1(i))
            b%ia2(nzb+k)  = iren(a%ia2(i))
          endif
        enddo
      else
        k = 0 
        do i=1,a%infoa(psb_nnz_)
          if ((a%ia1(i)>=irw).and.(a%ia1(i)<=lrw)) then 
            k = k + 1 
            if (k > nz) then
              nz = k 
              call psb_spreall(b,nzb+nz,iret)
            end if
            b%aspk(nzb+k) = a%aspk(i)
            b%ia1(nzb+k)  = (a%ia1(i))
            b%ia2(nzb+k)  = (a%ia2(i))
          endif
        enddo
        nzb=nzb+k
      end if
    end if

    b%infoa(psb_nnz_) = nzb
    b%m = b%m+lrw-irw+1
    b%k = max(b%k,a%k)
  end subroutine coo_dspgtrow

end subroutine psb_dspgtrow

