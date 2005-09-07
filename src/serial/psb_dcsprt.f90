! File:  psb_dcsprt.f90 
! Subroutine: 
! Parameters:

!*****************************************************************************
!*                                                                           *
!* Print out a matrix.                                                       *
!*  Should really align with the F77 version under the SERIAL dir, which     *
!*  does a nice printout in MatrixMarket format; this would be a quick job.  *
!*                                                                           *
!*  Handles both a shift in the row/col indices and a fuctional transform    *
!*  on the indices.                                                          *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*****************************************************************************
subroutine psb_dcsprt(iout,a,iv,eirs,eics,head,ivr,ivc)
  use psb_spmat_type
  implicit none 

  integer, intent(in)               :: iout
  type(psb_dspmat_type), intent(in) :: a
  integer, intent(in), optional     :: iv(:)
  integer, intent(in), optional     :: eirs,eics
  character(len=*), optional        :: head
  integer, intent(in), optional     :: ivr(:), ivc(:)

  character(len=*), parameter       :: frmtr='(2(i6,1x),e16.8,2(i6,1x))'
  integer  :: irs,ics,i,j

  if (present(eirs)) then 
    irs = eirs
  else
    irs = 0
  endif
  if (present(eics)) then 
    ics = eics
  else
    ics = 0
  endif

  if (present(head)) then 
    write(iout,'(a)') head 
  endif

  if (a%fida=='CSR') then 

    write(iout,*) a%m,a%k,a%ia2(a%m+1)-1
    
    if (present(iv)) then 
      do i=1, a%m
        do j=a%ia2(i),a%ia2(i+1)-1
          write(iout,frmtr) iv(irs+i),iv(ics+a%ia1(j)),a%aspk(j)
        enddo
      enddo
    else
      if (present(ivr).and..not.present(ivc)) then 
        do i=1, a%m
          do j=a%ia2(i),a%ia2(i+1)-1
            write(iout,frmtr) ivr(irs+i),(ics+a%ia1(j)),a%aspk(j)
          enddo
        enddo
      else if (present(ivr).and.present(ivc)) then 
        do i=1, a%m
          do j=a%ia2(i),a%ia2(i+1)-1
            write(iout,frmtr) ivr(irs+i),ivc(ics+a%ia1(j)),a%aspk(j)
          enddo
        enddo
      else if (.not.present(ivr).and.present(ivc)) then 
        do i=1, a%m
          do j=a%ia2(i),a%ia2(i+1)-1
            write(iout,frmtr) (irs+i),ivc(ics+a%ia1(j)),a%aspk(j)
          enddo
        enddo
      else if (.not.present(ivr).and..not.present(ivc)) then 
        do i=1, a%m
          do j=a%ia2(i),a%ia2(i+1)-1
            write(iout,frmtr) (irs+i),(ics+a%ia1(j)),a%aspk(j)
          enddo
        enddo
      endif
    endif

  else  if (a%fida=='COO') then 

    if (present(ivr).and..not.present(ivc)) then 
      write(iout,*) a%m,a%k,a%infoa(nnz_)
      do j=1,a%infoa(nnz_)
        write(iout,frmtr) ivr(a%ia1(j)),a%ia2(j),a%aspk(j)
      enddo
    else if (present(ivr).and.present(ivc)) then 
      write(iout,*) a%m,a%k,a%infoa(nnz_)
      do j=1,a%infoa(nnz_)
        write(iout,frmtr) ivr(a%ia1(j)),ivc(a%ia2(j)),a%aspk(j)
      enddo
    else if (.not.present(ivr).and.present(ivc)) then 
      write(iout,*) a%m,a%k,a%infoa(nnz_)
      do j=1,a%infoa(nnz_)
        write(iout,frmtr) a%ia1(j),ivc(a%ia2(j)),a%aspk(j)
      enddo
    else if (.not.present(ivr).and..not.present(ivc)) then 
      write(iout,*) a%m,a%k,a%infoa(nnz_)
      do j=1,a%infoa(nnz_)
        write(iout,frmtr) a%ia1(j),a%ia2(j),a%aspk(j)
      enddo
    endif
  else
    write(0,*) 'Feeling lazy today, format not implemented: "',a%fida,'"'
  endif
end subroutine psb_dcsprt
