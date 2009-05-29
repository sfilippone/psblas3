!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
! File:  psb_zcsprt.f90 
! Subroutine: 
! Arguments:

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
subroutine psb_zcsprt(iout,a,iv,eirs,eics,head,ivr,ivc)
  use psb_spmat_type
  use psb_string_mod
  implicit none 

  integer, intent(in)               :: iout
  type(psb_zspmat_type), intent(in) :: a
  integer, intent(in), optional     :: iv(:)
  integer, intent(in), optional     :: eirs,eics
  character(len=*), optional        :: head
  integer, intent(in), optional     :: ivr(:), ivc(:)

  character(len=80)                 :: frmtv 
  integer  :: irs,ics,i,j, nmx, ni

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
    write(iout,'(a)') '%%MatrixMarket matrix coordinate complex general'
    write(iout,'(a,a)') '% ',head 
    write(iout,'(a)') '%'    
    write(iout,'(a,a)') '% ',psb_toupper(a%fida)
  endif

  nmx = max(a%m,a%k,1)
  ni  = floor(log10(1.0*nmx)) + 1
  
  write(frmtv,'(a,i3.3,a,i3.3,a)') '(2(i',ni,',1x),2(es26.18,1x),2(i',ni,',1x))'

 
  select case(psb_toupper(a%fida)) 

  case ('CSR')

    write(iout,*) a%m,a%k,a%ia2(a%m+1)-1

    if (present(iv)) then 
      do i=1, a%m
        do j=a%ia2(i),a%ia2(i+1)-1
          write(iout,frmtv) iv(irs+i),iv(ics+a%ia1(j)),a%aspk(j)
        enddo
      enddo
    else
      if (present(ivr).and..not.present(ivc)) then 
        do i=1, a%m
          do j=a%ia2(i),a%ia2(i+1)-1
            write(iout,frmtv) ivr(irs+i),(ics+a%ia1(j)),a%aspk(j)
          enddo
        enddo
      else if (present(ivr).and.present(ivc)) then 
        do i=1, a%m
          do j=a%ia2(i),a%ia2(i+1)-1
            write(iout,frmtv) ivr(irs+i),ivc(ics+a%ia1(j)),a%aspk(j)
          enddo
        enddo
      else if (.not.present(ivr).and.present(ivc)) then 
        do i=1, a%m
          do j=a%ia2(i),a%ia2(i+1)-1
            write(iout,frmtv) (irs+i),ivc(ics+a%ia1(j)),a%aspk(j)
          enddo
        enddo
      else if (.not.present(ivr).and..not.present(ivc)) then 
        do i=1, a%m
          do j=a%ia2(i),a%ia2(i+1)-1
            write(iout,frmtv) (irs+i),(ics+a%ia1(j)),a%aspk(j)
          enddo
        enddo
      endif
    endif

  case ('CSC')

    write(iout,*) a%m,a%k,a%ia2(a%k+1)-1

    if (present(iv)) then 
      do i=1, a%k
        do j=a%ia2(i),a%ia2(i+1)-1
          write(iout,frmtv) iv(irs+a%ia1(j)),iv(ics+i),a%aspk(j)
        enddo
      enddo
    else
      if (present(ivr).and..not.present(ivc)) then 
        do i=1, a%k
          do j=a%ia2(i),a%ia2(i+1)-1
            write(iout,frmtv) ivr(irs+a%ia1(j)),(ics+i),a%aspk(j)
          enddo
        enddo
      else if (present(ivr).and.present(ivc)) then 
        do i=1, a%k
          do j=a%ia2(i),a%ia2(i+1)-1
            write(iout,frmtv) ivr(irs+a%ia1(j)),ivc(ics+i),a%aspk(j)
          enddo
        enddo
      else if (.not.present(ivr).and.present(ivc)) then 
        do i=1, a%m
          do j=a%ia2(i),a%ia2(i+1)-1
            write(iout,frmtv) (irs+a%ia1(j)),ivc(ics+i),a%aspk(j)
          enddo
        enddo
      else if (.not.present(ivr).and..not.present(ivc)) then 
        do i=1, a%k
          do j=a%ia2(i),a%ia2(i+1)-1
            write(iout,frmtv) (irs+a%ia1(j)),(ics+i),a%aspk(j)
          enddo
        enddo
      endif
    endif

  case ('COO') 
    if(present(iv)) then 
      write(iout,*) a%m,a%k,a%infoa(psb_nnz_)
      do j=1,a%infoa(psb_nnz_)
        write(iout,frmtv) iv(a%ia1(j)),iv(a%ia2(j)),a%aspk(j)
      enddo
    else      
      if (present(ivr).and..not.present(ivc)) then 
        write(iout,*) a%m,a%k,a%infoa(psb_nnz_)
        do j=1,a%infoa(psb_nnz_)
          write(iout,frmtv) ivr(a%ia1(j)),a%ia2(j),a%aspk(j)
        enddo
      else if (present(ivr).and.present(ivc)) then 
        write(iout,*) a%m,a%k,a%infoa(psb_nnz_)
        do j=1,a%infoa(psb_nnz_)
          write(iout,frmtv) ivr(a%ia1(j)),ivc(a%ia2(j)),a%aspk(j)
        enddo
      else if (.not.present(ivr).and.present(ivc)) then 
        write(iout,*) a%m,a%k,a%infoa(psb_nnz_)
        do j=1,a%infoa(psb_nnz_)
          write(iout,frmtv) a%ia1(j),ivc(a%ia2(j)),a%aspk(j)
        enddo
      else if (.not.present(ivr).and..not.present(ivc)) then 
        write(iout,*) a%m,a%k,a%infoa(psb_nnz_)
        do j=1,a%infoa(psb_nnz_)
          write(iout,frmtv) a%ia1(j),a%ia2(j),a%aspk(j)
        enddo
      endif
    endif
  case default
    write(0,*) 'Feeling lazy today, format not implemented: "',a%fida,'"'
  end select

end subroutine psb_zcsprt
