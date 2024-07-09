!                Parallel Sparse BLAS   GPU plugin 
!      (C) Copyright 2013
!  
!                         Salvatore Filippone
!                         Alessandro Fanfarillo
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
  

subroutine psb_z_dia_print(iout,a,iv,head,ivr,ivc)
  
  use psb_base_mod
  use psb_z_dia_mat_mod, psb_protect_name => psb_z_dia_print
  implicit none 

  integer(psb_ipk_), intent(in)           :: iout
  class(psb_z_dia_sparse_mat), intent(in) :: a   
  integer(psb_lpk_), intent(in), optional :: iv(:)
  character(len=*), optional              :: head
  integer(psb_lpk_), intent(in), optional :: ivr(:), ivc(:)

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='z_dia_print'
  logical, parameter :: debug=.false.

  class(psb_z_coo_sparse_mat),allocatable :: acoo

  character(len=80)  :: frmt 
  integer(psb_ipk_)  :: irs,ics,i,j, nmx, ni, nr, nc, nz, jc, ir1, ir2

  write(iout,'(a)') '%%MatrixMarket matrix coordinate complex general'
  if (present(head)) write(iout,'(a,a)') '% ',head 
  write(iout,'(a)') '%'    
  write(iout,'(a,a)') '% COO'

  if (a%is_dev()) call a%sync()

  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nz  = a%get_nzeros()
  frmt = psb_z_get_print_frmt(nr,nc,nz,iv,ivr,ivc)
  write(iout,*) nr, nc, nz 

  nc=size(a%data,2)

      
    
  if(present(iv)) then 
    do j=1,nc
      jc = a%offset(j)
      if (jc > 0) then 
        ir1 = 1
        ir2 = nr - jc
      else
        ir1 = 1 - jc
        ir2 = nr
      end if
      do i=ir1, ir2
        write(iout,frmt) iv(i),iv(i+jc),a%data(i,j)
      enddo
    enddo

  else  if (present(ivr).and..not.present(ivc)) then 
    do j=1,nc
      jc = a%offset(j)
      if (jc > 0) then 
        ir1 = 1
        ir2 = nr - jc
      else
        ir1 = 1 - jc
        ir2 = nr
      end if
      do i=ir1, ir2
        write(iout,frmt) ivr(i),(i+jc),a%data(i,j)
      enddo
    enddo

  else if (present(ivr).and.present(ivc)) then 
    do j=1,nc
      jc = a%offset(j)
      if (jc > 0) then 
        ir1 = 1
        ir2 = nr - jc
      else
        ir1 = 1 - jc
        ir2 = nr
      end if
      do i=ir1, ir2
        write(iout,frmt) ivr(i),ivc(i+jc),a%data(i,j)
      enddo
    enddo

  else if (.not.present(ivr).and.present(ivc)) then 
    do j=1,nc
      jc = a%offset(j)
      if (jc > 0) then 
        ir1 = 1
        ir2 = nr - jc
      else
        ir1 = 1 - jc
        ir2 = nr
      end if
      do i=ir1, ir2
        write(iout,frmt) (i),ivc(i+jc),a%data(i,j)
      enddo
    enddo

  else if (.not.present(ivr).and..not.present(ivc)) then 
    do j=1,nc
      jc = a%offset(j)
      if (jc > 0) then 
        ir1 = 1
        ir2 = nr - jc
      else
        ir1 = 1 - jc
        ir2 = nr
      end if
      do i=ir1, ir2
        write(iout,frmt) (i),(i+jc),a%data(i,j)
      enddo
    enddo
    
  endif
  
end subroutine psb_z_dia_print
