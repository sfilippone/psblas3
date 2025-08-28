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
  

subroutine psb_c_ell_print(iout,a,iv,head,ivr,ivc)
  
  use psb_base_mod
  use psb_c_ell_mat_mod, psb_protect_name => psb_c_ell_print
  implicit none 

  integer(psb_ipk_), intent(in)           :: iout
  class(psb_c_ell_sparse_mat), intent(in) :: a   
  integer(psb_lpk_), intent(in), optional :: iv(:)
  character(len=*), optional              :: head
  integer(psb_lpk_), intent(in), optional :: ivr(:), ivc(:)

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='c_ell_print'
  logical, parameter :: debug=.false.

  character(len=80)  :: frmt 
  integer(psb_ipk_)  :: irs,ics,i,j, nmx, ni, nr, nc, nz

  
  write(iout,'(a)') '%%MatrixMarket matrix coordinate complex general'
  if (present(head)) write(iout,'(a,a)') '% ',head 
  write(iout,'(a)') '%'    
  write(iout,'(a,a)') '% ELL'

  if (a%is_dev()) call a%sync()
  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nz  = a%get_nzeros()
  frmt = psb_c_get_print_frmt(nr,nc,nz,iv,ivr,ivc)

  write(iout,*) nr, nc, nz 
  if(present(iv)) then 
    do i=1, nr
      do j=1,a%irn(i)
        write(iout,frmt) iv(i),iv(a%ja(i,j)),a%val(i,j)
      end do
    enddo
  else      
    if (present(ivr).and..not.present(ivc)) then 
      do i=1, nr
        do j=1,a%irn(i)
          write(iout,frmt) ivr(i),(a%ja(i,j)),a%val(i,j)
        end do
      enddo
    else if (present(ivr).and.present(ivc)) then 
      do i=1, nr
        do j=1,a%irn(i)
          write(iout,frmt) ivr(i),ivc(a%ja(i,j)),a%val(i,j)
        end do
      enddo
    else if (.not.present(ivr).and.present(ivc)) then 
      do i=1, nr
        do j=1,a%irn(i)
          write(iout,frmt) (i),ivc(a%ja(i,j)),a%val(i,j)
        end do
      enddo
    else if (.not.present(ivr).and..not.present(ivc)) then 
      do i=1, nr
        do j=1,a%irn(i)
          write(iout,frmt) (i),(a%ja(i,j)),a%val(i,j)
        end do
      enddo
    endif
  endif

end subroutine psb_c_ell_print
