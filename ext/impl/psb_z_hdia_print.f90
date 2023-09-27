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
subroutine psb_z_hdia_print(iout,a,iv,head,ivr,ivc)
  
  use psb_base_mod
  use psb_z_hdia_mat_mod, psb_protect_name => psb_z_hdia_print
  use psi_ext_util_mod
  implicit none 

  integer(psb_ipk_), intent(in)           :: iout
  class(psb_z_hdia_sparse_mat), intent(in) :: a   
  integer(psb_lpk_), intent(in), optional :: iv(:)
  character(len=*), optional              :: head
  integer(psb_lpk_), intent(in), optional :: ivr(:), ivc(:)

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='hdia_print'
  logical, parameter :: debug=.false.

  class(psb_z_coo_sparse_mat),allocatable :: acoo

  character(len=80)  :: frmt 
  integer(psb_ipk_)  :: irs,ics,i,j, nmx, ni, nr, nc, nz
  integer(psb_ipk_)  :: nhacks, hacksize,maxnzhack, k, ncd,ib, nzhack, info,&
       & hackfirst, hacknext
  integer(psb_ipk_), allocatable :: ia(:), ja(:)
  complex(psb_dpk_), allocatable    :: val(:) 
  
  
  write(iout,'(a)') '%%MatrixMarket matrix coordinate complex general'
  if (present(head)) write(iout,'(a,a)') '% ',head 
  write(iout,'(a)') '%'    
  write(iout,'(a,a)') '% HDIA'

  if (a%is_dev()) call a%sync()
  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nz  = a%get_nzeros()
  frmt = psb_z_get_print_frmt(nr,nc,nz,iv,ivr,ivc)


  nhacks   = a%nhacks
  hacksize = a%hacksize
  maxnzhack = 0
  do k=1, nhacks
    maxnzhack = max(maxnzhack,(a%hackoffsets(k+1)-a%hackoffsets(k)))
  end do
  maxnzhack = hacksize*maxnzhack
  allocate(ia(maxnzhack),ja(maxnzhack),val(maxnzhack),stat=info)
  if (info /= 0) return 

  write(iout,*) nr, nc, nz 
  do k=1, nhacks
    i = (k-1)*hacksize + 1
    ib = min(hacksize,nr-i+1) 
    hackfirst = a%hackoffsets(k)
    hacknext  = a%hackoffsets(k+1)
    ncd = hacknext-hackfirst
    
    call psi_z_xtr_coo_from_dia(nr,nc,&
           & ia, ja, val, nzhack,&
           & hacksize,ncd,&
           & a%val((hacksize*hackfirst)+1:hacksize*hacknext),&
           & a%diaOffsets(hackfirst+1:hacknext),info,rdisp=(i-1))
    !nzhack = sum(ib - abs(a%diaOffsets(hackfirst+1:hacknext)))
    
    if(present(iv)) then 
      do j=1,nzhack
        write(iout,frmt) iv(ia(j)),iv(ja(j)),val(j)
      enddo
    else      
      if (present(ivr).and..not.present(ivc)) then 
        do j=1,nzhack
          write(iout,frmt) ivr(ia(j)),ja(j),val(j)
        enddo
      else if (present(ivr).and.present(ivc)) then 
        do j=1,nzhack
          write(iout,frmt) ivr(ia(j)),ivc(ja(j)),val(j)
        enddo
      else if (.not.present(ivr).and.present(ivc)) then 
        do j=1,nzhack
          write(iout,frmt) ia(j),ivc(ja(j)),val(j)
        enddo
      else if (.not.present(ivr).and..not.present(ivc)) then 
        do j=1,nzhack
          write(iout,frmt) ia(j),ja(j),val(j)
        enddo
      endif
    end if
      
  end do
  
end subroutine psb_z_hdia_print
