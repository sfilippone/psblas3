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
  

subroutine psb_z_cp_hdia_from_coo(a,b,info) 
  
  use psb_base_mod
  use psb_z_hdia_mat_mod, psb_protect_name => psb_z_cp_hdia_from_coo
  implicit none 

  class(psb_z_hdia_sparse_mat), intent(inout) :: a
  class(psb_z_coo_sparse_mat), intent(in)    :: b
  integer(psb_ipk_), intent(out)             :: info

  !locals
  type(psb_z_coo_sparse_mat) :: tmp

  info = psb_success_
  if (b%is_dev()) call b%sync()
  if (b%is_by_rows()) then 
    call inner_cp_hdia_from_coo(a,b,info)
    if (info /= psb_success_)  goto 9999
  else
    call b%cp_to_coo(tmp,info)
    if (info /= psb_success_)  goto 9999
    if (.not.tmp%is_by_rows()) call tmp%fix(info)
    if (info /= psb_success_)  goto 9999
    call inner_cp_hdia_from_coo(a,tmp,info)
    if (info /= psb_success_)  goto 9999
    call tmp%free()
  end if
  call a%set_host()

  return

9999 continue

  info = psb_err_alloc_dealloc_
  return

contains

  subroutine inner_cp_hdia_from_coo(a,tmp,info)
    use psb_base_mod
    use psi_ext_util_mod
  
    implicit none 
    class(psb_z_hdia_sparse_mat), intent(inout) :: a
    class(psb_z_coo_sparse_mat), intent(in)     :: tmp
    integer(psb_ipk_), intent(out)              :: info

    !locals
    integer(psb_ipk_)              :: ndiag,mi,mj,dm,bi,w
    integer(psb_ipk_),allocatable  :: d(:), offset(:), irsz(:) 
    integer(psb_ipk_)              :: k,i,j,nc,nr,nza, nzd,nd,hacksize,nhacks,iszd,&
         &  ib, ir, kfirst, klast1, hackfirst, hacknext, nzout
    integer(psb_ipk_)              :: debug_level, debug_unit
    character(len=20)              :: name
    logical, parameter :: debug=.false. 
    nr  = tmp%get_nrows()
    nc  = tmp%get_ncols()
    nza = tmp%get_nzeros()
    ! If it is sorted then we can lessen memory impact 
    a%psb_z_base_sparse_mat = tmp%psb_z_base_sparse_mat

    hacksize = a%hacksize
    a%nhacks = (nr+hacksize-1)/hacksize
    nhacks   = a%nhacks

    ndiag = nr+nc-1
    if (info == psb_success_) call psb_realloc(nr,irsz,info)
    if (info == psb_success_) call psb_realloc(ndiag,d,info)
    if (info == psb_success_) call psb_realloc(ndiag,offset,info)
    if (info == psb_success_) call psb_realloc(nhacks+1,a%hackoffsets,info)
    if (info /= psb_success_) return

    irsz = 0
    do k=1,nza
      ir = tmp%ia(k)
      irsz(ir) = irsz(ir)+1
    end do

    a%nzeros = 0
    d    = 0
    iszd = 0
    a%hackOffsets(1)=0
    klast1 = 1 
    do k=1, nhacks
      i = (k-1)*hacksize + 1
      ib = min(hacksize,nr-i+1) 
      kfirst = klast1 
      klast1 = kfirst + sum(irsz(i:i+ib-1))
      ! klast1 points to last element of chunk plus 1
      if (debug) then 
        write(*,*) 'Loop iteration ',k,nhacks,i,ib,nr
        write(*,*) 'RW:',tmp%ia(kfirst),tmp%ia(klast1-1)
        write(*,*) 'CL:',tmp%ja(kfirst),tmp%ja(klast1-1)
      end if
      call psi_dia_offset_from_coo(nr,nc,(klast1-kfirst),&
           & tmp%ia(kfirst:klast1-1), tmp%ja(kfirst:klast1-1),&
           & nd, d, offset, info, initd=.false., cleard=.true.)
      iszd = iszd + nd 
      a%hackOffsets(k+1)=iszd
      if (debug) write(*,*) 'From chunk ',k,i,ib,sum(irsz(i:i+ib-1)),': ',nd, iszd
      if (debug) write(*,*) 'offset ', offset(1:nd)
    end do
    if (debug) then 
      write(*,*) 'Hackcount ',nhacks,' Allocation height ',iszd
      write(*,*) 'Hackoffsets ',a%hackOffsets(:)
    end if
    if (info == psb_success_) call psb_realloc(hacksize*iszd,a%diaOffsets,info)
    if (info == psb_success_) call psb_realloc(hacksize*iszd,a%val,info)
    if (info /= psb_success_) return
    klast1 = 1 
    ! 
    ! Second run: copy elements 
    !
    do k=1, nhacks
      i = (k-1)*hacksize + 1
      ib = min(hacksize,nr-i+1) 
      kfirst = klast1 
      klast1 = kfirst + sum(irsz(i:i+ib-1))
      ! klast1 points to last element of chunk plus 1
      hackfirst = a%hackoffsets(k)
      hacknext  = a%hackoffsets(k+1)
      call psi_dia_offset_from_coo(nr,nc,(klast1-kfirst),&
           & tmp%ia(kfirst:klast1-1), tmp%ja(kfirst:klast1-1),&
           & nd, d, a%diaOffsets(hackfirst+1:hacknext), info, &
           & initd=.false., cleard=.false.)
      if (debug) write(*,*) 'Out from dia_offset: ', a%diaOffsets(hackfirst+1:hacknext)
      call psi_z_xtr_dia_from_coo(nr,nc,(klast1-kfirst),&
           & tmp%ia(kfirst:klast1-1), tmp%ja(kfirst:klast1-1),&
           & tmp%val(kfirst:klast1-1), &
           & d,hacksize,(hacknext-hackfirst),&
           & a%val((hacksize*hackfirst)+1:hacksize*hacknext),info,&
           & initdata=.true.,rdisp=(i-1))
          
      call countnz(nr,nc,(i-1),hacksize,(hacknext-hackfirst),&
           & a%diaOffsets(hackfirst+1:hacknext),nzout)
      a%nzeros = a%nzeros + nzout
      call cleand(nr,(hacknext-hackfirst),d,a%diaOffsets(hackfirst+1:hacknext))
      
    end do
    if (debug) then 
      write(*,*) 'NZEROS: ',a%nzeros, nza
      write(*,*) 'diaoffsets: ',a%diaOffsets(1:iszd)
      write(*,*) 'values: '
      j=0
      do k=1,nhacks
        write(*,*) 'Hack No. ',k
        do i=1,hacksize*(iszd/nhacks)
          j = j + 1
          write(*,*) j, a%val(j)
        end do
      end do
    end if
  end subroutine inner_cp_hdia_from_coo
  
  subroutine countnz(nr,nc,rdisp,nrd,ncd,offsets,nz)
    implicit none 
    integer(psb_ipk_), intent(in)  :: nr,nc,nrd,ncd,rdisp,offsets(:)
    integer(psb_ipk_), intent(out) :: nz
    !
    integer(psb_ipk_) :: i,j,k, ir, jc, m4, ir1, ir2, nrcmdisp, rdisp1    
    nz = 0 
    nrcmdisp = min(nr-rdisp,nc-rdisp) 
    rdisp1   = 1-rdisp
    do j=1, ncd
      if (offsets(j)>=0) then 
        ir1 = 1
        ! ir2 = min(nrd,nr - offsets(j) - rdisp_,nc-offsets(j)-rdisp_)
        ir2 = min(nrd, nrcmdisp - offsets(j))
      else
        ! ir1 = max(1,1-offsets(j)-rdisp_) 
        ir1 = max(1, rdisp1 - offsets(j))
        ir2 = min(nrd, nrcmdisp) 
      end if
      nz = nz + (ir2-ir1+1)
    end do
  end subroutine countnz
  
  subroutine cleand(nr,nd,d,offset)
    implicit none 
    integer(psb_ipk_), intent(in)    :: nr,nd,offset(:)
    integer(psb_ipk_), intent(inout) :: d(:)
    integer(psb_ipk_) :: i,id

    do i=1,nd
      id = offset(i) + nr
      d(id) = 0
    end do
  end subroutine cleand

end subroutine psb_z_cp_hdia_from_coo
