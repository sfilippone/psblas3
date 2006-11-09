!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
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
  use psb_string_mod
  implicit none

  type(psb_dspmat_type), intent(in), target :: a
  integer, intent(in)               :: ireq
  integer, intent(out)              :: ires, info
  integer, intent(in), optional     :: iaux

  integer :: i,j,k,ip,jp,nr,irw,nz, err_act, row, ipx, pia, pja, rb,idx, nc
  integer, pointer :: ia1(:), ia2(:), ia3(:), ja(:)
  character(len=20)                 :: name, ch_err

  name='psb_dspinfo'
  info  = 0
  call psb_erractionsave(err_act)


  if (ireq == psb_nztotreq_) then 
     ! The number of nonzeroes
     if (toupper(a%fida) == 'CSR') then 
        nr   = a%m
        ires = a%ia2(nr+1)-1
     else if ((toupper(a%fida) == 'COO').or.(toupper(a%fida) == 'COI')) then 
        ires = a%infoa(psb_nnz_)
     else if (toupper(a%fida) == 'JAD') then 
        ires = a%infoa(psb_nnz_)
     else if (toupper(a%fida) == 'CSC') then 
        nc   = a%k
        ires = a%ia2(nc+1)-1
      else
        ires=-1
        info=136
        ch_err=a%fida(1:3)
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

  else if (ireq == psb_nzrowreq_) then 
     ! The number of nonzeroes in row iaux
     if (.not.present(iaux)) then 
        write(0,*) 'Need IAUX when ireq=nzrowreq'
        ires=-1
        return
     endif
     irw = iaux
     if (toupper(a%fida) == 'CSR') then 
        ires = a%ia2(irw+1)-a%ia2(irw)
     else if ((toupper(a%fida) == 'COO').or.(toupper(a%fida) == 'COI')) then 

        if (a%infoa(psb_srtd_) == psb_isrtdcoo_) then 
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
     else if (toupper(a%fida) == 'JAD') then 
        pia = a%ia2(2) ! points to the beginning of ia(3,png)
        pja = a%ia2(3) ! points to the beginning of ja(:)
        ja  => a%ia2(pja:)             ! the array containing the pointers to ka and aspk
        ia1 => a%ia2(pia:pja-1:3)      ! the array containing the first row index of each block
        ia2 => a%ia2(pia+1:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first jad column
        ia3 => a%ia2(pia+2:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first csr column

        idx=a%pl(irw)
        j=0
        nz=0
        blkfnd: do
           j=j+1
           if(ia1(j).eq.idx) then
              nz=nz+ia3(j)-ia2(j)
              ipx = ia1(j)         ! the first row index of the block
              rb  = idx-ipx        ! the row offset within the block
              row = ia3(j)+rb
              nz  = nz+ja(row+1)-ja(row)
              exit blkfnd
           else if(ia1(j).gt.idx) then
              nz=nz+ia3(j-1)-ia2(j-1)
              ipx = ia1(j-1)         ! the first row index of the block
              rb  = idx-ipx          ! the row offset within the block
              row = ia3(j-1)+rb
              nz  = nz+ja(row+1)-ja(row)
              exit blkfnd
           end if
        end do blkfnd
        ires=nz
     else
        ires=-1
        info=136
        ch_err=a%fida(1:3)
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

  else  if (ireq == psb_nzsizereq_) then 
     if (toupper(a%fida) == 'CSR') then 
        ires = size(a%aspk)
     else if ((toupper(a%fida) == 'COO').or.(toupper(a%fida) == 'COI')) then 
        ires = size(a%aspk)
     else if (toupper(a%fida) == 'JAD') then 
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
