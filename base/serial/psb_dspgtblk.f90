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
! File:  psb_dspgtblk.f90 
! Subroutine: psb_dspgtblk
!    Gets one or more rows from a sparse matrix. 
! Parameters:

!*****************************************************************************
!*                                                                           *
!* Takes a specified row from matrix A and copies into matrix B (possibly    *
!*  appending to B). Output is always COO. Input might be anything, once     *
!* we get to actually write the code.....                                    *
!*                                                                           *
!*****************************************************************************
subroutine psb_dspgtblk(irw,a,b,info,append,iren,lrw)
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

  name='psb_spgtblk'
  info  = 0
  call psb_erractionsave(err_act)

  irw_ = irw 
  if (present(lrw)) then
    lrw_ = lrw
  else
    lrw_ = irw
  endif
  if (lrw_ < irw) then
    write(0,*) 'SPGTBLK input error: fixing lrw',irw,lrw_
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
     call csr_dspgtblk(irw_,a,b,append_,iren_,lrw_)
     
  else if (a%fida == 'COO') then 
     call coo_dspgtblk(irw_,a,b,append_,iren_,lrw_)
     
  else if (a%fida == 'JAD') then 
     call jad_dspgtblk(irw_,a,b,append_,iren_,lrw_)

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
  if (err_act.eq.psb_act_abort_) then
     call psb_error()
     return
  end if
  return
  
contains
  
  subroutine csr_dspgtblk(irw,a,b,append,iren,lrw)

    use psb_spmat_type
    use psb_const_mod
    implicit none

    type(psb_dspmat_type), intent(in)     :: a
    integer                               :: irw
    type(psb_dspmat_type), intent(inout)  :: b
    logical, intent(in)                   :: append
    integer, pointer                      :: iren(:)
    integer                               :: lrw

    integer :: idx,i,j ,nr,nz,nzb, row_idx
    integer, pointer :: indices(:)

    if (append) then 
      nzb = b%infoa(psb_nnz_)
    else
      nzb = 0
    endif

    if (a%pl(1) /= 0) then

      nr = lrw - irw + 1 
      allocate(indices(nr))
      nz = 0
      do i=1,nr
        indices(i)=a%pl(irw+i-1)
        nz=nz+a%ia2(indices(i)+1)-a%ia2(indices(i))
      end do

      if (min(size(b%ia1),size(b%ia2),size(b%aspk)) < nzb+nz) then 
        call psb_sp_reall(b,nzb+nz,iret)
      endif

      k=0
      if(associated(iren)) then
        do i=1,nr
          row_idx=indices(i)
          do j=a%ia2(row_idx),a%ia2(row_idx+1)-1
            k             = k + 1
            b%aspk(nzb+k) = a%aspk(j)
            b%ia1(nzb+k)  = iren(row_idx)
            b%ia2(nzb+k)  = iren(a%ia1(j))
          end do
        end do
      else
        do i=1,nr
          row_idx=indices(i)
          do j=a%ia2(row_idx),a%ia2(row_idx+1)-1
            k             = k + 1
            b%aspk(nzb+k) = a%aspk(j)
            b%ia1(nzb+k)  = row_idx
            b%ia2(nzb+k)  = a%ia1(j)
          end do
        end do
      end if

      b%infoa(psb_nnz_) = nzb+k
      b%m = b%m+nr
      b%k = max(b%k,a%k)

    else
      idx = irw

      if (idx<0) then 
        write(0,*) ' spgtblk Error : idx no good ',idx
        return
      end if
      nr = lrw - irw + 1 
      nz = a%ia2(idx+nr) - a%ia2(idx)

      if (min(size(b%ia1),size(b%ia2),size(b%aspk)) < nzb+nz) then 
        call psb_sp_reall(b,nzb+nz,iret)
      endif
      b%fida='COO'

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

        do i=irw,lrw
          do j=a%ia2(i),a%ia2(i+1)-1
            k             = k + 1
            b%ia1(nzb+k)  = i
            b%ia2(nzb+k)  = a%ia1(j)
            b%aspk(nzb+k) = a%aspk(j)
!!$          write(0,*) 'csr_gtblk: in:',a%aspk(j),i,a%ia1(j)
          end do
        enddo
      end if
      b%infoa(psb_nnz_) = nzb+nz
      if (a%pr(1) /= 0) then
        write(0,*) 'Feeling lazy today, Right Permutation will have to wait'
      endif
      b%m = b%m+lrw-irw+1
      b%k = max(b%k,a%k)

    endif

  end subroutine csr_dspgtblk

  subroutine coo_dspgtblk(irw,a,b,append,iren,lrw)

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
      write(0,*) 'Fatal error in SPGTBLK: do not feed a permuted mat so far!'
      idx = -1 
    else
      idx = irw
    endif
    if (idx<0) then 
      write(0,*) ' spgtblk Error : idx no good ',idx
      return
    end if

    if (a%infoa(psb_srtd_) == psb_isrtdcoo_) then 
      ! In this case we can do a binary search. 
      do
        call ibsrch(ip,irw,nza,a%ia1)
        if (ip /= -1) exit
        irw = irw + 1
        if (irw > lrw) then
          write(0,*) 'Warning : did not find any rows. Is this an error? ',irw,lrw,idx
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
          call psb_sp_reall(b,nzb+nz,iret)
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
        call psb_sp_reall(b,nzb+nz,iret)
      endif
      
      if (associated(iren)) then 
        k = 0 
        do i=1,a%infoa(psb_nnz_)
          if ((a%ia1(i)>=irw).and.(a%ia1(i)<=lrw)) then 
            k = k + 1 
            if (k > nz) then
              nz = k 
              call psb_sp_reall(b,nzb+nz,iret)
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
              call psb_sp_reall(b,nzb+nz,iret)
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
  end subroutine coo_dspgtblk




  subroutine jad_dspgtblk(irw,a,b,append,iren,lrw)

    type(psb_dspmat_type), intent(in), target :: a
    integer                               :: irw
    type(psb_dspmat_type), intent(inout)  :: b
    logical, intent(in)                   :: append
    integer, pointer                      :: iren(:)
    integer                               :: lrw

    integer, pointer                      :: ia1(:), ia2(:), ia3(:),&
         & ja(:), ka(:), indices(:), blks(:)
    integer  :: png, pia, pja, ipx, blk, rb, row, k_pt, npg, col, ng


    png = a%ia2(1) ! points to the number of blocks
    pia = a%ia2(2) ! points to the beginning of ia(3,png)
    pja = a%ia2(3) ! points to the beginning of ja(:)
    
    ng  =  a%ia2(png)              ! the number of blocks
    ja  => a%ia2(pja:)             ! the array containing the pointers to ka and aspk
    ka  => a%ia1(:)                ! the array containing the column indices
    ia1 => a%ia2(pia:pja-1:3)      ! the array containing the first row index of each block
    ia2 => a%ia2(pia+1:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first jad column
    ia3 => a%ia2(pia+2:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first csr column

    if (append) then 
       nzb = b%infoa(psb_nnz_)
    else
       nzb = 0
    endif

    if (a%pl(1) /= 0) then

       nr = lrw - irw + 1 
       allocate(indices(nr),blks(nr))
       nz = 0

       do i=1,nr
          indices(i)=a%pl(irw+i-1)
          j=0
          blkfnd: do
             j=j+1
             if(ia1(j).eq.indices(i)) then
                blks(i)=j
                nz=nz+ia3(j)-ia2(j)
                ipx = ia1(j)         ! the first row index of the block
                rb  = indices(i)-ipx   ! the row offset within the block
                row = ia3(j)+rb
                nz  = nz+ja(row+1)-ja(row)
                exit blkfnd
             else if(ia1(j).gt.indices(i)) then
                blks(i)=j-1
                nz=nz+ia3(j-1)-ia2(j-1)
                ipx = ia1(j-1)         ! the first row index of the block
                rb  = indices(i)-ipx   ! the row offset within the block
                row = ia3(j-1)+rb
                nz  = nz+ja(row+1)-ja(row)
                exit blkfnd
             end if
          end do blkfnd
       end do
       
       if (size(b%ia1) < nzb+nz) then 
          call psb_sp_reall(b,nzb+nz,iret)
       endif

       k=0
       ! cycle over rows
       do i=1,nr

          ! find which block the row belongs to
          blk = blks(i)

          ! extract first part of the row from the jad block
          ipx = ia1(blk)             ! the first row index of the block
          k_pt= ia2(blk)             ! the pointer to the beginning of a column in ja
          rb  = indices(i)-ipx       ! the row offset within the block
          npg = ja(k_pt+1)-ja(k_pt)  ! the number of rows in the block

          if(associated(iren))then
             do  col = ia2(blk), ia3(blk)-1 
                k=k+1
                b%aspk(nzb+k) = a%aspk(ja(col)+rb)
                b%ia1(nzb+k)  = iren(irw+i-1)
                b%ia2(nzb+k)  = iren(ka(ja(col)+rb))
             end do
          else
             do  col = ia2(blk), ia3(blk)-1 
                k=k+1
                b%aspk(nzb+k) = a%aspk(ja(col)+rb)
                b%ia1(nzb+k)  = irw+i-1
                b%ia2(nzb+k)  = ka(ja(col)+rb)
             end do
          end if
          ! extract second part of the row from the csr tail
          row=ia3(blk)+rb
          if(associated(iren))then
             do j=ja(row), ja(row+1)-1
                k=k+1
                b%aspk(nzb+k) = a%aspk(j)
                b%ia1(nzb+k)  = iren(irw+i-1)
                b%ia2(nzb+k)  = iren(ka(j))
             end do
          else
             do j=ja(row), ja(row+1)-1
                k=k+1
                b%aspk(nzb+k) = a%aspk(j)
                b%ia1(nzb+k)  = irw+i-1
                b%ia2(nzb+k)  = ka(j)
             end do
          end if
       end do

       b%infoa(psb_nnz_) = nzb+k
       b%m = b%m+lrw-irw+1
       b%k = max(b%k,a%k)
       b%fida='COO'
       
    else
       ! There might be some problems
       info=134
       ch_err=a%fida(1:3)
       call psb_errpush(info,name,a_err=ch_err)
    end if

    

  end subroutine jad_dspgtblk




end subroutine psb_dspgtblk

