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
module psb_getrow_mod

  interface csr_getrow
    module procedure csr_cspgtrow, csr_zspgtrow
  end interface
  interface coo_getrow
    module procedure coo_cspgtrow, coo_zspgtrow
  end interface
  interface jad_getrow
    module procedure jad_cspgtrow, jad_zspgtrow
  end interface

contains
  
!!$  subroutine csr_sspgtrow(irw,a,nz,ia,ja,val,nzin,append,lrw,info,iren)
!!$    
!!$    use psb_sort_mod
!!$    use psb_spmat_type
!!$    use psb_const_mod
!!$    implicit none
!!$
!!$    type(psb_sspmat_type), intent(in)    :: a
!!$    integer                              :: irw
!!$    integer, intent(out)                 :: nz
!!$    integer, allocatable, intent(inout)  :: ia(:), ja(:)
!!$    real(psb_spk_), allocatable,  intent(inout)    :: val(:)
!!$    integer                              :: nzin
!!$    logical, intent(in)                  :: append
!!$    integer                              :: lrw,info
!!$    integer, optional                    :: iren(:)
!!$
!!$    integer :: idx,i,j, k, nr, row_idx, nzin_
!!$    integer, allocatable  :: indices(:)
!!$
!!$    if (append) then 
!!$      nzin_ = nzin
!!$    else
!!$      nzin_ = 0
!!$    endif
!!$
!!$    if (a%pl(1) /= 0) then
!!$
!!$      nr = lrw - irw + 1 
!!$      allocate(indices(nr))
!!$      nz = 0
!!$      do i=1,nr
!!$        indices(i)=a%pl(irw+i-1)
!!$        nz=nz+a%ia2(indices(i)+1)-a%ia2(indices(i))
!!$      end do
!!$      
!!$      call psb_ensure_size(nzin_+nz,ia,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nz,val,info)
!!$      if (info /= 0) return
!!$
!!$      k=0
!!$      if(present(iren)) then
!!$        do i=1,nr
!!$          row_idx=indices(i)
!!$          do j=a%ia2(row_idx),a%ia2(row_idx+1)-1
!!$            k            = k + 1
!!$            val(nzin_+k) = a%aspk(j)
!!$            ia(nzin_+k)  = iren(row_idx)
!!$            ja(nzin_+k)  = iren(a%ia1(j))
!!$          end do
!!$        end do
!!$      else
!!$        do i=1,nr
!!$          row_idx=indices(i)
!!$          do j=a%ia2(row_idx),a%ia2(row_idx+1)-1
!!$            k            = k + 1
!!$            val(nzin_+k) = a%aspk(j)
!!$            ia(nzin_+k)  = row_idx
!!$            ja(nzin_+k)  = a%ia1(j)
!!$          end do
!!$        end do
!!$      end if
!!$
!!$    else
!!$
!!$      idx = irw
!!$
!!$      if (idx<0) then 
!!$        write(0,*) ' spgtrow Error : idx no good ',idx
!!$        info = 2
!!$        return
!!$      end if
!!$      nr = lrw - irw + 1 
!!$      nz = a%ia2(idx+nr) - a%ia2(idx)
!!$
!!$      call psb_ensure_size(nzin_+nz,ia,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nz,val,info)
!!$      if (info /= 0) return
!!$
!!$
!!$      if (present(iren)) then 
!!$        k=0
!!$        do i=irw,lrw
!!$          do j=a%ia2(i),a%ia2(i+1)-1
!!$            k            = k + 1
!!$            val(nzin_+k) = a%aspk(j)
!!$            ia(nzin_+k)  = iren(i)
!!$            ja(nzin_+k)  = iren(a%ia1(j))
!!$          end do
!!$        enddo
!!$      else
!!$        k=0
!!$
!!$        do i=irw,lrw
!!$          do j=a%ia2(i),a%ia2(i+1)-1
!!$            k            = k + 1
!!$            ia(nzin_+k)  = i
!!$            ja(nzin_+k)  = a%ia1(j)
!!$            val(nzin_+k) = a%aspk(j)
!!$          end do
!!$        enddo
!!$      end if
!!$! !$      if (nz /= k) then 
!!$! !$        write(0,*) 'csr getrow Size mismatch ',nz,k
!!$! !$      endif
!!$      if (a%pr(1) /= 0) then
!!$        write(0,*) 'Feeling lazy today, Right Permutation will have to wait'
!!$      endif
!!$
!!$    endif
!!$
!!$  end subroutine csr_sspgtrow
!!$
!!$  subroutine coo_sspgtrow(irw,a,nz,ia,ja,val,nzin,append,lrw,info,iren)
!!$
!!$    use psb_sort_mod
!!$    use psb_spmat_type
!!$    use psb_const_mod
!!$    use psb_error_mod
!!$    implicit none
!!$
!!$    type(psb_sspmat_type), intent(in)    :: a
!!$    integer                              :: irw
!!$    integer, intent(out)                 :: nz
!!$    integer, allocatable, intent(inout)  :: ia(:), ja(:)
!!$    real(psb_spk_), allocatable,  intent(inout)    :: val(:)
!!$    integer                              :: nzin
!!$    logical, intent(in)                  :: append
!!$    integer                              :: lrw,info
!!$    integer, optional                    :: iren(:)
!!$    integer  :: nzin_, nza, idx,ip,jp,i,k, nzt
!!$    integer  :: debug_level, debug_unit
!!$    character(len=20) :: name='coo_sspgtrow'
!!$
!!$    debug_unit  = psb_get_debug_unit()
!!$    debug_level = psb_get_debug_level()
!!$
!!$    nza = a%infoa(psb_nnz_)
!!$    if (a%pl(1) /= 0) then
!!$      write(debug_unit,*) 'Fatal error in SPGTROW: do not feed a permuted mat so far!'
!!$      idx = -1 
!!$    else
!!$      idx = irw
!!$    endif
!!$    if (idx<0) then 
!!$      write(debug_unit,*) ' spgtrow Error : idx no good ',idx
!!$      info = 2
!!$      return
!!$    end if
!!$
!!$    if (append) then 
!!$      nzin_ = nzin
!!$    else
!!$      nzin_ = 0
!!$    endif
!!$
!!$    if (a%infoa(psb_srtd_) == psb_isrtdcoo_) then 
!!$      ! In this case we can do a binary search. 
!!$      if (debug_level >= psb_debug_serial_)&
!!$           & write(debug_unit,*) trim(name), ': srtdcoo '
!!$      do
!!$        ip = psb_ibsrch(irw,nza,a%ia1)
!!$        if (ip /= -1) exit
!!$        irw = irw + 1
!!$        if (irw > lrw) then
!!$          write(debug_unit,*)  trim(name),&
!!$               & 'Warning : did not find any rows. Is this an error? ',&
!!$               & irw,lrw,idx
!!$          exit
!!$        end if
!!$      end do
!!$      
!!$      if (ip /= -1) then 
!!$        ! expand [ip,jp] to contain all row entries.
!!$        do 
!!$          if (ip < 2) exit
!!$          if (a%ia1(ip-1) == irw) then  
!!$            ip = ip -1 
!!$          else 
!!$            exit
!!$          end if
!!$        end do
!!$
!!$      end if
!!$      
!!$      do
!!$        jp = psb_ibsrch(lrw,nza,a%ia1)
!!$        if (jp /= -1) exit
!!$        lrw = lrw - 1
!!$        if (irw > lrw) then
!!$          write(debug_unit,*) trim(name),&
!!$               & 'Warning : did not find any rows. Is this an error?'
!!$          exit
!!$        end if
!!$      end do
!!$      
!!$      if (jp /= -1) then 
!!$        ! expand [ip,jp] to contain all row entries.
!!$        do 
!!$          if (jp == nza) exit
!!$          if (a%ia1(jp+1) == lrw) then  
!!$            jp = jp + 1
!!$          else 
!!$            exit
!!$          end if
!!$        end do
!!$      end if
!!$      if (debug_level >= psb_debug_serial_) &
!!$           & write(debug_unit,*)  trim(name),': ip jp',ip,jp,nza
!!$      if ((ip /= -1) .and.(jp /= -1)) then 
!!$        ! Now do the copy.
!!$        nz = jp - ip +1 
!!$
!!$        call psb_ensure_size(nzin_+nz,ia,info)
!!$        if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
!!$        if (info==0) call psb_ensure_size(nzin_+nz,val,info)
!!$        if (info /= 0) return
!!$
!!$        if (present(iren)) then 
!!$          do i=ip,jp
!!$            nzin_ = nzin_ + 1
!!$            val(nzin_) = a%aspk(i)
!!$            ia(nzin_)  = iren(a%ia1(i))
!!$            ja(nzin_)  = iren(a%ia2(i))
!!$          enddo
!!$        else
!!$          do i=ip,jp
!!$            nzin_ = nzin_ + 1
!!$            val(nzin_) = a%aspk(i)
!!$            ia(nzin_)  = a%ia1(i)
!!$            ja(nzin_)  = a%ia2(i)
!!$          enddo
!!$        end if
!!$      else 
!!$        nz = 0 
!!$      end if
!!$
!!$    else
!!$      if (debug_level >= psb_debug_serial_) &
!!$           & write(debug_unit,*)  trim(name),': unsorted '
!!$
!!$      nzt = (nza*(lrw-irw+1))/max(a%m,1)
!!$      call psb_ensure_size(nzin_+nzt,ia,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nzt,ja,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nzt,val,info)
!!$      if (info /= 0) return
!!$      
!!$      if (present(iren)) then 
!!$        k = 0 
!!$        do i=1, a%infoa(psb_nnz_)
!!$          if ((a%ia1(i)>=irw).and.(a%ia1(i)<=lrw)) then 
!!$            k = k + 1 
!!$            if (k > nzt) then
!!$              nzt = k 
!!$              call psb_ensure_size(nzin_+nzt,ia,info)
!!$              if (info==0) call psb_ensure_size(nzin_+nzt,ja,info)
!!$              if (info==0) call psb_ensure_size(nzin_+nzt,val,info)
!!$              if (info /= 0) return
!!$            end if
!!$            val(nzin_+k) = a%aspk(i)
!!$            ia(nzin_+k)  = iren(a%ia1(i))
!!$            ja(nzin_+k)  = iren(a%ia2(i))
!!$          endif
!!$        enddo
!!$      else
!!$        k = 0 
!!$        do i=1,a%infoa(psb_nnz_)
!!$          if ((a%ia1(i)>=irw).and.(a%ia1(i)<=lrw)) then 
!!$            k = k + 1 
!!$            if (k > nzt) then
!!$              nzt = k 
!!$              call psb_ensure_size(nzin_+nzt,ia,info)
!!$              if (info==0) call psb_ensure_size(nzin_+nzt,ja,info)
!!$              if (info==0) call psb_ensure_size(nzin_+nzt,val,info)
!!$              if (info /= 0) return
!!$
!!$            end if
!!$            val(nzin_+k) = a%aspk(i)
!!$            ia(nzin_+k)  = (a%ia1(i))
!!$            ja(nzin_+k)  = (a%ia2(i))
!!$          endif
!!$        enddo
!!$        nzin_=nzin_+k
!!$      end if
!!$      nz = k 
!!$    end if
!!$
!!$  end subroutine coo_sspgtrow
!!$
!!$
!!$  subroutine jad_sspgtrow(irw,a,nz,ia,ja,val,nzin,append,lrw,info,iren)
!!$
!!$    use psb_sort_mod
!!$    use psb_spmat_type
!!$    use psb_const_mod
!!$    
!!$    implicit none
!!$
!!$    type(psb_sspmat_type), intent(in), target :: a
!!$    integer                               :: irw
!!$    integer, intent(out)                 :: nz
!!$    integer, allocatable, intent(inout)  :: ia(:), ja(:)
!!$    real(psb_spk_), allocatable,  intent(inout)    :: val(:)
!!$    integer                               :: nzin
!!$    logical, intent(in)                   :: append
!!$    integer, optional                     :: iren(:)
!!$    integer                               :: lrw,info
!!$
!!$    integer, pointer                      :: ia1(:), ia2(:), ia3(:),&
!!$         & ja_(:), ka_(:), indices(:), blks(:)
!!$    integer  :: png, pia, pja, ipx, blk, rb, row, k_pt, npg, col, ng, nzin_,&
!!$         & i,j,k,nr
!!$
!!$
!!$    png = a%ia2(1) ! points to the number of blocks
!!$    pia = a%ia2(2) ! points to the beginning of ia(3,png)
!!$    pja = a%ia2(3) ! points to the beginning of ja(:)
!!$
!!$    ng  =  a%ia2(png)              ! the number of blocks
!!$    ja_  => a%ia2(pja:)             ! the array containing the pointers to ka and aspk
!!$    ka_  => a%ia1(:)                ! the array containing the column indices
!!$    ia1 => a%ia2(pia:pja-1:3)      ! the array containing the first row index of each block
!!$    ia2 => a%ia2(pia+1:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first jad column
!!$    ia3 => a%ia2(pia+2:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first csr column
!!$
!!$    if (append) then 
!!$      nzin_ = nzin
!!$    else
!!$      nzin_ = 0
!!$    endif
!!$
!!$    if (a%pl(1) /= 0) then
!!$
!!$      nr = lrw - irw + 1 
!!$      allocate(indices(nr),blks(nr))
!!$      nz = 0
!!$
!!$      do i=1,nr
!!$        indices(i)=a%pl(irw+i-1)
!!$        j=0
!!$        blkfnd: do
!!$          j=j+1
!!$          if(ia1(j) == indices(i)) then
!!$            blks(i)=j
!!$            nz=nz+ia3(j)-ia2(j)
!!$            ipx = ia1(j)         ! the first row index of the block
!!$            rb  = indices(i)-ipx   ! the row offset within the block
!!$            row = ia3(j)+rb
!!$            nz  = nz+ja_(row+1)-ja_(row)
!!$            exit blkfnd
!!$          else if(ia1(j) > indices(i)) then
!!$            blks(i)=j-1
!!$            nz=nz+ia3(j-1)-ia2(j-1)
!!$            ipx = ia1(j-1)         ! the first row index of the block
!!$            rb  = indices(i)-ipx   ! the row offset within the block
!!$            row = ia3(j-1)+rb
!!$            nz  = nz+ja_(row+1)-ja_(row)
!!$            exit blkfnd
!!$          end if
!!$        end do blkfnd
!!$      end do
!!$
!!$
!!$      call psb_ensure_size(nzin_+nz,ia,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nz,val,info)
!!$      if (info /= 0) return
!!$
!!$      k=0
!!$      ! cycle over rows
!!$      do i=1,nr
!!$
!!$        ! find which block the row belongs to
!!$        blk = blks(i)
!!$
!!$        ! extract first part of the row from the jad block
!!$        ipx = ia1(blk)             ! the first row index of the block
!!$        k_pt= ia2(blk)             ! the pointer to the beginning of a column in ja
!!$        rb  = indices(i)-ipx       ! the row offset within the block
!!$        npg = ja_(k_pt+1)-ja_(k_pt)  ! the number of rows in the block
!!$
!!$        if(present(iren))then
!!$          do  col = ia2(blk), ia3(blk)-1 
!!$            k=k+1
!!$            val(nzin_+k) = a%aspk(ja_(col)+rb)
!!$            ia(nzin_+k)  = iren(irw+i-1)
!!$            ja(nzin_+k)  = iren(ka_(ja_(col)+rb))
!!$          end do
!!$        else
!!$          do  col = ia2(blk), ia3(blk)-1 
!!$            k=k+1
!!$            val(nzin_+k) = a%aspk(ja_(col)+rb)
!!$            ia(nzin_+k)  = irw+i-1
!!$            ja(nzin_+k)  = ka_(ja_(col)+rb)
!!$          end do
!!$        end if
!!$        ! extract second part of the row from the csr tail
!!$        row=ia3(blk)+rb
!!$        if(present(iren))then
!!$          do j=ja_(row), ja_(row+1)-1
!!$            k=k+1
!!$            val(nzin_+k) = a%aspk(j)
!!$            ia(nzin_+k)  = iren(irw+i-1)
!!$            ja(nzin_+k)  = iren(ka_(j))
!!$          end do
!!$        else
!!$          do j=ja_(row), ja_(row+1)-1
!!$            k=k+1
!!$            val(nzin_+k) = a%aspk(j)
!!$            ia(nzin_+k)  = irw+i-1
!!$            ja(nzin_+k)  = ka_(j)
!!$          end do
!!$        end if
!!$      end do
!!$
!!$    else
!!$      ! There might be some problems
!!$      info=134
!!$    end if
!!$
!!$  end subroutine jad_sspgtrow

  
!!$  subroutine csr_dspgtrow(irw,a,nz,ia,ja,val,nzin,append,lrw,info,iren)
!!$
!!$    use psb_sort_mod
!!$    use psb_spmat_type
!!$    use psb_const_mod
!!$    implicit none
!!$
!!$    type(psb_dspmat_type), intent(in)    :: a
!!$    integer                              :: irw
!!$    integer, intent(out)                 :: nz
!!$    integer, allocatable, intent(inout)  :: ia(:), ja(:)
!!$    real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
!!$    integer                              :: nzin
!!$    logical, intent(in)                  :: append
!!$    integer                              :: lrw,info
!!$    integer, optional                    :: iren(:)
!!$
!!$    integer :: idx,i,j, k, nr, row_idx, nzin_
!!$    integer, allocatable  :: indices(:)
!!$
!!$    if (append) then 
!!$      nzin_ = nzin
!!$    else
!!$      nzin_ = 0
!!$    endif
!!$
!!$    if (a%pl(1) /= 0) then
!!$
!!$      nr = lrw - irw + 1 
!!$      allocate(indices(nr))
!!$      nz = 0
!!$      do i=1,nr
!!$        indices(i)=a%pl(irw+i-1)
!!$        nz=nz+a%ia2(indices(i)+1)-a%ia2(indices(i))
!!$      end do
!!$      
!!$      call psb_ensure_size(nzin_+nz,ia,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nz,val,info)
!!$      if (info /= 0) return
!!$
!!$      k=0
!!$      if(present(iren)) then
!!$        do i=1,nr
!!$          row_idx=indices(i)
!!$          do j=a%ia2(row_idx),a%ia2(row_idx+1)-1
!!$            k            = k + 1
!!$            val(nzin_+k) = a%aspk(j)
!!$            ia(nzin_+k)  = iren(row_idx)
!!$            ja(nzin_+k)  = iren(a%ia1(j))
!!$          end do
!!$        end do
!!$      else
!!$        do i=1,nr
!!$          row_idx=indices(i)
!!$          do j=a%ia2(row_idx),a%ia2(row_idx+1)-1
!!$            k            = k + 1
!!$            val(nzin_+k) = a%aspk(j)
!!$            ia(nzin_+k)  = row_idx
!!$            ja(nzin_+k)  = a%ia1(j)
!!$          end do
!!$        end do
!!$      end if
!!$
!!$    else
!!$
!!$      idx = irw
!!$
!!$      if (idx<0) then 
!!$        write(0,*) ' spgtrow Error : idx no good ',idx
!!$        info = 2
!!$        return
!!$      end if
!!$      nr = lrw - irw + 1 
!!$      nz = a%ia2(idx+nr) - a%ia2(idx)
!!$
!!$      call psb_ensure_size(nzin_+nz,ia,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nz,val,info)
!!$      if (info /= 0) return
!!$
!!$
!!$      if (present(iren)) then 
!!$        k=0
!!$        do i=irw,lrw
!!$          do j=a%ia2(i),a%ia2(i+1)-1
!!$            k            = k + 1
!!$            val(nzin_+k) = a%aspk(j)
!!$            ia(nzin_+k)  = iren(i)
!!$            ja(nzin_+k)  = iren(a%ia1(j))
!!$          end do
!!$        enddo
!!$      else
!!$        k=0
!!$
!!$        do i=irw,lrw
!!$          do j=a%ia2(i),a%ia2(i+1)-1
!!$            k            = k + 1
!!$            ia(nzin_+k)  = i
!!$            ja(nzin_+k)  = a%ia1(j)
!!$            val(nzin_+k) = a%aspk(j)
!!$          end do
!!$        enddo
!!$      end if
! !$      if (nz /= k) then 
! !$        write(0,*) 'csr getrow Size mismatch ',nz,k
! !$      endif
!!$      if (a%pr(1) /= 0) then
!!$        write(0,*) 'Feeling lazy today, Right Permutation will have to wait'
!!$      endif
!!$
!!$    endif
!!$
!!$  end subroutine csr_dspgtrow

!!$  subroutine coo_dspgtrow(irw,a,nz,ia,ja,val,nzin,append,lrw,info,iren)
!!$
!!$    use psb_sort_mod
!!$    use psb_spmat_type
!!$    use psb_const_mod
!!$    use psb_error_mod
!!$    implicit none
!!$
!!$    type(psb_dspmat_type), intent(in)    :: a
!!$    integer                              :: irw
!!$    integer, intent(out)                 :: nz
!!$    integer, allocatable, intent(inout)  :: ia(:), ja(:)
!!$    real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
!!$    integer                              :: nzin
!!$    logical, intent(in)                  :: append
!!$    integer                              :: lrw,info
!!$    integer, optional                    :: iren(:)
!!$    integer  :: nzin_, nza, idx,ip,jp,i,k, nzt
!!$    integer  :: debug_level, debug_unit
!!$    character(len=20) :: name='coo_dspgtrow'
!!$
!!$    debug_unit  = psb_get_debug_unit()
!!$    debug_level = psb_get_debug_level()
!!$
!!$    nza = a%infoa(psb_nnz_)
!!$    if (a%pl(1) /= 0) then
!!$      write(debug_unit,*) 'Fatal error in SPGTROW: do not feed a permuted mat so far!'
!!$      idx = -1 
!!$    else
!!$      idx = irw
!!$    endif
!!$    if (idx<0) then 
!!$      write(debug_unit,*) ' spgtrow Error : idx no good ',idx
!!$      info = 2
!!$      return
!!$    end if
!!$
!!$    if (append) then 
!!$      nzin_ = nzin
!!$    else
!!$      nzin_ = 0
!!$    endif
!!$
!!$    if (a%infoa(psb_srtd_) == psb_isrtdcoo_) then 
!!$      ! In this case we can do a binary search. 
!!$      if (debug_level >= psb_debug_serial_)&
!!$           & write(debug_unit,*) trim(name), ': srtdcoo '
!!$      do
!!$        ip = psb_ibsrch(irw,nza,a%ia1)
!!$        if (ip /= -1) exit
!!$        irw = irw + 1
!!$        if (irw > lrw) then
!!$          write(debug_unit,*)  trim(name),&
!!$               & 'Warning : did not find any rows. Is this an error? ',&
!!$               & irw,lrw,idx
!!$          exit
!!$        end if
!!$      end do
!!$      
!!$      if (ip /= -1) then 
!!$        ! expand [ip,jp] to contain all row entries.
!!$        do 
!!$          if (ip < 2) exit
!!$          if (a%ia1(ip-1) == irw) then  
!!$            ip = ip -1 
!!$          else 
!!$            exit
!!$          end if
!!$        end do
!!$
!!$      end if
!!$      
!!$      do
!!$        jp = psb_ibsrch(lrw,nza,a%ia1)
!!$        if (jp /= -1) exit
!!$        lrw = lrw - 1
!!$        if (irw > lrw) then
!!$          write(debug_unit,*) trim(name),&
!!$               & 'Warning : did not find any rows. Is this an error?'
!!$          exit
!!$        end if
!!$      end do
!!$      
!!$      if (jp /= -1) then 
!!$        ! expand [ip,jp] to contain all row entries.
!!$        do 
!!$          if (jp == nza) exit
!!$          if (a%ia1(jp+1) == lrw) then  
!!$            jp = jp + 1
!!$          else 
!!$            exit
!!$          end if
!!$        end do
!!$      end if
!!$      if (debug_level >= psb_debug_serial_) &
!!$           & write(debug_unit,*)  trim(name),': ip jp',ip,jp,nza
!!$      if ((ip /= -1) .and.(jp /= -1)) then 
!!$        ! Now do the copy.
!!$        nz = jp - ip +1 
!!$
!!$        call psb_ensure_size(nzin_+nz,ia,info)
!!$        if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
!!$        if (info==0) call psb_ensure_size(nzin_+nz,val,info)
!!$        if (info /= 0) return
!!$
!!$        if (present(iren)) then 
!!$          do i=ip,jp
!!$            nzin_ = nzin_ + 1
!!$            val(nzin_) = a%aspk(i)
!!$            ia(nzin_)  = iren(a%ia1(i))
!!$            ja(nzin_)  = iren(a%ia2(i))
!!$          enddo
!!$        else
!!$          do i=ip,jp
!!$            nzin_ = nzin_ + 1
!!$            val(nzin_) = a%aspk(i)
!!$            ia(nzin_)  = a%ia1(i)
!!$            ja(nzin_)  = a%ia2(i)
!!$          enddo
!!$        end if
!!$      else 
!!$        nz = 0 
!!$      end if
!!$
!!$    else
!!$      if (debug_level >= psb_debug_serial_) &
!!$           & write(debug_unit,*)  trim(name),': unsorted '
!!$
!!$      nzt = (nza*(lrw-irw+1))/max(a%m,1)
!!$      call psb_ensure_size(nzin_+nzt,ia,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nzt,ja,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nzt,val,info)
!!$      if (info /= 0) return
!!$      
!!$      if (present(iren)) then 
!!$        k = 0 
!!$        do i=1, a%infoa(psb_nnz_)
!!$          if ((a%ia1(i)>=irw).and.(a%ia1(i)<=lrw)) then 
!!$            k = k + 1 
!!$            if (k > nzt) then
!!$              nzt = k 
!!$              call psb_ensure_size(nzin_+nzt,ia,info)
!!$              if (info==0) call psb_ensure_size(nzin_+nzt,ja,info)
!!$              if (info==0) call psb_ensure_size(nzin_+nzt,val,info)
!!$              if (info /= 0) return
!!$            end if
!!$            val(nzin_+k) = a%aspk(i)
!!$            ia(nzin_+k)  = iren(a%ia1(i))
!!$            ja(nzin_+k)  = iren(a%ia2(i))
!!$          endif
!!$        enddo
!!$      else
!!$        k = 0 
!!$        do i=1,a%infoa(psb_nnz_)
!!$          if ((a%ia1(i)>=irw).and.(a%ia1(i)<=lrw)) then 
!!$            k = k + 1 
!!$            if (k > nzt) then
!!$              nzt = k 
!!$              call psb_ensure_size(nzin_+nzt,ia,info)
!!$              if (info==0) call psb_ensure_size(nzin_+nzt,ja,info)
!!$              if (info==0) call psb_ensure_size(nzin_+nzt,val,info)
!!$              if (info /= 0) return
!!$
!!$            end if
!!$            val(nzin_+k) = a%aspk(i)
!!$            ia(nzin_+k)  = (a%ia1(i))
!!$            ja(nzin_+k)  = (a%ia2(i))
!!$          endif
!!$        enddo
!!$        nzin_=nzin_+k
!!$      end if
!!$      nz = k 
!!$    end if
!!$
!!$  end subroutine coo_dspgtrow
!!$
!!$
!!$  subroutine jad_dspgtrow(irw,a,nz,ia,ja,val,nzin,append,lrw,info,iren)
!!$
!!$    use psb_sort_mod
!!$    use psb_spmat_type
!!$    use psb_const_mod
!!$    
!!$    implicit none
!!$
!!$    type(psb_dspmat_type), intent(in), target :: a
!!$    integer                               :: irw
!!$    integer, intent(out)                 :: nz
!!$    integer, allocatable, intent(inout)  :: ia(:), ja(:)
!!$    real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
!!$    integer                               :: nzin
!!$    logical, intent(in)                   :: append
!!$    integer, optional                     :: iren(:)
!!$    integer                               :: lrw,info
!!$
!!$    integer, pointer                      :: ia1(:), ia2(:), ia3(:),&
!!$         & ja_(:), ka_(:), indices(:), blks(:)
!!$    integer  :: png, pia, pja, ipx, blk, rb, row, k_pt, npg, col, ng, nzin_,&
!!$         & i,j,k,nr
!!$
!!$
!!$    png = a%ia2(1) ! points to the number of blocks
!!$    pia = a%ia2(2) ! points to the beginning of ia(3,png)
!!$    pja = a%ia2(3) ! points to the beginning of ja(:)
!!$
!!$    ng  =  a%ia2(png)              ! the number of blocks
!!$    ja_  => a%ia2(pja:)             ! the array containing the pointers to ka and aspk
!!$    ka_  => a%ia1(:)                ! the array containing the column indices
!!$    ia1 => a%ia2(pia:pja-1:3)      ! the array containing the first row index of each block
!!$    ia2 => a%ia2(pia+1:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first jad column
!!$    ia3 => a%ia2(pia+2:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first csr column
!!$
!!$    if (append) then 
!!$      nzin_ = nzin
!!$    else
!!$      nzin_ = 0
!!$    endif
!!$
!!$    if (a%pl(1) /= 0) then
!!$
!!$      nr = lrw - irw + 1 
!!$      allocate(indices(nr),blks(nr))
!!$      nz = 0
!!$
!!$      do i=1,nr
!!$        indices(i)=a%pl(irw+i-1)
!!$        j=0
!!$        blkfnd: do
!!$          j=j+1
!!$          if(ia1(j) == indices(i)) then
!!$            blks(i)=j
!!$            nz=nz+ia3(j)-ia2(j)
!!$            ipx = ia1(j)         ! the first row index of the block
!!$            rb  = indices(i)-ipx   ! the row offset within the block
!!$            row = ia3(j)+rb
!!$            nz  = nz+ja_(row+1)-ja_(row)
!!$            exit blkfnd
!!$          else if(ia1(j) > indices(i)) then
!!$            blks(i)=j-1
!!$            nz=nz+ia3(j-1)-ia2(j-1)
!!$            ipx = ia1(j-1)         ! the first row index of the block
!!$            rb  = indices(i)-ipx   ! the row offset within the block
!!$            row = ia3(j-1)+rb
!!$            nz  = nz+ja_(row+1)-ja_(row)
!!$            exit blkfnd
!!$          end if
!!$        end do blkfnd
!!$      end do
!!$
!!$
!!$      call psb_ensure_size(nzin_+nz,ia,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
!!$      if (info==0) call psb_ensure_size(nzin_+nz,val,info)
!!$      if (info /= 0) return
!!$
!!$      k=0
!!$      ! cycle over rows
!!$      do i=1,nr
!!$
!!$        ! find which block the row belongs to
!!$        blk = blks(i)
!!$
!!$        ! extract first part of the row from the jad block
!!$        ipx = ia1(blk)             ! the first row index of the block
!!$        k_pt= ia2(blk)             ! the pointer to the beginning of a column in ja
!!$        rb  = indices(i)-ipx       ! the row offset within the block
!!$        npg = ja_(k_pt+1)-ja_(k_pt)  ! the number of rows in the block
!!$
!!$        if(present(iren))then
!!$          do  col = ia2(blk), ia3(blk)-1 
!!$            k=k+1
!!$            val(nzin_+k) = a%aspk(ja_(col)+rb)
!!$            ia(nzin_+k)  = iren(irw+i-1)
!!$            ja(nzin_+k)  = iren(ka_(ja_(col)+rb))
!!$          end do
!!$        else
!!$          do  col = ia2(blk), ia3(blk)-1 
!!$            k=k+1
!!$            val(nzin_+k) = a%aspk(ja_(col)+rb)
!!$            ia(nzin_+k)  = irw+i-1
!!$            ja(nzin_+k)  = ka_(ja_(col)+rb)
!!$          end do
!!$        end if
!!$        ! extract second part of the row from the csr tail
!!$        row=ia3(blk)+rb
!!$        if(present(iren))then
!!$          do j=ja_(row), ja_(row+1)-1
!!$            k=k+1
!!$            val(nzin_+k) = a%aspk(j)
!!$            ia(nzin_+k)  = iren(irw+i-1)
!!$            ja(nzin_+k)  = iren(ka_(j))
!!$          end do
!!$        else
!!$          do j=ja_(row), ja_(row+1)-1
!!$            k=k+1
!!$            val(nzin_+k) = a%aspk(j)
!!$            ia(nzin_+k)  = irw+i-1
!!$            ja(nzin_+k)  = ka_(j)
!!$          end do
!!$        end if
!!$      end do
!!$
!!$    else
!!$      ! There might be some problems
!!$      info=134
!!$    end if
!!$
!!$  end subroutine jad_dspgtrow

  
  subroutine csr_cspgtrow(irw,a,nz,ia,ja,val,nzin,append,lrw,info,iren)

    use psb_sort_mod
    use psb_spmat_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    type(psb_cspmat_type), intent(in)    :: a
    integer                              :: irw
    integer, intent(out)                 :: nz
    integer, allocatable, intent(inout)  :: ia(:), ja(:)
    complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
    integer                              :: nzin
    logical, intent(in)                  :: append
    integer                              :: lrw,info
    integer, optional                    :: iren(:)

    integer :: idx,i,j, k, nr, row_idx, nzin_
    integer, allocatable  :: indices(:)
    integer             :: debug_level, debug_unit

    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    if (append) then 
      nzin_ = nzin
    else
      nzin_ = 0
    endif

    if (a%pl(1) /= 0) then

      nr = lrw - irw + 1 
      allocate(indices(nr))
      nz = 0
      do i=1,nr
        indices(i)=a%pl(irw+i-1)
        nz=nz+a%ia2(indices(i)+1)-a%ia2(indices(i))
      end do
      
      call psb_ensure_size(nzin_+nz,ia,info)
      if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
      if (info==0) call psb_ensure_size(nzin_+nz,val,info)
      if (info /= 0) return

      k=0
      if(present(iren)) then
        do i=1,nr
          row_idx=indices(i)
          do j=a%ia2(row_idx),a%ia2(row_idx+1)-1
            k            = k + 1
            val(nzin_+k) = a%aspk(j)
            ia(nzin_+k)  = iren(row_idx)
            ja(nzin_+k)  = iren(a%ia1(j))
          end do
        end do
      else
        do i=1,nr
          row_idx=indices(i)
          do j=a%ia2(row_idx),a%ia2(row_idx+1)-1
            k            = k + 1
            val(nzin_+k) = a%aspk(j)
            ia(nzin_+k)  = row_idx
            ja(nzin_+k)  = a%ia1(j)
          end do
        end do
      end if

    else
      idx = irw

      if (idx<0) then 
        write(0,*) ' spgtrow Error : idx no good ',idx
        return
      end if
      nr = lrw - irw + 1 
      nz = a%ia2(idx+nr) - a%ia2(idx)

      call psb_ensure_size(nzin_+nz,ia,info)
      if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
      if (info==0) call psb_ensure_size(nzin_+nz,val,info)
      if (info /= 0) return


      if (present(iren)) then 
        k=0
        do i=irw,lrw
          do j=a%ia2(i),a%ia2(i+1)-1
            k            = k + 1
            val(nzin_+k) = a%aspk(j)
            ia(nzin_+k)  = iren(i)
            ja(nzin_+k)  = iren(a%ia1(j))
          end do
        enddo
      else
        k=0

        do i=irw,lrw
          do j=a%ia2(i),a%ia2(i+1)-1
            k            = k + 1
            ia(nzin_+k)  = i
            ja(nzin_+k)  = a%ia1(j)
            val(nzin_+k) = a%aspk(j)
          end do
        enddo
      end if
      if (a%pr(1) /= 0) then
        write(debug_unit,*)&
             & 'Feeling lazy today, Right Permutation will have to wait'
      endif

    endif

  end subroutine csr_cspgtrow


  subroutine coo_cspgtrow(irw,a,nz,ia,ja,val,nzin,append,lrw,info,iren)

    use psb_sort_mod
    use psb_spmat_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    type(psb_cspmat_type), intent(in)    :: a
    integer                              :: irw
    integer, intent(out)                 :: nz
    integer, allocatable, intent(inout)  :: ia(:), ja(:)
    complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
    integer                              :: nzin
    logical, intent(in)                  :: append
    integer                              :: lrw,info
    integer, optional                    :: iren(:)
    integer  :: nzin_, nza, idx,ip,jp,i,k, nzt
    integer           :: debug_level, debug_unit
    character(len=20) :: name='coo_zspgtrow'

    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    nza = a%infoa(psb_nnz_)
    if (a%pl(1) /= 0) then
      write(debug_unit,*) trim(name),&
           & 'Fatal error: do not feed a permuted mat so far!'
      idx = -1 
    else
      idx = irw
    endif
    if (idx<0) then 
      write(debug_unit,*)  trim(name),&
           &' Error : idx no good ',idx
      info = 2
      return
    end if

    if (append) then 
      nzin_ = nzin
    else
      nzin_ = 0
    endif

    if (a%infoa(psb_srtd_) == psb_isrtdcoo_) then 
      ! In this case we can do a binary search. 
      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*)  trim(name),': srtdcoo '
      do
        ip = psb_ibsrch(irw,nza,a%ia1)
        if (ip /= -1) exit
        irw = irw + 1
        if (irw > lrw) then
          write(debug_unit,*) trim(name),&
               & 'Warning : did not find any rows. Is this an error? ',&
               & irw,lrw,idx
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
        jp = psb_ibsrch(lrw,nza,a%ia1)
        if (jp /= -1) exit
        lrw = lrw - 1
        if (irw > lrw) then
          write(debug_unit,*)  trim(name),&
               & 'Warning : did not find any rows. Is this an error?'
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
      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*)  trim(name),': ip jp',ip,jp,nza
      if ((ip /= -1) .and.(jp /= -1)) then 
        ! Now do the copy.
        nz = jp - ip +1 

        call psb_ensure_size(nzin_+nz,ia,info)
        if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
        if (info==0) call psb_ensure_size(nzin_+nz,val,info)
        if (info /= 0) return

        if (present(iren)) then 
          do i=ip,jp
            nzin_ = nzin_ + 1
            val(nzin_) = a%aspk(i)
            ia(nzin_)  = iren(a%ia1(i))
            ja(nzin_)  = iren(a%ia2(i))
          enddo
        else
          do i=ip,jp
            nzin_ = nzin_ + 1
            val(nzin_) = a%aspk(i)
            ia(nzin_)  = a%ia1(i)
            ja(nzin_)  = a%ia2(i)
          enddo
        end if
      else 
        nz = 0 
      end if

    else
      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*) trim(name),': unsorted '
      nzt = (nza*(lrw-irw+1))/max(a%m,1)
      
      call psb_ensure_size(nzin_+nzt,ia,info)
      if (info==0) call psb_ensure_size(nzin_+nzt,ja,info)
      if (info==0) call psb_ensure_size(nzin_+nzt,val,info)
      if (info /= 0) return
      
      if (present(iren)) then 
        k = 0 
        do i=1, a%infoa(psb_nnz_)
          if ((a%ia1(i)>=irw).and.(a%ia1(i)<=lrw)) then 
            k = k + 1 
            if (k > nzt) then
              nzt = k 
              call psb_ensure_size(nzin_+nzt,ia,info)
              if (info==0) call psb_ensure_size(nzin_+nzt,ja,info)
              if (info==0) call psb_ensure_size(nzin_+nzt,val,info)
              if (info /= 0) return
            end if
            val(nzin_+k) = a%aspk(i)
            ia(nzin_+k)  = iren(a%ia1(i))
            ja(nzin_+k)  = iren(a%ia2(i))
          endif
        enddo
      else
        k = 0 
        do i=1,a%infoa(psb_nnz_)
          if ((a%ia1(i)>=irw).and.(a%ia1(i)<=lrw)) then 
            k = k + 1 
            if (k > nzt) then
              nzt = k 
              call psb_ensure_size(nzin_+nzt,ia,info)
              if (info==0) call psb_ensure_size(nzin_+nzt,ja,info)
              if (info==0) call psb_ensure_size(nzin_+nzt,val,info)
              if (info /= 0) return

            end if
            val(nzin_+k) = a%aspk(i)
            ia(nzin_+k)  = (a%ia1(i))
            ja(nzin_+k)  = (a%ia2(i))
          endif
        enddo
        nzin_=nzin_+k
      end if
      nz = k 
    end if

  end subroutine coo_cspgtrow


  subroutine jad_cspgtrow(irw,a,nz,ia,ja,val,nzin,append,lrw,info,iren)

    use psb_sort_mod
    use psb_spmat_type
    use psb_const_mod
    implicit none

    type(psb_cspmat_type), intent(in), target :: a
    integer                               :: irw
    integer, intent(out)                 :: nz
    integer, allocatable, intent(inout)  :: ia(:), ja(:)
    complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
    integer                               :: nzin
    logical, intent(in)                   :: append
    integer, optional                     :: iren(:)
    integer                               :: lrw,info

    integer, pointer                      :: ia1(:), ia2(:), ia3(:),&
         & ja_(:), ka_(:), indices(:), blks(:)
    integer  :: png, pia, pja, ipx, blk, rb, row, k_pt, npg, col, ng, nzin_,&
         & i,j,k,nr


    png = a%ia2(1) ! points to the number of blocks
    pia = a%ia2(2) ! points to the beginning of ia(3,png)
    pja = a%ia2(3) ! points to the beginning of ja(:)

    ng  =  a%ia2(png)              ! the number of blocks
    ja_  => a%ia2(pja:)             ! the array containing the pointers to ka and aspk
    ka_  => a%ia1(:)                ! the array containing the column indices
    ia1 => a%ia2(pia:pja-1:3)      ! the array containing the first row index of each block
    ia2 => a%ia2(pia+1:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first jad column
    ia3 => a%ia2(pia+2:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first csr column

    if (append) then 
      nzin_ = nzin
    else
      nzin_ = 0
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
          if(ia1(j) == indices(i)) then
            blks(i)=j
            nz=nz+ia3(j)-ia2(j)
            ipx = ia1(j)         ! the first row index of the block
            rb  = indices(i)-ipx   ! the row offset within the block
            row = ia3(j)+rb
            nz  = nz+ja_(row+1)-ja_(row)
            exit blkfnd
          else if(ia1(j) > indices(i)) then
            blks(i)=j-1
            nz=nz+ia3(j-1)-ia2(j-1)
            ipx = ia1(j-1)         ! the first row index of the block
            rb  = indices(i)-ipx   ! the row offset within the block
            row = ia3(j-1)+rb
            nz  = nz+ja_(row+1)-ja_(row)
            exit blkfnd
          end if
        end do blkfnd
      end do


      call psb_ensure_size(nzin_+nz,ia,info)
      if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
      if (info==0) call psb_ensure_size(nzin_+nz,val,info)
      if (info /= 0) return

      k=0
      ! cycle over rows
      do i=1,nr

        ! find which block the row belongs to
        blk = blks(i)

        ! extract first part of the row from the jad block
        ipx = ia1(blk)             ! the first row index of the block
        k_pt= ia2(blk)             ! the pointer to the beginning of a column in ja
        rb  = indices(i)-ipx       ! the row offset within the block
        npg = ja_(k_pt+1)-ja_(k_pt)  ! the number of rows in the block

        if(present(iren))then
          do  col = ia2(blk), ia3(blk)-1 
            k=k+1
            val(nzin_+k) = a%aspk(ja_(col)+rb)
            ia(nzin_+k)  = iren(irw+i-1)
            ja(nzin_+k)  = iren(ka_(ja_(col)+rb))
          end do
        else
          do  col = ia2(blk), ia3(blk)-1 
            k=k+1
            val(nzin_+k) = a%aspk(ja_(col)+rb)
            ia(nzin_+k)  = irw+i-1
            ja(nzin_+k)  = ka_(ja_(col)+rb)
          end do
        end if
        ! extract second part of the row from the csr tail
        row=ia3(blk)+rb
        if(present(iren))then
          do j=ja_(row), ja_(row+1)-1
            k=k+1
            val(nzin_+k) = a%aspk(j)
            ia(nzin_+k)  = iren(irw+i-1)
            ja(nzin_+k)  = iren(ka_(j))
          end do
        else
          do j=ja_(row), ja_(row+1)-1
            k=k+1
            val(nzin_+k) = a%aspk(j)
            ia(nzin_+k)  = irw+i-1
            ja(nzin_+k)  = ka_(j)
          end do
        end if
      end do

    else
      ! There might be some problems
      info=134
    end if

  end subroutine jad_cspgtrow

  
  subroutine csr_zspgtrow(irw,a,nz,ia,ja,val,nzin,append,lrw,info,iren)

    use psb_sort_mod
    use psb_spmat_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    type(psb_zspmat_type), intent(in)    :: a
    integer                              :: irw
    integer, intent(out)                 :: nz
    integer, allocatable, intent(inout)  :: ia(:), ja(:)
    complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
    integer                              :: nzin
    logical, intent(in)                  :: append
    integer                              :: lrw,info
    integer, optional                    :: iren(:)

    integer :: idx,i,j, k, nr, row_idx, nzin_
    integer, allocatable  :: indices(:)
    integer             :: debug_level, debug_unit

    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    if (append) then 
      nzin_ = nzin
    else
      nzin_ = 0
    endif

    if (a%pl(1) /= 0) then

      nr = lrw - irw + 1 
      allocate(indices(nr))
      nz = 0
      do i=1,nr
        indices(i)=a%pl(irw+i-1)
        nz=nz+a%ia2(indices(i)+1)-a%ia2(indices(i))
      end do
      
      call psb_ensure_size(nzin_+nz,ia,info)
      if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
      if (info==0) call psb_ensure_size(nzin_+nz,val,info)
      if (info /= 0) return

      k=0
      if(present(iren)) then
        do i=1,nr
          row_idx=indices(i)
          do j=a%ia2(row_idx),a%ia2(row_idx+1)-1
            k            = k + 1
            val(nzin_+k) = a%aspk(j)
            ia(nzin_+k)  = iren(row_idx)
            ja(nzin_+k)  = iren(a%ia1(j))
          end do
        end do
      else
        do i=1,nr
          row_idx=indices(i)
          do j=a%ia2(row_idx),a%ia2(row_idx+1)-1
            k            = k + 1
            val(nzin_+k) = a%aspk(j)
            ia(nzin_+k)  = row_idx
            ja(nzin_+k)  = a%ia1(j)
          end do
        end do
      end if

    else
      idx = irw

      if (idx<0) then 
        write(0,*) ' spgtrow Error : idx no good ',idx
        return
      end if
      nr = lrw - irw + 1 
      nz = a%ia2(idx+nr) - a%ia2(idx)

      call psb_ensure_size(nzin_+nz,ia,info)
      if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
      if (info==0) call psb_ensure_size(nzin_+nz,val,info)
      if (info /= 0) return


      if (present(iren)) then 
        k=0
        do i=irw,lrw
          do j=a%ia2(i),a%ia2(i+1)-1
            k            = k + 1
            val(nzin_+k) = a%aspk(j)
            ia(nzin_+k)  = iren(i)
            ja(nzin_+k)  = iren(a%ia1(j))
          end do
        enddo
      else
        k=0

        do i=irw,lrw
          do j=a%ia2(i),a%ia2(i+1)-1
            k            = k + 1
            ia(nzin_+k)  = i
            ja(nzin_+k)  = a%ia1(j)
            val(nzin_+k) = a%aspk(j)
          end do
        enddo
      end if
      if (a%pr(1) /= 0) then
        write(debug_unit,*)&
             & 'Feeling lazy today, Right Permutation will have to wait'
      endif

    endif

  end subroutine csr_zspgtrow


  subroutine coo_zspgtrow(irw,a,nz,ia,ja,val,nzin,append,lrw,info,iren)

    use psb_sort_mod
    use psb_spmat_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    type(psb_zspmat_type), intent(in)    :: a
    integer                              :: irw
    integer, intent(out)                 :: nz
    integer, allocatable, intent(inout)  :: ia(:), ja(:)
    complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
    integer                              :: nzin
    logical, intent(in)                  :: append
    integer                              :: lrw,info
    integer, optional                    :: iren(:)
    integer  :: nzin_, nza, idx,ip,jp,i,k, nzt
    integer           :: debug_level, debug_unit
    character(len=20) :: name='coo_zspgtrow'

    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    nza = a%infoa(psb_nnz_)
    if (a%pl(1) /= 0) then
      write(debug_unit,*) trim(name),&
           & 'Fatal error: do not feed a permuted mat so far!'
      idx = -1 
    else
      idx = irw
    endif
    if (idx<0) then 
      write(debug_unit,*)  trim(name),&
           &' Error : idx no good ',idx
      info = 2
      return
    end if

    if (append) then 
      nzin_ = nzin
    else
      nzin_ = 0
    endif

    if (a%infoa(psb_srtd_) == psb_isrtdcoo_) then 
      ! In this case we can do a binary search. 
      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*)  trim(name),': srtdcoo '
      do
        ip = psb_ibsrch(irw,nza,a%ia1)
        if (ip /= -1) exit
        irw = irw + 1
        if (irw > lrw) then
          write(debug_unit,*) trim(name),&
               & 'Warning : did not find any rows. Is this an error? ',&
               & irw,lrw,idx
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
        jp = psb_ibsrch(lrw,nza,a%ia1)
        if (jp /= -1) exit
        lrw = lrw - 1
        if (irw > lrw) then
          write(debug_unit,*)  trim(name),&
               & 'Warning : did not find any rows. Is this an error?'
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
      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*)  trim(name),': ip jp',ip,jp,nza
      if ((ip /= -1) .and.(jp /= -1)) then 
        ! Now do the copy.
        nz = jp - ip +1 

        call psb_ensure_size(nzin_+nz,ia,info)
        if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
        if (info==0) call psb_ensure_size(nzin_+nz,val,info)
        if (info /= 0) return

        if (present(iren)) then 
          do i=ip,jp
            nzin_ = nzin_ + 1
            val(nzin_) = a%aspk(i)
            ia(nzin_)  = iren(a%ia1(i))
            ja(nzin_)  = iren(a%ia2(i))
          enddo
        else
          do i=ip,jp
            nzin_ = nzin_ + 1
            val(nzin_) = a%aspk(i)
            ia(nzin_)  = a%ia1(i)
            ja(nzin_)  = a%ia2(i)
          enddo
        end if
      else 
        nz = 0 
      end if

    else
      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*) trim(name),': unsorted '
      nzt = (nza*(lrw-irw+1))/max(a%m,1)
      
      call psb_ensure_size(nzin_+nzt,ia,info)
      if (info==0) call psb_ensure_size(nzin_+nzt,ja,info)
      if (info==0) call psb_ensure_size(nzin_+nzt,val,info)
      if (info /= 0) return
      
      if (present(iren)) then 
        k = 0 
        do i=1, a%infoa(psb_nnz_)
          if ((a%ia1(i)>=irw).and.(a%ia1(i)<=lrw)) then 
            k = k + 1 
            if (k > nzt) then
              nzt = k 
              call psb_ensure_size(nzin_+nzt,ia,info)
              if (info==0) call psb_ensure_size(nzin_+nzt,ja,info)
              if (info==0) call psb_ensure_size(nzin_+nzt,val,info)
              if (info /= 0) return
            end if
            val(nzin_+k) = a%aspk(i)
            ia(nzin_+k)  = iren(a%ia1(i))
            ja(nzin_+k)  = iren(a%ia2(i))
          endif
        enddo
      else
        k = 0 
        do i=1,a%infoa(psb_nnz_)
          if ((a%ia1(i)>=irw).and.(a%ia1(i)<=lrw)) then 
            k = k + 1 
            if (k > nzt) then
              nzt = k 
              call psb_ensure_size(nzin_+nzt,ia,info)
              if (info==0) call psb_ensure_size(nzin_+nzt,ja,info)
              if (info==0) call psb_ensure_size(nzin_+nzt,val,info)
              if (info /= 0) return

            end if
            val(nzin_+k) = a%aspk(i)
            ia(nzin_+k)  = (a%ia1(i))
            ja(nzin_+k)  = (a%ia2(i))
          endif
        enddo
        nzin_=nzin_+k
      end if
      nz = k 
    end if

  end subroutine coo_zspgtrow


  subroutine jad_zspgtrow(irw,a,nz,ia,ja,val,nzin,append,lrw,info,iren)

    use psb_sort_mod
    use psb_spmat_type
    use psb_const_mod
    implicit none

    type(psb_zspmat_type), intent(in), target :: a
    integer                               :: irw
    integer, intent(out)                 :: nz
    integer, allocatable, intent(inout)  :: ia(:), ja(:)
    complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
    integer                               :: nzin
    logical, intent(in)                   :: append
    integer, optional                     :: iren(:)
    integer                               :: lrw,info

    integer, pointer                      :: ia1(:), ia2(:), ia3(:),&
         & ja_(:), ka_(:), indices(:), blks(:)
    integer  :: png, pia, pja, ipx, blk, rb, row, k_pt, npg, col, ng, nzin_,&
         & i,j,k,nr


    png = a%ia2(1) ! points to the number of blocks
    pia = a%ia2(2) ! points to the beginning of ia(3,png)
    pja = a%ia2(3) ! points to the beginning of ja(:)

    ng  =  a%ia2(png)              ! the number of blocks
    ja_  => a%ia2(pja:)             ! the array containing the pointers to ka and aspk
    ka_  => a%ia1(:)                ! the array containing the column indices
    ia1 => a%ia2(pia:pja-1:3)      ! the array containing the first row index of each block
    ia2 => a%ia2(pia+1:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first jad column
    ia3 => a%ia2(pia+2:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first csr column

    if (append) then 
      nzin_ = nzin
    else
      nzin_ = 0
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
          if(ia1(j) == indices(i)) then
            blks(i)=j
            nz=nz+ia3(j)-ia2(j)
            ipx = ia1(j)         ! the first row index of the block
            rb  = indices(i)-ipx   ! the row offset within the block
            row = ia3(j)+rb
            nz  = nz+ja_(row+1)-ja_(row)
            exit blkfnd
          else if(ia1(j) > indices(i)) then
            blks(i)=j-1
            nz=nz+ia3(j-1)-ia2(j-1)
            ipx = ia1(j-1)         ! the first row index of the block
            rb  = indices(i)-ipx   ! the row offset within the block
            row = ia3(j-1)+rb
            nz  = nz+ja_(row+1)-ja_(row)
            exit blkfnd
          end if
        end do blkfnd
      end do


      call psb_ensure_size(nzin_+nz,ia,info)
      if (info==0) call psb_ensure_size(nzin_+nz,ja,info)
      if (info==0) call psb_ensure_size(nzin_+nz,val,info)
      if (info /= 0) return

      k=0
      ! cycle over rows
      do i=1,nr

        ! find which block the row belongs to
        blk = blks(i)

        ! extract first part of the row from the jad block
        ipx = ia1(blk)             ! the first row index of the block
        k_pt= ia2(blk)             ! the pointer to the beginning of a column in ja
        rb  = indices(i)-ipx       ! the row offset within the block
        npg = ja_(k_pt+1)-ja_(k_pt)  ! the number of rows in the block

        if(present(iren))then
          do  col = ia2(blk), ia3(blk)-1 
            k=k+1
            val(nzin_+k) = a%aspk(ja_(col)+rb)
            ia(nzin_+k)  = iren(irw+i-1)
            ja(nzin_+k)  = iren(ka_(ja_(col)+rb))
          end do
        else
          do  col = ia2(blk), ia3(blk)-1 
            k=k+1
            val(nzin_+k) = a%aspk(ja_(col)+rb)
            ia(nzin_+k)  = irw+i-1
            ja(nzin_+k)  = ka_(ja_(col)+rb)
          end do
        end if
        ! extract second part of the row from the csr tail
        row=ia3(blk)+rb
        if(present(iren))then
          do j=ja_(row), ja_(row+1)-1
            k=k+1
            val(nzin_+k) = a%aspk(j)
            ia(nzin_+k)  = iren(irw+i-1)
            ja(nzin_+k)  = iren(ka_(j))
          end do
        else
          do j=ja_(row), ja_(row+1)-1
            k=k+1
            val(nzin_+k) = a%aspk(j)
            ia(nzin_+k)  = irw+i-1
            ja(nzin_+k)  = ka_(j)
          end do
        end if
      end do

    else
      ! There might be some problems
      info=134
    end if

  end subroutine jad_zspgtrow


end module psb_getrow_mod
