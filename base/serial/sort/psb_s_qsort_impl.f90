!!$ 
!!$              Parallel Sparse BLAS  version 3.4
!!$    (C) Copyright 2006, 2010, 2015
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
!
!  The  quicksort routines 
!  References:
!  D. Knuth
!  The Art of Computer Programming, vol. 3
!  Addison-Wesley
!  
!  Aho, Hopcroft, Ullman
!  Data Structures and Algorithms
!  Addison-Wesley
!
subroutine psb_sqsort(x,ix,dir,flag)
  use psb_s_sort_mod, psb_protect_name => psb_sqsort
  use psb_error_mod
  implicit none 
  real(psb_spk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, err_act, i

  integer(psb_ipk_)  :: ierr(5)
  character(len=20)  :: name

  name='psb_sqsort'
  call psb_erractionsave(err_act)

  if (present(flag)) then 
    flag_ = flag
  else 
    flag_ = psb_sort_ovw_idx_
  end if
  select case(flag_) 
  case( psb_sort_ovw_idx_, psb_sort_keep_idx_)
    ! OK keep going
  case default
    ierr(1) = 4; ierr(2) = flag_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
    goto 9999
  end select

  if (present(dir)) then 
    dir_ = dir
  else
    dir_= psb_sort_up_
  end if

  n = size(x)

  if (present(ix)) then
    if (size(ix) < n) then 
      ierr(1) = 2; ierr(2) = size(ix); 
      call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
      goto 9999
    end if
    if (flag_==psb_sort_ovw_idx_) then
      do i=1,n
        ix(i) = i
      end do
    end if

    select case(dir_) 
    case (psb_sort_up_)
      call psi_sqsrx_up(n,x,ix)
    case (psb_sort_down_)
      call psi_sqsrx_dw(n,x,ix)
    case (psb_asort_up_)
        call psi_saqsrx_up(n,x,ix)
    case (psb_asort_down_)
        call psi_saqsrx_dw(n,x,ix)
    case default
      ierr(1) = 3; ierr(2) = dir_; 
      call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
      goto 9999
    end select
  else 
    select case(dir_) 
    case (psb_sort_up_)
      call psi_sqsr_up(n,x)
    case (psb_sort_down_)
      call psi_sqsr_dw(n,x)
    case (psb_asort_up_)
        call psi_saqsr_up(n,x)
    case (psb_asort_down_)
        call psi_saqsr_dw(n,x)
    case default
      ierr(1) = 3; ierr(2) = dir_; 
      call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
      goto 9999
    end select

  end if

  return

9999 call psb_error_handler(err_act)

  return
end subroutine psb_sqsort

subroutine psi_sqsrx_up(n,x,idx)
  use psb_s_sort_mod, psb_protect_name => psi_sqsrx_up
  use psb_error_mod
  implicit none 

  real(psb_spk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), intent(inout) :: idx(:)
  integer(psb_ipk_), intent(in)   :: n
  !     .. Local Scalars ..
  real(psb_spk_) :: piv, xk, xt
  integer(psb_ipk_) :: i, j, ilx, iux, istp, lpiv
  integer(psb_ipk_) :: ixt, n1, n2
  integer(psb_ipk_), parameter :: maxstack=64,nparms=3,ithrs=72
  integer(psb_ipk_) :: istack(nparms,maxstack)

  if (n > ithrs) then          
    !
    !     Init stack pointer
    !
    istp = 1
    istack(1,istp) = 1
    istack(2,istp) = n

    do 
      if (istp <= 0) exit
      ilx  = istack(1,istp)
      iux  = istack(2,istp)
      istp = istp - 1
      !
      !       Choose a pivot with median-of-three heuristics, leave it 
      !       in the LPIV location
      !            
      i = ilx
      j = iux 
      lpiv = (i+j)/2
      piv  = x(lpiv)
      if (piv < x(i)) then
        xt        = x(i)
        ixt       = idx(i)
        x(i)      = x(lpiv)
        idx(i)    = idx(lpiv)
        x(lpiv)   = xt
        idx(lpiv) = ixt
        piv       = x(lpiv)
      endif
      if (piv > x(j)) then
        xt        = x(j)
        ixt       = idx(j)
        x(j)      = x(lpiv)
        idx(j)    = idx(lpiv)
        x(lpiv)   = xt
        idx(lpiv) = ixt
        piv       = x(lpiv)
      endif
      if (piv < x(i)) then
        xt        = x(i)
        ixt       = idx(i)
        x(i)      = x(lpiv)
        idx(i)    = idx(lpiv)
        x(lpiv)   = xt
        idx(lpiv) = ixt
        piv       = x(lpiv)
      endif
      !
      !     now piv is correct;  place it into first location
      xt        = x(i)
      ixt       = idx(i)
      x(i)      = x(lpiv)
      idx(i)    = idx(lpiv)
      x(lpiv)   = xt
      idx(lpiv) = ixt
      piv       = x(lpiv)

      i = ilx - 1 
      j = iux + 1 

      outer_up: do
        in_up1: do
          i = i + 1
          xk = x(i)
          if (xk >= piv) exit in_up1
        end do in_up1
        !
        !     Ensure finite termination for next loop
        !
        xt  = xk
        x(i) = piv
        in_up2:do 
          j = j - 1
          xk = x(j)
          if (xk <= piv) exit in_up2
        end do in_up2
        x(i) = xt  

        if (j > i) then
          xt     = x(i)
          ixt    = idx(i)
          x(i)   = x(j)
          idx(i) = idx(j)
          x(j)   = xt 
          idx(j) = ixt  
        else
          exit outer_up
        end if
      end do outer_up
      if (i == ilx) then 
        if (x(i) /= piv) then
          call psb_errpush(psb_err_internal_error_,&
               & r_name='psi_sqsrx',a_err='impossible pivot condition')
          call psb_error()
        endif
        i = i + 1 
      endif

      n1 = (i-1)-ilx+1
      n2 = iux-(i)+1
      if (n1 > n2) then
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_sisrx_up(n1,x(ilx:i-1),idx(ilx:i-1))
        endif
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_sisrx_up(n2,x(i:iux),idx(i:iux))
        endif
      else
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_sisrx_up(n2,x(i:iux),idx(i:iux))
        endif
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_sisrx_up(n1,x(ilx:i-1),idx(ilx:i-1))
        endif
      endif
    enddo
  else
    call psi_sisrx_up(n,x,idx)
  endif
end subroutine psi_sqsrx_up

subroutine psi_sqsrx_dw(n,x,idx)
  use psb_s_sort_mod, psb_protect_name => psi_sqsrx_dw
  use psb_error_mod
  implicit none 

  real(psb_spk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), intent(inout) :: idx(:)
  integer(psb_ipk_), intent(in)   :: n
  !     .. Local Scalars ..
  real(psb_spk_) :: piv, xk, xt
  integer(psb_ipk_) :: i, j, ilx, iux, istp, lpiv
  integer(psb_ipk_) :: ixt, n1, n2

  integer(psb_ipk_), parameter :: maxstack=64,nparms=3,ithrs=72
  integer(psb_ipk_) :: istack(nparms,maxstack)

  if (n > ithrs) then          
    !
    !     Init stack pointer
    !
    istp = 1
    istack(1,istp) = 1
    istack(2,istp) = n

    do 
      if (istp <= 0) exit
      ilx  = istack(1,istp)
      iux  = istack(2,istp)
      istp = istp - 1
      !
      !       Choose a pivot with median-of-three heuristics, leave it 
      !       in the LPIV location
      !            
      i = ilx
      j = iux 
      lpiv = (i+j)/2
      piv  = x(lpiv)
      if (piv > x(i)) then
        xt        = x(i)
        ixt       = idx(i)
        x(i)      = x(lpiv)
        idx(i)    = idx(lpiv)
        x(lpiv)   = xt
        idx(lpiv) = ixt
        piv       = x(lpiv)
      endif
      if (piv < x(j)) then
        xt        = x(j)
        ixt       = idx(j)
        x(j)      = x(lpiv)
        idx(j)    = idx(lpiv)
        x(lpiv)   = xt
        idx(lpiv) = ixt
        piv       = x(lpiv)
      endif
      if (piv > x(i)) then
        xt        = x(i)
        ixt       = idx(i)
        x(i)      = x(lpiv)
        idx(i)    = idx(lpiv)
        x(lpiv)   = xt
        idx(lpiv) = ixt
        piv       = x(lpiv)
      endif
      !
      !     now piv is correct;  place it into first location
      xt        = x(i)
      ixt       = idx(i)
      x(i)      = x(lpiv)
      idx(i)    = idx(lpiv)
      x(lpiv)   = xt
      idx(lpiv) = ixt
      piv       = x(lpiv)

      i = ilx - 1 
      j = iux + 1 

      outer_dw: do
        in_dw1: do
          i = i + 1
          xk = x(i)
          if (xk <= piv) exit in_dw1
        end do in_dw1
        !
        !     Ensure finite termination for next loop
        !
        xt  = xk
        x(i) = piv
        in_dw2:do 
          j = j - 1
          xk = x(j)
          if (xk >= piv) exit in_dw2
        end do in_dw2
        x(i) = xt  

        if (j > i) then
          xt     = x(i)
          ixt    = idx(i)
          x(i)   = x(j)
          idx(i) = idx(j)
          x(j)   = xt  
          idx(j) = ixt  
        else
          exit outer_dw
        end if
      end do outer_dw
      if (i == ilx) then 
        if (x(i) /= piv) then
          call psb_errpush(psb_err_internal_error_,& 
               & r_name='psi_sqsrx',a_err='impossible pivot condition')
          call psb_error()
        endif
        i = i + 1 
      endif

      n1 = (i-1)-ilx+1
      n2 = iux-(i)+1
      if (n1 > n2) then
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_sisrx_dw(n1,x(ilx:i-1),idx(ilx:i-1))
        endif
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_sisrx_dw(n2,x(i:iux),idx(i:iux))
        endif
      else
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_sisrx_dw(n2,x(i:iux),idx(i:iux))
        endif
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_sisrx_dw(n1,x(ilx:i-1),idx(ilx:i-1))
        endif
      endif
    enddo
  else
    call psi_sisrx_dw(n,x,idx)
  endif

end subroutine psi_sqsrx_dw

subroutine psi_sqsr_up(n,x)
  use psb_s_sort_mod, psb_protect_name => psi_sqsr_up
  use psb_error_mod
  implicit none 

  real(psb_spk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), intent(in)   :: n
  !     ..
  !     .. Local Scalars ..
  real(psb_spk_) :: piv, xt, xk
  integer(psb_ipk_) :: i, j, ilx, iux, istp, lpiv
  integer(psb_ipk_) :: n1, n2

  integer(psb_ipk_), parameter :: maxstack=64,nparms=3,ithrs=72
  integer(psb_ipk_) :: istack(nparms,maxstack)


  if (n > ithrs) then          
    !
    !     Init stack pointer
    !
    istp = 1
    istack(1,istp) = 1
    istack(2,istp) = n

    do 
      if (istp <= 0) exit
      ilx  = istack(1,istp)
      iux  = istack(2,istp)
      istp = istp - 1
      !
      !       Choose a pivot with median-of-three heuristics, leave it 
      !       in the LPIV location
      !            
      i = ilx
      j = iux 
      lpiv = (i+j)/2
      piv  = x(lpiv)
      if (piv < x(i)) then
        xt = x(i)
        x(i) = x(lpiv)
        x(lpiv) = xt
        piv = x(lpiv)
      endif
      if (piv > x(j)) then
        xt = x(j)
        x(j) = x(lpiv)
        x(lpiv) = xt
        piv = x(lpiv)
      endif
      if (piv < x(i)) then
        xt = x(i)
        x(i) = x(lpiv)
        x(lpiv) = xt
        piv = x(lpiv)
      endif
      !
      !     now piv is correct;  place it into first location

      xt = x(i)
      x(i) = x(lpiv)
      x(lpiv) = xt

      i = ilx - 1 
      j = iux + 1 

      outer_up: do
        in_up1: do
          i = i + 1
          xk = x(i)
          if (xk >= piv) exit in_up1
        end do in_up1
        !
        !     Ensure finite termination for next loop
        !
        xt  = xk
        x(i) = piv
        in_up2:do 
          j = j - 1
          xk = x(j)
          if (xk <= piv) exit in_up2
        end do in_up2
        x(i) = xt  

        if (j > i) then
          xt  = x(i)
          x(i) = x(j)
          x(j) = xt 
        else
          exit outer_up
        end if
      end do outer_up
      if (i == ilx) then 
        if (x(i) /= piv) then
          call psb_errpush(psb_err_internal_error_,&
               & r_name='psi_sqsr',a_err='impossible pivot condition')
          call psb_error()
        endif
        i = i + 1 
      endif

      n1 = (i-1)-ilx+1
      n2 = iux-(i)+1
      if (n1 > n2) then
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_sisr_up(n1,x(ilx:i-1))
        endif
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_sisr_up(n2,x(i:iux))
        endif
      else
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_sisr_up(n2,x(i:iux))
        endif
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_sisr_up(n1,x(ilx:i-1))
        endif
      endif
    enddo
  else
    call psi_sisr_up(n,x)
  endif

end subroutine psi_sqsr_up

subroutine psi_sqsr_dw(n,x)
  use psb_s_sort_mod, psb_protect_name => psi_sqsr_dw
  use psb_error_mod
  implicit none 

  real(psb_spk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), intent(in)   :: n
  !     ..
  !     .. Local Scalars ..
  real(psb_spk_) :: piv, xt, xk
  integer(psb_ipk_) :: i, j, ilx, iux, istp, lpiv
  integer(psb_ipk_) :: n1, n2

  integer(psb_ipk_), parameter :: maxstack=64,nparms=3,ithrs=72
  integer(psb_ipk_) :: istack(nparms,maxstack)


  if (n > ithrs) then          
    !
    !     Init stack pointer
    !
    istp = 1
    istack(1,istp) = 1
    istack(2,istp) = n

    do 
      if (istp <= 0) exit
      ilx  = istack(1,istp)
      iux  = istack(2,istp)
      istp = istp - 1
      !
      !       Choose a pivot with median-of-three heuristics, leave it 
      !       in the LPIV location
      !            
      i = ilx
      j = iux 
      lpiv = (i+j)/2
      piv  = x(lpiv)
      if (piv > x(i)) then
        xt = x(i)
        x(i) = x(lpiv)
        x(lpiv) = xt
        piv = x(lpiv)
      endif
      if (piv < x(j)) then
        xt = x(j)
        x(j) = x(lpiv)
        x(lpiv) = xt
        piv = x(lpiv)
      endif
      if (piv > x(i)) then
        xt = x(i)
        x(i) = x(lpiv)
        x(lpiv) = xt
        piv = x(lpiv)
      endif
      !
      !     now piv is correct;  place it into first location

      xt = x(i)
      x(i) = x(lpiv)
      x(lpiv) = xt

      i = ilx - 1 
      j = iux + 1 

      outer_dw: do
        in_dw1: do
          i = i + 1
          xk = x(i)
          if (xk <= piv) exit in_dw1
        end do in_dw1
        !
        !     Ensure finite termination for next loop
        !
        xt  = xk
        x(i) = piv
        in_dw2:do 
          j = j - 1
          xk = x(j)
          if (xk >= piv) exit in_dw2
        end do in_dw2
        x(i) = xt  

        if (j > i) then
          xt  = x(i)
          x(i) = x(j)
          x(j) = xt  
        else
          exit outer_dw
        end if
      end do outer_dw
      if (i == ilx) then 
        if (x(i) /= piv) then
          call psb_errpush(psb_err_internal_error_, &
               & r_name='psi_sqsr',a_err='impossible pivot condition')
          call psb_error()
        endif
        i = i + 1 
      endif

      n1 = (i-1)-ilx+1
      n2 = iux-(i)+1
      if (n1 > n2) then
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_sisr_dw(n1,x(ilx:i-1))
        endif
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_sisr_dw(n2,x(i:iux))
        endif
      else
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_sisr_dw(n2,x(i:iux))
        endif
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_sisr_dw(n1,x(ilx:i-1))
        endif
      endif
    enddo
  else
    call psi_sisr_dw(n,x)
  endif

end subroutine psi_sqsr_dw

subroutine psi_saqsrx_up(n,x,idx)
  use psb_s_sort_mod, psb_protect_name => psi_saqsrx_up
  use psb_error_mod
  implicit none 

  real(psb_spk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), intent(inout) :: idx(:)
  integer(psb_ipk_), intent(in)   :: n
  !     .. Local Scalars ..
  real(psb_spk_) :: piv, xk
  real(psb_spk_) :: xt
  integer(psb_ipk_) :: i, j, ilx, iux, istp, lpiv
  integer(psb_ipk_) :: ixt, n1, n2

  integer(psb_ipk_), parameter :: maxstack=64,nparms=3,ithrs=72
  integer(psb_ipk_) :: istack(nparms,maxstack)

  if (n > ithrs) then          
    !
    !     Init stack pointer
    !
    istp = 1
    istack(1,istp) = 1
    istack(2,istp) = n

    do 
      if (istp <= 0) exit
      ilx  = istack(1,istp)
      iux  = istack(2,istp)
      istp = istp - 1
      !
      !       Choose a pivot with median-of-three heuristics, leave it 
      !       in the LPIV location
      !            
      i = ilx
      j = iux 
      lpiv = (i+j)/2
      piv  = abs(x(lpiv))
      if (piv < abs(x(i))) then
        xt   = x(i)
        ixt  = idx(i)
        x(i) = x(lpiv)
        idx(i) = idx(lpiv)
        x(lpiv) = xt
        idx(lpiv) = ixt
        piv = abs(x(lpiv))
      endif
      if (piv > abs(x(j))) then
        xt        = x(j)
        ixt       = idx(j)
        x(j)      = x(lpiv)
        idx(j)    = idx(lpiv)
        x(lpiv)   = xt
        idx(lpiv) = ixt
        piv = abs(x(lpiv))
      endif
      if (piv < abs(x(i))) then
        xt        = x(i)
        ixt       = idx(i)
        x(i)      = x(lpiv)
        idx(i)    = idx(lpiv)
        x(lpiv)   = xt
        idx(lpiv) = ixt
        piv = abs(x(lpiv))
      endif
      !
      !     now piv is correct;  place it into first location
      xt        = x(i)
      ixt       = idx(i)
      x(i)      = x(lpiv)
      idx(i)    = idx(lpiv)
      x(lpiv)   = xt
      idx(lpiv) = ixt

      i = ilx - 1 
      j = iux + 1 

      outer_up: do
        in_up1: do
          i = i + 1
          xk = abs(x(i))
          if (xk >= piv) exit in_up1
        end do in_up1
        !
        !     Ensure finite termination for next loop
        !
        xt  = x(i)
        x(i) = piv
        in_up2:do 
          j = j - 1
          xk = abs(x(j))
          if (xk <= piv) exit in_up2
        end do in_up2
        x(i) = xt  

        if (j > i) then
          xt     = x(i)
          ixt    = idx(i)
          x(i)   = x(j)
          idx(i) = idx(j)
          x(j)   = xt 
          idx(j) = ixt  
        else
          exit outer_up
        end if
      end do outer_up
      if (i == ilx) then 
        if (x(i) /= piv) then
          call psb_errpush(psb_err_internal_error_, &
               & r_name='psi_saqsrx',a_err='impossible pivot condition')
          call psb_error()
        endif
        i = i + 1 
      endif

      n1 = (i-1)-ilx+1
      n2 = iux-(i)+1
      if (n1 > n2) then
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_saisrx_up(n1,x(ilx:i-1),idx(ilx:i-1))
        endif
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_saisrx_up(n2,x(i:iux),idx(i:iux))
        endif
      else
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_saisrx_up(n2,x(i:iux),idx(i:iux))
        endif
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_saisrx_up(n1,x(ilx:i-1),idx(ilx:i-1))
        endif
      endif
    enddo
  else
    call psi_saisrx_up(n,x,idx)
  endif


end subroutine psi_saqsrx_up

subroutine psi_saqsrx_dw(n,x,idx)
  use psb_s_sort_mod, psb_protect_name => psi_saqsrx_dw
  use psb_error_mod
  implicit none 

  real(psb_spk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), intent(inout) :: idx(:)
  integer(psb_ipk_), intent(in)   :: n
  !     .. Local Scalars ..
  real(psb_spk_) :: piv, xk
  real(psb_spk_) :: xt
  integer(psb_ipk_) :: i, j, ilx, iux, istp, lpiv
  integer(psb_ipk_) :: ixt, n1, n2

  integer(psb_ipk_), parameter :: maxstack=64,nparms=3,ithrs=72
  integer(psb_ipk_) :: istack(nparms,maxstack)
  if (n > ithrs) then          
    !
    !     Init stack pointer
    !
    istp = 1
    istack(1,istp) = 1
    istack(2,istp) = n

    do 
      if (istp <= 0) exit
      ilx  = istack(1,istp)
      iux  = istack(2,istp)
      istp = istp - 1
      !
      !       Choose a pivot with median-of-three heuristics, leave it 
      !       in the LPIV location
      !            
      i = ilx
      j = iux 
      lpiv = (i+j)/2
      piv  = abs(x(lpiv))
      if (piv > abs(x(i))) then
        xt        = x(i)
        ixt       = idx(i)
        x(i)      = x(lpiv)
        idx(i)    = idx(lpiv)
        x(lpiv)   = xt
        idx(lpiv) = ixt
        piv = abs(x(lpiv))
      endif
      if (piv < abs(x(j))) then
        xt        = x(j)
        ixt       = idx(j)
        x(j)      = x(lpiv)
        idx(j)    = idx(lpiv)
        x(lpiv)   = xt
        idx(lpiv) = ixt
        piv = abs(x(lpiv))
      endif
      if (piv > abs(x(i))) then
        xt        = x(i)
        ixt       = idx(i)
        x(i)      = x(lpiv)
        idx(i)    = idx(lpiv)
        x(lpiv)   = xt
        idx(lpiv) = ixt
        piv = abs(x(lpiv))
      endif
      !
      !     now piv is correct;  place it into first location
      xt        = x(i)
      ixt       = idx(i)
      x(i)      = x(lpiv)
      idx(i)    = idx(lpiv)
      x(lpiv)   = xt
      idx(lpiv) = ixt

      i = ilx - 1 
      j = iux + 1 

      outer_dw: do
        in_dw1: do
          i = i + 1
          xk = abs(x(i))
          if (xk <= piv) exit in_dw1
        end do in_dw1
        !
        !     Ensure finite termination for next loop
        !
        xt  = x(i)
        x(i) = piv
        in_dw2:do 
          j = j - 1
          xk = abs(x(j))
          if (xk >= piv) exit in_dw2
        end do in_dw2
        x(i) = xt  

        if (j > i) then
          xt     = x(i)
          ixt    = idx(i)
          x(i)   = x(j)
          idx(i) = idx(j)
          x(j)   = xt  
          idx(j) = ixt  
        else
          exit outer_dw
        end if
      end do outer_dw
      if (i == ilx) then 
        if (x(i) /= piv) then
          call psb_errpush(psb_err_internal_error_,& 
               & r_name='psi_saqsrx',a_err='impossible pivot condition')
          call psb_error()
        endif
        i = i + 1 
      endif

      n1 = (i-1)-ilx+1
      n2 = iux-(i)+1
      if (n1 > n2) then
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_saisrx_dw(n1,x(ilx:i-1),idx(ilx:i-1))
        endif
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_saisrx_dw(n2,x(i:iux),idx(i:iux))
        endif
      else
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_saisrx_dw(n2,x(i:iux),idx(i:iux))
        endif
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_saisrx_dw(n1,x(ilx:i-1),idx(ilx:i-1))
        endif
      endif
    enddo
  else
    call psi_saisrx_dw(n,x,idx)
  endif

end subroutine psi_saqsrx_dw

subroutine psi_saqsr_up(n,x)
  use psb_s_sort_mod, psb_protect_name => psi_saqsr_up
  use psb_error_mod
  implicit none 

  real(psb_spk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), intent(in)   :: n
  !     .. Local Scalars ..
  real(psb_spk_) :: piv, xk
  real(psb_spk_) :: xt
  integer(psb_ipk_) :: i, j, ilx, iux, istp, lpiv
  integer(psb_ipk_) :: ixt, n1, n2

  integer(psb_ipk_), parameter :: maxstack=64,nparms=3,ithrs=72
  integer(psb_ipk_) :: istack(nparms,maxstack)

  if (n > ithrs) then          
    !
    !     Init stack pointer
    !
    istp = 1
    istack(1,istp) = 1
    istack(2,istp) = n

    do 
      if (istp <= 0) exit
      ilx  = istack(1,istp)
      iux  = istack(2,istp)
      istp = istp - 1
      !
      !       Choose a pivot with median-of-three heuristics, leave it 
      !       in the LPIV location
      !            
      i = ilx
      j = iux 
      lpiv = (i+j)/2
      piv  = abs(x(lpiv))
      if (piv < abs(x(i))) then
        xt = x(i)
        x(i) = x(lpiv)
        x(lpiv) = xt
        piv = abs(x(lpiv))
      endif
      if (piv > abs(x(j))) then
        xt = x(j)
        x(j) = x(lpiv)
        x(lpiv) = xt
        piv = abs(x(lpiv))
      endif
      if (piv < abs(x(i))) then
        xt = x(i)
        x(i) = x(lpiv)
        x(lpiv) = xt
        piv = abs(x(lpiv))
      endif
      !
      !     now piv is correct;  place it into first location

      xt = x(i)
      x(i) = x(lpiv)
      x(lpiv) = xt

      i = ilx - 1 
      j = iux + 1 

      outer_up: do
        in_up1: do
          i = i + 1
          xk = abs(x(i))
          if (xk >= piv) exit in_up1
        end do in_up1
        !
        !     Ensure finite termination for next loop
        !
        xt  = x(i)
        x(i) = piv
        in_up2:do 
          j = j - 1
          xk = abs(x(j))
          if (xk <= piv) exit in_up2
        end do in_up2
        x(i) = xt  

        if (j > i) then
          xt  = x(i)
          x(i) = x(j)
          x(j) = xt 
        else
          exit outer_up
        end if
      end do outer_up
      if (i == ilx) then 
        if (x(i) /= piv) then
          call psb_errpush(psb_err_internal_error_, & 
               & r_name='psi_sqasr',a_err='impossible pivot condition')
          call psb_error()
        endif
        i = i + 1 
      endif

      n1 = (i-1)-ilx+1
      n2 = iux-(i)+1
      if (n1 > n2) then
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_saisr_up(n1,x(ilx:i-1))
        endif
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_saisr_up(n2,x(i:iux))
        endif
      else
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_saisr_up(n2,x(i:iux))
        endif
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_saisr_up(n1,x(ilx:i-1))
        endif
      endif
    enddo
  else
    call psi_saisr_up(n,x)
  endif

end subroutine psi_saqsr_up

subroutine psi_saqsr_dw(n,x)
  use psb_s_sort_mod, psb_protect_name => psi_saqsr_dw
  use psb_error_mod
  implicit none 

  real(psb_spk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), intent(in)   :: n
  !     .. Local Scalars ..
  real(psb_spk_) :: piv, xk
  real(psb_spk_) :: xt
  integer(psb_ipk_) :: i, j, ilx, iux, istp, lpiv
  integer(psb_ipk_) :: ixt, n1, n2

  integer(psb_ipk_), parameter :: maxstack=64,nparms=3,ithrs=72
  integer(psb_ipk_) :: istack(nparms,maxstack)

  if (n > ithrs) then          
    !
    !     Init stack pointer
    !
    istp = 1
    istack(1,istp) = 1
    istack(2,istp) = n

    do 
      if (istp <= 0) exit
      ilx  = istack(1,istp)
      iux  = istack(2,istp)
      istp = istp - 1
      !
      !       Choose a pivot with median-of-three heuristics, leave it 
      !       in the LPIV location
      !            
      i = ilx
      j = iux 
      lpiv = (i+j)/2
      piv  = abs(x(lpiv))
      if (piv > abs(x(i))) then
        xt = x(i)
        x(i) = x(lpiv)
        x(lpiv) = xt
        piv = abs(x(lpiv))
      endif
      if (piv < abs(x(j))) then
        xt = x(j)
        x(j) = x(lpiv)
        x(lpiv) = xt
        piv = abs(x(lpiv))
      endif
      if (piv > abs(x(i))) then
        xt = x(i)
        x(i) = x(lpiv)
        x(lpiv) = xt
        piv = abs(x(lpiv))
      endif
      !
      !     now piv is correct;  place it into first location

      xt = x(i)
      x(i) = x(lpiv)
      x(lpiv) = xt

      i = ilx - 1 
      j = iux + 1 

      outer_dw: do
        in_dw1: do
          i = i + 1
          xk = abs(x(i))
          if (xk <= piv) exit in_dw1
        end do in_dw1
        !
        !     Ensure finite termination for next loop
        !
        xt  = x(i)
        x(i) = piv
        in_dw2:do 
          j = j - 1
          xk = abs(x(j))
          if (xk >= piv) exit in_dw2
        end do in_dw2
        x(i) = xt  

        if (j > i) then
          xt  = x(i)
          x(i) = x(j)
          x(j) = xt  
        else
          exit outer_dw
        end if
      end do outer_dw
      if (i == ilx) then 
        if (x(i) /= piv) then
          call psb_errpush(psb_err_internal_error_,& 
               & r_name='psi_sqasr',a_err='impossible pivot condition')
          call psb_error()
        endif
        i = i + 1 
      endif

      n1 = (i-1)-ilx+1
      n2 = iux-(i)+1
      if (n1 > n2) then
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_saisr_dw(n1,x(ilx:i-1))
        endif
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_saisr_dw(n2,x(i:iux))
        endif
      else
        if (n2 > ithrs) then
          istp = istp + 1
          istack(1,istp) = i
          istack(2,istp) = iux
        else
          call psi_saisr_dw(n2,x(i:iux))
        endif
        if (n1 > ithrs) then 
          istp = istp + 1
          istack(1,istp) = ilx
          istack(2,istp) = i-1
        else
          call psi_saisr_dw(n1,x(ilx:i-1))
        endif
      endif
    enddo
  else
    call psi_saisr_dw(n,x)
  endif

end subroutine psi_saqsr_dw


