!   
!                Parallel Sparse BLAS  version 3.4
!      (C) Copyright 2006, 2010, 2015
!                         Salvatore Filippone    University of Rome Tor Vergata
!                         Alfredo Buttari        CNRS-IRIT, Toulouse
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
!    
  !
  !  The merge-sort routines
  !  References:
  !  D. Knuth
  !  The Art of Computer Programming, vol. 3
  !  Addison-Wesley
  !  
  !  Aho, Hopcroft, Ullman
  !  Data Structures and Algorithms
  !  Addison-Wesley
  !

  subroutine psb_smsort(x,ix,dir,flag)
    use psb_s_sort_mod, psb_protect_name => psb_smsort
    use psb_error_mod
    use psb_ip_reord_mod
    implicit none 
    real(psb_spk_), intent(inout)           :: x(:) 
    integer(psb_ipk_), optional, intent(in)    :: dir, flag
    integer(psb_ipk_), optional, intent(inout) :: ix(:)

    integer(psb_ipk_) :: dir_, flag_, n, err_act

    integer(psb_ipk_), allocatable :: iaux(:)
    integer(psb_ipk_) :: iret, info, i 
    integer(psb_ipk_)  :: ierr(5)
    character(len=20)  :: name

    name='psb_smsort'
    call psb_erractionsave(err_act)

    if (present(dir)) then 
      dir_ = dir
    else
      dir_= psb_sort_up_
    end if
    select case(dir_) 
    case( psb_sort_up_, psb_sort_down_, psb_asort_up_, psb_asort_down_)
      ! OK keep going
    case default
      ierr(1) = 3; ierr(2) = dir_; 
      call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
      goto 9999
    end select

    n = size(x)

    if (present(ix)) then 
      if (size(ix) < n) then 
        ierr(1) = 2; ierr(2) = size(ix); 
        call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
        goto 9999
      end if
      if (present(flag)) then 
        flag_ = flag
      else 
        flag_ = psb_sort_ovw_idx_
      end if
      select case(flag_) 
      case(psb_sort_ovw_idx_)
        do i=1,n 
          ix(i) = i
        end do
      case (psb_sort_keep_idx_)
        ! OK keep going
      case default
        ierr(1) = 4; ierr(2) = flag_; 
        call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
        goto 9999
      end select

    end if

    allocate(iaux(0:n+1),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_alloc_dealloc_,r_name='psb_s_msort')
      goto 9999
    endif

    select case(dir_)
    case (psb_sort_up_)
      call psi_s_msort_up(n,x,iaux,iret)
    case (psb_sort_down_)
      call psi_s_msort_dw(n,x,iaux,iret)
    case (psb_asort_up_)
      call psi_s_amsort_up(n,x,iaux,iret)
    case (psb_asort_down_)
      call psi_s_amsort_dw(n,x,iaux,iret)
    end select
    !
    ! Do the actual reordering, since the inner routines
    ! only provide linked pointers. 
    !
    if (iret == 0 ) then 
      if (present(ix)) then
        call psb_ip_reord(n,x,ix,iaux)
      else      
        call psb_ip_reord(n,x,iaux)
      end if
    end if


    return

9999 call psb_error_handler(err_act)

    return


  end subroutine psb_smsort

  subroutine psi_s_msort_up(n,k,l,iret)
    use psb_const_mod
    implicit none
    integer(psb_ipk_) :: n, iret
    real(psb_spk_)  ::  k(n)
    integer(psb_ipk_) :: l(0:n+1)
    !
    integer(psb_ipk_) :: p,q,s,t
    !     ..
    iret = 0
    !  first step: we are preparing ordered sublists, exploiting
    !  what order was already in the input data; negative links
    !  mark the end of the sublists
    l(0) = 1
    t = n + 1
    do  p = 1,n - 1
      if (k(p) <= k(p+1)) then
        l(p) = p + 1
      else
        l(t) = - (p+1)
        t = p
      end if
    end do
    l(t) = 0
    l(n) = 0
    ! see if the input was already sorted
    if (l(n+1) == 0) then
      iret = 1
      return 
    else
      l(n+1) = abs(l(n+1))
    end if

    mergepass: do 
      ! otherwise, begin a pass through the list.
      ! throughout all the subroutine we have:
      !  p, q: pointing to the sublists being merged
      !  s: pointing to the most recently processed record
      !  t: pointing to the end of previously completed sublist
      s = 0
      t = n + 1
      p = l(s)
      q = l(t)
      if (q == 0) exit mergepass

      outer: do 

        if (k(p) > k(q)) then 

          l(s) = sign(q,l(s))
          s = q
          q = l(q)
          if (q > 0) then 
            do 
              if (k(p) <= k(q)) cycle outer
              s = q
              q = l(q)
              if (q <= 0) exit
            end do
          end if
          l(s) = p
          s = t
          do 
            t = p
            p = l(p)
            if (p <= 0) exit
          end do

        else 

          l(s) = sign(p,l(s))
          s = p
          p = l(p)
          if (p>0) then 
            do 
              if (k(p) > k(q)) cycle outer 
              s = p
              p = l(p)
              if (p <= 0) exit
            end do
          end if
          !  otherwise, one sublist ended, and we append to it the rest
          !  of the other one.
          l(s) = q
          s = t
          do 
            t = q
            q = l(q)
            if (q <= 0) exit
          end do
        end if

        p = -p
        q = -q
        if (q == 0) then
          l(s) = sign(p,l(s))
          l(t) = 0
          exit outer 
        end if
      end do outer
    end do mergepass

  end subroutine psi_s_msort_up

  subroutine psi_s_msort_dw(n,k,l,iret)
    use psb_const_mod
    implicit none
    integer(psb_ipk_) :: n, iret
    real(psb_spk_)   :: k(n)
    integer(psb_ipk_) :: l(0:n+1)
    !
    integer(psb_ipk_) :: p,q,s,t
    !     ..
    iret = 0
    !  first step: we are preparing ordered sublists, exploiting
    !  what order was already in the input data; negative links
    !  mark the end of the sublists
    l(0) = 1
    t = n + 1
    do  p = 1,n - 1
      if (k(p) >= k(p+1)) then
        l(p) = p + 1
      else
        l(t) = - (p+1)
        t = p
      end if
    end do
    l(t) = 0
    l(n) = 0
    ! see if the input was already sorted
    if (l(n+1) == 0) then
      iret = 1
      return 
    else
      l(n+1) = abs(l(n+1))
    end if

    mergepass: do 
      ! otherwise, begin a pass through the list.
      ! throughout all the subroutine we have:
      !  p, q: pointing to the sublists being merged
      !  s: pointing to the most recently processed record
      !  t: pointing to the end of previously completed sublist
      s = 0
      t = n + 1
      p = l(s)
      q = l(t)
      if (q == 0) exit mergepass

      outer: do 

        if (k(p) < k(q)) then 

          l(s) = sign(q,l(s))
          s = q
          q = l(q)
          if (q > 0) then 
            do 
              if (k(p) >= k(q)) cycle outer
              s = q
              q = l(q)
              if (q <= 0) exit
            end do
          end if
          l(s) = p
          s = t
          do 
            t = p
            p = l(p)
            if (p <= 0) exit
          end do

        else 

          l(s) = sign(p,l(s))
          s = p
          p = l(p)
          if (p>0) then 
            do 
              if (k(p) < k(q)) cycle outer 
              s = p
              p = l(p)
              if (p <= 0) exit
            end do
          end if
          !  otherwise, one sublist ended, and we append to it the rest
          !  of the other one.
          l(s) = q
          s = t
          do 
            t = q
            q = l(q)
            if (q <= 0) exit
          end do
        end if

        p = -p
        q = -q
        if (q == 0) then
          l(s) = sign(p,l(s))
          l(t) = 0
          exit outer 
        end if
      end do outer
    end do mergepass

  end subroutine psi_s_msort_dw

  subroutine psi_s_amsort_up(n,k,l,iret)
    use psb_const_mod
    implicit none
    integer(psb_ipk_) :: n, iret
    real(psb_spk_)   :: k(n)
    integer(psb_ipk_) :: l(0:n+1)
    !
    integer(psb_ipk_) :: p,q,s,t
    !     ..
    iret = 0
    !  first step: we are preparing ordered sublists, exploiting
    !  what order was already in the input data; negative links
    !  mark the end of the sublists
    l(0) = 1
    t = n + 1
    do  p = 1,n - 1
      if (abs(k(p)) <= abs(k(p+1))) then
        l(p) = p + 1
      else
        l(t) = - (p+1)
        t = p
      end if
    end do
    l(t) = 0
    l(n) = 0
    ! see if the input was already sorted
    if (l(n+1) == 0) then
      iret = 1
      return 
    else
      l(n+1) = abs(l(n+1))
    end if

    mergepass: do 
      ! otherwise, begin a pass through the list.
      ! throughout all the subroutine we have:
      !  p, q: pointing to the sublists being merged
      !  s: pointing to the most recently processed record
      !  t: pointing to the end of previously completed sublist
      s = 0
      t = n + 1
      p = l(s)
      q = l(t)
      if (q == 0) exit mergepass

      outer: do 

        if (abs(k(p)) > abs(k(q))) then 

          l(s) = sign(q,l(s))
          s = q
          q = l(q)
          if (q > 0) then 
            do 
              if (abs(k(p)) <= abs(k(q))) cycle outer
              s = q
              q = l(q)
              if (q <= 0) exit
            end do
          end if
          l(s) = p
          s = t
          do 
            t = p
            p = l(p)
            if (p <= 0) exit
          end do

        else 

          l(s) = sign(p,l(s))
          s = p
          p = l(p)
          if (p>0) then 
            do 
              if (abs(k(p)) > abs(k(q))) cycle outer 
              s = p
              p = l(p)
              if (p <= 0) exit
            end do
          end if
          !  otherwise, one sublist ended, and we append to it the rest
          !  of the other one.
          l(s) = q
          s = t
          do 
            t = q
            q = l(q)
            if (q <= 0) exit
          end do
        end if

        p = -p
        q = -q
        if (q == 0) then
          l(s) = sign(p,l(s))
          l(t) = 0
          exit outer 
        end if
      end do outer
    end do mergepass

  end subroutine psi_s_amsort_up

  subroutine psi_s_amsort_dw(n,k,l,iret)
    use psb_const_mod
    implicit none
    integer(psb_ipk_) :: n, iret
    real(psb_spk_)   :: k(n)
    integer(psb_ipk_) :: l(0:n+1)
    !
    integer(psb_ipk_) :: p,q,s,t
    !     ..
    iret = 0
    !  first step: we are preparing ordered sublists, exploiting
    !  what order was already in the input data; negative links
    !  mark the end of the sublists
    l(0) = 1
    t = n + 1
    do  p = 1,n - 1
      if (abs(k(p)) >= abs(k(p+1))) then
        l(p) = p + 1
      else
        l(t) = - (p+1)
        t = p
      end if
    end do
    l(t) = 0
    l(n) = 0
    ! see if the input was already sorted
    if (l(n+1) == 0) then
      iret = 1
      return 
    else
      l(n+1) = abs(l(n+1))
    end if

    mergepass: do 
      ! otherwise, begin a pass through the list.
      ! throughout all the subroutine we have:
      !  p, q: pointing to the sublists being merged
      !  s: pointing to the most recently processed record
      !  t: pointing to the end of previously completed sublist
      s = 0
      t = n + 1
      p = l(s)
      q = l(t)
      if (q == 0) exit mergepass

      outer: do 

        if (abs(k(p)) < abs(k(q))) then 

          l(s) = sign(q,l(s))
          s = q
          q = l(q)
          if (q > 0) then 
            do 
              if (abs(k(p)) >= abs(k(q))) cycle outer
              s = q
              q = l(q)
              if (q <= 0) exit
            end do
          end if
          l(s) = p
          s = t
          do 
            t = p
            p = l(p)
            if (p <= 0) exit
          end do

        else 

          l(s) = sign(p,l(s))
          s = p
          p = l(p)
          if (p>0) then 
            do 
              if (abs(k(p)) < abs(k(q))) cycle outer 
              s = p
              p = l(p)
              if (p <= 0) exit
            end do
          end if
          !  otherwise, one sublist ended, and we append to it the rest
          !  of the other one.
          l(s) = q
          s = t
          do 
            t = q
            q = l(q)
            if (q <= 0) exit
          end do
        end if

        p = -p
        q = -q
        if (q == 0) then
          l(s) = sign(p,l(s))
          l(t) = 0
          exit outer 
        end if
      end do outer
    end do mergepass

  end subroutine psi_s_amsort_dw








