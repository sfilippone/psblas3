!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
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
! Original code adapted from:
! == =====================================================================
! Sparse Matrix Multiplication Package
!
! Randolph E. Bank and Craig C. Douglas
!
! na.bank@na-net.ornl.gov and na.cdouglas@na-net.ornl.gov
!
! Compile this with the following command (or a similar one):
!
!     f77 -c -O smmp.f
!
! == =====================================================================
subroutine lsymbmm(n, m, l, ia, ja, diaga, ib, jb, diagb,&
     &  ic, jc, diagc,  index)
  use psb_const_mod
  use psb_realloc_mod
  use psb_sort_mod, only: psb_msort
  !
  integer(psb_lpk_) ::       ia(*), ja(*), diaga, &
       &  ib(*), jb(*), diagb,  diagc,  index(*)
  integer(psb_lpk_), allocatable :: ic(:),jc(:)
  integer(psb_lpk_) :: nze
  integer(psb_ipk_) :: info

  !
  !       symbolic matrix multiply c=a*b
  !
  if (size(ic) < n+1) then 
    write(psb_err_unit,*)&
         &  'Called realloc in SYMBMM '
    call psb_realloc(n+1,ic,info)
    if (info /= psb_success_) then 
      write(psb_err_unit,*)&
           &  'realloc failed in SYMBMM ',info
    end if
  endif
  maxlmn = max(l,m,n)
  do i=1,maxlmn
    index(i)=0
  end do
  if (diagc.eq.0) then
    ic(1)=1
  else
    ic(1)=n+2
  endif
  minlm = min(l,m)
  minmn = min(m,n)
  !
  !    main loop
  !
  do i=1,n
    istart=-1
    length=0
    !
    !    merge row lists
    !
    rowi: do jj=ia(i),ia(i+1)
      !    a = d + ...
      if (jj.eq.ia(i+1)) then
        if (diaga.eq.0 .or. i.gt.minmn) cycle rowi
        j = i
      else
        j=ja(jj)
      endif
      !    b = d + ...
      if (index(j).eq.0 .and. diagb.eq.1 .and. j.le.minlm)then
        index(j)=istart
        istart=j
        length=length+1
      endif
      if ((j<1).or.(j>m)) then 
        write(psb_err_unit,*)&
             &  ' SymbMM: Problem with A ',i,jj,j,m
      endif
      do k=ib(j),ib(j+1)-1 
        if ((jb(k)<1).or.(jb(k)>maxlmn)) then 
          write(psb_err_unit,*)&
               & 'Problem in SYMBMM 1:',j,k,jb(k),maxlmn
        else
          if(index(jb(k)).eq.0) then
            index(jb(k))=istart
            istart=jb(k)
            length=length+1
          endif
        endif
      end do
    end do rowi

    !
    !   row i of jc
    !
    if (diagc.eq.1 .and. index(i).ne.0) length = length - 1
    ic(i+1)=ic(i)+length

    if (ic(i+1) > size(jc)) then 
      if (n > (2*i)) then 
        nze = max(ic(i+1), ic(i)*((n+i-1)/i))
      else
        nze = max(ic(i+1), nint((dble(ic(i))*(dble(n)/i)))   )
      endif
      call psb_realloc(nze,jc,info)
    end if

    do j= ic(i),ic(i+1)-1
      if (diagc.eq.1 .and. istart.eq.i) then
        istart = index(istart)
        index(i) = 0
      endif
      jc(j)=istart
      istart=index(istart)
      index(jc(j))=0
    end do
    call psb_msort(jc(ic(i):ic(i)+length -1))
    index(i) = 0
  end do
  return
end subroutine lsymbmm
! == =====================================================================
! Sparse Matrix Multiplication Package
!
! Randolph E. Bank and Craig C. Douglas
!
! na.bank@na-net.ornl.gov and na.cdouglas@na-net.ornl.gov
!
! Compile this with the following command (or a similar one):
!
!     f77 -c -O smmp.f
!
! == =====================================================================
subroutine lcnumbmm(n, m, l, ia, ja, diaga, a, ib, jb, diagb, b,&
     & ic, jc, diagc, c, temp)
  !
  use psb_const_mod
  integer(psb_lpk_) :: ia(*), ja(*), diaga,&
       &  ib(*), jb(*), diagb,  ic(*), jc(*), diagc 
  !
  complex(psb_spk_) :: a(*), b(*), c(*), temp(*),ajj
  !
  !       numeric matrix multiply c=a*b
  !
  maxlmn = max(l,m,n)
  do i = 1,maxlmn
    temp(i) = 0.
  end do
  minlm = min(l,m)
  minln = min(l,n)
  minmn = min(m,n)
  !
  !   c = a*b
  !
  do i = 1,n
    rowi: do  jj = ia(i),ia(i+1)
      !    a = d + ...
      if (jj.eq.ia(i+1)) then
        if (diaga.eq.0 .or. i.gt.minmn) cycle rowi
        j = i
        ajj = a(i)
      else
        j=ja(jj)
        ajj = a(jj)
      endif
      !    b = d + ...
      if (diagb.eq.1 .and. j.le.minlm) &
           &    temp(j) = temp(j) + ajj * b(j)
      if ((j<1).or.(j>m)) then 
        write(psb_err_unit,*)&
             & ' NUMBMM: Problem with A ',i,jj,j,m
      endif

      do  k = ib(j),ib(j+1)-1
        if((jb(k)<1).or. (jb(k) > maxlmn))  then 
          write(psb_err_unit,*)&
               & ' NUMBMM: jb  problem',j,k,jb(k),maxlmn
        else
          temp(jb(k)) = temp(jb(k)) + ajj * b(k)
        endif
      end do
    end do rowi

    !    c = d + ...
    if (diagc.eq.1 .and. i.le.minln) then
      c(i) = temp(i)
      temp(i) = 0.
    endif
    !$$$        if (mod(i,100) == 1)
    !$$$     +    write(psb_err_unit,*)
    !$$$     ' NUMBMM: Fixing row ',i,ic(i),ic(i+1)-1
    do j = ic(i),ic(i+1)-1
      if((jc(j)<1).or. (jc(j) > maxlmn))  then 
        write(psb_err_unit,*)&
             & ' NUMBMM: output problem',i,j,jc(j),maxlmn
      else
        c(j) = temp(jc(j))
        temp(jc(j)) = 0.
      endif
    end do
  end do

  return
end subroutine lcnumbmm
! == =====================================================================
! Sparse Matrix Multiplication Package
!
! Randolph E. Bank and Craig C. Douglas
!
! na.bank@na-net.ornl.gov and na.cdouglas@na-net.ornl.gov
!
! Compile this with the following command (or a similar one):
!
!     f77 -c -O smmp.f
!
! == =====================================================================
subroutine ldnumbmm(n, m, l,  ia, ja, diaga, a, ib, jb, diagb, b,&
     &  ic, jc, diagc, c, temp)
  use psb_const_mod
  !
  integer(psb_lpk_) :: ia(*), ja(*), diaga, ib(*), jb(*), diagb,&
       &  ic(*), jc(*), diagc 
  !
  real(psb_dpk_) :: a(*), b(*), c(*), temp(*),ajj
  !
  !       numeric matrix multiply c=a*b
  !
  maxlmn = max(l,m,n)
  do i = 1,maxlmn
    temp(i) = 0.
  end do
  minlm = min(l,m)
  minln = min(l,n)
  minmn = min(m,n)
  !
  !   c = a*b
  !
  do i = 1,n
    rowi: do jj = ia(i),ia(i+1)
      !    a = d + ...
      if (jj.eq.ia(i+1)) then
        if (diaga.eq.0 .or. i.gt.minmn) cycle rowi
        j = i
        ajj = a(i)
      else
        j=ja(jj)
        ajj = a(jj)
      endif
      !    b = d + ...
      if (diagb.eq.1 .and. j.le.minlm) &
           &  temp(j) = temp(j) + ajj * b(j)
      if ((j<1).or.(j>m)) then 
        write(psb_err_unit,*)&
             & ' NUMBMM: Problem with A ',i,jj,j,m
      endif

      do k = ib(j),ib(j+1)-1
        if((jb(k)<1).or. (jb(k) > maxlmn))  then 
          write(psb_err_unit,*)&
               & ' NUMBMM: jb  problem',j,k,jb(k),maxlmn
        else
          temp(jb(k)) = temp(jb(k)) + ajj * b(k)
        endif
      end do
    end do rowi

    !    c = d + ...
    if (diagc.eq.1 .and. i.le.minln) then
      c(i) = temp(i)
      temp(i) = 0.
    endif
    !$$$        if (mod(i,100) == 1)
    !$$$     +    write(psb_err_unit,*)
    !$$$     ' NUMBMM: Fixing row ',i,ic(i),ic(i+1)-1
    do j = ic(i),ic(i+1)-1
      if((jc(j)<1).or. (jc(j) > maxlmn))  then 
        write(psb_err_unit,*)&
             & ' NUMBMM: output problem',i,j,jc(j),maxlmn
      else
        c(j) = temp(jc(j))
        temp(jc(j)) = 0.
      endif
    end do
  end do

  return
end subroutine ldnumbmm
! == =====================================================================
! Sparse Matrix Multiplication Package
!
! Randolph E. Bank and Craig C. Douglas
!
! na.bank@na-net.ornl.gov and na.cdouglas@na-net.ornl.gov
!
! Compile this with the following command (or a similar one):
!
!     f77 -c -O smmp.f
!
! == =====================================================================
subroutine lsnumbmm(n, m, l,  ia, ja, diaga, a, ib, jb, diagb, b,&
     &  ic, jc, diagc, c, temp)
  use psb_const_mod
  !
  integer(psb_lpk_) :: ia(*), ja(*), diaga, ib(*), jb(*), diagb,&
       &  ic(*), jc(*), diagc 
  !
  real(psb_spk_) :: a(*), b(*), c(*), temp(*),ajj
  !
  !       numeric matrix multiply c=a*b
  !
  maxlmn = max(l,m,n)
  do i = 1,maxlmn
    temp(i) = 0.
  end do
  minlm = min(l,m)
  minln = min(l,n)
  minmn = min(m,n)
  !
  !   c = a*b
  !
  do  i = 1,n
    rowi: do jj = ia(i),ia(i+1)
      !    a = d + ...
      if (jj.eq.ia(i+1)) then
        if (diaga.eq.0 .or. i.gt.minmn) cycle rowi
        j = i
        ajj = a(i)
      else
        j=ja(jj)
        ajj = a(jj)
      endif
      !    b = d + ...
      if (diagb.eq.1 .and. j.le.minlm) &
           &  temp(j) = temp(j) + ajj * b(j)
      if ((j<1).or.(j>m)) then 
        write(psb_err_unit,*)&
             & ' NUMBMM: Problem with A ',i,jj,j,m
      endif

      do  k = ib(j),ib(j+1)-1
        if((jb(k)<1).or. (jb(k) > maxlmn))  then 
          write(psb_err_unit,*)&
               & ' NUMBMM: jb  problem',j,k,jb(k),maxlmn
        else
          temp(jb(k)) = temp(jb(k)) + ajj * b(k)
        endif
      end do
    end do rowi

    !    c = d + ...
    if (diagc.eq.1 .and. i.le.minln) then
      c(i) = temp(i)
      temp(i) = 0.
    endif
    !$$$        if (mod(i,100) == 1)
    !$$$     +    write(psb_err_unit,*)
    !$$$     ' NUMBMM: Fixing row ',i,ic(i),ic(i+1)-1
    do j = ic(i),ic(i+1)-1
      if((jc(j)<1).or. (jc(j) > maxlmn))  then 
        write(psb_err_unit,*)&
             &  ' NUMBMM: output problem',i,j,jc(j),maxlmn
      else
        c(j) = temp(jc(j))
        temp(jc(j)) = 0.
      endif
    end do
  end do

  return
end subroutine lsnumbmm
! == =====================================================================
! Sparse Matrix Multiplication Package
!
! Randolph E. Bank and Craig C. Douglas
!
! na.bank@na-net.ornl.gov and na.cdouglas@na-net.ornl.gov
!
! Compile this with the following command (or a similar one):
!
!     f77 -c -O smmp.f
!
! == =====================================================================
subroutine lznumbmm(n, m, l, ia, ja, diaga, a, ib, jb, diagb, b,&
     &  ic, jc, diagc, c, temp)
  !
  use psb_const_mod
  integer(psb_lpk_) :: ia(*), ja(*), diaga, ib(*), jb(*), diagb,&
       &  ic(*), jc(*), diagc 
  !
  complex(psb_dpk_) :: a(*), b(*), c(*), temp(*),ajj
  !
  !       numeric matrix multiply c=a*b
  !
  maxlmn = max(l,m,n)
  do i = 1,maxlmn
    temp(i) = 0.
  end do
  minlm = min(l,m)
  minln = min(l,n)
  minmn = min(m,n)
  !
  !   c = a*b
  !
  do i = 1,n
    rowi: do jj = ia(i),ia(i+1)
      !    a = d + ...
      if (jj.eq.ia(i+1)) then
        if (diaga.eq.0 .or. i.gt.minmn) cycle rowi
        j = i
        ajj = a(i)
      else
        j=ja(jj)
        ajj = a(jj)
      endif
      !    b = d + ...
      if (diagb.eq.1 .and. j.le.minlm) &
           & temp(j) = temp(j) + ajj * b(j)
      if ((j<1).or.(j>m)) then 
        write(psb_err_unit,*)&
             &  ' NUMBMM: Problem with A ',i,jj,j,m
      endif

      do k = ib(j),ib(j+1)-1
        if((jb(k)<1).or. (jb(k) > maxlmn))  then 
          write(psb_err_unit,*)&
               & ' NUMBMM: jb  problem',j,k,jb(k),maxlmn
        else
          temp(jb(k)) = temp(jb(k)) + ajj * b(k)
        endif
      end do
    end do rowi

    !    c = d + ...
    if (diagc.eq.1 .and. i.le.minln) then
      c(i) = temp(i)
      temp(i) = 0.
    endif
    do j = ic(i),ic(i+1)-1
      if((jc(j)<1).or. (jc(j) > maxlmn))  then 
        write(psb_err_unit,*)&
             & ' NUMBMM: output problem',i,j,jc(j),maxlmn
      else
        c(j) = temp(jc(j))
        temp(jc(j)) = 0.
      endif
    end do
  end do

  return
end subroutine lznumbmm
