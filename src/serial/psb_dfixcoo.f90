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
! File:  psb_dfixcoo.f90 
! Subroutine: 
! Parameters:

Subroutine psb_dfixcoo(A,INFO,idir)
  use psb_spmat_type
  use psb_const_mod
  use psb_string_mod
  implicit none

  !....Parameters...
  Type(psb_dspmat_type), intent(inout) :: A
  Integer, intent(out)                 :: info
  integer, intent(in), optional :: idir

  integer, allocatable :: iaux(:)
  !locals
  Integer              :: nza, nzl,iret,idir_
  integer              :: i,j, irw, icl
  logical, parameter   :: debug=.false.

  info  = 0
  if(debug) write(0,*)'fixcoo: ',size(a%ia1),size(a%ia2)
  if (toupper(a%fida) /= 'COO') then 
    write(0,*) 'Fixcoo Invalid input ',a%fida
    info = -1
    return
  end if
  if (present(idir)) then 
    idir_ = idir
  else
    idir_ = 0
  endif

  nza = a%infoa(psb_nnz_)
  if (nza < 2) return

  allocate(iaux(nza+2),stat=info) 
  if (info /= 0) return

  select case(idir_) 

  case(0) !  Row major order

    call mrgsrt(nza,a%ia1(1),iaux(1),iret)
    if (iret.eq.0) call reordvn(nza,a%aspk(1),a%ia1(1),a%ia2(1),iaux(1))
    i    = 1
    j    = i
    do while (i.le.nza)
      do while ((a%ia1(j).eq.a%ia1(i)))
        j = j+1
        if (j > nza) exit
      enddo
      nzl = j - i
      call mrgsrt(nzl,a%ia2(i),iaux(1),iret)
      if (iret.eq.0) &
           & call reordvn(nzl,a%aspk(i),a%ia1(i),a%ia2(i),iaux(1))
      i = j
    enddo

    i = 1
    irw = a%ia1(i)
    icl = a%ia2(i)
    j = 1
    do 
      j = j + 1
      if (j > nza) exit
      if ((a%ia1(j) == irw).and.(a%ia2(j) == icl)) then 
        a%aspk(i) = a%aspk(i) + a%aspk(j)
      else
        i = i+1
        a%aspk(i) = a%aspk(j)
        a%ia1(i) = a%ia1(j)
        a%ia2(i) = a%ia2(j)
        irw = a%ia1(i) 
        icl = a%ia2(i) 
      endif
    enddo
    a%infoa(psb_nnz_) = i    
    a%infoa(psb_srtd_) = psb_isrtdcoo_

    if(debug) write(0,*)'FIXCOO: end second loop'

  case(1) !  Col major order

    call mrgsrt(nza,a%ia2(1),iaux(1),iret)
    if (iret.eq.0) call reordvn(nza,a%aspk(1),a%ia1(1),a%ia2(1),iaux(1))
    i    = 1
    j    = i
    do while (i.le.nza)
      do while ((a%ia2(j).eq.a%ia2(i)))
        j = j+1
        if (j > nza) exit
      enddo
      nzl = j - i
      call mrgsrt(nzl,a%ia1(i),iaux(1),iret)
      if (iret.eq.0) &
           & call reordvn(nzl,a%aspk(i),a%ia1(i),a%ia2(i),iaux(1))
      i = j
    enddo

    i = 1
    irw = a%ia1(i)
    icl = a%ia2(i)
    j = 1
    do 
      j = j + 1
      if (j > nza) exit
      if ((a%ia1(j) == irw).and.(a%ia2(j) == icl)) then 
        a%aspk(i) = a%aspk(i) + a%aspk(j)
      else
        i = i+1
        a%aspk(i) = a%aspk(j)
        a%ia1(i) = a%ia1(j)
        a%ia2(i) = a%ia2(j)
        irw = a%ia1(i) 
        icl = a%ia2(i) 
      endif
    enddo
    a%infoa(psb_nnz_) = i    
    a%infoa(psb_srtd_) = psb_isrtdcoo_

    if(debug) write(0,*)'FIXCOO: end second loop'
  case default
    write(0,*) 'Fixcoo: unknown direction ',idir_
  end select
  a%infoa(psb_upd_) = psb_upd_srch_
  deallocate(iaux)
  return
end Subroutine psb_dfixcoo
