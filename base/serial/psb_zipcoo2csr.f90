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
! File:  psb_zipcoo2csr.f90 
! Subroutine: 
! Parameters:

subroutine psb_zipcoo2csr(a,info,rwshr)
  use psb_spmat_type
  use psb_const_mod
  use psb_realloc_mod
  use psb_serial_mod, only : psb_fixcoo
  use psb_error_mod
  use psb_string_mod
  implicit none

  !....Parameters...
  Type(psb_zspmat_type), intent(inout) :: A
  Integer, intent(out)                 :: info
  logical, optional                    :: rwshr

  integer, allocatable :: iaux(:), itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act
  Integer, Parameter  :: maxtry=8
  logical, parameter  :: debug=.false.
  character(len=20)   :: name

  name='psb_ipcoo2csr'
  info  = 0
  call psb_erractionsave(err_act)

  if(debug) write(0,*)'Inside ipcoo2csr ',a%fida,a%m
  if (toupper(a%fida) /= 'COO') then 
    write(0,*) 'ipcoo2csr Invalid input ',a%fida
    info = -1
    call psb_errpush(info,name)
    goto 9999
  end if
  if (present(rwshr)) then 
    rwshr_ = rwshr
  else
    rwshr_ = .false.
  end if

  call psb_fixcoo(a,info)
  nr  = a%m 
  nza = a%infoa(psb_nnz_)
  allocate(iaux(nr+1),stat=info)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if

  if(debug) write(0,*)'DIPCOO2CSR: out of fixcoo',nza,nr,size(a%ia2),size(iaux)

  call psb_transfer(a%ia1,itemp,info)
  call psb_transfer(a%ia2,a%ia1,info)
  call psb_transfer(iaux,a%ia2,info)

  !
  ! This routine can be used in two modes:
  ! 1. Normal: just look at the row indices and trust them. This
  !    implies putting in empty rows where needed. In this case you
  !    can get in trouble if A%M < A%IA1(NZA)
  ! 2. Shrink mode: disregard the actual value of the row indices,
  !    just treat them as ident markers. In this case you can get in
  !    trouble when the number of distinct row indices is greater 
  !    than A%M
  !
  !

  a%ia2(1) = 1

  if (nza <= 0) then 
    do i=1,nr
      a%ia2(i+1) = a%ia2(i)
    end do
  else

    if (rwshr_) then 
      
      
      j = 1
      i = 1
      irw = itemp(j) 

      do j=1, nza
        if (itemp(j) /= irw) then 
          a%ia2(i+1) = j
          irw = itemp(j) 
          i = i + 1
          if (i>nr) then 
            write(0,*) 'IPCOO2CSR: RWSHR=.true. : ',&
             & i, nr,' Expect trouble!'
            exit
          end if
        endif
      enddo 
!      write(0,*) 'Exit from loop',j,nza,i
      do 
        if (i>=nr+1) exit
        a%ia2(i+1) = j
        i = i + 1
      end do

    else
      
      if (nr < itemp(nza)) then 
        write(0,*) 'IPCOO2CSR: RWSHR=.false. : ',&
             &nr,itemp(nza),' Expect trouble!'
      end if
             

      j = 1 
      i = 1
      irw = itemp(j) 

      outer: do 
        inner: do 
          if (i >= irw) exit inner
          if (i>nr) then 
            write(0,*) 'Strange situation: i>nr ',i,nr,j,nza,irw,idl
            exit outer
          end if
          a%ia2(i+1) = a%ia2(i) 
          i = i + 1
        end do inner
        j = j + 1
        if (j > nza) exit
        if (itemp(j) /= irw) then 
          a%ia2(i+1) = j
          irw = itemp(j) 
          i = i + 1
        endif
        if (i>nr) exit
      enddo outer
      !
      ! Cleanup empty rows at the end
      !
      if (j /= (nza+1)) then 
        write(0,*) 'IPCOO2CSR : Problem from loop :',j,nza
      endif
      do 
        if (i>nr) exit
        a%ia2(i+1) = j
        i = i + 1
      end do

    endif

  end if

!!$  write(0,*) 'IPcoo2csr end loop ',i,nr,a%ia2(nr+1),nza
  a%fida='CSR'
  a%infoa(psb_upd_) = psb_upd_srch_

  deallocate(itemp)
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
     call psb_error()
     return
  end if
  return

end Subroutine psb_zipcoo2csr
