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
! File:  psb_cfixcoo.f90 
! Subroutine: 
! Arguments:

Subroutine psb_cfixcoo(a,info,idir)
  use psb_spmat_type
  use psb_const_mod
  use psb_string_mod
  use psb_serial_mod, psb_protect_name => psb_cfixcoo
  use psb_error_mod
  implicit none

  !....Parameters...
  Type(psb_cspmat_type), intent(inout) :: A
  Integer, intent(out)                 :: info
  integer, intent(in), optional :: idir

  integer, allocatable :: iaux(:)
  !locals
  Integer              :: nza, nzl,iret,idir_, dupl_
  integer              :: i,j, irw, icl, err_act
  integer              :: debug_level, debug_unit
  character(len=20)    :: name = 'psb_fixcoo'

  info  = 0
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  if(debug_level >= psb_debug_serial_) &
       & write(debug_unit,*)  trim(name),': start ',&
       & size(a%ia1),size(a%ia2)
  if (psb_toupper(a%fida) /= 'COO') then 
    write(debug_unit,*) 'Fixcoo Invalid input ',a%fida
    info = -1
    return
  end if
  if (present(idir)) then 
    idir_ = idir
  else
    idir_ = 0
  endif

  nza = psb_sp_getifld(psb_nnz_,a,info)
  if (nza < 2) return

  dupl_ = psb_sp_getifld(psb_dupl_,a,info)

  allocate(iaux(nza+2),stat=info) 
  if (info /= 0) return

  select case(idir_) 

  case(0) !  Row major order

    call msort_up(nza,a%ia1(1),iaux(1),iret)
    if (iret == 0) call creordvn(nza,a%aspk(1),a%ia1(1),a%ia2(1),iaux(1))
    i    = 1
    j    = i
    do while (i <= nza)
      do while ((a%ia1(j) == a%ia1(i)))
        j = j+1
        if (j > nza) exit
      enddo
      nzl = j - i
      call msort_up(nzl,a%ia2(i),iaux(1),iret)
      if (iret == 0) &
           & call creordvn(nzl,a%aspk(i),a%ia1(i),a%ia2(i),iaux(1))
      i = j
    enddo

    i = 1
    irw = a%ia1(i)
    icl = a%ia2(i)
    j = 1

    select case(dupl_)
    case(psb_dupl_ovwrt_)

      do 
        j = j + 1
        if (j > nza) exit
        if ((a%ia1(j) == irw).and.(a%ia2(j) == icl)) then 
          a%aspk(i) = a%aspk(j)
        else
          i = i+1
          a%aspk(i) = a%aspk(j)
          a%ia1(i) = a%ia1(j)
          a%ia2(i) = a%ia2(j)
          irw = a%ia1(i) 
          icl = a%ia2(i) 
        endif
      enddo

    case(psb_dupl_add_)

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

    case(psb_dupl_err_)
      do 
        j = j + 1
        if (j > nza) exit
        if ((a%ia1(j) == irw).and.(a%ia2(j) == icl)) then 
          call psb_errpush(130,name)          
          goto 9999
        else
          i = i+1
          a%aspk(i) = a%aspk(j)
          a%ia1(i) = a%ia1(j)
          a%ia2(i) = a%ia2(j)
          irw = a%ia1(i) 
          icl = a%ia2(i) 
        endif
      enddo

    end select


    if(debug_level >= psb_debug_serial_)&
         & write(debug_unit,*)  trim(name),': end second loop'

  case(1) !  Col major order

    call msort_up(nza,a%ia2(1),iaux(1),iret)
    if (iret == 0) call creordvn(nza,a%aspk(1),a%ia1(1),a%ia2(1),iaux(1))
    i    = 1
    j    = i
    do while (i <= nza)
      do while ((a%ia2(j) == a%ia2(i)))
        j = j+1
        if (j > nza) exit
      enddo
      nzl = j - i
      call msort_up(nzl,a%ia1(i),iaux(1),iret)
      if (iret == 0) &
           & call creordvn(nzl,a%aspk(i),a%ia1(i),a%ia2(i),iaux(1))
      i = j
    enddo

    i = 1
    irw = a%ia1(i)
    icl = a%ia2(i)
    j = 1


    select case(dupl_)
    case(psb_dupl_ovwrt_)
      do 
        j = j + 1
        if (j > nza) exit
        if ((a%ia1(j) == irw).and.(a%ia2(j) == icl)) then 
          a%aspk(i) = a%aspk(j)
        else
          i = i+1
          a%aspk(i) = a%aspk(j)
          a%ia1(i) = a%ia1(j)
          a%ia2(i) = a%ia2(j)
          irw = a%ia1(i) 
          icl = a%ia2(i) 
        endif
      enddo

    case(psb_dupl_add_)
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

    case(psb_dupl_err_)
      do 
        j = j + 1
        if (j > nza) exit
        if ((a%ia1(j) == irw).and.(a%ia2(j) == icl)) then 
          call psb_errpush(130,name)
          goto 9999
        else
          i = i+1
          a%aspk(i) = a%aspk(j)
          a%ia1(i) = a%ia1(j)
          a%ia2(i) = a%ia2(j)
          irw = a%ia1(i) 
          icl = a%ia2(i) 
        endif
      enddo
    end select
    if (debug_level >= psb_debug_serial_)&
         & write(debug_unit,*)  trim(name),': end second loop'
  case default
    write(debug_unit,*) trim(name),': unknown direction ',idir_
  end select

  call psb_sp_setifld(psb_isrtdcoo_,psb_srtd_,a,info)
  call psb_sp_setifld(i,psb_nnz_,a,info)
  call psb_sp_setifld(psb_spmat_asb_,psb_state_,a,info)
  call psb_sp_setifld(psb_upd_srch_,psb_upd_,a,info)

  deallocate(iaux)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end Subroutine psb_cfixcoo
