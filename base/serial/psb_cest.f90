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
subroutine psb_cest(afmt, m,n,nnz, lia1, lia2, lar, iup, info)

  use psb_error_mod
  use psb_const_mod
  use psb_string_mod
  use psb_spmat_type
  implicit none

  !     .. scalar arguments ..
  integer, intent(in) ::  m,n,nnz,iup
  integer, intent(out) :: lia1, lia2, lar, info
  !     .. array arguments..
  character(len=5)  ::  afmt
  integer           ::  int_val(5), err_act
  character(len=20) ::  name

  name = 'psb_cest'      
  call psb_erractionsave(err_act)

  if (afmt == '???') then 
    afmt = psb_fidef_
  endif

  select case(iup)
  case (psb_upd_perm_)
    if (toupper(afmt) == 'JAD') then 
      lia1 = 2*(nnz + nnz/5) +1000
      lia2 = 2*(nnz + nnz/5) +1000 +m
      lar = nnz + nnz/5
    else if (toupper(afmt) == 'COO') then 
      lia1 = nnz
      lia2 = 2*nnz + 1000
      lar = nnz
    else if(toupper(afmt) == 'CSR') then
      lia1 = nnz
      lia2 = 2*nnz + 1000 + m + 1
      lar = nnz
    else
      info = 136
      call psb_errpush(info,name,a_err=toupper(afmt))
      goto 9999
    endif

  case (psb_upd_dflt_, psb_upd_srch_)

    if (toupper(afmt) == 'JAD') then 
      lia1 = nnz + nnz/5
      lia2 = nnz + nnz/5
      lar = nnz + nnz/5
    else if (toupper(afmt) == 'COO') then 
      lia1 = nnz
      lia2 = nnz
      lar = nnz
    else if(toupper(afmt) == 'CSR') then
      lia1 = nnz
      lia2 = max(nnz,m+1)
      lar = nnz
    else
      info = 136
      call psb_errpush(info,name,a_err=toupper(afmt))
      goto 9999
    endif

  case default

    info = 3012
    call psb_errpush(info,name,int_val)
    goto 9999

  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if ( err_act .ne. 0 ) then 
     call psb_error()
     return
  endif

  return

end subroutine psb_cest

