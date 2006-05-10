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
! File: psb_zasb.f90
!
! Subroutine: psb_zasb
!    Assembles a dense matrix for PSBLAS routines
! 
! Parameters: 
!    x       - real,pointer(dim=2).    The matrix to be assembled.
!    desc_a  - type(<psb_desc_type>).  The communication descriptor.
!    info    - integer.                Eventually returns an error code
subroutine psb_zasb(x, desc_a, info)
  !....assembly dense matrix x .....
  use psb_descriptor_type
  use psb_const_mod
  use psb_comm_mod
  use psb_error_mod
  use psb_realloc_mod
  implicit none

  type(psb_desc_type), intent(in) ::  desc_a
  complex(kind(1.d0)), pointer       ::  x(:,:)
  integer, intent(out)            ::  info

  ! local variables
  integer :: err, icontxt,nprow,npcol,me,mypcol,temp,lwork,nrow,ncol, err_act
  integer :: int_err(5), i1sz, i2sz, dectype, i,j
  double precision :: real_err(5)
  logical, parameter :: debug=.false.
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  name='psb_zasb'
  call psb_erractionsave(err_act)

  icontxt=desc_a%matrix_data(psb_ctxt_)
  dectype=desc_a%matrix_data(psb_dec_type_)

  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)

  if ((.not.associated(desc_a%matrix_data))) then
    info=3110
    call psb_errpush(info,name)
    goto 9999
  endif

  if (debug) write(*,*) 'asb start: ',nprow,npcol,me,&
       &desc_a%matrix_data(psb_dec_type_)
  !     ....verify blacs grid correctness..
  if (nprow.eq.-1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (npcol.ne.1) then
    info = 2030
    int_err(1) = npcol
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  else if (.not.psb_is_asb_dec(dectype)) then
    if (debug) write(*,*) 'asb error ',&
         &dectype
    info = 3110
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  ! check size
  icontxt=desc_a%matrix_data(psb_ctxt_)
  nrow=desc_a%matrix_data(psb_n_row_)
  ncol=desc_a%matrix_data(psb_n_col_)
  i1sz = size(x,dim=1)
  i2sz = size(x,dim=2)
  if (debug) write(*,*) 'asb: ',i1sz,i2sz,nrow,ncol
  if (i1sz.lt.ncol) then
    call psb_realloc(ncol,i2sz,x,info)
    if (info.ne.0) then
      info=2025
      int_err(1)=ncol
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    endif
  endif

  ! ..update halo elements..
  call psb_halo(x,desc_a,info)
  if(info /= 0) then
     info=4010
     ch_err='psb_halo'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  
  if (err_act.eq.act_ret) then
     return
  else
     call psb_error(icontxt)
  end if
  return
  
end subroutine psb_zasb


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
! Subroutine: psb_zasb
!    Assembles a dense matrix for PSBLAS routines
! 
! Parameters: 
!    x       - real,pointer(dim=1).    The matrix to be assembled.
!    desc_a  - type(<psb_desc_type>).  The communication descriptor.
!    info    - integer.                Eventually returns an error code
subroutine psb_zasbv(x, desc_a, info)
  !....assembly dense matrix x .....
  use psb_descriptor_type
  use psb_const_mod
  use psb_comm_mod
  use psb_error_mod
  use psb_realloc_mod
  implicit none

  type(psb_desc_type), intent(in) ::  desc_a
  complex(kind(1.d0)), pointer   ::  x(:)
  integer, intent(out)        ::  info

  ! local variables
  integer :: err, icontxt,nprow,npcol,me,mypcol,temp,lwork
  integer :: int_err(5), i1sz,nrow,ncol, dectype, i, err_act
  double precision :: real_err(5)

  logical, parameter :: debug=.false.
  character(len=20)             :: name,ch_err

  info = 0
  int_err(1) = 0
  name = 'psb_zasbv'
  
  icontxt=desc_a%matrix_data(psb_ctxt_)
  dectype=desc_a%matrix_data(psb_dec_type_)

  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)

  !     ....verify blacs grid correctness..
  if (nprow.eq.-1) then
     info = 2010
     call psb_errpush(info,name)
     goto 9999
  else if (npcol.ne.1) then
     info = 2030
     int_err(1) = npcol
     call psb_errpush(info,name,i_err=int_err)
     goto 9999
  else if (.not.psb_is_asb_dec(dectype)) then
     info = 3110
     call psb_errpush(info,name)
     goto 9999
  endif

  nrow=desc_a%matrix_data(psb_n_row_)
  ncol=desc_a%matrix_data(psb_n_col_)
  if (debug) write(*,*) name,' sizes: ',nrow,ncol
  i1sz = size(x)
  if (debug) write(*,*) 'dasb: sizes ',i1sz,ncol
  if (i1sz.lt.ncol) then
    call psb_realloc(ncol,x,info)
    if (info.ne.0) then           
      info=2025
      int_err(1)=ncol
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    endif

  endif  
  
  ! ..update halo elements..
  call psb_halo(x,desc_a,info)
  if(info /= 0) then
     info=4010
     ch_err='f90_pshalo'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  
  if (err_act.eq.act_ret) then
     return
  else
     call psb_error(icontxt)
  end if
  return
  
end subroutine psb_zasbv

