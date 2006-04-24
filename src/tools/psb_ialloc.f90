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
! File: psb_ialloc.f90
!
! Function: psb_ialloc
!    Allocates dense integer matrix for PSBLAS routines
! 
! Parameters: 
!    x      - the matrix to be allocated.
!    desc_a - the communication descriptor.
!    info   - possibly returns an error code
!    n      - optional number of columns.
subroutine psb_ialloc(x, desc_a, info, n)
  !....allocate dense  matrix for psblas routines.....
  use psb_descriptor_type
  use psb_const_mod
  use psb_error_mod
  implicit none
  
  !....parameters...
  integer, pointer                :: x(:,:)
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(out)            :: info
  integer, optional, intent(in)         :: n


  !locals
  integer             :: nprow,npcol,myrow,mypcol,err,n_col,n_row,i,j,jj,err_act
  integer             :: icontxt,dectype,n_
  integer             :: int_err(5),temp(1),exch(3)
  real(kind(1.d0))    :: real_err(5)
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  name='psb_ialloc'
  call psb_erractionsave(err_act)
  
  icontxt=desc_a%matrix_data(psb_ctxt_)
  
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mypcol)
  !     ....verify blacs grid correctness..
  if (nprow.eq.-1) then
     info = 2010
     call psb_errpush(info,name)
     goto 9999
  else if (npcol.ne.1) then
     info = 2030
     int_err(1) = npcol
     call psb_errpush(info,name,int_err)
     goto 9999
  endif

  dectype=desc_a%matrix_data(psb_dec_type_)
  !... check m and n parameters....
  if (.not.psb_is_ok_dec(dectype)) then 
    info = 3110
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(n)) then 
    n_ = n
  else
    n_ = 1
  endif
  !global check on n parameters
  if (myrow.eq.psb_root_) then
    exch(1)=n_
    call igebs2d(icontxt,psb_all_,psb_topdef_, ione,ione, exch, ione)
  else
    call igebr2d(icontxt,psb_all_,psb_topdef_, ione,ione, exch, ione, psb_root_, 0)
    if (exch(1).ne.n_) then
      info=550
      int_err(1)=1
      call psb_errpush(info,name,int_err)
      goto 9999
    endif
  endif

  !....allocate x .....
  if (psb_is_asb_dec(dectype).or.psb_is_upd_dec(dectype)) then
    n_col = max(1,desc_a%matrix_data(psb_n_col_))
    allocate(x(n_col,n_),stat=info)
    if (info.ne.0) then
      info=4010
      ch_err='allocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    endif
    do j=1,n_
      do i=1,n_col
        x(i,j) = 0
      end do
    end do
  else if (psb_is_bld_dec(dectype)) then
    n_row = max(1,desc_a%matrix_data(psb_n_row_))
    allocate(x(n_row,n_),stat=info)
    if (info.ne.0) then
      info=4010
      ch_err='allocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    endif
    do j = 1, n_
      do i=1,n_row
        x(i,j) = 0
      end do
    end do
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error(icontxt)
    return
  end if
  return

end subroutine psb_ialloc



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
! Function: psb_iallocv
!    Allocates dense matrix for PSBLAS routines
! 
! Parameters: 
!    m      - integer.                  The number of rows.
!    x      - integer,dimension(:).     The matrix to be allocated.
!    desc_a - type(<psb_desc_type>).    The communication descriptor.
!    info   - integer.                  Eventually returns an error code
subroutine psb_iallocv(x, desc_a, info,n)
  !....allocate sparse matrix structure for psblas routines.....
  use psb_descriptor_type
  use psb_const_mod
  use psb_realloc_mod
  use psb_error_mod

  implicit none
  
  !....parameters...
  integer, pointer                :: x(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(out)            :: info
  integer, optional, intent(in)   :: n

  !locals
  integer             :: nprow,npcol,myrow,mycol,err,n_col,n_row,dectype,i,err_act
  integer             :: icontxt, n_
  integer             :: int_err(5),temp(1),exch
  real(kind(1.d0))    :: real_err(5)
  logical, parameter  :: debug=.false. 
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  name='psb_iallocv'
  call psb_erractionsave(err_act)
  
  icontxt=desc_a%matrix_data(psb_ctxt_)

  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
  !     ....verify blacs grid correctness..
  if (nprow.eq.-1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (npcol.ne.1) then
    info = 2030
    int_err(1) = npcol
    call psb_errpush(info,name,int_err)
    goto 9999
  endif

  dectype=desc_a%matrix_data(psb_dec_type_)
  if (debug) write(0,*) 'dall: dectype',dectype
  if (debug) write(0,*) 'dall: is_ok? dectype',psb_is_ok_dec(dectype)
  !... check m and n parameters....
  if (.not.psb_is_ok_dec(dectype)) then 
    info = 3110
    call psb_errpush(info,name)
    goto 9999
  endif

  ! As this is a rank-1 array, optional parameter N is actually ignored.

  !....allocate x .....
  if (psb_is_asb_dec(dectype).or.psb_is_upd_dec(dectype)) then
     n_col = max(1,desc_a%matrix_data(psb_n_col_))
     allocate(x(n_col),stat=info)
     if (info.ne.0) then
        info=2025
        int_err(1)=n_col
        call psb_errpush(info,name,int_err)
        goto 9999
     endif
  else if (psb_is_bld_dec(dectype)) then
     n_row = max(1,desc_a%matrix_data(psb_n_row_))
     allocate(x(n_row),stat=info)
     if (info.ne.0) then
        info=2025
        int_err(1)=n_row
        call psb_errpush(info,name,int_err)
        goto 9999
     endif
  endif
     
  x = 0

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error(icontxt)
    return
  end if
  return

end subroutine psb_iallocv

