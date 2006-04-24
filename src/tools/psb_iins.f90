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
! File: psb_iins.f90
!
! Subroutine: psb_iins
!    Insert dense integer submatrix to dense integer matrix.
! 
! Parameters: 
!    m       - integer.                          Rows number of submatrix belonging to blck to be inserted.
!    n       - integer.                          Cols number of submatrix belonging to blck to be inserted.
!    x       - integer, pointer, dimension(:,:). The destination dense matrix.  
!    ix      - integer.                          x global-row corresponding to position at which blck submatrix must be inserted.
!    jx      - integer.                          x global-col corresponding to position at which blck submatrix must be inserted.
!    blck    - integer, pointer, dimension(:,:). The source dense submatrix.  
!    desc_a  - type(<psb_desc_type>).            The communication descriptor.
!    info    - integer.                          Eventually returns an error code
!    iblck   - integer(optional).                First row of submatrix belonging to blck to be inserted.
!    jblck   - integer(optional).                First col of submatrix belonging to blck to be inserted.
subroutine psb_iins(m, n, blck, x, ix, jx, desc_a, info,&
     & iblck, jblck,dupl)
  !....insert dense submatrix to dense matrix .....
  use psb_descriptor_type
  use psb_const_mod
  use psb_error_mod
  implicit none

  !....parameters...
  integer, intent(in)              ::  m,n
  type(psb_desc_type), intent(in)  ::  desc_a
  integer, pointer                 ::  x(:,:)
  integer, intent(in)              ::  ix,jx
  integer, intent(in)              ::  blck(:,:)
  integer, intent(out)             ::  info
  integer, optional, intent(in)    ::  iblck,jblck
  integer, optional, intent(in)      ::  dupl

  !locals.....

  integer                :: icontxt,i,loc_row,glob_row,&
       & loc_cols,col,iblock, jblock, mglob,dupl_
  integer                :: nprow,npcol, myrow ,mycol, int_err(5),err_act
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)
  name = 'psb_iins'


  if ((.not.associated(desc_a%matrix_data))) then
    info=3110
    call psb_errpush(info, name)
    return
  end if

  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
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

  if (.not.associated(desc_a%glob_to_loc)) then
    info=3110
    call psb_errpush(info,name,int_err)
    goto 9999
  end if

  !... check parameters....
  if (m.lt.0) then
    info = 10
    int_err(1) = 1
    int_err(2) = m
    call psb_errpush(info,name,int_err)
    goto 9999
  else if (n.lt.0) then
    info = 10
    int_err(1) = 2
    int_err(2) = n
    call psb_errpush(info,name,int_err)
    goto 9999
  else if (ix.lt.1) then
    info = 20
    int_err(1) = 6
    int_err(2) = ix
    call psb_errpush(info,name,int_err)
    goto 9999
  else if (jx.lt.1) then
    info = 20
    int_err(1) = 7
    int_err(2) = jx
    call psb_errpush(info,name,int_err)
    goto 9999
  else if (.not.psb_is_ok_dec(desc_a%matrix_data(psb_dec_type_))) then
    info = 3110
    int_err(1) = desc_a%matrix_data(psb_dec_type_)
    call psb_errpush(info,name,int_err)
    goto 9999
  else if (size(x, dim=1).lt.desc_a%matrix_data(psb_n_row_)) then
    info = 310
    int_err(1) = 5
    int_err(2) = 4
    call psb_errpush(info,name,int_err)
    goto 9999
  else if (size(x, dim=2).lt.n) then
    ! check if dimension of x is greater than dimension of submatrix
    ! to insert
    info = 320
    int_err(1) = 2
    int_err(2) = size(x, dim=2)
    int_err(3) = n
    call psb_errpush(info,name,int_err)
    goto 9999
  endif

  loc_cols = desc_a%matrix_data(psb_n_col_)
  mglob    = desc_a%matrix_data(psb_m_)
  if (present(iblck)) then
    iblock = iblck
  else
    iblock = 1
  endif

  if (present(jblck)) then
    jblock = jblck
  else
    jblock = 1
  endif
  if (present(dupl)) then 
    dupl_ = dupl
  else
    dupl_ = psb_dupl_ovwrt_
  endif

  select case(dupl_) 
  case(psb_dupl_ovwrt_) 
    do i = 1, m
      !loop over all blck's rows

      ! row actual block row 
      glob_row=ix+i-1
      if (glob_row > mglob) exit
      loc_row=desc_a%glob_to_loc(glob_row)
      if (loc_row.ge.1) then
        ! this row belongs to me
        ! copy i-th row of block blck in x
        do col = 1, n
          x(loc_row,jx+col-1) = blck(iblock+i-1,jblock+col-1)
        enddo
      end if
    enddo
  case(psb_dupl_add_) 
    do i = 1, m
      !loop over all blck's rows

      ! row actual block row 
      glob_row=ix+i-1
      if (glob_row > mglob) exit
      loc_row=desc_a%glob_to_loc(glob_row)
      if (loc_row.ge.1) then
        ! this row belongs to me
        ! copy i-th row of block blck in x
        do col = 1, n
          x(loc_row,jx+col-1) = x(loc_row,jx+col-1) + blck(iblock+i-1,jblock+col-1)
        enddo
      end if
    enddo
  case default
    info = 321
    call psb_errpush(info,name)
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error(icontxt)
    return
  end if
  return

end subroutine psb_iins





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
! Subroutine: psb_iinsvm
!    Insert dense integer submatrix to dense integer matrix.
! 
! Parameters: 
!    m       - integer.                          Rows number of submatrix belonging to blck to be inserted.
!    x       - integer, pointer, dimension(:,:). The destination dense matrix.  
!    ix      - integer.                          x global-row corresponding to position at which blck submatrix must be inserted.
!    jx      - integer.                          x global-col corresponding to position at which blck submatrix must be inserted.
!    blck    - integer, pointer, dimension(:,:). The source dense submatrix.  
!    desc_a  - type(<psb_desc_type>).            The communication descriptor.
!    info    - integer.                          Eventually returns an error code
!    iblck   - integer(optional).                First row of submatrix belonging to blck to be inserted.
subroutine psb_iinsvm(m, blck, x, ix, jx, desc_a, info,&
     & iblck,dupl)
  !....insert dense submatrix to dense matrix .....
  use psb_descriptor_type
  use psb_const_mod
  use psb_error_mod
  implicit none

  ! m rows number of submatrix belonging to blck to be inserted

  ! iblck first row of submatrix belonging to blck to be inserted

  ! ix  x global-row corresponding to position at which blck submatrix
  !     must be inserted

  ! jx  x global-col corresponding to position at which blck submatrix
  !     must be inserted

  !....parameters...
  integer, intent(in)              ::  m
  type(psb_desc_type), intent(in)  ::  desc_a
  integer, pointer                 ::  x(:,:)
  integer, intent(in)              ::  ix,jx
  integer, intent(in)              ::  blck(:)
  integer, intent(out)             ::  info
  integer, optional, intent(in)    ::  iblck
  integer, optional, intent(in)      ::  dupl

  !locals.....
  integer                :: icontxt,i,loc_row,glob_row,&
       & loc_cols,iblock, jblock,mglob, err_act, int_err(5)
  integer                :: nprow,npcol, myrow ,mycol,dupl_
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  name = 'psb_iinsvm'
  call psb_erractionsave(err_act)


  loc_cols=desc_a%matrix_data(psb_n_col_)
  mglob    = desc_a%matrix_data(psb_m_)

  if (present(iblck)) then
    iblock = iblck
  else
    iblock = 1
  endif

  if (present(dupl)) then 
    dupl_ = dupl
  else
    dupl_ = psb_dupl_ovwrt_
  endif

  select case(dupl_) 
  case(psb_dupl_ovwrt_) 
    do i = 1, m
      !loop over all blck's rows

      ! row actual block row 
      glob_row=ix+i-1
      if (glob_row > mglob) exit

      loc_row=desc_a%glob_to_loc(glob_row)
      if (loc_row.ge.1) then
        ! this row belongs to me
        ! copy i-th row of block blck in x
        x(loc_row,jx) = blck(iblock+i-1)
      end if
    enddo
  case(psb_dupl_add_) 
    do i = 1, m
      !loop over all blck's rows

      ! row actual block row 
      glob_row=ix+i-1
      if (glob_row > mglob) exit

      loc_row=desc_a%glob_to_loc(glob_row)
      if (loc_row.ge.1) then
        ! this row belongs to me
        ! copy i-th row of block blck in x
        x(loc_row,jx) = x(loc_row,jx) + blck(iblock+i-1)
      end if
    enddo
  case default
    info = 321
    call psb_errpush(info,name)
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error(icontxt)
    return
  end if
  return

end subroutine psb_iinsvm


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
! Subroutine: psb_iinsvv
!    Insert dense integer submatrix to dense integer matrix.
! 
! Parameters: 
!    m       - integer.                          Rows number of submatrix belonging to blck to be inserted.
!    x       - integer, pointer, dimension(:,:). The destination dense matrix.  
!    ix      - integer.                          x global-row corresponding to position at which blck submatrix must be inserted.
!    blck    - integer, pointer, dimension(:,:). The source dense submatrix.  
!    desc_a  - type(<psb_desc_type>).            The communication descriptor.
!    info    - integer.                          Eventually returns an error code
!    iblck   - integer(optional).                First row of submatrix belonging to blck to be inserted.
subroutine psb_iinsvv(m, blck, x, ix, desc_a, info,&
     & iblck,dupl)
  !....insert dense submatrix to dense matrix .....
  use psb_descriptor_type
  use psb_const_mod
  use psb_error_mod
  implicit none

  ! m rows number of submatrix belonging to blck to be inserted

  ! iblck first row of submatrix belonging to blck to be inserted

  ! ix  x global-row corresponding to position at which blck submatrix
  !     must be inserted

  !....parameters...
  integer, intent(in)              ::  m
  type(psb_desc_type), intent(in)  ::  desc_a
  integer,  pointer                ::  x(:)
  integer, intent(in)              ::  ix
  integer,  intent(in)             ::  blck(:)
  integer, intent(out)             ::  info
  integer, optional, intent(in)    ::  iblck
  integer, optional, intent(in)      ::  dupl

  !locals.....
  integer                :: icontxt,i,loc_row,glob_row,k,&
       & loc_rows,loc_cols,col,iblock, jblock, mglob, err_act, int_err(5)
  integer                :: nprow,npcol, myrow ,mycol,dupl_
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  name = 'psb_iinsvv'
  call psb_erractionsave(err_act)

  if ((.not.associated(desc_a%matrix_data))) then
    info=3110
    call psb_errpush(info,name)
    return
  end if
  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
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

  if (.not.associated(desc_a%glob_to_loc)) then
    info=3110
    call psb_errpush(info,name)
    goto 9999
  end if

  !... check parameters....
  if (m.lt.0) then
    info = 10
    int_err(1) = 1
    int_err(2) = m
    call psb_errpush(info,name,int_err)
    goto 9999
  else if (ix.lt.1) then
    info = 20
    int_err(1) = 6
    int_err(2) = ix
    call psb_errpush(info,name,int_err)
    goto 9999
  else if (.not.psb_is_ok_dec(desc_a%matrix_data(psb_dec_type_))) then
    info = 3110
    int_err(1) = desc_a%matrix_data(psb_dec_type_)
    call psb_errpush(info,name,int_err)
    goto 9999
  else if (size(x, dim=1).lt.desc_a%matrix_data(psb_n_row_)) then
    info = 310
    int_err(1) = 5
    int_err(2) = 4
    call psb_errpush(info,name,int_err)
    goto 9999
  endif

  loc_rows=desc_a%matrix_data(psb_n_row_)
  loc_cols=desc_a%matrix_data(psb_n_col_)
  mglob    = desc_a%matrix_data(psb_m_)

  if (present(iblck)) then
    iblock = iblck
  else
    iblock = 1
  endif
  if (present(dupl)) then 
    dupl_ = dupl
  else
    dupl_ = psb_dupl_ovwrt_
  endif


  select case(dupl_) 
  case(psb_dupl_ovwrt_) 
    do i = 1, m
      !loop over all blck's rows

      ! row actual block row 
      glob_row=ix+i-1
      if (glob_row > mglob) exit

      loc_row=desc_a%glob_to_loc(glob_row)
      if (loc_row.ge.1) then
        ! this row belongs to me
        ! copy i-th row of block blck in x
	x(loc_row) = blck(iblock+i-1)
      end if
    enddo
  case(psb_dupl_add_) 
    do i = 1, m
      !loop over all blck's rows

      ! row actual block row 
      glob_row=ix+i-1
      if (glob_row > mglob) exit

      loc_row=desc_a%glob_to_loc(glob_row)
      if (loc_row.ge.1) then
        ! this row belongs to me
        ! copy i-th row of block blck in x
	x(loc_row) = x(loc_row) + blck(iblock+i-1)
      end if
    enddo
  case default
    info = 321
    call psb_errpush(info,name)
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error(icontxt)
    return
  end if
  return

end subroutine psb_iinsvv


