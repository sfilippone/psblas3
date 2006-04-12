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
! File: psb_dins.f90
!
! Subroutine: psb_dins
!    Insert dense submatrix to dense matrix.
! 
! Parameters: 
!    m       - integer.                       Rows number of submatrix belonging to blck to be inserted.
!    n       - integer.                       Cols number of submatrix belonging to blck to be inserted.
!    x       - real, pointer, dimension(:,:). The destination dense matrix.  
!    ix      - integer.                       x global-row corresponding to position at which blck submatrix must be inserted.
!    jx      - integer.                       x global-col corresponding to position at which blck submatrix must be inserted.
!    blck    - real, pointer, dimension(:,:). The source dense submatrix.  
!    desc_a  - type(<psb_desc_type>).         The communication descriptor.
!    info    - integer.                       Eventually returns an error code
!    iblck   - integer(optional).             First row of submatrix belonging to blck to be inserted.
!    jblck   - integer(optional).             First col of submatrix belonging to blck to be inserted.
subroutine psb_dins(m, n, x, ix, jx, blck, desc_a, info,&
     & iblck, jblck)
  !....insert dense submatrix to dense matrix .....
  use psb_descriptor_type
  use psb_const_mod
  use psb_error_mod
  implicit none


  !....parameters...
  integer, intent(in)                ::  m,n
  type(psb_desc_type), intent(in)    ::  desc_a
  real(kind(1.d0)),pointer           ::  x(:,:)
  integer, intent(in)                ::  ix,jx
  real(kind(1.d0)), intent(in)       ::  blck(:,:)
  integer,intent(out)                ::  info
  integer, optional, intent(in)      ::  iblck,jblck

  !locals.....

  integer                :: icontxt,i,loc_row,glob_row,row,k,err_act,&
       & nprocs,mode, loc_cols,col,iblock, jblock, mglob, int_err(5), err
  integer                :: nprow,npcol, me ,mypcol
  character              :: temp_descra*11,temp_fida*5
  character(len=20)   :: name, char_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)
  name = 'psb_dins'

  if (.not.associated(desc_a%glob_to_loc)) then
    info=3110
    call psb_errpush(info,name)
    return
  end if
  if ((.not.associated(desc_a%matrix_data))) then
    int_err(1)=3110
    call psb_errpush(info,name)
    return
  end if

  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
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

end subroutine psb_dins





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
! Subroutine: psb_dinsvm
!    Insert dense submatrix to dense matrix.
! 
! Parameters: 
!    m       - integer.                       Rows number of submatrix belonging to blck to be inserted.
!    x       - real, pointer, dimension(:,:). The destination dense matrix.  
!    ix      - integer.                       x global-row corresponding to position at which blck submatrix must be inserted.
!    jx      - integer.                       x global-col corresponding to position at which blck submatrix must be inserted.
!    blck    - real, pointer, dimension(:,:). The source dense submatrix.  
!    desc_a  - type(<psb_desc_type>).         The communication descriptor.
!    info    - integer.                       Eventually returns an error code
!    iblck   - integer(optional).             First row of submatrix belonging to blck to be inserted.
subroutine psb_dinsvm(m, x, ix, jx, blck, desc_a,info,&
     & iblck)
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
  integer, intent(in)                ::  m
  type(psb_desc_type), intent(in)    ::  desc_a
  real(kind(1.d0)),pointer           ::  x(:,:)
  integer, intent(in)                ::  ix,jx
  real(kind(1.d0)), intent(in)       ::  blck(:)
  integer, intent(out)               ::  info
  integer, optional, intent(in)      ::  iblck

  !locals.....
  integer                :: icontxt,i,loc_row,glob_row,loc_cols,mglob,err_act, int_err(5),err
  integer                :: nprow,npcol, me ,mypcol, iblock
  character(len=20)   :: name, char_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)
  name = 'psb_dinsvm'

  if (.not.associated(desc_a%glob_to_loc)) then
    info=3110
    call psb_errpush(info,name)
    return
  end if
  if ((.not.associated(desc_a%matrix_data))) then
    int_err(1)=3110
    call psb_errpush(info,name)
    return
  end if

  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
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
  else if (size(x, dim=2).lt.1) then
    ! check if dimension of x is greater than dimension of submatrix
    ! to insert
    info = 320
    int_err(1) = 2
    int_err(2) = size(x, dim=2)
    int_err(3) = 1
    call psb_errpush(info,name,int_err)
    goto 9999
  endif

  loc_cols=desc_a%matrix_data(psb_n_col_)
  mglob    = desc_a%matrix_data(psb_m_)

  if (present(iblck)) then
    iblock = iblck
  else
    iblock = 1
  endif

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

end subroutine psb_dinsvm



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
! Subroutine: psb_dinsvv
!    Insert dense submatrix to dense matrix.
! 
! Parameters: 
!    m       - integer.                       Rows number of submatrix belonging to blck to be inserted.
!    x       - real, pointer, dimension(:).   The destination dense matrix.  
!    ix      - integer.                       x global-row corresponding to position at which blck submatrix must be inserted.
!    blck    - real, pointer, dimension(:).   The source dense submatrix.  
!    desc_a  - type(<psb_desc_type>).         The communication descriptor.
!    info    - integer.                       Eventually returns an error code
!    iblck   - integer(optional).             First row of submatrix belonging to blck to be inserted.
!    insflag - integer(optional).             ???                                                                             
subroutine psb_dinsvv(m, x, ix, blck, desc_a, info,&
     & iblck,insflag)
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
  integer, intent(in)                ::  m
  type(psb_desc_type), intent(in)    ::  desc_a
  real(kind(1.d0)),pointer           ::  x(:)
  integer, intent(in)                ::  ix
  real(kind(1.d0)), intent(in)       ::  blck(:)
  integer, intent(out)               ::  info
  integer, optional, intent(in)      ::  iblck
  integer, optional, intent(in)      ::  insflag

  !locals.....
  integer                :: icontxt,i,loc_row,glob_row,row,k,&
       & loc_rows,loc_cols,iblock, liflag,mglob,err_act, int_err(5), err
  integer                :: nprow,npcol, me ,mypcol
  character(len=20)   :: name, char_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)
  name = 'psb_dinsvv'

  if (.not.associated(desc_a%glob_to_loc)) then
    info=3110
    call psb_errpush(info,name)
    return
  end if
  if ((.not.associated(desc_a%matrix_data))) then
    int_err(1)=3110
    call psb_errpush(info,name)
    return
  end if

  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
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
  if (present(insflag)) then
    liflag = insflag
  else
    liflag = psb_upd_glbnum_
  end if

  if (liflag == psb_upd_glbnum_) then 
    do i = 1, m
      !loop over all blck's rows

      ! row actual block row 
      glob_row=ix+i-1
      if (glob_row > mglob) exit

      loc_row=desc_a%glob_to_loc(glob_row)
      if (loc_row.ge.1) then
        ! this row belongs to me
        ! copy i-th row of block blck in x
        x(loc_row) = x(loc_row) +  blck(iblock+i-1)
      end if
    enddo
  else if (liflag == psb_upd_locnum_) then 
    k = min(ix+m-1,loc_rows)
    do i=ix,k
      x(i) = x(i) + blck(i-ix+1)
    enddo
  else
    info=-1
    call psb_errpush(info,name)
    goto 9999
  endif

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

end subroutine psb_dinsvv


