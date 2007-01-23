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
! File: psb_zspalloc.f90
!
! Subroutine: psb_zspalloc
!    Allocate sparse matrix structure for psblas routines.
! 
! Parameters: 
!    a        - type(<psb_zspmat_type>).       The sparse matrix to be allocated.      
!    desc_a   - type(<psb_desc_type>).         The communication descriptor to be updated.
!    info     - integer.                       Possibly returns an error code.
!    nnz      - integer(optional).             The number of nonzeroes in the matrix.
!
subroutine psb_zspalloc(a, desc_a, info, nnz)

  use psb_descriptor_type
  use psb_spmat_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  !....parameters...
  type(psb_desc_type), intent(inout) :: desc_a
  type(psb_zspmat_type), intent(out) :: a
  integer, intent(out)               :: info
  integer, optional, intent(in)      :: nnz

  !locals
  integer             :: ictxt
  integer             :: np,me,loc_row,&
       &  length_ia1,length_ia2, err_act,m,n
  integer             :: int_err(5)
  logical, parameter  :: debug=.false.
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=0
  call psb_erractionsave(err_act)
  name = 'psb_zspalloc'

  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  !     ....verify blacs grid correctness..
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  !
  ! hmm, not a good idea, not all compilers can rely on any given
  ! value for non initialized pointers. let's avoid this,
  ! and just rely on documentation. 
  ! check if psdalloc is already called for this matrix

  ! set fields in desc_a%matrix_data....
  loc_row = psb_cd_get_local_rows(desc_a)
  m       = psb_cd_get_global_rows(desc_a)
  n       = psb_cd_get_global_cols(desc_a)

  !...allocate matrix data...
  if (present(nnz))then 
    if (nnz < 0) then
      info=45
      int_err(1)=7
      int_err(2)=nnz
      call psb_errpush(info,name,int_err)
      goto 9999
    endif
    length_ia1=nnz
    length_ia2=nnz
  else 
    length_ia1=max(1,4*loc_row)
    length_ia2=max(1,4*loc_row)
  endif

  if (debug) write(*,*) 'allocating size:',length_ia1

  !....allocate aspk, ia1, ia2.....
  call psb_sp_all(loc_row,loc_row,a,length_ia1,info)
  if(info /= 0) then
    info=4010
    ch_err='sp_all'
    call psb_errpush(info,name,int_err)
    goto 9999
  end if

  ! set permutation matrices
  a%pl(1)=0
  a%pr(1)=0
  ! set infoa fields
  a%fida   = 'COO'
  a%descra = 'GUN'
  a%infoa(psb_nnz_)  = 0
  a%infoa(psb_srtd_) = 0
  a%infoa(psb_state_) = psb_spmat_bld_

  if (debug) write(0,*) 'spall: ',  &
       & psb_cd_get_dectype(desc_a),psb_desc_bld_

  return

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error(ictxt)
  end if
  return

end subroutine psb_zspalloc
