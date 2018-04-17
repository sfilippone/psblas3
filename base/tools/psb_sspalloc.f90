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
! File: psb_sspalloc.f90
!
! Subroutine: psb_sspalloc
!    Allocate sparse matrix structure for psblas routines.
! 
! Arguments: 
!    a        - type(psb_sspmat_type).       The sparse matrix to be allocated.      
!    desc_a   - type(psb_desc_type).         The communication descriptor to be updated.
!    info     - integer.                       Return code.
!    nnz      - integer(optional).             The number of nonzeroes in the matrix.
!                                              (local, user estimate)
!
subroutine psb_sspalloc(a, desc_a, info, nnz)
  use psb_base_mod, psb_protect_name => psb_sspalloc
  implicit none

  !....parameters...
  type(psb_desc_type), intent(in) :: desc_a
  type(psb_sspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(out)               :: info
  integer(psb_ipk_), optional, intent(in)      :: nnz

  !locals
  integer(psb_ipk_) :: ictxt, dectype
  integer(psb_ipk_) :: np,me,loc_row,loc_col,&
       &  length_ia1,length_ia2, err_act,m,n
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_sspall'
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt   = desc_a%get_context()
  dectype = desc_a%get_dectype()

  call psb_info(ictxt, me, np)
  !     ....verify blacs grid correctness..
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif


  loc_row = desc_a%get_local_rows()
  loc_col = desc_a%get_local_cols()
  m       = desc_a%get_global_rows()
  n       = desc_a%get_global_cols()

  !...allocate matrix data...
  if (present(nnz))then 
    if (nnz < 0) then
      info=45
      call psb_errpush(info,name,i_err=(/7_psb_ipk_,nnz/))
      goto 9999
    endif
    length_ia1=nnz
    length_ia2=nnz
  else 
    length_ia1=max(1,5*loc_row)
  endif

  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),':allocating size:',length_ia1
  call a%free()
  !....allocate aspk, ia1, ia2.....
  call a%csall(loc_row,loc_col,info,nz=length_ia1)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='sp_all')
    goto 9999
  end if

  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': ',  &
       & desc_a%get_dectype(),psb_desc_bld_

  return

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_sspalloc
