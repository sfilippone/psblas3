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
! File: psb_dspasb.f90
!
! Subroutine: psb_dspasb
!    Assemble sparse matrix
!
! Arguments: 
!    a        - type(psb_dspmat_type).     The sparse matrix to be allocated.      
!    desc_a   - type(psb_desc_type).       The communication descriptor.
!    info     - integer.                     return code.
!    afmt     - character(optional)          The desired output storage format.
!    upd      - character(optional).         How will the matrix be updated? 
!                                            psb_upd_srch_    Simple strategy  
!                                            psb_upd_perm_    Permutation(more memory)
!    dupl     - integer(optional).           Duplicate coefficient handling:
!                                            psb_dupl_ovwrt_     overwrite
!                                            psb_dupl_add_       add 
!                                            psb_dupl_err_       raise an error. 
! 
!
subroutine psb_dspasb(a,desc_a, info, afmt, upd, dupl, mold)
  use psb_descriptor_type
  use psb_spmat_type
  use psb_serial_mod
  use psb_const_mod
  use psi_mod
  use psb_error_mod
  use psb_string_mod
  use psb_penv_mod
  use psbn_d_mat_mod
  implicit none


  !...Parameters....
  type(psbn_d_sparse_mat), intent (inout)  :: a
  type(psb_desc_type), intent(in)         :: desc_a
  integer, intent(out)                    :: info
  integer,optional, intent(in)            :: dupl, upd
  character(len=*), optional, intent(in)         :: afmt
  class(psbn_d_base_sparse_mat), intent(in), optional :: mold
  !....Locals....
  integer               :: int_err(5)
  integer               :: np,me,n_col, err_act
  integer               :: spstate
  integer               :: ictxt,n_row
  integer              :: debug_level, debug_unit
  character(len=20)     :: name, ch_err

  info = 0
  int_err(1)=0
  name = 'psb_spasb'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt    = psb_cd_get_context(desc_a)
  n_row    = psb_cd_get_local_rows(desc_a)
  n_col    = psb_cd_get_local_cols(desc_a)

  ! check on BLACS grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.psb_is_asb_desc(desc_a)) then 
    info = 600
    int_err(1) = psb_cd_get_dectype(desc_a)
    call psb_errpush(info,name)
    goto 9999
  endif

  if (debug_level >= psb_debug_ext_)&
       & write(debug_unit, *) me,' ',trim(name),&
       & '   Begin matrix assembly...'

  !check on errors encountered in psdspins


  if (a%is_bld()) then 
    !
    ! First case: we come from a fresh build. 
    ! 

    n_row = psb_cd_get_local_rows(desc_a)
    n_col = psb_cd_get_local_cols(desc_a)
    call a%set_nrows(n_row)
    call a%set_ncols(n_col)
  end if

  call a%cscnv(info,type=afmt,dupl=dupl, mold=mold)

  
  IF (debug_level >= psb_debug_ext_) then 
    ch_err=a%get_fmt()
    write(debug_unit, *) me,' ',trim(name),':  From SPCNV',&
         & info,' ',ch_err
  end IF
  
  if (info /= psb_no_err_) then    
    info=4010
    ch_err='psb_spcnv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_dspasb
