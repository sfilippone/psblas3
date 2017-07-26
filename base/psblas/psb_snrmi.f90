!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006, 2010, 2015, 2017
!        Salvatore Filippone    Cranfield University
!        Alfredo Buttari        CNRS-IRIT, Toulouse
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
! File: psb_snrmi.f90
!
! Function: psb_snrmi
!    Forms the infinity norm of a sparse matrix,                           
!
!    nrmi := max_i(abs(sum(A(i,:))))                                              
!
! Arguments:
!    a      -  type(psb_sspmat_type).   The sparse matrix containing A.
!    desc_a -  type(psb_desc_type).     The communication descriptor.
!    info   -  integer.                   Return code
!
function psb_snrmi(a,desc_a,info)  result(res)
  use psb_base_mod, psb_protect_name => psb_snrmi
  implicit none

  type(psb_sspmat_type), intent(in) :: a
  integer(psb_ipk_), intent(out)    :: info
  type(psb_desc_type), intent(in)   :: desc_a
  real(psb_spk_)                    :: res

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, n, iia, jja, ia, ja, mdim, ndim, m
  character(len=20)        :: name, ch_err

  name='psb_snrmi'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ia = 1
  ja = 1
  m = desc_a%get_global_rows()
  n = desc_a%get_global_cols()

  call psb_chkmat(m,n,ia,ja,desc_a,info,iia,jja)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkmat'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if ((iia /= 1).or.(jja /= 1)) then
    info=psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info,name)
    goto 9999
  end if

  if ((m /= 0).and.(n /= 0)) then
    res = a%spnmi()
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_snrmi'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  else
    res = szero
  end if

  ! compute global max
  call psb_amx(ictxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return
end function psb_snrmi
