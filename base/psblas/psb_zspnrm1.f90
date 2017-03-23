!   
!                Parallel Sparse BLAS  version 3.4
!      (C) Copyright 2006, 2010, 2015
!                         Salvatore Filippone    University of Rome Tor Vergata
!                         Alfredo Buttari        CNRS-IRIT, Toulouse
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
! File: psb_znrm1.f90
!
! Function: psb_znrm1
!    Forms the  norm1 of a sparse matrix,       
!
!    norm1 := max_j(sum(abs(A(:,j))))                                                 
!
! Arguments:
!    a      -  type(psb_zspmat_type).   The sparse matrix containing A.
!    desc_a -  type(psb_desc_type).     The communication descriptor.
!    info   -  integer.                   Return code
!
function psb_zspnrm1(a,desc_a,info)  result(res)
  use psb_base_mod, psb_protect_name => psb_zspnrm1
  implicit none

  type(psb_zspmat_type), intent(in) :: a
  integer(psb_ipk_), intent(out)      :: info
  type(psb_desc_type), intent(in)     :: desc_a
  real(psb_dpk_)                      :: res

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, nr,nc,&
       & err_act, n, iia, jja, ia, ja, mdim, ndim, m
  character(len=20)      :: name, ch_err
  real(psb_dpk_), allocatable :: v(:)

  name='psb_zspnrm1'
  if (psb_errstatus_fatal()) return 
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
  nr = desc_a%get_local_rows()
  nc = desc_a%get_local_cols()

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

!!$  call psb_geall(v,desc_a,info)
!!$  if(info == psb_success_) then 
!!$    v = zzero
!!$    call psb_geasb(v,desc_a,info)
!!$  end if
!!$  if(info /= psb_success_) then
!!$    info=psb_err_from_subroutine_
!!$    ch_err='geall/asb'
!!$    call psb_errpush(info,name,a_err=ch_err)
!!$    goto 9999
!!$  end if

  if ((m /= 0).and.(n /= 0)) then
    v = a%aclsum(info)
    if (info == psb_success_) &
         & call psb_realloc(desc_a%get_local_cols(),v,info,pad=dzero)
    if (info == psb_success_) call psb_halo(v,desc_a,info,tran='T')
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_halo'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    res = maxval(v(1:nr))
  else
    res = dzero 
  end if
  ! compute global max
  call psb_amx(ictxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return
end function psb_zspnrm1
