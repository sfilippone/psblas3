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
!         software without specific prior written permission.
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
!
! File: psb_richardson_mod.f90
!  Interfaces for Richardson iterative methods.
!
!
! Subroutine: psb_zrichardson
! 
!    Front-end for the Richardson iterations, complexversion
!    
! Arguments:
!
!    a      -  type(psb_zspmat_type)      Input: sparse matrix containing A.
!    prec   -  class(psb_zprec_type)       Input: preconditioner
!    b      -  complex,dimension(:)         Input: vector containing the
!                                           right hand side B
!    x      -  complex,dimension(:)         Input/Output: vector containing the
!                                           initial guess and final solution X.
!    eps    -  real                         Input: Stopping tolerance; the iteration is
!                                           stopped when the error
!                                           estimate |err| <= eps
!                                           
!    desc_a -  type(psb_desc_type).       Input: The communication descriptor.
!    info   -  integer.                     Output: Return code
!
!    itmax  -  integer(optional)            Input: maximum number of iterations to be
!                                           performed.
!    iter   -  integer(optional)            Output: how many iterations have been
!                                           performed.
!    err    -  real   (optional)            Output: error estimate on exit
!    itrace -  integer(optional)            Input: print an informational message
!                                           with the error estimate every itrace
!                                           iterations
!    istop  -  integer(optional)            Input: stopping criterion, or how
!                                           to estimate the error. 
!                                           1: err =  |r|/(|a||x|+|b|)
!                                           2: err =  |r|/|b|
!                                           where r is the (preconditioned, recursive
!                                           estimate of) residual 
! 
Subroutine psb_zrichardson_vect(a,prec,b,x,eps,desc_a,info,&
     & itmax,iter,err,itrace,istop)

  use psb_base_mod
  use psb_prec_mod
  use psb_z_linsolve_conv_mod
  use psb_linsolve_mod, psb_protect_name => psb_zrichardson_vect

  Type(psb_zspmat_type), Intent(in)    :: a
  Type(psb_desc_type), Intent(in)      :: desc_a
  class(psb_zprec_type), intent(inout) :: prec 
  type(psb_z_vect_type), Intent(inout) :: b
  type(psb_z_vect_type), Intent(inout) :: x
  Real(psb_dpk_), Intent(in)           :: eps
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace, istop
  integer(psb_ipk_), Optional, Intent(out)       :: iter
  Real(psb_dpk_), Optional, Intent(out) :: err


  logical           :: do_alloc_wrk
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: me,np,err_act
  complex(psb_dpk_), allocatable, target   :: aux(:)
  type(psb_z_vect_type), allocatable, target :: wwrk(:)
  type(psb_z_vect_type), pointer  :: q, p, r, z, w
  real(psb_dpk_) :: derr
  integer(psb_ipk_) :: itmax_, istop_, naux, it, itx, itrace_,&
       &  n_col, n_row,ieg,nspl, istebz
  integer(psb_lpk_) :: mglob
  integer(psb_ipk_) :: debug_level, debug_unit
  type(psb_itconv_type)       :: stopdat
  character(len=20)           :: name
  character(len=*), parameter :: methdname='RICHARDSON'

  info = psb_success_
  name = 'psb_zrichardson'
  call psb_erractionsave(err_act)

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)

  if (present(itrace)) then
    itrace_ = itrace
  else
    itrace_ = -1
  end if

  if (present(istop)) then 
    istop_ = istop 
  else
    istop_ = 2
  endif
  if (present(itmax)) then 
    itmax_ = itmax
  else
    itmax_ = 1000
  endif

  do_alloc_wrk = .not.prec%is_allocated_wrk()
  if (do_alloc_wrk) call prec%allocate_wrk(info,vmold=x%v,desc=desc_a)

  if (.not.allocated(b%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(x%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  mglob = desc_a%get_global_rows()
  n_row = desc_a%get_local_rows()
  n_col = desc_a%get_local_cols()

  call psb_chkvect(mglob,lone,x%get_nrows(),lone,lone,desc_a,info)
  if (info == psb_success_)&
       & call psb_chkvect(mglob,lone,b%get_nrows(),lone,lone,desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_    
    call psb_errpush(info,name,a_err='psb_chkvect on X/B')
    goto 9999
  end if

  naux=4*n_col
  allocate(aux(naux), stat=info)
  if (info == psb_success_) call psb_geall(wwrk,desc_a,info,n=5_psb_ipk_)
  if (info == psb_success_) call psb_geasb(wwrk,desc_a,info,mold=x%v,scratch=.true.)  
  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_non_
    call psb_errpush(info,name)
    goto 9999
  end if
  p  => wwrk(1)
  q  => wwrk(2)
  r  => wwrk(3)
  z  => wwrk(4) 
  w  => wwrk(5)

  call psb_geaxpby(zone,b,zzero,r,desc_a,info)
  if (info == psb_success_) call psb_spmm(-zone,a,x,zone,r,desc_a,info,work=aux)
  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_non_
    call psb_errpush(info,name)
    goto 9999
  end if
  
  
  call psb_init_conv(methdname,istop_,itrace_,itmax_,a,x,b,eps,desc_a,stopdat,info)
  if (info /= psb_success_) Then 
    call psb_errpush(psb_err_from_subroutine_non_,name)
    goto 9999
  End If
  
  loop: do itx=1,itmax_
    call prec%apply(r,z,desc_a,info,work=aux)
    call psb_geaxpby(zone,z,zone,x,desc_a,info)
    call psb_geaxpby(zone,b,zzero,r,desc_a,info)
    if (info == psb_success_) call psb_spmm(-zone,a,x,zone,r,desc_a,info,work=aux)
    if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit loop
  end do loop
  call psb_end_conv(methdname,itx,desc_a,stopdat,info,derr,iter)
  if (present(err)) err = derr

  if (info == psb_success_) call psb_gefree(wwrk,desc_a,info)
  if (info == psb_success_) deallocate(aux,stat=info)
  if ((info==psb_success_).and.do_alloc_wrk) call prec%free_wrk(info)
  
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err=trim(methdname))
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_zrichardson_vect

