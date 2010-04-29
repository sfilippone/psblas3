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
!
! File: psb_krylov_mod.f90
!  Interfaces for Krylov subspace iterative methods.
!
Module psb_base_inner_krylov_mod

  use psb_const_mod

  interface psb_end_conv
    module procedure psb_d_end_conv
  end interface

  integer, parameter :: psb_ik_bni_=1, psb_ik_rni_=2, psb_ik_ani_=3
  integer, parameter :: psb_ik_xni_=4, psb_ik_bn2_=5, psb_ik_xn2_=6
  integer, parameter :: psb_ik_errnum_=7, psb_ik_errden_=8, psb_ik_eps_=9, psb_ik_rn2_=10
  integer, parameter :: psb_ik_stopc_=1, psb_ik_trace_=2, psb_ik_itmax_=3
  integer, parameter :: psb_ik_ivsz_=16
  type psb_itconv_type
    integer        :: controls(psb_ik_ivsz_)
    real(psb_dpk_) :: values(psb_ik_ivsz_)
  end type psb_itconv_type

contains

  subroutine log_header(methdname)
    !use psb_base_mod
    implicit none 
    character(len=*), intent(in)   :: methdname
    character(len=*), parameter    :: fmt='(a18,1x,a4,3(2x,a15))'
    integer, parameter             :: outlen=18 
    character(len=len(methdname))  :: mname
    character(len=outlen)          :: outname

    mname = adjustl(trim(methdname))
    write(outname,'(a)') mname(1:min(len_trim(mname),outlen-1))//':'
    write(*,fmt) adjustl(outname),'Iteration','Error Estimate','Tolerance'

  end subroutine log_header


  subroutine log_conv(methdname,me,itx,itrace,errnum,errden,eps)
    !use psb_base_mod
    implicit none 
    character(len=*), intent(in)  :: methdname
    integer, intent(in)           :: me, itx, itrace
    real(psb_dpk_), intent(in)    :: errnum, errden, eps
    character(len=*), parameter   :: fmt='(a18,1x,i4,3(2x,es15.9))'
    integer, parameter            :: outlen=18 
    character(len=len(methdname)) :: mname
    character(len=outlen)         :: outname

    if ((mod(itx,itrace) == 0).and.(me == 0)) then 
      mname = adjustl(trim(methdname))
      write(outname,'(a)') mname(1:min(len_trim(mname),outlen-1))//':'
      if (errden > dzero ) then 
        write(*,fmt) adjustl(outname),itx,errnum/errden,eps
      else
        write(*,fmt) adjustl(outname),itx,errnum,eps
      end if
    endif

  end subroutine log_conv

  subroutine log_end(methdname,me,it,errnum,errden,eps,err,iter)
    !use psb_base_mod
    implicit none 
    character(len=*), intent(in) :: methdname
    integer, intent(in)          :: me, it
    real(psb_dpk_), intent(in) :: errnum, errden, eps
    real(psb_dpk_), optional, intent(out) :: err
    integer, optional, intent(out)  :: iter

    character(len=*), parameter  :: fmt='(a,2x,es15.9,1x,a,1x,i4,1x,a)'
    character(len=*), parameter  :: fmt1='(a,3(2x,es15.9))'

    if (errden == dzero) then 
      if (errnum > eps) then         
        if (me == 0) then 
          write(*,fmt) trim(methdname)//' failed to converge to ',eps,&
               & ' in ',it,' iterations. '
          write(*,fmt1) 'Last iteration error estimate: ',&
               & errnum
        end if
      end if
      if (present(err)) err=errnum
    else
      if (errnum/errden > eps) then         
        if (me == 0) then 
          write(*,fmt) trim(methdname)//' failed to converge to ',eps,&
               & ' in ',it,' iterations. '
          write(*,fmt1) 'Last iteration error estimate: ',&
               & errnum/errden
        end if
      endif
      if (present(err)) err=errnum/errden
    end if
    if (present(iter)) iter = it

  end subroutine log_end

  
  subroutine psb_d_end_conv(methdname,it,desc_a,stopdat,info,err,iter)
    use psb_sparse_mod
    implicit none 
    character(len=*), intent(in)    :: methdname
    integer, intent(in)             :: it
    type(psb_desc_type), intent(in) :: desc_a
    type(psb_itconv_type)           :: stopdat
    integer, intent(out)            :: info
    real(psb_dpk_), optional, intent(out) :: err
    integer, optional, intent(out)  :: iter

    integer                         :: ictxt, me, np, err_act
    real(psb_dpk_)                  :: errnum, errden, eps
    character(len=20)               :: name

    info = psb_success_
    name = 'psb_end_conv'

    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt,me,np)

    errnum = stopdat%values(psb_ik_errnum_) 
    errden = stopdat%values(psb_ik_errden_) 
    eps    = stopdat%values(psb_ik_eps_) 
    call log_end(methdname,me,it,errnum,errden,eps,err,iter)

  end subroutine psb_d_end_conv

end module psb_base_inner_krylov_mod
