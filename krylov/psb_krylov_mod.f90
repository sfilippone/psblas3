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
Module psb_krylov_mod

  use psb_base_mod
  use psb_prec_mod


  interface psb_krylov
    module procedure psb_dkrylov, psb_zkrylov
  end interface

  interface psb_cg
     subroutine psb_dcg(a,prec,b,x,eps,&
	  & desc_a,info,itmax,iter,err,itrace,istop)
       use psb_base_mod
       use psb_prec_mod
       type(psb_dspmat_type), intent(in)  :: a
       type(psb_desc_type), intent(in)    :: desc_a
       real(kind(1.d0)), intent(in)       :: b(:)
       real(kind(1.d0)), intent(inout)    :: x(:)
       real(kind(1.d0)), intent(in)       :: eps
       type(psb_dprec_type), intent(in)   :: prec
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: itmax, itrace,istop
       integer, optional, intent(out)     :: iter
       real(kind(1.d0)), optional, intent(out) :: err
     end subroutine psb_dcg
  end interface

  interface psb_bicg
     subroutine psb_dbicg(a,prec,b,x,eps,&
	  & desc_a,info,itmax,iter,err,itrace,istop)
       use psb_base_mod
       use psb_prec_mod
       type(psb_dspmat_type), intent(in)  :: a
       type(psb_desc_type), intent(in)    :: desc_a
       real(kind(1.d0)), intent(in)       :: b(:)
       real(kind(1.d0)), intent(inout)    :: x(:)
       real(kind(1.d0)), intent(in)       :: eps
       type(psb_dprec_type), intent(in)   :: prec
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: itmax, itrace,istop
       integer, optional, intent(out)     :: iter
       real(kind(1.d0)), optional, intent(out) :: err
     end subroutine psb_dbicg
  end interface

  interface psb_bicgstab
     subroutine psb_dcgstab(a,prec,b,x,eps,&
	  & desc_a,info,itmax,iter,err,itrace,istop)
       use psb_base_mod
       use psb_prec_mod
       type(psb_dspmat_type), intent(in)  :: a
       type(psb_desc_type), intent(in)    :: desc_a
       real(kind(1.d0)), intent(in)       :: b(:)
       real(kind(1.d0)), intent(inout)    :: x(:)
       real(kind(1.d0)), intent(in)       :: eps
       type(psb_dprec_type), intent(in)   :: prec
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: itmax, itrace,istop
       integer, optional, intent(out)     :: iter
       real(kind(1.d0)), optional, intent(out) :: err
     end subroutine psb_dcgstab
     subroutine psb_zcgstab(a,prec,b,x,eps,&
	  & desc_a,info,itmax,iter,err,itrace,istop)
       use psb_base_mod
       use psb_prec_mod
       type(psb_zspmat_type), intent(in)  :: a
       type(psb_desc_type), intent(in)    :: desc_a
       complex(kind(1.d0)), intent(in)       :: b(:)
       complex(kind(1.d0)), intent(inout)    :: x(:)
       real(kind(1.d0)), intent(in)       :: eps
       type(psb_zprec_type), intent(in)   :: prec
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: itmax, itrace,istop
       integer, optional, intent(out)     :: iter
       real(kind(1.d0)), optional, intent(out) :: err
     end subroutine psb_zcgstab
  end interface

  interface psb_bicgstabl
    Subroutine psb_dcgstabl(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err, itrace,irst,istop)
      use psb_base_mod
      use psb_prec_mod
      Type(psb_dspmat_type), Intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      type(psb_dprec_type), intent(in)   :: prec
      Real(Kind(1.d0)), Intent(in)       :: b(:)
      Real(Kind(1.d0)), Intent(inout)    :: x(:)
      Real(Kind(1.d0)), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(Kind(1.d0)), Optional, Intent(out) :: err
    end subroutine psb_dcgstabl
  end interface

  interface psb_rgmres
    Subroutine psb_dgmresr(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
      use psb_base_mod
      use psb_prec_mod
!!$  parameters 
      Type(psb_dspmat_type), Intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      type(psb_dprec_type), intent(in)   :: prec 
      Real(Kind(1.d0)), Intent(in)       :: b(:)
      Real(Kind(1.d0)), Intent(inout)    :: x(:)
      Real(Kind(1.d0)), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(Kind(1.d0)), Optional, Intent(out) :: err
    end subroutine psb_dgmresr
  end interface

  interface psb_cgs
    subroutine psb_dcgs(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
      use psb_base_mod
      use psb_prec_mod
!!$  parameters 
      type(psb_dspmat_type), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a 
      type(psb_dprec_type), intent(in)   :: prec 
      real(kind(1.d0)), intent(in)       :: b(:)
      real(kind(1.d0)), intent(inout)    :: x(:)
      real(kind(1.d0)), intent(in)       :: eps
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(kind(1.d0)), optional, intent(out) :: err
    end subroutine psb_dcgs
     subroutine psb_zcgs(a,prec,b,x,eps,&
	  & desc_a,info,itmax,iter,err,itrace,istop)
       use psb_base_mod
       use psb_prec_mod
       type(psb_zspmat_type), intent(in)  :: a
       type(psb_desc_type), intent(in)    :: desc_a
       complex(kind(1.d0)), intent(in)       :: b(:)
       complex(kind(1.d0)), intent(inout)    :: x(:)
       real(kind(1.d0)), intent(in)       :: eps
       type(psb_zprec_type), intent(in)   :: prec
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: itmax, itrace,istop
       integer, optional, intent(out)     :: iter
       real(kind(1.d0)), optional, intent(out) :: err
     end subroutine psb_zcgs
  end interface
  
contains


  Subroutine psb_dkrylov(method,a,prec,b,x,eps,desc_a,info,&
       &itmax,iter,err,itrace,irst,istop)
    character(len=*)                   :: method
    Type(psb_dspmat_type), Intent(in)  :: a
    Type(psb_desc_type), Intent(in)    :: desc_a
    type(psb_dprec_type), intent(in)   :: prec 
    Real(Kind(1.d0)), Intent(in)       :: b(:)
    Real(Kind(1.d0)), Intent(inout)    :: x(:)
    Real(Kind(1.d0)), Intent(in)       :: eps
    integer, intent(out)               :: info
    Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
    Integer, Optional, Intent(out)     :: iter
    Real(Kind(1.d0)), Optional, Intent(out) :: err
    
    integer             :: itmax_, itrace_, irst_, istop_, iter_
    real(kind(1.d0))    :: err_
    
    if (present(itmax)) then 
      itmax_ = itmax
    else
      itmax_ = 1000
    end if
    
    if (present(itrace)) then 
      itrace_ = itrace
    else
      itrace_ = -1
    end if
    
    if (present(irst)) then 
      irst_ = irst 
    else 
      irst_ = 1
    end if
    
    if (present(istop)) then 
      istop_ = istop 
    else
      istop_ = 1
    end if


    select case(toupper(method))
    case('CG') 
      call  psb_cg(a,prec,b,x,eps,desc_a,info,&
         &itmax_,iter_,err_,itrace_,istop_)
    case('CGS') 
      call  psb_cgs(a,prec,b,x,eps,desc_a,info,&
         &itmax_,iter_,err_,itrace_,istop_)
    case('BICG') 
      call  psb_bicg(a,prec,b,x,eps,desc_a,info,&
         &itmax_,iter_,err_,itrace_,istop_)
    case('BICGSTAB') 
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
         &itmax_,iter_,err_,itrace_,istop_)
    case('RGMRES')
      call  psb_rgmres(a,prec,b,x,eps,desc_a,info,&
         &itmax_,iter_,err_,itrace_,irst_,istop_)
    case('BICGSTABL')
      call  psb_bicgstabl(a,prec,b,x,eps,desc_a,info,&
         &itmax_,iter_,err_,itrace_,irst_,istop_)
    case default
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
         &itmax_,iter_,err_,itrace_,istop_)
    end select
    
    if (present(err)) then 
      err = err_
    endif
    
    if (present(iter)) then 
      iter = iter_
    endif

  end subroutine psb_dkrylov


  Subroutine psb_zkrylov(method,a,prec,b,x,eps,desc_a,info,&
       &itmax,iter,err,itrace,irst,istop)
    character(len=*)                   :: method
    Type(psb_zspmat_type), Intent(in)  :: a
    Type(psb_desc_type), Intent(in)    :: desc_a
    type(psb_zprec_type), intent(in)   :: prec 
    complex(Kind(1.d0)), Intent(in)    :: b(:)
    complex(Kind(1.d0)), Intent(inout) :: x(:)
    Real(Kind(1.d0)), Intent(in)       :: eps
    integer, intent(out)               :: info
    Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
    Integer, Optional, Intent(out)     :: iter
    Real(Kind(1.d0)), Optional, Intent(out) :: err
    
    integer             :: itmax_, itrace_, irst_, istop_, iter_
    real(kind(1.d0))    :: err_
    
    if (present(itmax)) then 
      itmax_ = itmax
    else
      itmax_ = 1000
    end if
    
    if (present(itrace)) then 
      itrace_ = itrace
    else
      itrace_ = -1
    end if
    
    if (present(irst)) then 
      irst_ = irst 
    else 
      irst_ = 1
    end if
    
    if (present(istop)) then 
      istop_ = istop 
    else
      istop_ = 1
    end if


    select case(toupper(method))
!!$    case('CG') 
!!$      call  psb_cg(a,prec,b,x,eps,desc_a,info,&
!!$         &itmax_,iter_,err_,itrace_,istop_)
    case('CGS') 
      call  psb_cgs(a,prec,b,x,eps,desc_a,info,&
         &itmax_,iter_,err_,itrace_,istop_)
!!$    case('BICG') 
!!$      call  psb_bicg(a,prec,b,x,eps,desc_a,info,&
!!$         &itmax_,iter_,err_,itrace_,istop_)
    case('BICGSTAB') 
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
         &itmax_,iter_,err_,itrace_,istop_)
!!$    case('RGMRES')
!!$      call  psb_rgmres(a,prec,b,x,eps,desc_a,info,&
!!$         &itmax_,iter_,err_,itrace_,irst_,istop_)
!!$    case('BICGSTABL')
!!$      call  psb_bicgstabl(a,prec,b,x,eps,desc_a,info,&
!!$         &itmax_,iter_,err_,itrace_,irst_,istop_)
    case default
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
         &itmax_,iter_,err_,itrace_,istop_)
    end select
    
    if (present(err)) then 
      err = err_
    endif
    
    if (present(iter)) then 
      iter = iter_
    endif

  end subroutine psb_zkrylov



end module psb_krylov_mod


  
