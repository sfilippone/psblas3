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
Module psb_methd_mod

  interface psb_cg
     subroutine psb_dcg(a,prec,b,x,eps,&
	  & desc_a,info,itmax,iter,err,itrace,istop)
       use psb_serial_mod
       use psb_descriptor_type
       use psb_prec_type
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
       use psb_serial_mod
       use psb_descriptor_type
       use psb_prec_type
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
       use psb_serial_mod
       use psb_descriptor_type
       use psb_prec_type
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
  end interface

  interface psb_bicgstabl
    Subroutine psb_dcgstabl(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err, itrace,irst,istop)
      use psb_serial_mod
      use psb_descriptor_type
      Use psb_prec_type
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
      use psb_serial_mod
      use psb_descriptor_type
      Use psb_prec_type
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
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
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
  end interface

end module psb_methd_mod


  
