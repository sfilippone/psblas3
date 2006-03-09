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
module psb_psblas_mod
  use psb_comm_mod

  interface psb_gedot
    function psb_ddotv(x, y, desc_a,info) 
      use psb_descriptor_type
      real(kind(1.d0))                   :: psb_ddotv
      real(kind(1.d0)), intent(in)       :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end function psb_ddotv
    function psb_ddot(x, y, desc_a, info, jx, jy) 
      use psb_descriptor_type
      real(kind(1.d0))                   :: psb_ddot
      real(kind(1.d0)), intent(in)       :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, optional, intent(in)      :: jx, jy
      integer, intent(out)               :: info
    end function psb_ddot
 end interface
  
  interface psb_gedots
    subroutine  psb_ddotvs(res,x, y, desc_a, info) 
      use psb_descriptor_type
      real(kind(1.d0)), intent(out)      :: res
      real(kind(1.d0)), intent(in)       :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_ddotvs
    subroutine  psb_dmdots(res,x, y, desc_a,info) 
      use psb_descriptor_type
      real(kind(1.d0)), intent(out)      :: res(:)
      real(kind(1.d0)), intent(in)       :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_dmdots
  end interface

  interface psb_geaxpby
     subroutine psb_daxpbyv(alpha, x, beta, y,&
	  & desc_a, info)
       use psb_descriptor_type
       real(kind(1.d0)), intent (in)       ::  x(:)
       real(kind(1.d0)), intent (inout)    ::  y(:)
       real(kind(1.d0)), intent (in)       :: alpha, beta
       type(psb_desc_type), intent (in)    :: desc_a
       integer, intent(out)                :: info
     end subroutine psb_daxpbyv
     subroutine psb_daxpby(alpha, x, beta, y,&
	  & desc_a, info, n, jx, jy)
       use psb_descriptor_type
       real(kind(1.d0)), intent (in)       ::  x(:,:)
       real(kind(1.d0)), intent (inout)    ::  y(:,:)
       real(kind(1.d0)), intent (in)       ::  alpha, beta
       type(psb_desc_type), intent (in)    :: desc_a
       integer, optional :: n, jx, jy
       integer, intent(out)                :: info
     end subroutine psb_daxpby
  end interface

  interface psb_geamax
     function psb_damax(x, desc_a, info, jx)
       use psb_descriptor_type
       real(kind(1.d0))   psb_damax
       real(kind(1.d0)), intent (in)       :: x(:,:)
       type(psb_desc_type), intent (in)    :: desc_a
       integer, optional, intent (in)      :: jx
       integer, intent(out)                :: info
     end function psb_damax
     function psb_damaxv(x, desc_a,info)
       use psb_descriptor_type
       real(kind(1.d0)) psb_damaxv
       real(kind(1.d0)), intent (in)       :: x(:)
       type(psb_desc_type), intent (in)    :: desc_a
       integer, intent(out)                :: info
     end function psb_damaxv
  end interface

  interface psb_geamaxs
     subroutine  psb_damaxvs(res,x,desc_a,info)
       use psb_descriptor_type
       real(kind(1.d0)), intent (out)      :: res
       real(kind(1.d0)), intent (in)       :: x(:)
       type(psb_desc_type), intent (in) :: desc_a
       integer, intent(out)               :: info
     end subroutine psb_damaxvs
     subroutine  psb_dmamax(res,x,desc_a,info,jx)
       use psb_descriptor_type
       real(kind(1.d0)), intent (out)      :: res(:)
       real(kind(1.d0)), intent (in)       :: x(:,:)
       type(psb_desc_type), intent (in)    :: desc_a
       integer, intent(out)                :: info
       integer, optional                   :: jx
     end subroutine psb_dmamax
  end interface

  interface psb_geasum
     function psb_dasum(x, desc_a, info, jx)
       use psb_descriptor_type
       real(kind(1.d0))   psb_dasum
       real(kind(1.d0)), intent (in)       :: x(:,:)
       type(psb_desc_type), intent (in)    :: desc_a
       integer, optional, intent (in)      :: jx
       integer, intent(out)                :: info
     end function psb_dasum
     function psb_dasumv(x, desc_a, info)
       use psb_descriptor_type
       real(kind(1.d0)) psb_dasumv
       real(kind(1.d0)), intent (in)       :: x(:)
       type(psb_desc_type), intent (in)    :: desc_a
       integer, intent(out)                :: info
     end function psb_dasumv
   end interface

  interface psb_geasums
     subroutine  psb_dasumvs(res,x,desc_a,info)
       use psb_descriptor_type
       real(kind(1.d0)), intent (out)      :: res
       real(kind(1.d0)), intent (in)       :: x(:)
       type(psb_desc_type), intent (in)    :: desc_a
       integer, intent(out)                :: info
     end subroutine psb_dasumvs
     subroutine  psb_dmasum(res,x,desc_a,info)
       use psb_descriptor_type
       real(kind(1.d0)), intent (out)      :: res(:)
       real(kind(1.d0)), intent (in)       :: x(:,:)
       type(psb_desc_type), intent (in)    :: desc_a
       integer, intent(out)                :: info
     end subroutine psb_dmasum
  end interface


  interface psb_genrm2
     function psb_dnrm2(x, desc_a, info, jx)
       use psb_descriptor_type
       real(kind(1.d0))   psb_dnrm2
       real(kind(1.d0)), intent (in)       :: x(:,:)
       type(psb_desc_type), intent (in)    :: desc_a
       integer, optional, intent (in)      :: jx
       integer, intent(out)                :: info
     end function psb_dnrm2
     function psb_dnrm2v(x, desc_a, info)
       use psb_descriptor_type
       real(kind(1.d0)) psb_dnrm2v
       real(kind(1.d0)), intent (in)       :: x(:)
       type(psb_desc_type), intent (in)    :: desc_a
       integer, intent(out)                :: info
     end function psb_dnrm2v
  end interface

  interface psb_genrm2s
     subroutine  psb_dnrm2vs(res,x,desc_a,info)
       use psb_descriptor_type
       real(kind(1.d0)), intent (out)      :: res
       real(kind(1.d0)), intent (in)       :: x(:)
       type(psb_desc_type), intent (in)    :: desc_a
       integer, intent(out)                :: info
     end subroutine psb_dnrm2vs
  end interface
  

  interface psb_spnrmi
     function psb_dnrmi(a, desc_a,info)
       use psb_serial_mod
       use psb_descriptor_type
       real(kind(1.d0))                    :: psb_dnrmi
       type (psb_dspmat_type), intent (in) :: a
       type (psb_desc_type), intent (in)   :: desc_a
       integer, intent(out)                :: info
     end function psb_dnrmi
  end interface

  interface psb_spmm
     subroutine psb_dspmm(alpha, a, x, beta, y, desc_a, info,&
          &trans, k, jx, jy,work,doswap)
       use psb_serial_mod
       use psb_descriptor_type
       type (psb_dspmat_type), intent(in)   :: a
       real(kind(1.d0)), intent(inout)      :: x(:,:)
       real(kind(1.d0)), intent(inout)      :: y(:,:)
       real(kind(1.d0)), intent(in)         :: alpha, beta
       type(psb_desc_type), intent(in)      :: desc_a
       character, optional, intent(in)      :: trans
       real(kind(1.d0)), optional, intent(inout) :: work(:)
       integer, optional, intent(in)        :: k, jx, jy,doswap
       integer, intent(out)                 :: info
     end subroutine psb_dspmm
     subroutine psb_dspmv(alpha, a, x, beta, y,&
	  & desc_a, info, trans, work,doswap)
       use psb_serial_mod
       use psb_descriptor_type
       type (psb_dspmat_type), intent(in)   :: a
       real(kind(1.d0)), intent(inout)      :: x(:)
       real(kind(1.d0)), intent(inout)      :: y(:)
       real(kind(1.d0)), intent(in)         :: alpha, beta
       type(psb_desc_type), intent(in)      :: desc_a
       character, optional, intent(in)      :: trans
       real(kind(1.d0)), optional, intent(inout) :: work(:)
       integer, optional, intent(in)        :: doswap
       integer, intent(out)                 :: info
     end subroutine psb_dspmv
  end interface

  interface psb_spsm
     subroutine psb_dspsm(alpha, t, x, beta, y,&
	  & desc_a, info, trans, unit, choice,& 
	  & diag, n, jx, jy, work)
       use psb_serial_mod
       use psb_descriptor_type
       type (psb_dspmat_type), intent(in)     :: t
       real(kind(1.d0)), intent(in)           :: x(:,:)
       real(kind(1.d0)), intent(inout)        :: y(:,:)
       real(kind(1.d0)), intent(in)           :: alpha, beta
       type(psb_desc_type), intent(in)     :: desc_a
       character, optional, intent(in)        :: trans, unit
       integer, optional, intent(in)          :: n, jx, jy
       integer, optional, intent(in)          :: choice
       real(kind(1.d0)), optional, intent(in) :: work(:), diag(:)
       integer, intent(out)               :: info
     end subroutine psb_dspsm
     subroutine psb_dspsv(alpha, t, x, beta, y,&
	  & desc_a, info, trans, unit, choice,& 
	  & diag, work)
       use psb_serial_mod
       use psb_descriptor_type
       type (psb_dspmat_type), intent(in)     :: t
       real(kind(1.d0)), intent(in)           :: x(:)
       real(kind(1.d0)), intent(inout)        :: y(:)
       real(kind(1.d0)), intent(in)           :: alpha, beta
       type(psb_desc_type), intent(in)        :: desc_a
       character, optional, intent(in)        :: trans, unit
       integer, optional, intent(in)          :: choice
       real(kind(1.d0)), optional, intent(in) :: work(:), diag(:)
       integer, intent(out)                   :: info
     end subroutine psb_dspsv
  end interface


!   interface psb_gelp
!      subroutine psb_dgelp(trans,iperm,x,desc_a,info)
!        use psb_descriptor_type
!        type(psb_desc_type), intent(in)      ::  desc_a
!        real(kind(1.d0)), intent(inout)      ::  x(:,:)
!        integer, intent(inout)               ::  iperm(:),info
!        character, intent(in)                ::  trans
!      end subroutine psb_dgelp
!      subroutine psb_dgelpv(trans,iperm,x,desc_a,info)
!        use psb_descriptor_type
!        type(psb_desc_type), intent(in)      ::  desc_a
!        real(kind(1.d0)), intent(inout)      ::  x(:)
!        integer, intent(inout)               ::  iperm(:),info
!        character, intent(in)                ::  trans
!      end subroutine psb_dgelpv
!   end interface
     
end module psb_psblas_mod
