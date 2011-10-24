!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
module psb_i_comm_mod

  interface psb_ovrl
    subroutine  psb_iovrlm(x,desc_a,info,jx,ik,work,update,mode)
      use psb_descriptor_type
      integer,          intent(inout), target :: x(:,:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      integer, intent(inout), optional, target  :: work(:)
      integer, intent(in), optional           :: update,jx,ik,mode
    end subroutine psb_iovrlm
    subroutine  psb_iovrlv(x,desc_a,info,work,update,mode)
      use psb_descriptor_type
      integer, intent(inout), target          :: x(:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      integer, intent(inout), optional, target :: work(:)
      integer, intent(in), optional           :: update,mode
    end subroutine psb_iovrlv
  end interface

  interface psb_halo
    subroutine  psb_ihalom(x,desc_a,info,alpha,jx,ik,work,tran,mode,data)
      use psb_descriptor_type
      integer, intent(inout), target         :: x(:,:)
      type(psb_desc_type), intent(in)        :: desc_a
      integer, intent(out)                   :: info
      real(psb_dpk_), intent(in), optional   :: alpha
      integer, intent(inout), optional, target  :: work(:)
      integer, intent(in), optional          :: mode,jx,ik,data
      character, intent(in), optional        :: tran
    end subroutine psb_ihalom
    subroutine  psb_ihalov(x,desc_a,info,alpha,work,tran,mode,data)
      use psb_descriptor_type
      integer, intent(inout)                 :: x(:)
      type(psb_desc_type), intent(in)        :: desc_a
      integer, intent(out)                   :: info
      real(psb_dpk_), intent(in), optional   :: alpha
      integer, intent(inout), optional, target :: work(:)
      integer, intent(in), optional          :: mode,data
      character, intent(in), optional        :: tran
    end subroutine psb_ihalov
  end interface


  interface psb_scatter
    subroutine  psb_iscatterm(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      integer, intent(out)             :: locx(:,:)
      integer, intent(in)              :: globx(:,:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer, intent(out)             :: info
      integer, intent(in), optional    :: root
    end subroutine psb_iscatterm
    subroutine  psb_iscatterv(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      integer, intent(out)             :: locx(:)
      integer, intent(in)              :: globx(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer, intent(out)             :: info
      integer, intent(in), optional    :: root
    end subroutine psb_iscatterv
  end interface

  interface psb_gather
    subroutine  psb_igatherm(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      integer, intent(in)             :: locx(:,:)
      integer, intent(out)            :: globx(:,:)
      type(psb_desc_type), intent(in) :: desc_a
      integer, intent(out)            :: info
      integer, intent(in), optional   :: root
    end subroutine psb_igatherm
    subroutine  psb_igatherv(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      integer, intent(in)             :: locx(:)
      integer, intent(out)            :: globx(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer, intent(out)            :: info
      integer, intent(in), optional   :: root
    end subroutine psb_igatherv
  end interface
  
end module psb_i_comm_mod
