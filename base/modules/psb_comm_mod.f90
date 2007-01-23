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
module psb_comm_mod

  interface psb_ovrl
     subroutine  psb_dovrlm(x,desc_a,info,jx,ik,work,update)
       use psb_descriptor_type
       real(kind(1.d0)), intent(inout)           :: x(:,:)
       type(psb_desc_type), intent(in)           :: desc_a
       integer, intent(out)                      :: info
       real(kind(1.d0)), intent(inout), optional :: work(:)
       integer, intent(in), optional             :: update,jx,ik
     end subroutine psb_dovrlm
     subroutine  psb_dovrlv(x,desc_a,info,work,update)
       use psb_descriptor_type
       real(kind(1.d0)), intent(inout)           :: x(:)
       type(psb_desc_type), intent(in)           :: desc_a
       integer, intent(out)                      :: info
       real(kind(1.d0)), intent(inout), optional :: work(:)
       integer, intent(in), optional             :: update
     end subroutine psb_dovrlv
     subroutine  psb_zovrlm(x,desc_a,info,jx,ik,work,update)
       use psb_descriptor_type
       complex(kind(1.d0)), intent(inout)           :: x(:,:)
       type(psb_desc_type), intent(in)           :: desc_a
       integer, intent(out)                      :: info
       complex(kind(1.d0)), intent(inout), optional :: work(:)
       integer, intent(in), optional             :: update,jx,ik
     end subroutine psb_zovrlm
     subroutine  psb_zovrlv(x,desc_a,info,work,update)
       use psb_descriptor_type
       complex(kind(1.d0)), intent(inout)           :: x(:)
       type(psb_desc_type), intent(in)           :: desc_a
       integer, intent(out)                      :: info
       complex(kind(1.d0)), intent(inout), optional :: work(:)
       integer, intent(in), optional             :: update
     end subroutine psb_zovrlv
  end interface

  interface psb_halo
     subroutine  psb_dhalom(x,desc_a,info,alpha,jx,ik,work,tran,mode,data)
       use psb_descriptor_type
       real(kind(1.d0)), intent(inout)           :: x(:,:)
       type(psb_desc_type), intent(in)           :: desc_a
       integer, intent(out)                      :: info
       real(kind(1.d0)), intent(in), optional    :: alpha
       real(kind(1.d0)), target, optional        :: work(:)
       integer, intent(in), optional             :: mode,jx,ik,data
       character, intent(in), optional           :: tran
     end subroutine psb_dhalom
     subroutine  psb_dhalov(x,desc_a,info,alpha,work,tran,mode,data)
       use psb_descriptor_type
       real(kind(1.d0)), intent(inout)           :: x(:)
       type(psb_desc_type), intent(in)           :: desc_a
       integer, intent(out)                      :: info
       real(kind(1.d0)), intent(in), optional    :: alpha
       real(kind(1.d0)), target, optional :: work(:)
       integer, intent(in), optional             :: mode,data
       character, intent(in), optional           :: tran
     end subroutine psb_dhalov
     subroutine  psb_ihalom(x,desc_a,info,alpha,jx,ik,work,tran,mode)
       use psb_descriptor_type
       integer, intent(inout) :: x(:,:)
       type(psb_desc_type), intent(in)        :: desc_a
       integer, intent(out)                   :: info
       real(kind(1.d0)), intent(in), optional :: alpha
       integer, intent(inout), optional       :: work(:)
       integer, intent(in), optional          :: mode,jx,ik
       character, intent(in), optional        :: tran
     end subroutine psb_ihalom
     subroutine  psb_ihalov(x,desc_a,info,alpha,work,tran,mode)
       use psb_descriptor_type
       integer, intent(inout)                 :: x(:)
       type(psb_desc_type), intent(in)        :: desc_a
       integer, intent(out)                   :: info
       real(kind(1.d0)), intent(in), optional :: alpha
       integer, intent(inout), optional       :: work(:)
       integer, intent(in), optional          :: mode
       character, intent(in), optional        :: tran
     end subroutine psb_ihalov
     subroutine  psb_zhalom(x,desc_a,info,alpha,jx,ik,work,tran,mode)
       use psb_descriptor_type
       complex(kind(1.d0)), intent(inout)           :: x(:,:)
       type(psb_desc_type), intent(in)           :: desc_a
       integer, intent(out)                      :: info
       complex(kind(1.d0)), intent(in), optional    :: alpha
       complex(kind(1.d0)), target, optional        :: work(:)
       integer, intent(in), optional             :: mode,jx,ik
       character, intent(in), optional           :: tran
     end subroutine psb_zhalom
     subroutine  psb_zhalov(x,desc_a,info,alpha,work,tran,mode)
       use psb_descriptor_type
       complex(kind(1.d0)), intent(inout)           :: x(:)
       type(psb_desc_type), intent(in)           :: desc_a
       integer, intent(out)                      :: info
       complex(kind(1.d0)), intent(in), optional    :: alpha
       complex(kind(1.d0)), target, optional :: work(:)
       integer, intent(in), optional             :: mode
       character, intent(in), optional           :: tran
     end subroutine psb_zhalov
  end interface


  interface psb_dscatter
     subroutine  psb_dscatterm(globx, locx, desc_a, info, iroot,&
          & iiglobx, ijglobx, iilocx,ijlocx,ik)
       use psb_descriptor_type
       real(kind(1.d0)), intent(out)    :: locx(:,:)
       real(kind(1.d0)), intent(in)     :: globx(:,:)
       type(psb_desc_type), intent(in)  :: desc_a
       integer, intent(out)             :: info
       integer, intent(in), optional    :: iroot,iiglobx,&
            & ijglobx,iilocx,ijlocx,ik
     end subroutine psb_dscatterm
     subroutine  psb_dscatterv(globx, locx, desc_a, info, iroot)
       use psb_descriptor_type
       real(kind(1.d0)), intent(out)    :: locx(:)
       real(kind(1.d0)), intent(in)     :: globx(:)
       type(psb_desc_type), intent(in)  :: desc_a
       integer, intent(out)             :: info
       integer, intent(in), optional    :: iroot
     end subroutine psb_dscatterv
     subroutine  psb_zscatterm(globx, locx, desc_a, info, iroot,&
          & iiglobx, ijglobx, iilocx,ijlocx,ik)
       use psb_descriptor_type
       complex(kind(1.d0)), intent(out)    :: locx(:,:)
       complex(kind(1.d0)), intent(in)     :: globx(:,:)
       type(psb_desc_type), intent(in)  :: desc_a
       integer, intent(out)             :: info
       integer, intent(in), optional    :: iroot,iiglobx,&
            & ijglobx,iilocx,ijlocx,ik
     end subroutine psb_zscatterm
     subroutine  psb_zscatterv(globx, locx, desc_a, info, iroot)
       use psb_descriptor_type
       complex(kind(1.d0)), intent(out)    :: locx(:)
       complex(kind(1.d0)), intent(in)     :: globx(:)
       type(psb_desc_type), intent(in)  :: desc_a
       integer, intent(out)             :: info
       integer, intent(in), optional    :: iroot
     end subroutine psb_zscatterv
  end interface

  interface psb_gather
     subroutine  psb_dgatherm(globx, locx, desc_a, info, iroot,&
          & iiglobx, ijglobx, iilocx,ijlocx,ik)
       use psb_descriptor_type
       real(kind(1.d0)), intent(in)    :: locx(:,:)
       real(kind(1.d0)), intent(out)   :: globx(:,:)
       type(psb_desc_type), intent(in) :: desc_a
       integer, intent(out)            :: info
       integer, intent(in), optional   :: iroot, iiglobx, ijglobx, iilocx, ijlocx, ik
     end subroutine psb_dgatherm
     subroutine  psb_dgatherv(globx, locx, desc_a, info, iroot,&
          & iiglobx, iilocx)
       use psb_descriptor_type
       real(kind(1.d0)), intent(in)    :: locx(:)
       real(kind(1.d0)), intent(out)   :: globx(:)
       type(psb_desc_type), intent(in) :: desc_a
       integer, intent(out)            :: info
       integer, intent(in), optional   :: iroot, iiglobx, iilocx
     end subroutine psb_dgatherv
     subroutine  psb_zgatherm(globx, locx, desc_a, info, iroot,&
          & iiglobx, ijglobx, iilocx,ijlocx,ik)
       use psb_descriptor_type
       complex(kind(1.d0)), intent(in)    :: locx(:,:)
       complex(kind(1.d0)), intent(out)   :: globx(:,:)
       type(psb_desc_type), intent(in) :: desc_a
       integer, intent(out)            :: info
       integer, intent(in), optional   :: iroot, iiglobx, ijglobx, iilocx, ijlocx, ik
     end subroutine psb_zgatherm
     subroutine  psb_zgatherv(globx, locx, desc_a, info, iroot,&
          & iiglobx, iilocx)
       use psb_descriptor_type
       complex(kind(1.d0)), intent(in)    :: locx(:)
       complex(kind(1.d0)), intent(out)   :: globx(:)
       type(psb_desc_type), intent(in) :: desc_a
       integer, intent(out)            :: info
       integer, intent(in), optional   :: iroot, iiglobx, iilocx
     end subroutine psb_zgatherv
  end interface
  
end module psb_comm_mod
