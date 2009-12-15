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
module psb_comm_mod

  interface psb_ovrl
    subroutine  psb_sovrlm(x,desc_a,info,jx,ik,work,update,mode)
      use psb_descriptor_type
      real(psb_spk_), intent(inout)           :: x(:,:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      real(psb_spk_), intent(inout), optional, target :: work(:)
      integer, intent(in), optional           :: update,jx,ik,mode
    end subroutine psb_sovrlm
    subroutine  psb_sovrlv(x,desc_a,info,work,update,mode)
      use psb_descriptor_type
      real(psb_spk_), intent(inout)           :: x(:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      real(psb_spk_), intent(inout), optional, target :: work(:)
      integer, intent(in), optional           :: update,mode
    end subroutine psb_sovrlv
    subroutine  psb_dovrlm(x,desc_a,info,jx,ik,work,update,mode)
      use psb_descriptor_type
      real(psb_dpk_), intent(inout)           :: x(:,:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      real(psb_dpk_), intent(inout), optional, target :: work(:)
      integer, intent(in), optional           :: update,jx,ik,mode
    end subroutine psb_dovrlm
    subroutine  psb_dovrlv(x,desc_a,info,work,update,mode)
      use psb_descriptor_type
      real(psb_dpk_), intent(inout)           :: x(:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      real(psb_dpk_), intent(inout), optional, target :: work(:)
      integer, intent(in), optional           :: update,mode
    end subroutine psb_dovrlv
    subroutine  psb_iovrlm(x,desc_a,info,jx,ik,work,update,mode)
      use psb_descriptor_type
      integer,          intent(inout)         :: x(:,:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      integer, intent(inout), optional, target  :: work(:)
      integer, intent(in), optional           :: update,jx,ik,mode
    end subroutine psb_iovrlm
    subroutine  psb_iovrlv(x,desc_a,info,work,update,mode)
      use psb_descriptor_type
      integer, intent(inout)                  :: x(:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      integer, intent(inout), optional, target :: work(:)
      integer, intent(in), optional           :: update,mode
    end subroutine psb_iovrlv
    subroutine  psb_covrlm(x,desc_a,info,jx,ik,work,update,mode)
      use psb_descriptor_type
      complex(psb_spk_), intent(inout)           :: x(:,:)
      type(psb_desc_type), intent(in)            :: desc_a
      integer, intent(out)                       :: info
      complex(psb_spk_), intent(inout), optional, target :: work(:)
      integer, intent(in), optional              :: update,jx,ik,mode
    end subroutine psb_covrlm
    subroutine  psb_covrlv(x,desc_a,info,work,update,mode)
      use psb_descriptor_type
      complex(psb_spk_), intent(inout)           :: x(:)
      type(psb_desc_type), intent(in)            :: desc_a
      integer, intent(out)                       :: info
      complex(psb_spk_), intent(inout), optional, target :: work(:)
      integer, intent(in), optional              :: update,mode
    end subroutine psb_covrlv
    subroutine  psb_zovrlm(x,desc_a,info,jx,ik,work,update,mode)
      use psb_descriptor_type
      complex(psb_dpk_), intent(inout)           :: x(:,:)
      type(psb_desc_type), intent(in)            :: desc_a
      integer, intent(out)                       :: info
      complex(psb_dpk_), intent(inout), optional, target :: work(:)
      integer, intent(in), optional              :: update,jx,ik,mode
    end subroutine psb_zovrlm
    subroutine  psb_zovrlv(x,desc_a,info,work,update,mode)
      use psb_descriptor_type
      complex(psb_dpk_), intent(inout)           :: x(:)
      type(psb_desc_type), intent(in)            :: desc_a
      integer, intent(out)                       :: info
      complex(psb_dpk_), intent(inout), optional, target :: work(:)
      integer, intent(in), optional              :: update,mode
    end subroutine psb_zovrlv
  end interface

  interface psb_halo
    subroutine  psb_shalom(x,desc_a,info,alpha,jx,ik,work,tran,mode,data)
      use psb_descriptor_type
      real(psb_spk_), intent(inout)           :: x(:,:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      real(psb_spk_), intent(in), optional    :: alpha
      real(psb_spk_), target, optional, intent(inout) :: work(:)
      integer, intent(in), optional           :: mode,jx,ik,data
      character, intent(in), optional         :: tran
    end subroutine psb_shalom
    subroutine  psb_shalov(x,desc_a,info,alpha,work,tran,mode,data)
      use psb_descriptor_type
      real(psb_spk_), intent(inout)           :: x(:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      real(psb_spk_), intent(in), optional    :: alpha
      real(psb_spk_), target, optional, intent(inout) :: work(:)
      integer, intent(in), optional           :: mode,data
      character, intent(in), optional         :: tran
    end subroutine psb_shalov
    subroutine  psb_dhalom(x,desc_a,info,alpha,jx,ik,work,tran,mode,data)
      use psb_descriptor_type
      real(psb_dpk_), intent(inout)           :: x(:,:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      real(psb_dpk_), intent(in), optional    :: alpha
      real(psb_dpk_), target, optional, intent(inout) :: work(:)
      integer, intent(in), optional           :: mode,jx,ik,data
      character, intent(in), optional         :: tran
    end subroutine psb_dhalom
    subroutine  psb_dhalov(x,desc_a,info,alpha,work,tran,mode,data)
      use psb_descriptor_type
      real(psb_dpk_), intent(inout)           :: x(:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      real(psb_dpk_), intent(in), optional    :: alpha
      real(psb_dpk_), target, optional, intent(inout) :: work(:)
      integer, intent(in), optional           :: mode,data
      character, intent(in), optional         :: tran
    end subroutine psb_dhalov
    subroutine  psb_ihalom(x,desc_a,info,alpha,jx,ik,work,tran,mode,data)
      use psb_descriptor_type
      integer, intent(inout) :: x(:,:)
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
    subroutine  psb_chalom(x,desc_a,info,alpha,jx,ik,work,tran,mode,data)
      use psb_descriptor_type
      complex(psb_spk_), intent(inout)        :: x(:,:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      complex(psb_spk_), intent(in), optional :: alpha
      complex(psb_spk_), target, optional, intent(inout) :: work(:)
      integer, intent(in), optional           :: mode,jx,ik,data
      character, intent(in), optional         :: tran
    end subroutine psb_chalom
    subroutine  psb_chalov(x,desc_a,info,alpha,work,tran,mode,data)
      use psb_descriptor_type
      complex(psb_spk_), intent(inout)        :: x(:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      complex(psb_spk_), intent(in), optional :: alpha
      complex(psb_spk_), target, optional, intent(inout) :: work(:)
      integer, intent(in), optional           :: mode,data
      character, intent(in), optional         :: tran
    end subroutine psb_chalov
    subroutine  psb_zhalom(x,desc_a,info,alpha,jx,ik,work,tran,mode,data)
      use psb_descriptor_type
      complex(psb_dpk_), intent(inout)        :: x(:,:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      complex(psb_dpk_), intent(in), optional :: alpha
      complex(psb_dpk_), target, optional, intent(inout) :: work(:)
      integer, intent(in), optional           :: mode,jx,ik,data
      character, intent(in), optional         :: tran
    end subroutine psb_zhalom
    subroutine  psb_zhalov(x,desc_a,info,alpha,work,tran,mode,data)
      use psb_descriptor_type
      complex(psb_dpk_), intent(inout)        :: x(:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      complex(psb_dpk_), intent(in), optional :: alpha
      complex(psb_dpk_), target, optional, intent(inout) :: work(:)
      integer, intent(in), optional           :: mode,data
      character, intent(in), optional         :: tran
    end subroutine psb_zhalov
  end interface


  interface psb_scatter
    subroutine  psb_dscatterm(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      real(psb_dpk_), intent(out)    :: locx(:,:)
      real(psb_dpk_), intent(in)     :: globx(:,:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer, intent(out)             :: info
      integer, intent(in), optional    :: root
    end subroutine psb_dscatterm
    subroutine  psb_dscatterv(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      real(psb_dpk_), intent(out)    :: locx(:)
      real(psb_dpk_), intent(in)     :: globx(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer, intent(out)             :: info
      integer, intent(in), optional    :: root
    end subroutine psb_dscatterv
    subroutine  psb_zscatterm(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      complex(psb_dpk_), intent(out) :: locx(:,:)
      complex(psb_dpk_), intent(in)  :: globx(:,:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer, intent(out)             :: info
      integer, intent(in), optional    :: root
    end subroutine psb_zscatterm
    subroutine  psb_zscatterv(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      complex(psb_dpk_), intent(out) :: locx(:)
      complex(psb_dpk_), intent(in)  :: globx(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer, intent(out)             :: info
      integer, intent(in), optional    :: root
    end subroutine psb_zscatterv
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
    subroutine  psb_sscatterm(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      real(psb_spk_), intent(out)    :: locx(:,:)
      real(psb_spk_), intent(in)     :: globx(:,:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer, intent(out)             :: info
      integer, intent(in), optional    :: root
    end subroutine psb_sscatterm
    subroutine  psb_sscatterv(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      real(psb_spk_), intent(out)    :: locx(:)
      real(psb_spk_), intent(in)     :: globx(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer, intent(out)             :: info
      integer, intent(in), optional    :: root
    end subroutine psb_sscatterv
    subroutine  psb_cscatterm(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      complex(psb_spk_), intent(out) :: locx(:,:)
      complex(psb_spk_), intent(in)  :: globx(:,:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer, intent(out)             :: info
      integer, intent(in), optional    :: root
    end subroutine psb_cscatterm
    subroutine  psb_cscatterv(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      complex(psb_spk_), intent(out) :: locx(:)
      complex(psb_spk_), intent(in)  :: globx(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer, intent(out)             :: info
      integer, intent(in), optional    :: root
    end subroutine psb_cscatterv
  end interface

  interface psb_gather
    subroutine  psb_dsp_allgather(globa, loca, desc_a, info, root, dupl,keepnum,keeploc)
      use psb_descriptor_type
      use psb_mat_mod
      implicit none
      type(psb_d_sparse_mat), intent(inout) :: loca
      type(psb_d_sparse_mat), intent(out)   :: globa
      type(psb_desc_type), intent(in) :: desc_a
      integer, intent(out)            :: info
      integer, intent(in), optional   :: root,dupl
      logical, intent(in), optional   :: keepnum,keeploc
    end subroutine psb_dsp_allgather
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
    subroutine  psb_sgatherm(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      real(psb_spk_), intent(in)    :: locx(:,:)
      real(psb_spk_), intent(out)   :: globx(:,:)
      type(psb_desc_type), intent(in) :: desc_a
      integer, intent(out)            :: info
      integer, intent(in), optional   :: root
    end subroutine psb_sgatherm
    subroutine  psb_sgatherv(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      real(psb_spk_), intent(in)    :: locx(:)
      real(psb_spk_), intent(out)   :: globx(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer, intent(out)            :: info
      integer, intent(in), optional   :: root
    end subroutine psb_sgatherv
    subroutine  psb_dgatherm(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      real(psb_dpk_), intent(in)    :: locx(:,:)
      real(psb_dpk_), intent(out)   :: globx(:,:)
      type(psb_desc_type), intent(in) :: desc_a
      integer, intent(out)            :: info
      integer, intent(in), optional   :: root
    end subroutine psb_dgatherm
    subroutine  psb_dgatherv(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      real(psb_dpk_), intent(in)    :: locx(:)
      real(psb_dpk_), intent(out)   :: globx(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer, intent(out)            :: info
      integer, intent(in), optional   :: root
    end subroutine psb_dgatherv
    subroutine  psb_cgatherm(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      complex(psb_spk_), intent(in)  :: locx(:,:)
      complex(psb_spk_), intent(out) :: globx(:,:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer, intent(out)             :: info
      integer, intent(in), optional    :: root
    end subroutine psb_cgatherm
    subroutine  psb_cgatherv(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      complex(psb_spk_), intent(in)  :: locx(:)
      complex(psb_spk_), intent(out) :: globx(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer, intent(out)             :: info
      integer, intent(in), optional    :: root
    end subroutine psb_cgatherv
    subroutine  psb_zgatherm(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      complex(psb_dpk_), intent(in)  :: locx(:,:)
      complex(psb_dpk_), intent(out) :: globx(:,:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer, intent(out)             :: info
      integer, intent(in), optional    :: root
    end subroutine psb_zgatherm
    subroutine  psb_zgatherv(globx, locx, desc_a, info, root)
      use psb_descriptor_type
      complex(psb_dpk_), intent(in)  :: locx(:)
      complex(psb_dpk_), intent(out) :: globx(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer, intent(out)             :: info
      integer, intent(in), optional    :: root
    end subroutine psb_zgatherv
  end interface

end module psb_comm_mod
