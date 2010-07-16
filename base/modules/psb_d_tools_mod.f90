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

Module psb_d_tools_mod

  interface  psb_geall
    subroutine psb_dalloc(x, desc_a, info, n, lb)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      implicit none
      real(psb_dpk_), allocatable, intent(out) :: x(:,:)
      type(psb_desc_type), intent(in) :: desc_a
      integer,intent(out)             :: info
      integer, optional, intent(in)   :: n, lb
    end subroutine psb_dalloc
    subroutine psb_dallocv(x, desc_a,info,n)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_dpk_), allocatable, intent(out)       :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer,intent(out)             :: info
      integer, optional, intent(in)   :: n
    end subroutine psb_dallocv
  end interface


  interface psb_geasb
    subroutine psb_dasb(x, desc_a, info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      type(psb_desc_type), intent(in) ::  desc_a
      real(psb_dpk_), allocatable, intent(inout)       ::  x(:,:)
      integer, intent(out)            ::  info
    end subroutine psb_dasb
    subroutine psb_dasbv(x, desc_a, info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      type(psb_desc_type), intent(in) ::  desc_a
      real(psb_dpk_), allocatable, intent(inout)   ::  x(:)
      integer, intent(out)        ::  info
    end subroutine psb_dasbv
  end interface

  interface psb_sphalo
    Subroutine psb_dsphalo(a,desc_a,blk,info,rowcnv,colcnv,&
         & rowscale,colscale,outfmt,data)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      use psb_mat_mod, only : psb_d_sparse_mat
      Type(psb_d_sparse_mat),Intent(in)    :: a
      Type(psb_d_sparse_mat),Intent(inout) :: blk
      Type(psb_desc_type),Intent(in),target :: desc_a
      integer, intent(out)                :: info
      logical, optional, intent(in)       :: rowcnv,colcnv,rowscale,colscale
      character(len=5), optional          :: outfmt 
      integer, intent(in), optional       :: data
    end Subroutine psb_dsphalo
  end interface

  interface psb_gefree
    subroutine psb_dfree(x, desc_a, info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_dpk_),allocatable, intent(inout)        :: x(:,:)
      type(psb_desc_type), intent(in) :: desc_a
      integer, intent(out)            :: info
    end subroutine psb_dfree
    subroutine psb_dfreev(x, desc_a, info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_dpk_),allocatable, intent(inout)        :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer, intent(out)       :: info
    end subroutine psb_dfreev
  end interface

  interface psb_geins
    subroutine psb_dinsi(m,irw,val, x,desc_a,info,dupl)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      integer, intent(in)                ::  m
      type(psb_desc_type), intent(in)    ::  desc_a
      real(psb_dpk_),intent(inout)           ::  x(:,:)
      integer, intent(in)                ::  irw(:)
      real(psb_dpk_), intent(in)       ::  val(:,:)
      integer, intent(out)               ::  info
      integer, optional, intent(in)      ::  dupl
    end subroutine psb_dinsi
    subroutine psb_dinsvi(m,irw,val,x,desc_a,info,dupl)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      integer, intent(in)              ::  m
      type(psb_desc_type), intent(in)  ::  desc_a
      real(psb_dpk_),intent(inout)     ::  x(:)
      integer, intent(in)              ::  irw(:)
      real(psb_dpk_), intent(in)       ::  val(:)
      integer, intent(out)             ::  info
      integer, optional, intent(in)    ::  dupl
    end subroutine psb_dinsvi
  end interface


  interface psb_cdbldext
    Subroutine psb_dcdbldext(a,desc_a,novr,desc_ov,info,extype)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      Use psb_mat_mod, only : psb_d_sparse_mat
      integer, intent(in)                     :: novr
      Type(psb_d_sparse_mat), Intent(in)       :: a
      Type(psb_desc_type), Intent(in), target :: desc_a
      Type(psb_desc_type), Intent(out)        :: desc_ov
      integer, intent(out)                    :: info
      integer, intent(in),optional            :: extype
    end Subroutine psb_dcdbldext
  end interface

  interface psb_spall
    subroutine psb_dspalloc(a, desc_a, info, nnz)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      use psb_mat_mod, only : psb_d_sparse_mat
      type(psb_desc_type), intent(in)       :: desc_a
      type(psb_d_sparse_mat), intent(inout) :: a
      integer, intent(out)                :: info
      integer, optional, intent(in)       :: nnz
    end subroutine psb_dspalloc
  end interface

  interface psb_spasb
    subroutine psb_dspasb(a,desc_a, info, afmt, upd, dupl,mold)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      use psb_mat_mod, only : psb_d_sparse_mat, psb_d_base_sparse_mat
      type(psb_d_sparse_mat), intent (inout)  :: a
      type(psb_desc_type), intent(in)         :: desc_a
      integer, intent(out)                    :: info
      integer,optional, intent(in)            :: dupl, upd
      character(len=*), optional, intent(in)  :: afmt
      class(psb_d_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_dspasb
  end interface

  interface psb_spfree
    subroutine psb_dspfree(a, desc_a,info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      use psb_mat_mod, only : psb_d_sparse_mat
      type(psb_desc_type), intent(in)       :: desc_a
      type(psb_d_sparse_mat), intent(inout) :: a
      integer, intent(out)                  :: info
    end subroutine psb_dspfree
  end interface


  interface psb_spins
    subroutine psb_dspins(nz,ia,ja,val,a,desc_a,info,rebuild)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      use psb_mat_mod, only : psb_d_sparse_mat
      type(psb_desc_type), intent(inout)   :: desc_a
      type(psb_d_sparse_mat), intent(inout) :: a
      integer, intent(in)                  :: nz,ia(:),ja(:)
      real(psb_dpk_), intent(in)         :: val(:)
      integer, intent(out)                 :: info
      logical, intent(in), optional        :: rebuild
    end subroutine psb_dspins
    subroutine psb_dspins_2desc(nz,ia,ja,val,a,desc_ar,desc_ac,info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      use psb_mat_mod, only : psb_d_sparse_mat
      type(psb_d_sparse_mat), intent(inout) :: a
      type(psb_desc_type), intent(in)      :: desc_ar
      type(psb_desc_type), intent(inout)   :: desc_ac
      integer, intent(in)                  :: nz,ia(:),ja(:)
      real(psb_dpk_), intent(in)           :: val(:)
      integer, intent(out)                 :: info
    end subroutine psb_dspins_2desc
  end interface


  interface psb_sprn
    subroutine psb_dsprn(a, desc_a,info,clear)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      use psb_mat_mod, only : psb_d_sparse_mat
      type(psb_desc_type), intent(in)       :: desc_a
      type(psb_d_sparse_mat), intent(inout) :: a
      integer, intent(out)                  :: info
      logical, intent(in), optional         :: clear
    end subroutine psb_dsprn
  end interface

!!$
!!$  interface psb_linmap_init
!!$    module procedure psb_dlinmap_init
!!$  end interface
!!$
!!$  interface psb_linmap_ins
!!$    module procedure psb_dlinmap_ins
!!$  end interface
!!$
!!$  interface psb_linmap_asb
!!$    module procedure psb_dlinmap_asb
!!$  end interface
!!$
!!$contains
!!$
!!$  subroutine psb_dlinmap_init(a_map,cd_xt,descin,descout)
!!$    use psb_descriptor_type
!!$    use psb_serial_mod
!!$    use psb_penv_mod
!!$    use psb_error_mod
!!$    use psb_base_tools_mod
!!$    use psb_d_mat_mod
!!$    implicit none 
!!$    type(psb_d_sparse_mat), intent(out) :: a_map
!!$    type(psb_desc_type), intent(out)   :: cd_xt
!!$    type(psb_desc_type), intent(in)    :: descin, descout 
!!$
!!$    integer :: nrow_in, nrow_out, ncol_in, info, ictxt
!!$
!!$    ictxt = psb_cd_get_context(descin)
!!$    call psb_cdcpy(descin,cd_xt,info)
!!$    if (info == psb_success_) call psb_cd_reinit(cd_xt,info)
!!$    if (info /= psb_success_) then 
!!$      write(psb_err_unit,*) 'Error on reinitialising the extension map'
!!$      call psb_error(ictxt)
!!$      call psb_abort(ictxt)
!!$      stop
!!$    end if
!!$
!!$    nrow_in  = psb_cd_get_local_rows(cd_xt)
!!$    ncol_in  = psb_cd_get_local_cols(cd_xt)
!!$    nrow_out = psb_cd_get_local_rows(descout)
!!$
!!$    call a_map%csall(nrow_out,ncol_in,info)
!!$
!!$  end subroutine psb_dlinmap_init
!!$
!!$  subroutine psb_dlinmap_ins(nz,ir,ic,val,a_map,cd_xt,descin,descout)
!!$    use psb_d_mat_mod
!!$    use psb_descriptor_type
!!$    implicit none 
!!$    integer, intent(in)                  :: nz
!!$    integer, intent(in)                  :: ir(:),ic(:)
!!$    real(psb_dpk_), intent(in)         :: val(:)
!!$    type(psb_d_sparse_mat), intent(inout) :: a_map
!!$    type(psb_desc_type), intent(inout)   :: cd_xt
!!$    type(psb_desc_type), intent(in)      :: descin, descout 
!!$    integer :: info
!!$    call psb_spins(nz,ir,ic,val,a_map,descout,cd_xt,info)
!!$
!!$  end subroutine psb_dlinmap_ins
!!$
!!$  subroutine psb_dlinmap_asb(a_map,cd_xt,descin,descout,afmt)
!!$    use psb_base_tools_mod
!!$    use psb_d_mat_mod
!!$    use psb_descriptor_type
!!$    use psb_serial_mod
!!$    implicit none 
!!$    type(psb_d_sparse_mat), intent(inout)   :: a_map
!!$    type(psb_desc_type), intent(inout)     :: cd_xt
!!$    type(psb_desc_type), intent(in)        :: descin, descout 
!!$    character(len=*), optional, intent(in) :: afmt
!!$
!!$
!!$    integer :: nrow_in, nrow_out, ncol_in, info, ictxt
!!$
!!$    ictxt = psb_cd_get_context(descin)
!!$
!!$    call psb_cdasb(cd_xt,info)
!!$    call a_map%set_ncols(psb_cd_get_local_cols(cd_xt))
!!$    call a_map%cscnv(info,type=afmt)
!!$
!!$  end subroutine psb_dlinmap_asb
!!$
end module psb_d_tools_mod
