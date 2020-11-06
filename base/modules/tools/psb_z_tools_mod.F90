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
Module psb_z_tools_mod
  use psb_desc_mod, only : psb_desc_type, psb_dpk_, psb_ipk_, psb_lpk_
  use psb_z_vect_mod, only : psb_z_base_vect_type, psb_z_vect_type
  use psb_z_mat_mod, only : psb_zspmat_type, psb_lzspmat_type, psb_z_base_sparse_mat, &
       & psb_lz_csr_sparse_mat, psb_lz_coo_sparse_mat, &
       & psb_z_csr_sparse_mat, psb_z_coo_sparse_mat
  use psb_l_vect_mod, only : psb_l_vect_type
  use psb_z_multivect_mod, only : psb_z_base_multivect_type, psb_z_multivect_type
  use psi_mod, only : psb_snd, psb_rcv ! Needed only for psb_getelem

  interface  psb_geall
    subroutine psb_zalloc_vect(x, desc_a,info)
      import
      implicit none
      type(psb_z_vect_type), intent(out)  :: x
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_),intent(out)             :: info
    end subroutine psb_zalloc_vect
    subroutine psb_zalloc_vect_r2(x, desc_a,info,n,lb)
      import
      implicit none
      type(psb_z_vect_type), allocatable, intent(out)  :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_),intent(out)             :: info
      integer(psb_ipk_), optional, intent(in)   :: n, lb
    end subroutine psb_zalloc_vect_r2
    subroutine psb_zalloc_multivect(x, desc_a,info,n)
      import
      implicit none
      type(psb_z_multivect_type), intent(out)  :: x
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_),intent(out)             :: info
      integer(psb_ipk_), optional, intent(in)   :: n
    end subroutine psb_zalloc_multivect
  end interface


  interface psb_geasb
    subroutine psb_zasb_vect(x, desc_a, info,mold, scratch)
      import
      implicit none
      type(psb_desc_type), intent(in)      ::  desc_a
      type(psb_z_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)                 ::  info
      class(psb_z_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine psb_zasb_vect
    subroutine psb_zasb_vect_r2(x, desc_a, info,mold, scratch)
      import
      implicit none
      type(psb_desc_type), intent(in)      ::  desc_a
      type(psb_z_vect_type), intent(inout) :: x(:)
      integer(psb_ipk_), intent(out)                 ::  info
      class(psb_z_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine psb_zasb_vect_r2
    subroutine psb_zasb_multivect(x, desc_a, info,mold, scratch, n)
      import
      implicit none
      type(psb_desc_type), intent(in)      ::  desc_a
      type(psb_z_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)                 ::  info
      class(psb_z_base_multivect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
      integer(psb_ipk_), optional, intent(in)   :: n
    end subroutine psb_zasb_multivect
  end interface

  interface psb_gefree
    subroutine psb_zfree_vect(x, desc_a, info)
      import
      implicit none
      type(psb_desc_type), intent(in)  ::  desc_a
      type(psb_z_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             ::  info
    end subroutine psb_zfree_vect
    subroutine psb_zfree_vect_r2(x, desc_a, info)
      import
      implicit none
      type(psb_desc_type), intent(in)  ::  desc_a
      type(psb_z_vect_type), allocatable, intent(inout) :: x(:)
      integer(psb_ipk_), intent(out)             ::  info
    end subroutine psb_zfree_vect_r2
    subroutine psb_zfree_multivect(x, desc_a, info)
      import
      implicit none
      type(psb_desc_type), intent(in)  ::  desc_a
      type(psb_z_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             ::  info
    end subroutine psb_zfree_multivect
  end interface


  interface psb_geins
    subroutine psb_zins_vect(m,irw,val,x,desc_a,info,dupl,local)
      import
      implicit none
      integer(psb_ipk_), intent(in)              :: m
      type(psb_desc_type), intent(in)  :: desc_a
      type(psb_z_vect_type), intent(inout) :: x
      integer(psb_lpk_), intent(in)              :: irw(:)
      complex(psb_dpk_), intent(in)    :: val(:)
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), optional, intent(in)    :: dupl
      logical, intent(in), optional        :: local
    end subroutine psb_zins_vect
    subroutine psb_zins_vect_v(m,irw,val,x,desc_a,info,dupl,local)
      import
      implicit none
      integer(psb_ipk_), intent(in)              :: m
      type(psb_desc_type), intent(in)  :: desc_a
      type(psb_z_vect_type), intent(inout) :: x
      type(psb_l_vect_type), intent(inout)       :: irw
      type(psb_z_vect_type), intent(inout)    :: val
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), optional, intent(in)    :: dupl
      logical, intent(in), optional        :: local
    end subroutine psb_zins_vect_v
    subroutine psb_zins_vect_r2(m,irw,val,x,desc_a,info,dupl,local)
      import
      implicit none
      integer(psb_ipk_), intent(in)              :: m
      type(psb_desc_type), intent(in)  :: desc_a
      type(psb_z_vect_type), intent(inout) :: x(:)
      integer(psb_lpk_), intent(in)              :: irw(:)
      complex(psb_dpk_), intent(in)    :: val(:,:)
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), optional, intent(in)    :: dupl
      logical, intent(in), optional        :: local
    end subroutine psb_zins_vect_r2
    subroutine psb_zins_multivect(m,irw,val,x,desc_a,info,dupl,local)
      import
      implicit none
      integer(psb_ipk_), intent(in)              :: m
      type(psb_desc_type), intent(in)  :: desc_a
      type(psb_z_multivect_type), intent(inout) :: x
      integer(psb_lpk_), intent(in)              :: irw(:)
      complex(psb_dpk_), intent(in)    :: val(:,:)
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), optional, intent(in)    :: dupl
      logical, intent(in), optional        :: local
    end subroutine psb_zins_multivect
  end interface

  interface psb_cdbldext
    Subroutine psb_zcdbldext(a,desc_a,novr,desc_ov,info,extype)
      import
      implicit none
      integer(psb_ipk_), intent(in)                     :: novr
      Type(psb_zspmat_type), Intent(in)        :: a
      Type(psb_desc_type), Intent(inout), target :: desc_a
      Type(psb_desc_type), Intent(out)           :: desc_ov
      integer(psb_ipk_), intent(out)                    :: info
      integer(psb_ipk_), intent(in),optional            :: extype
    end Subroutine psb_zcdbldext
  end interface

  interface psb_sphalo
    Subroutine psb_zsphalo(a,desc_a,blk,info,rowcnv,colcnv,&
         & rowscale,colscale,outfmt,data)
      import
      implicit none
      Type(psb_zspmat_type),Intent(in)       :: a
      Type(psb_zspmat_type),Intent(inout)    :: blk
      Type(psb_desc_type),Intent(in), target :: desc_a
      integer(psb_ipk_), intent(out)                   :: info
      logical, optional, intent(in)          :: rowcnv,colcnv,rowscale,colscale
      character(len=5), optional             :: outfmt
      integer(psb_ipk_), intent(in), optional          :: data
    end Subroutine psb_zsphalo
    Subroutine psb_lzsphalo(a,desc_a,blk,info,rowcnv,colcnv,&
         & rowscale,colscale,outfmt,data)
      import
      implicit none
      Type(psb_lzspmat_type),Intent(in)       :: a
      Type(psb_lzspmat_type),Intent(inout)    :: blk
      Type(psb_desc_type),Intent(in), target :: desc_a
      integer(psb_ipk_), intent(out)                   :: info
      logical, optional, intent(in)          :: rowcnv,colcnv,rowscale,colscale
      character(len=5), optional             :: outfmt
      integer(psb_ipk_), intent(in), optional          :: data
    end Subroutine psb_lzsphalo
    Subroutine psb_lz_csr_halo(a,desc_a,blk,info,rowcnv,colcnv,&
         & rowscale,colscale,data,outcol_glob,col_desc)
      import
      implicit none
      Type(psb_lz_csr_sparse_mat),Intent(in)       :: a
      Type(psb_lz_csr_sparse_mat),Intent(inout)    :: blk
      Type(psb_desc_type),Intent(in), target :: desc_a
      integer(psb_ipk_), intent(out)                   :: info
      logical, optional, intent(in)          :: rowcnv,colcnv,rowscale,colscale,outcol_glob
      integer(psb_ipk_), intent(in), optional          :: data
      type(psb_desc_type),Intent(in), optional, target :: col_desc
    end Subroutine psb_lz_csr_halo
    Subroutine psb_z_lz_csr_halo(a,desc_a,blk,info,rowcnv,colcnv,&
         &  rowscale,colscale,data,outcol_glob,col_desc)
      import
      implicit none
      type(psb_z_csr_sparse_mat),Intent(in)    :: a
      type(psb_lz_csr_sparse_mat),Intent(inout) :: blk
      type(psb_desc_type),intent(in), target :: desc_a
      integer(psb_ipk_), intent(out)                :: info
      logical, optional, intent(in)       :: rowcnv,colcnv,rowscale,colscale,outcol_glob
      integer(psb_ipk_), intent(in), optional       :: data
      type(psb_desc_type),Intent(in), optional, target :: col_desc
    end Subroutine psb_z_lz_csr_halo
  end interface


  interface psb_spall
    subroutine psb_zspalloc(a, desc_a, info, nnz)
      import
      implicit none
      type(psb_desc_type), intent(in) :: desc_a
      type(psb_zspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_ipk_), optional, intent(in)      :: nnz
    end subroutine psb_zspalloc
  end interface

  interface psb_spasb
    subroutine psb_zspasb(a,desc_a, info, afmt, upd, dupl,mold)
      import
      implicit none
      type(psb_zspmat_type), intent (inout)   :: a
      type(psb_desc_type), intent(in)         :: desc_a
      integer(psb_ipk_), intent(out)                    :: info
      integer(psb_ipk_),optional, intent(in)            :: dupl, upd
      character(len=*), optional, intent(in)  :: afmt
      class(psb_z_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_zspasb
  end interface

  interface psb_spfree
    subroutine psb_zspfree(a, desc_a,info)
      import
      implicit none
      type(psb_desc_type), intent(in) :: desc_a
      type(psb_zspmat_type), intent(inout)       ::a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_zspfree
  end interface


  interface psb_spins
    subroutine psb_zspins(nz,ia,ja,val,a,desc_a,info,rebuild,local)
      import
      implicit none
      type(psb_desc_type), intent(inout)   :: desc_a
      type(psb_zspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in)          :: nz
      integer(psb_lpk_), intent(in)          :: ia(:),ja(:)
      complex(psb_dpk_), intent(in)      :: val(:)
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional        :: rebuild
      logical, intent(in), optional        :: local
    end subroutine psb_zspins
    subroutine psb_zspins_csr_lirp(nr,irp,ja,val,irw,a,desc_a,info,rebuild,local)
      import
      implicit none
      type(psb_desc_type), intent(inout)     :: desc_a
      type(psb_zspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in)          :: nr
      integer(psb_lpk_), intent(in)          :: irw,irp(:),ja(:)
      complex(psb_dpk_), intent(in)            :: val(:)
      integer(psb_ipk_), intent(out)         :: info
      logical, intent(in), optional         :: rebuild, local
    end subroutine psb_zspins_csr_lirp
#if defined(IPK4) && defined(LPK8)
    subroutine psb_zspins_csr_iirp(nr,irw,irp,ja,val,a,desc_a,info,rebuild,local)
      import
      implicit none
      type(psb_desc_type), intent(inout)     :: desc_a
      type(psb_zspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in)          :: nr,irp(:)
      integer(psb_lpk_), intent(in)          :: irw,ja(:)
      complex(psb_dpk_), intent(in)            :: val(:)
      integer(psb_ipk_), intent(out)         :: info
      logical, intent(in), optional         :: rebuild, local
    end subroutine psb_zspins_csr_iirp
#endif
    subroutine psb_zspins_v(nz,ia,ja,val,a,desc_a,info,rebuild,local)
      use psb_i_vect_mod, only : psb_i_vect_type
      import
      implicit none
      type(psb_desc_type), intent(inout)   :: desc_a
      type(psb_zspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in)        :: nz
      type(psb_l_vect_type), intent(inout) :: ia,ja
      type(psb_z_vect_type), intent(inout) :: val
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: rebuild
      logical, intent(in), optional        :: local
    end subroutine psb_zspins_v
    subroutine psb_zspins_2desc(nz,ia,ja,val,a,desc_ar,desc_ac,info)
      import
      implicit none
      type(psb_desc_type), intent(in)      :: desc_ar
      type(psb_desc_type), intent(inout)   :: desc_ac
      type(psb_zspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in)       :: nz
      integer(psb_lpk_), intent(in)       :: ia(:),ja(:)
      complex(psb_dpk_), intent(in)        :: val(:)
      integer(psb_ipk_), intent(out)      :: info
    end subroutine psb_zspins_2desc
  end interface


  interface psb_sprn
    subroutine psb_zsprn(a, desc_a,info,clear)
      import
      implicit none
      type(psb_desc_type), intent(in)      :: desc_a
      type(psb_zspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(out)                 :: info
      logical, intent(in), optional        :: clear
    end subroutine psb_zsprn
  end interface

  interface psb_par_spspmm
    subroutine psb_z_par_csr_spspmm(acsr,desc_a,bcsr,ccsr,desc_c,info,data)
      import :: psb_z_csr_sparse_mat, psb_desc_type, psb_ipk_
      Implicit None
      type(psb_z_csr_sparse_mat),intent(in)    :: acsr
      type(psb_z_csr_sparse_mat),intent(inout) :: bcsr
      type(psb_z_csr_sparse_mat),intent(out)   :: ccsr
      type(psb_desc_type),intent(in)           :: desc_a
      type(psb_desc_type),intent(inout)        :: desc_c
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), intent(in), optional  :: data
    End Subroutine psb_z_par_csr_spspmm
    subroutine psb_lz_par_csr_spspmm(acsr,desc_a,bcsr,ccsr,desc_c,info,data)
      import :: psb_lz_csr_sparse_mat, psb_desc_type, psb_ipk_
      Implicit None
      type(psb_lz_csr_sparse_mat),intent(in)    :: acsr
      type(psb_lz_csr_sparse_mat),intent(inout) :: bcsr
      type(psb_lz_csr_sparse_mat),intent(out)   :: ccsr
      type(psb_desc_type),intent(in)           :: desc_a
      type(psb_desc_type),intent(inout)        :: desc_c
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), intent(in), optional  :: data
    End Subroutine psb_lz_par_csr_spspmm
  end interface psb_par_spspmm

  interface psb_glob_transpose
    subroutine psb_z_coo_glob_transpose(ain,desc_r,info,atrans,desc_c,desc_rx)
      import
      type(psb_z_coo_sparse_mat), intent(inout) :: ain
      type(psb_desc_type), intent(inout), target   :: desc_r
      type(psb_z_coo_sparse_mat), intent(out), optional :: atrans
      type(psb_desc_type), intent(inout), target, optional :: desc_c
      type(psb_desc_type), intent(out), optional   :: desc_rx
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_z_coo_glob_transpose
    subroutine psb_lz_coo_glob_transpose(ain,desc_r,info,atrans,desc_c,desc_rx)
      import
      type(psb_lz_coo_sparse_mat), intent(inout) :: ain
      type(psb_desc_type), intent(inout), target   :: desc_r
      type(psb_lz_coo_sparse_mat), intent(out), optional :: atrans
      type(psb_desc_type), intent(inout), target, optional :: desc_c
      type(psb_desc_type), intent(out), optional   :: desc_rx
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_lz_coo_glob_transpose
    subroutine psb_lz_simple_glob_transpose(ain,aout,desc_a,info)
      import
      type(psb_lzspmat_type), intent(in)  :: ain
      type(psb_lzspmat_type), intent(out) :: aout
      type(psb_desc_type)           :: desc_a
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_lz_simple_glob_transpose
    subroutine psb_lz_simple_glob_transpose_ip(ain,desc_a,info)
      import
      type(psb_lzspmat_type), intent(inout)  :: ain
      type(psb_desc_type)           :: desc_a
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_lz_simple_glob_transpose_ip
    subroutine psb_z_simple_glob_transpose(ain,aout,desc_a,info)
      import
      type(psb_zspmat_type), intent(in)  :: ain
      type(psb_zspmat_type), intent(out) :: aout
      type(psb_desc_type)           :: desc_a
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_z_simple_glob_transpose
    subroutine psb_z_simple_glob_transpose_ip(ain,desc_a,info)
      import
      type(psb_zspmat_type), intent(inout)  :: ain
      type(psb_desc_type)           :: desc_a
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_z_simple_glob_transpose_ip
  end interface psb_glob_transpose

  interface psb_getelem
    function psb_z_getelem(x,index,desc_a,info) result(res)
      import
      type(psb_z_vect_type), intent(inout) :: x
      integer(psb_lpk_), intent(in)          :: index
      type(psb_desc_type), intent(inout)     :: desc_a
      integer(psb_ipk_), intent(out)         :: info
      complex(psb_dpk_)                        :: res
    end function
  end interface

  interface psb_remap
    subroutine psb_z_remap(np_remap, desc_in, a_in, desc_out, a_out, info)      
      import      
      implicit none
      !....parameters...
      integer(psb_ipk_), intent(in)        :: np_remap
      type(psb_desc_type), intent(inout)   :: desc_in
      type(psb_zspmat_type), intent(inout) :: a_in
      type(psb_zspmat_type), intent(out)   :: a_out
      type(psb_desc_type), intent(out)     :: desc_out
      integer(psb_ipk_), intent(out)       :: info
    end subroutine psb_z_remap
  end interface psb_remap


end module psb_z_tools_mod
