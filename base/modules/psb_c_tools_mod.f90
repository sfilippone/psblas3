!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
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
Module psb_c_tools_mod
  use psb_desc_mod, only : psb_desc_type, psb_spk_, psb_ipk_
  use psb_c_vect_mod, only : psb_c_base_vect_type, psb_c_vect_type, psb_i_vect_type
  use psb_c_mat_mod, only : psb_cspmat_type, psb_c_base_sparse_mat

  interface  psb_geall
    subroutine psb_calloc(x, desc_a, info, n, lb)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      implicit none
      complex(psb_spk_), allocatable, intent(out)    :: x(:,:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), optional, intent(in)   :: n, lb
    end subroutine psb_calloc
    subroutine psb_callocv(x, desc_a,info,n)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      complex(psb_spk_), allocatable, intent(out)    :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), optional, intent(in)   :: n
    end subroutine psb_callocv
    subroutine psb_calloc_vect(x, desc_a,info,n)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      type(psb_c_vect_type), intent(out)  :: x
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_),intent(out)             :: info
      integer(psb_ipk_), optional, intent(in)   :: n
    end subroutine psb_calloc_vect
    subroutine psb_calloc_vect_r2(x, desc_a,info,n,lb)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      type(psb_c_vect_type), allocatable, intent(out)  :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_),intent(out)             :: info
      integer(psb_ipk_), optional, intent(in)   :: n, lb
    end subroutine psb_calloc_vect_r2
  end interface


  interface psb_geasb
    subroutine psb_casb(x, desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      type(psb_desc_type), intent(in) ::  desc_a
      complex(psb_spk_), allocatable, intent(inout)       ::  x(:,:)
      integer(psb_ipk_), intent(out)            ::  info
    end subroutine psb_casb
    subroutine psb_casbv(x, desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      type(psb_desc_type), intent(in) ::  desc_a
      complex(psb_spk_), allocatable, intent(inout)   ::  x(:)
      integer(psb_ipk_), intent(out)        ::  info
    end subroutine psb_casbv
    subroutine psb_casb_vect(x, desc_a, info,mold, scratch)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      type(psb_desc_type), intent(in)      ::  desc_a
      type(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)                 ::  info
      class(psb_c_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine psb_casb_vect
    subroutine psb_casb_vect_r2(x, desc_a, info,mold, scratch)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      type(psb_desc_type), intent(in)      ::  desc_a
      type(psb_c_vect_type), intent(inout) :: x(:)
      integer(psb_ipk_), intent(out)                 ::  info
      class(psb_c_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine psb_casb_vect_r2
  end interface

  interface psb_sphalo
    Subroutine psb_csphalo(a,desc_a,blk,info,rowcnv,colcnv,&
         & rowscale,colscale,outfmt,data)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      Type(psb_cspmat_type),Intent(in)       :: a
      Type(psb_cspmat_type),Intent(inout)    :: blk
      Type(psb_desc_type),Intent(in), target :: desc_a
      integer(psb_ipk_), intent(out)                   :: info
      logical, optional, intent(in)          :: rowcnv,colcnv,rowscale,colscale
      character(len=5), optional             :: outfmt 
      integer(psb_ipk_), intent(in), optional          :: data
    end Subroutine psb_csphalo
  end interface

  interface psb_gefree
    subroutine psb_cfree(x, desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      complex(psb_spk_),allocatable, intent(inout)        :: x(:,:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_cfree
    subroutine psb_cfreev(x, desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      complex(psb_spk_),allocatable, intent(inout)        :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_cfreev
    subroutine psb_cfree_vect(x, desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      type(psb_desc_type), intent(in)  ::  desc_a
      type(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             ::  info
    end subroutine psb_cfree_vect
    subroutine psb_cfree_vect_r2(x, desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      type(psb_desc_type), intent(in)  ::  desc_a
      type(psb_c_vect_type), allocatable, intent(inout) :: x(:)
      integer(psb_ipk_), intent(out)             ::  info
    end subroutine psb_cfree_vect_r2
  end interface


  interface psb_geins
    subroutine psb_cinsi(m,irw,val, x, desc_a,info,dupl,local)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      integer(psb_ipk_), intent(in)              ::  m
      type(psb_desc_type), intent(in)  ::  desc_a
      complex(psb_spk_),intent(inout)      ::  x(:,:)
      integer(psb_ipk_), intent(in)              ::  irw(:)
      complex(psb_spk_), intent(in)  ::  val(:,:)
      integer(psb_ipk_), intent(out)             ::  info
      integer(psb_ipk_), optional, intent(in)    ::  dupl
      logical, intent(in), optional        :: local
    end subroutine psb_cinsi
    subroutine psb_cinsvi(m, irw,val, x,desc_a,info,dupl,local)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      integer(psb_ipk_), intent(in)              ::  m
      type(psb_desc_type), intent(in)  ::  desc_a
      complex(psb_spk_),intent(inout)      ::  x(:)
      integer(psb_ipk_), intent(in)              ::  irw(:)
      complex(psb_spk_), intent(in)  ::  val(:)
      integer(psb_ipk_), intent(out)             ::  info
      integer(psb_ipk_), optional, intent(in)    ::  dupl
      logical, intent(in), optional        :: local
    end subroutine psb_cinsvi
    subroutine psb_cins_vect(m,irw,val,x,desc_a,info,dupl,local)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      integer(psb_ipk_), intent(in)              :: m
      type(psb_desc_type), intent(in)  :: desc_a
      type(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)              :: irw(:)
      complex(psb_spk_), intent(in)    :: val(:)
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), optional, intent(in)    :: dupl
      logical, intent(in), optional        :: local
    end subroutine psb_cins_vect
    subroutine psb_cins_vect_v(m,irw,val,x,desc_a,info,dupl,local)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, psb_i_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      integer(psb_ipk_), intent(in)              :: m
      type(psb_desc_type), intent(in)  :: desc_a
      type(psb_c_vect_type), intent(inout) :: x
      type(psb_i_vect_type), intent(inout)       :: irw
      type(psb_c_vect_type), intent(inout)    :: val
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), optional, intent(in)    :: dupl
      logical, intent(in), optional        :: local
    end subroutine psb_cins_vect_v
    subroutine psb_cins_vect_r2(m,irw,val,x,desc_a,info,dupl,local)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      integer(psb_ipk_), intent(in)              :: m
      type(psb_desc_type), intent(in)  :: desc_a
      type(psb_c_vect_type), intent(inout) :: x(:)
      integer(psb_ipk_), intent(in)              :: irw(:)
      complex(psb_spk_), intent(in)    :: val(:,:)
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), optional, intent(in)    :: dupl
      logical, intent(in), optional        :: local
    end subroutine psb_cins_vect_r2
  end interface

  interface psb_cdbldext
    Subroutine psb_ccdbldext(a,desc_a,novr,desc_ov,info,extype)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      integer(psb_ipk_), intent(in)                     :: novr
      Type(psb_cspmat_type), Intent(in)        :: a
      Type(psb_desc_type), Intent(inout), target :: desc_a
      Type(psb_desc_type), Intent(out)           :: desc_ov
      integer(psb_ipk_), intent(out)                    :: info
      integer(psb_ipk_), intent(in),optional            :: extype
    end Subroutine psb_ccdbldext
  end interface

  interface psb_spall
    subroutine psb_cspalloc(a, desc_a, info, nnz)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      type(psb_desc_type), intent(in) :: desc_a
      type(psb_cspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_ipk_), optional, intent(in)      :: nnz
    end subroutine psb_cspalloc
  end interface

  interface psb_spasb
    subroutine psb_cspasb(a,desc_a, info, afmt, upd, dupl,mold)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      type(psb_cspmat_type), intent (inout)   :: a
      type(psb_desc_type), intent(in)         :: desc_a
      integer(psb_ipk_), intent(out)                    :: info
      integer(psb_ipk_),optional, intent(in)            :: dupl, upd
      character(len=*), optional, intent(in)  :: afmt
      class(psb_c_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_cspasb
  end interface

  interface psb_spfree
    subroutine psb_cspfree(a, desc_a,info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      type(psb_desc_type), intent(in) :: desc_a
      type(psb_cspmat_type), intent(inout)       ::a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_cspfree
  end interface


  interface psb_spins
    subroutine psb_cspins(nz,ia,ja,val,a,desc_a,info,rebuild,local)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      type(psb_desc_type), intent(inout)   :: desc_a
      type(psb_cspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in)                  :: nz,ia(:),ja(:)
      complex(psb_spk_), intent(in)      :: val(:)
      integer(psb_ipk_), intent(out)                 :: info
      logical, intent(in), optional        :: rebuild
      logical, intent(in), optional        :: local
    end subroutine psb_cspins
    subroutine psb_cspins_v(nz,ia,ja,val,a,desc_a,info,rebuild,local)
      use psb_i_vect_mod, only : psb_i_vect_type
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type,&
           & psb_cspmat_type, psb_c_base_sparse_mat
      type(psb_desc_type), intent(inout)   :: desc_a
      type(psb_cspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in)        :: nz
      type(psb_i_vect_type), intent(inout) :: ia,ja
      type(psb_c_vect_type), intent(inout) :: val
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: rebuild
      logical, intent(in), optional        :: local
    end subroutine psb_cspins_v
    subroutine psb_cspins_2desc(nz,ia,ja,val,a,desc_ar,desc_ac,info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      type(psb_desc_type), intent(in)      :: desc_ar
      type(psb_desc_type), intent(inout)   :: desc_ac
      type(psb_cspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(in)                  :: nz,ia(:),ja(:)
      complex(psb_spk_), intent(in)        :: val(:)
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine psb_cspins_2desc
  end interface


  interface psb_sprn
    subroutine psb_csprn(a, desc_a,info,clear)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_base_vect_type, psb_c_vect_type, &
           & psb_cspmat_type, psb_c_base_sparse_mat
      type(psb_desc_type), intent(in)      :: desc_a
      type(psb_cspmat_type), intent(inout) :: a
      integer(psb_ipk_), intent(out)                 :: info
      logical, intent(in), optional        :: clear
    end subroutine psb_csprn
  end interface

end module psb_c_tools_mod
