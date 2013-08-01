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
module psb_cd_tools_mod

  use psb_const_mod
  use psb_desc_mod
  use psb_gen_block_map_mod
  use psb_list_map_mod
  use psb_glist_map_mod
  use psb_hash_map_mod
  use psb_repl_map_mod
    
  interface psb_cd_set_bld
    subroutine psb_cd_set_bld(desc,info)
      import :: psb_ipk_, psb_desc_type
      type(psb_desc_type), intent(inout) :: desc
      integer(psb_ipk_) :: info
    end subroutine psb_cd_set_bld
  end interface

  interface psb_cd_set_ovl_bld
    subroutine psb_cd_set_ovl_bld(desc,info)
      import :: psb_ipk_, psb_desc_type
      type(psb_desc_type), intent(inout) :: desc
      integer(psb_ipk_) :: info
    end subroutine psb_cd_set_ovl_bld
  end interface

  interface psb_cd_reinit
    Subroutine psb_cd_reinit(desc,info)
      import :: psb_ipk_, psb_desc_type
      Implicit None

      !     .. Array Arguments ..
      Type(psb_desc_type), Intent(inout) :: desc
      integer(psb_ipk_), intent(out)               :: info
    end Subroutine psb_cd_reinit
  end interface

  interface psb_cdcpy
    subroutine psb_cdcpy(desc_in, desc_out, info)
      import :: psb_ipk_, psb_desc_type

      implicit none
      !....parameters...

      type(psb_desc_type), intent(inout) :: desc_in
      type(psb_desc_type), intent(out)   :: desc_out
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psb_cdcpy
  end interface


  interface psb_cdprt
    subroutine psb_cdprt(iout,desc_p,glob,short)
      import :: psb_ipk_, psb_desc_type
      implicit none 
      type(psb_desc_type), intent(in)    :: desc_p
      integer(psb_ipk_), intent(in)                :: iout
      logical, intent(in), optional      :: glob,short
    end subroutine psb_cdprt
  end interface

  interface psb_cdins
    subroutine psb_cdinsrc(nz,ia,ja,desc_a,info,ila,jla)
      import :: psb_ipk_, psb_desc_type
      type(psb_desc_type), intent(inout) :: desc_a
      integer(psb_ipk_), intent(in)                :: nz,ia(:),ja(:)
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_ipk_), optional, intent(out)     :: ila(:), jla(:)
    end subroutine psb_cdinsrc
    subroutine psb_cdinsc(nz,ja,desc,info,jla,mask,lidx)
      import :: psb_ipk_, psb_desc_type
      type(psb_desc_type), intent(inout)         :: desc
      integer(psb_ipk_), intent(in)              :: nz,ja(:)
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), optional, intent(out)   :: jla(:)
      logical, optional, target, intent(in)      :: mask(:)
      integer(psb_ipk_), intent(in), optional    :: lidx(:)
    end subroutine psb_cdinsc
  end interface

  interface psb_cdbldext
    Subroutine psb_cd_lstext(desc_a,in_list,desc_ov,info, mask,extype)
      import :: psb_ipk_, psb_desc_type
      Implicit None
      Type(psb_desc_type), Intent(inout), target :: desc_a
      integer(psb_ipk_), intent(in)                     :: in_list(:)
      Type(psb_desc_type), Intent(out)        :: desc_ov
      integer(psb_ipk_), intent(out)                    :: info
      logical, intent(in), optional, target   :: mask(:)
      integer(psb_ipk_), intent(in),optional            :: extype
    end Subroutine psb_cd_lstext
  end interface


  interface psb_cdren
    subroutine psb_cdren(trans,iperm,desc_a,info)
      import :: psb_ipk_, psb_desc_type
      type(psb_desc_type), intent(inout)    :: desc_a
      integer(psb_ipk_), intent(inout)                :: iperm(:)
      character, intent(in)                 :: trans
      integer(psb_ipk_), intent(out)                  :: info
    end subroutine psb_cdren
  end interface

  interface psb_get_overlap
    subroutine psb_get_ovrlap(ovrel,desc,info)
      import :: psb_ipk_, psb_desc_type
      implicit none 
      integer(psb_ipk_), allocatable, intent(out) :: ovrel(:)
      type(psb_desc_type), intent(in) :: desc
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_get_ovrlap
  end interface

  interface psb_icdasb
    subroutine psb_icdasb(desc,info,ext_hv, mold)
      import :: psb_ipk_, psb_desc_type, psb_i_base_vect_type
      Type(psb_desc_type), intent(inout) :: desc
      integer(psb_ipk_), intent(out)               :: info
      logical, intent(in),optional       :: ext_hv
      type(psb_i_base_vect_type), optional, intent(in) :: mold
    end subroutine psb_icdasb
  end interface

  interface psb_cdall

    subroutine psb_cdall(ictxt, desc, info,mg,ng,parts,vg,vl,flag,nl,repl,&
         & globalcheck,lidx)
      import :: psb_ipk_, psb_desc_type, psb_parts
      implicit None
      procedure(psb_parts)                :: parts
      integer(psb_ipk_), intent(in)       :: mg,ng,ictxt, vg(:), vl(:),nl,lidx(:)
      integer(psb_ipk_), intent(in)       :: flag
      logical, intent(in)                 :: repl, globalcheck
      integer(psb_ipk_), intent(out)      :: info
      type(psb_desc_type), intent(out)    :: desc      
      optional :: mg,ng,parts,vg,vl,flag,nl,repl, globalcheck,lidx
    end subroutine psb_cdall
    
  end interface

  interface psb_cdasb
    module procedure psb_cdasb
  end interface

  interface psb_get_boundary
    module procedure psb_get_boundary
  end interface

  interface 
    subroutine psb_cd_switch_ovl_indxmap(desc,info) 
      import :: psb_ipk_, psb_desc_type
      implicit None
      include 'parts.fh'
      type(psb_desc_type), intent(inout) :: desc
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_cd_switch_ovl_indxmap
  end interface

contains

  subroutine psb_get_boundary(bndel,desc,info)
    use psi_mod, only : psi_crea_bnd_elem
    implicit none 
    integer(psb_ipk_), allocatable, intent(out) :: bndel(:)
    type(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(out)            :: info

    call psi_crea_bnd_elem(bndel,desc,info)

  end subroutine psb_get_boundary

  subroutine psb_cdasb(desc,info,mold)

    Type(psb_desc_type), intent(inout) :: desc
    integer(psb_ipk_), intent(out)     :: info
    type(psb_i_base_vect_type), optional, intent(in) :: mold

    call psb_icdasb(desc,info,ext_hv=.false.,mold=mold)
  end subroutine psb_cdasb

end module psb_cd_tools_mod


