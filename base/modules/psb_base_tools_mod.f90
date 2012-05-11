!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012
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
module psb_iv_tools_mod
  use psb_const_mod
  use psb_descriptor_type
  

  interface  psb_geall
    subroutine psb_ialloc(x, desc_a, info,n, lb)
      import :: psb_ipk_, psb_desc_type
      integer(psb_ipk_), allocatable, intent(out) :: x(:,:)
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      integer(psb_ipk_), optional, intent(in)     :: n, lb
    end subroutine psb_ialloc
    subroutine psb_iallocv(x, desc_a,info,n)
      import :: psb_ipk_, psb_desc_type
      integer(psb_ipk_), allocatable, intent(out) :: x(:)
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      integer(psb_ipk_), optional, intent(in)     :: n
    end subroutine psb_iallocv
  end interface


  interface psb_geasb
    subroutine psb_iasb(x, desc_a, info)
      import :: psb_ipk_, psb_desc_type
      type(psb_desc_type), intent(in) ::  desc_a
      integer(psb_ipk_), allocatable, intent(inout)  ::  x(:,:)
      integer(psb_ipk_), intent(out)            ::  info
    end subroutine psb_iasb
    subroutine psb_iasbv(x, desc_a, info)
      import :: psb_ipk_, psb_desc_type
      type(psb_desc_type), intent(in) ::  desc_a
      integer(psb_ipk_), allocatable, intent(inout)   ::  x(:)
      integer(psb_ipk_), intent(out)        ::  info
    end subroutine psb_iasbv
  end interface


  interface psb_gefree
    subroutine psb_ifree(x, desc_a, info)
      import :: psb_ipk_, psb_desc_type
      integer(psb_ipk_),allocatable, intent(inout) :: x(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_ifree
    subroutine psb_ifreev(x, desc_a, info)
      import :: psb_ipk_, psb_desc_type
      integer(psb_ipk_), allocatable, intent(inout) :: x(:)
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_ifreev
  end interface

  interface psb_geins
    subroutine psb_iinsi(m,irw,val, x,desc_a,info,dupl)
      import :: psb_ipk_, psb_desc_type
      integer(psb_ipk_), intent(in)              ::  m
      type(psb_desc_type), intent(in)  ::  desc_a
      integer(psb_ipk_),intent(inout)                  ::  x(:,:)
      integer(psb_ipk_), intent(in)              ::  irw(:)
      integer(psb_ipk_), intent(in)              ::  val(:,:)
      integer(psb_ipk_), intent(out)             ::  info
      integer(psb_ipk_), optional, intent(in)    ::  dupl
    end subroutine psb_iinsi
    subroutine psb_iinsvi(m, irw,val, x,desc_a,info,dupl)
      import :: psb_ipk_, psb_desc_type
      integer(psb_ipk_), intent(in)             ::  m
      type(psb_desc_type), intent(in) ::  desc_a
      integer(psb_ipk_),intent(inout)                 ::  x(:)
      integer(psb_ipk_), intent(in)             ::  irw(:)
      integer(psb_ipk_), intent(in)             ::  val(:)
      integer(psb_ipk_), intent(out)            ::  info
      integer(psb_ipk_), optional, intent(in)   ::  dupl
    end subroutine psb_iinsvi
  end interface


  interface psb_glob_to_loc
    subroutine psb_glob_to_loc2(x,y,desc_a,info,iact,owned)
      import :: psb_ipk_, psb_desc_type
      type(psb_desc_type), intent(in)    ::  desc_a
      integer(psb_ipk_),intent(in)                 ::  x(:)  
      integer(psb_ipk_),intent(out)                ::  y(:)  
      integer(psb_ipk_), intent(out)               ::  info
      logical, intent(in),  optional     :: owned
      character, intent(in), optional    ::  iact
    end subroutine psb_glob_to_loc2
    subroutine psb_glob_to_loc(x,desc_a,info,iact,owned)
      import :: psb_ipk_, psb_desc_type
      type(psb_desc_type), intent(in)    ::  desc_a
      integer(psb_ipk_),intent(inout)              ::  x(:)  
      integer(psb_ipk_), intent(out)               ::  info
      logical, intent(in),  optional     :: owned
      character, intent(in), optional    ::  iact
    end subroutine psb_glob_to_loc
    subroutine psb_glob_to_loc2s(x,y,desc_a,info,iact,owned)
      import :: psb_ipk_, psb_desc_type
      implicit none 
      type(psb_desc_type), intent(in)    ::  desc_a
      integer(psb_ipk_),intent(in)                 ::  x
      integer(psb_ipk_),intent(out)                ::  y  
      integer(psb_ipk_), intent(out)               ::  info
      character, intent(in), optional    ::  iact
      logical, intent(in),  optional     :: owned

    end subroutine psb_glob_to_loc2s

    subroutine psb_glob_to_locs(x,desc_a,info,iact,owned)
      import :: psb_ipk_, psb_desc_type
      implicit none 
      type(psb_desc_type), intent(in)    ::  desc_a
      integer(psb_ipk_),intent(inout)              ::  x  
      integer(psb_ipk_), intent(out)               ::  info
      character, intent(in), optional    ::  iact
      logical, intent(in),  optional     :: owned
    end subroutine psb_glob_to_locs
  end interface

  interface psb_loc_to_glob
    subroutine psb_loc_to_glob2(x,y,desc_a,info,iact)
      import :: psb_ipk_, psb_desc_type
      type(psb_desc_type), intent(in)    ::  desc_a
      integer(psb_ipk_),intent(in)                 ::  x(:)  
      integer(psb_ipk_),intent(out)                ::  y(:)  
      integer(psb_ipk_), intent(out)               ::  info
      character, intent(in), optional    ::  iact
    end subroutine psb_loc_to_glob2
    subroutine psb_loc_to_glob(x,desc_a,info,iact)
      import :: psb_ipk_, psb_desc_type
      type(psb_desc_type), intent(in)    ::  desc_a
      integer(psb_ipk_),intent(inout)              ::  x(:)  
      integer(psb_ipk_), intent(out)               ::  info
      character, intent(in), optional    ::  iact
    end subroutine psb_loc_to_glob
    subroutine psb_loc_to_glob2s(x,y,desc_a,info,iact)
      import :: psb_ipk_, psb_desc_type
      implicit none 
      type(psb_desc_type), intent(in)    ::  desc_a
      integer(psb_ipk_),intent(in)                 ::  x
      integer(psb_ipk_),intent(out)                ::  y  
      integer(psb_ipk_), intent(out)               ::  info
      character, intent(in), optional    ::  iact

    end subroutine psb_loc_to_glob2s
    subroutine psb_loc_to_globs(x,desc_a,info,iact)
      import :: psb_ipk_, psb_desc_type
      type(psb_desc_type), intent(in)    ::  desc_a
      integer(psb_ipk_),intent(inout)              ::  x  
      integer(psb_ipk_), intent(out)               ::  info
      character, intent(in), optional    ::  iact
    end subroutine psb_loc_to_globs

  end interface


  interface psb_is_owned
    module procedure psb_is_owned
  end interface

  interface psb_is_local
    module procedure psb_is_local
  end interface

  interface psb_owned_index
    module procedure psb_owned_index, psb_owned_index_v
  end interface

  interface psb_local_index
    module procedure psb_local_index, psb_local_index_v
  end interface

contains

  function psb_is_owned(idx,desc)
    implicit none 
    integer(psb_ipk_), intent(in) :: idx
    type(psb_desc_type), intent(in) :: desc
    logical :: psb_is_owned
    logical :: res
    integer(psb_ipk_) :: info

    call psb_owned_index(res,idx,desc,info)
    if (info /= psb_success_) res=.false.
    psb_is_owned = res
  end function psb_is_owned

  function psb_is_local(idx,desc)
    implicit none 
    integer(psb_ipk_), intent(in) :: idx
    type(psb_desc_type), intent(in) :: desc
    logical :: psb_is_local
    logical :: res
    integer(psb_ipk_) :: info

    call psb_local_index(res,idx,desc,info)
    if (info /= psb_success_) res=.false.
    psb_is_local = res
  end function psb_is_local

  subroutine psb_owned_index(res,idx,desc,info)
    implicit none 
    integer(psb_ipk_), intent(in) :: idx
    type(psb_desc_type), intent(in) :: desc
    logical, intent(out) :: res
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: lx

    call psb_glob_to_loc(idx,lx,desc,info,iact='I',owned=.true.)

    res = (lx>0)
  end subroutine psb_owned_index

  subroutine psb_owned_index_v(res,idx,desc,info)
    implicit none 
    integer(psb_ipk_), intent(in) :: idx(:)
    type(psb_desc_type), intent(in) :: desc
    logical, intent(out) :: res(:)
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), allocatable  :: lx(:)

    allocate(lx(size(idx)),stat=info)
    res=.false.
    if (info /= psb_success_) return
    call psb_glob_to_loc(idx,lx,desc,info,iact='I',owned=.true.)

    res = (lx>0)
  end subroutine psb_owned_index_v

  subroutine psb_local_index(res,idx,desc,info)
    implicit none 
    integer(psb_ipk_), intent(in) :: idx
    type(psb_desc_type), intent(in) :: desc
    logical, intent(out) :: res
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: lx

    call psb_glob_to_loc(idx,lx,desc,info,iact='I',owned=.false.)

    res = (lx>0)
  end subroutine psb_local_index

  subroutine psb_local_index_v(res,idx,desc,info)
    implicit none 
    integer(psb_ipk_), intent(in) :: idx(:)
    type(psb_desc_type), intent(in) :: desc
    logical, intent(out) :: res(:)
    integer(psb_ipk_), intent(out) :: info    
    integer(psb_ipk_), allocatable  :: lx(:)

    allocate(lx(size(idx)),stat=info)
    res=.false.
    if (info /= psb_success_) return
    call psb_glob_to_loc(idx,lx,desc,info,iact='I',owned=.false.)

    res = (lx>0)
  end subroutine psb_local_index_v

end module psb_iv_tools_mod

module psb_cd_if_tools_mod

  use psb_const_mod
  use psb_descriptor_type
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

      type(psb_desc_type), intent(in)  :: desc_in
      type(psb_desc_type), intent(out) :: desc_out
      integer(psb_ipk_), intent(out)             :: info
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
      Type(psb_desc_type), Intent(in), target :: desc_a
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
    subroutine psb_icdasb(desc,info,ext_hv)
      import :: psb_ipk_, psb_desc_type
      Type(psb_desc_type), intent(inout) :: desc
      integer(psb_ipk_), intent(out)               :: info
      logical, intent(in),optional       :: ext_hv
    end subroutine psb_icdasb
  end interface

end module psb_cd_if_tools_mod

module psb_cd_tools_mod

  use psb_const_mod
  use psb_cd_if_tools_mod

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

  subroutine psb_cdasb(desc,info)

    Type(psb_desc_type), intent(inout) :: desc
    integer(psb_ipk_), intent(out)               :: info

    call psb_icdasb(desc,info,ext_hv=.false.)
  end subroutine psb_cdasb

end module psb_cd_tools_mod


module psb_base_tools_mod
  use psb_iv_tools_mod
  use psb_cd_tools_mod
  
end module psb_base_tools_mod
 
