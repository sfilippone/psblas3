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
Module psb_tools_mod
  use psb_const_mod
  use psb_gps_mod
    
  interface  psb_geall
     ! 2-D double precision version
     subroutine psb_dalloc(x, desc_a, info, n)
       use psb_descriptor_type
       implicit none
       real(kind(1.d0)), allocatable, intent(out) :: x(:,:)
       type(psb_desc_type), intent(in) :: desc_a
       integer,intent(out)             :: info
       integer, optional, intent(in)   :: n
     end subroutine psb_dalloc
     ! 1-D double precision version
     subroutine psb_dallocv(x, desc_a,info,n)
       use psb_descriptor_type
       real(kind(1.d0)), allocatable, intent(out)       :: x(:)
       type(psb_desc_type), intent(in) :: desc_a
       integer,intent(out)             :: info
       integer, optional, intent(in)   :: n
     end subroutine psb_dallocv
     ! 2-D integer version
     subroutine psb_ialloc(x, desc_a, info,n)
       use psb_descriptor_type
       integer, allocatable, intent(out)                 :: x(:,:)
       type(psb_desc_type), intent(in)  :: desc_a
       integer, intent(out)             :: info
       integer, optional, intent(in)    :: n
     end subroutine psb_ialloc
     subroutine psb_iallocv(x, desc_a,info,n)
       use psb_descriptor_type
       integer, allocatable, intent(out)                :: x(:)
       type(psb_desc_type), intent(in) :: desc_a
       integer                         :: info
       integer, optional, intent(in)   :: n
     end subroutine psb_iallocv
     ! 2-D double precision version
     subroutine psb_zalloc(x, desc_a, info, n)
       use psb_descriptor_type
       implicit none
       complex(kind(1.d0)), allocatable, intent(out)    :: x(:,:)
       type(psb_desc_type), intent(in) :: desc_a
       integer                         :: info
       integer, optional, intent(in)   :: n
     end subroutine psb_zalloc
     ! 1-D double precision version
     subroutine psb_zallocv(x, desc_a,info,n)
       use psb_descriptor_type
       complex(kind(1.d0)), allocatable, intent(out)    :: x(:)
       type(psb_desc_type), intent(in) :: desc_a
       integer                         :: info
       integer, optional, intent(in)   :: n
     end subroutine psb_zallocv
  end interface


  interface psb_geasb
     ! 2-D double precision version
     subroutine psb_dasb(x, desc_a, info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in) ::  desc_a
       real(kind(1.d0)), allocatable, intent(inout)       ::  x(:,:)
       integer, intent(out)            ::  info
     end subroutine psb_dasb
     ! 1-D double precision version
     subroutine psb_dasbv(x, desc_a, info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in) ::  desc_a
       real(kind(1.d0)), allocatable, intent(inout)   ::  x(:)
       integer, intent(out)        ::  info
     end subroutine psb_dasbv
     ! 2-D integer version
     subroutine psb_iasb(x, desc_a, info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in) ::  desc_a
       integer, allocatable, intent(inout)                ::  x(:,:)
       integer, intent(out)            ::  info
     end subroutine psb_iasb
     ! 1-D integer version
     subroutine psb_iasbv(x, desc_a, info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in) ::  desc_a
       integer, allocatable, intent(inout)   ::  x(:)
       integer, intent(out)        ::  info
     end subroutine psb_iasbv
     ! 2-D double precision version
     subroutine psb_zasb(x, desc_a, info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in) ::  desc_a
       complex(kind(1.d0)), allocatable, intent(inout)       ::  x(:,:)
       integer, intent(out)            ::  info
     end subroutine psb_zasb
     ! 1-D double precision version
     subroutine psb_zasbv(x, desc_a, info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in) ::  desc_a
       complex(kind(1.d0)), allocatable, intent(inout)   ::  x(:)
       integer, intent(out)        ::  info
     end subroutine psb_zasbv
   end interface

  interface psb_sphalo
     Subroutine psb_dsphalo(a,desc_a,blk,info,rwcnv,clcnv,outfmt)
       use psb_descriptor_type
       use psb_spmat_type
       Type(psb_dspmat_type),Intent(in)    :: a
       Type(psb_dspmat_type),Intent(inout) :: blk
       Type(psb_desc_type),Intent(in)      :: desc_a
       integer, intent(out)                :: info
       logical, optional, intent(in)       :: rwcnv,clcnv
       character(len=5), optional          :: outfmt 
     end Subroutine psb_dsphalo
     Subroutine psb_zsphalo(a,desc_a,blk,info,rwcnv,clcnv,outfmt)
       use psb_descriptor_type
       use psb_spmat_type
       Type(psb_zspmat_type),Intent(in)    :: a
       Type(psb_zspmat_type),Intent(inout) :: blk
       Type(psb_desc_type),Intent(in)      :: desc_a
       integer, intent(out)                :: info
       logical, optional, intent(in)       :: rwcnv,clcnv
       character(len=5), optional          :: outfmt 
     end Subroutine psb_zsphalo
  end interface


  interface psb_csrp
     subroutine psb_dcsrp(trans,iperm,a, desc_a, info)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_dspmat_type), intent(inout)  ::  a
       type(psb_desc_type), intent(in)       ::  desc_a
       integer, intent(inout)                :: iperm(:), info
       character, intent(in)                 :: trans
     end subroutine psb_dcsrp
  end interface


  interface psb_cdprt
     subroutine psb_cdprt(iout,desc_p,glob,short)
       use psb_const_mod
       use psb_descriptor_type
       implicit none 
       type(psb_desc_type), intent(in)    :: desc_p
       integer, intent(in)                :: iout
       logical, intent(in), optional      :: glob,short
     end subroutine psb_cdprt
  end interface


  interface psb_gefree
     ! 2-D double precision version
     subroutine psb_dfree(x, desc_a, info)
       use psb_descriptor_type
       real(kind(1.d0)),allocatable, intent(inout)        :: x(:,:)
       type(psb_desc_type), intent(in) :: desc_a
       integer                         :: info
     end subroutine psb_dfree
     ! 1-D double precision version
     subroutine psb_dfreev(x, desc_a, info)
       use psb_descriptor_type
       real(kind(1.d0)),allocatable, intent(inout)        :: x(:)
       type(psb_desc_type), intent(in) :: desc_a
       integer                         :: info
     end subroutine psb_dfreev
     ! 2-D integer version
     subroutine psb_ifree(x, desc_a, info)
       use psb_descriptor_type
       integer,allocatable, intent(inout)                 :: x(:,:)
       type(psb_desc_type), intent(in) :: desc_a
       integer                         :: info
     end subroutine psb_ifree
     ! 1-D integer version
     subroutine psb_ifreev(x, desc_a, info)
       use psb_descriptor_type
       integer, allocatable, intent(inout)                :: x(:)
       type(psb_desc_type), intent(in) :: desc_a
       integer                         :: info
     end subroutine psb_ifreev
     ! 2-D double precision version
     subroutine psb_zfree(x, desc_a, info)
       use psb_descriptor_type
       complex(kind(1.d0)),allocatable, intent(inout)        :: x(:,:)
       type(psb_desc_type), intent(in) :: desc_a
       integer                         :: info
     end subroutine psb_zfree
     ! 1-D double precision version
     subroutine psb_zfreev(x, desc_a, info)
       use psb_descriptor_type
       complex(kind(1.d0)),allocatable, intent(inout)        :: x(:)
       type(psb_desc_type), intent(in) :: desc_a
       integer                         :: info
     end subroutine psb_zfreev
  end interface


  interface psb_gelp
     ! 2-D version
     subroutine psb_dgelp(trans,iperm,x,desc_a,info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in)      ::  desc_a
       real(kind(1.d0)), intent(inout)      ::  x(:,:)
       integer, intent(inout)               ::  iperm(:),info
       character, intent(in)                :: trans
     end subroutine psb_dgelp
     ! 1-D version
     subroutine psb_dgelpv(trans,iperm,x,desc_a,info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in) ::  desc_a
       real(kind(1.d0)), intent(inout)    ::  x(:)
       integer, intent(inout)             ::  iperm(:), info
       character, intent(in)              :: trans
     end subroutine psb_dgelpv
     ! 2-D version
     subroutine psb_zgelp(trans,iperm,x,desc_a,info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in)      ::  desc_a
       complex(kind(1.d0)), intent(inout)      ::  x(:,:)
       integer, intent(inout)               ::  iperm(:),info
       character, intent(in)                :: trans
     end subroutine psb_zgelp
     ! 1-D version
     subroutine psb_zgelpv(trans,iperm,x,desc_a,info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in) ::  desc_a
       complex(kind(1.d0)), intent(inout)    ::  x(:)
       integer, intent(inout)             ::  iperm(:), info
       character, intent(in)              :: trans
     end subroutine psb_zgelpv
  end interface


  interface psb_geins
     ! 2-D double precision version
     subroutine psb_dinsi(m,irw,val, x,desc_a,info,dupl)
       use psb_descriptor_type
       integer, intent(in)                ::  m
       type(psb_desc_type), intent(in)    ::  desc_a
       real(kind(1.d0)),intent(inout)           ::  x(:,:)
       integer, intent(in)                ::  irw(:)
       real(kind(1.d0)), intent(in)       ::  val(:,:)
       integer, intent(out)               ::  info
       integer, optional, intent(in)      ::  dupl
     end subroutine psb_dinsi
     ! 1-D double precision version
     subroutine psb_dinsvi(m,irw,val,x,desc_a,info,dupl)
       use psb_descriptor_type
       integer, intent(in)                ::  m
       type(psb_desc_type), intent(in)    ::  desc_a
       real(kind(1.d0)),intent(inout)           ::  x(:)
       integer, intent(in)                ::  irw(:)
       real(kind(1.d0)), intent(in)       ::  val(:)
       integer, intent(out)               ::  info
       integer, optional, intent(in)      ::  dupl
     end subroutine psb_dinsvi
     ! 2-D double precision version
     subroutine psb_iinsi(m,irw,val, x,desc_a,info,dupl)
       use psb_descriptor_type
       integer, intent(in)              ::  m
       type(psb_desc_type), intent(in)  ::  desc_a
       integer,intent(inout)                  ::  x(:,:)
       integer, intent(in)              ::  irw(:)
       integer, intent(in)              ::  val(:,:)
       integer, intent(out)             ::  info
       integer, optional, intent(in)    ::  dupl
     end subroutine psb_iinsi
     ! 1-D double precision version
     subroutine psb_iinsvi(m, irw,val, x,desc_a,info,dupl)
       use psb_descriptor_type
       integer, intent(in)             ::  m
       type(psb_desc_type), intent(in) ::  desc_a
       integer,intent(inout)                 ::  x(:)
       integer, intent(in)             ::  irw(:)
       integer, intent(in)             ::  val(:)
       integer, intent(out)            ::  info
       integer, optional, intent(in)   ::  dupl
     end subroutine psb_iinsvi
     ! 2-D double precision version
     subroutine psb_zinsi(m,irw,val, x, desc_a,info,dupl)
       use psb_descriptor_type
       integer, intent(in)              ::  m
       type(psb_desc_type), intent(in)  ::  desc_a
       complex(kind(1.d0)),intent(inout)      ::  x(:,:)
       integer, intent(in)              ::  irw(:)
       complex(kind(1.d0)), intent(in)  ::  val(:,:)
       integer, intent(out)             ::  info
       integer, optional, intent(in)    ::  dupl
     end subroutine psb_zinsi
     ! 1-D double precision version
     subroutine psb_zinsvi(m, irw,val, x,desc_a,info,dupl)
       use psb_descriptor_type
       integer, intent(in)              ::  m
       type(psb_desc_type), intent(in)  ::  desc_a
       complex(kind(1.d0)),intent(inout)      ::  x(:)
       integer, intent(in)              ::  irw(:)
       complex(kind(1.d0)), intent(in)  ::  val(:)
       integer, intent(out)             ::  info
       integer, optional, intent(in)    ::  dupl
     end subroutine psb_zinsvi
  end interface


  interface psb_cdall
    module procedure psb_cdall
  end interface

  interface psb_cdrep
     subroutine psb_cdrep(m, ictxt, desc_a,info)
       use psb_descriptor_type
       Integer, intent(in)               :: m,ictxt
       Type(psb_desc_type), intent(out)  :: desc_a
       integer, intent(out)              :: info
     end subroutine psb_cdrep
  end interface

  interface psb_cdasb
     subroutine psb_cdasb(desc_a,info)
       use psb_descriptor_type
       Type(psb_desc_type), intent(inout) :: desc_a
       integer, intent(out)               :: info
     end subroutine psb_cdasb
  end interface



  interface psb_cdcpy
     subroutine psb_cdcpy(desc_in, desc_out, info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in)  :: desc_in
       type(psb_desc_type), intent(out) :: desc_out
       integer, intent(out)             :: info
     end subroutine psb_cdcpy
  end interface

  interface psb_cdtransfer
     subroutine psb_cdtransfer(desc_in, desc_out, info)
       use psb_descriptor_type
       type(psb_desc_type), intent(inout) :: desc_in
       type(psb_desc_type), intent(inout)   :: desc_out
       integer, intent(out)               :: info
     end subroutine psb_cdtransfer
  end interface
  
 
  interface psb_cdfree
     subroutine psb_cdfree(desc_a,info)
       use psb_descriptor_type
       type(psb_desc_type), intent(inout) :: desc_a
       integer, intent(out)               :: info
     end subroutine psb_cdfree
  end interface
  
  interface psb_cdins
     subroutine psb_cdins(nz,ia,ja,desc_a,info,ila,jla)
       use psb_descriptor_type
       type(psb_desc_type), intent(inout) :: desc_a
       integer, intent(in)                :: nz,ia(:),ja(:)
       integer, intent(out)               :: info
       integer, optional, intent(out)     :: ila(:), jla(:)
     end subroutine psb_cdins
  end interface


  interface psb_cdbldovr
     Subroutine psb_dcdovr(a,desc_a,novr,desc_ov,info)
       use psb_descriptor_type
       Use psb_spmat_type
       integer, intent(in)                :: novr
       Type(psb_dspmat_type), Intent(in)  ::  a
       Type(psb_desc_type), Intent(in)    :: desc_a
       Type(psb_desc_type), Intent(inout) :: desc_ov
       integer, intent(out)               :: info
     end Subroutine psb_dcdovr
     Subroutine psb_zcdovr(a,desc_a,novr,desc_ov,info)
       use psb_descriptor_type
       Use psb_spmat_type
       integer, intent(in)                :: novr
       Type(psb_zspmat_type), Intent(in)  ::  a
       Type(psb_desc_type), Intent(in)    :: desc_a
       Type(psb_desc_type), Intent(inout) :: desc_ov
       integer, intent(out)               :: info
     end Subroutine psb_zcdovr
  end interface

  interface psb_cdren
     subroutine psb_cdren(trans,iperm,desc_a,info)
       use psb_descriptor_type
       type(psb_desc_type), intent(inout)    ::  desc_a
       integer, intent(inout)                ::  iperm(:)
       character, intent(in)                 :: trans
       integer, intent(out)                  :: info
     end subroutine psb_cdren
  end interface
  
  interface psb_spall
     subroutine psb_dspalloc(a, desc_a, info, nnz)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_desc_type), intent(inout) :: desc_a
       type(psb_dspmat_type), intent(out) :: a
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: nnz
     end subroutine psb_dspalloc
     subroutine psb_zspalloc(a, desc_a, info, nnz)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_desc_type), intent(inout) :: desc_a
       type(psb_zspmat_type), intent(out) :: a
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: nnz
     end subroutine psb_zspalloc
  end interface

  interface psb_spasb
     subroutine psb_dspasb(a,desc_a, info, afmt, upd, dupl)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_dspmat_type), intent (inout)   :: a
       type(psb_desc_type), intent(in)         :: desc_a
       integer, intent(out)                    :: info
       integer,optional, intent(in)            :: dupl, upd
       character, optional, intent(in)         :: afmt*5
     end subroutine psb_dspasb
     subroutine psb_zspasb(a,desc_a, info, afmt, upd, dupl)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_zspmat_type), intent (inout)   :: a
       type(psb_desc_type), intent(in)         :: desc_a
       integer, intent(out)                    :: info
       integer,optional, intent(in)            :: dupl, upd
       character, optional, intent(in)         :: afmt*5
     end subroutine psb_zspasb
  end interface


  interface psb_spcnv
     subroutine psb_dspcnv(a,b,desc_a,info)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_dspmat_type), intent(in)   :: a
       type(psb_dspmat_type), intent(out)  :: b
       type(psb_desc_type), intent(in)     :: desc_a
       integer, intent(out)                :: info
     end subroutine psb_dspcnv
     subroutine psb_zspcnv(a,b,desc_a,info)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_zspmat_type), intent(in)   :: a
       type(psb_zspmat_type), intent(out)  :: b
       type(psb_desc_type), intent(in)     :: desc_a
       integer, intent(out)                :: info
     end subroutine psb_zspcnv
  end interface


  interface psb_spfree
     subroutine psb_dspfree(a, desc_a,info)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_desc_type), intent(in) :: desc_a
       type(psb_dspmat_type), intent(inout)       ::a
       integer, intent(out)        :: info
     end subroutine psb_dspfree
     subroutine psb_zspfree(a, desc_a,info)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_desc_type), intent(in) :: desc_a
       type(psb_zspmat_type), intent(inout)       ::a
       integer, intent(out)        :: info
     end subroutine psb_zspfree
  end interface


  interface psb_spins
     subroutine psb_dspins(nz,ia,ja,val,a,desc_a,info,rebuild)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_desc_type), intent(inout)   :: desc_a
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(in)                  :: nz,ia(:),ja(:)
       real(kind(1.d0)), intent(in)         :: val(:)
       integer, intent(out)                 :: info
       logical, intent(in), optional        :: rebuild
     end subroutine psb_dspins
     subroutine psb_zspins(nz,ia,ja,val,a,desc_a,info,rebuild)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_desc_type), intent(inout)   :: desc_a
       type(psb_zspmat_type), intent(inout) :: a
       integer, intent(in)                  :: nz,ia(:),ja(:)
       complex(kind(1.d0)), intent(in)      :: val(:)
       integer, intent(out)                 :: info
       logical, intent(in), optional        :: rebuild
     end subroutine psb_zspins
  end interface


  interface psb_sprn
     subroutine psb_dsprn(a, desc_a,info,clear)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_desc_type), intent(in)      :: desc_a
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(out)                 :: info
       logical, intent(in), optional        :: clear
     end subroutine psb_dsprn
     subroutine psb_zsprn(a, desc_a,info,clear)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_desc_type), intent(in)      :: desc_a
       type(psb_zspmat_type), intent(inout) :: a
       integer, intent(out)                 :: info
       logical, intent(in), optional        :: clear
     end subroutine psb_zsprn
  end interface


  interface psb_glob_to_loc
     subroutine psb_glob_to_loc2(x,y,desc_a,info,iact)
       use psb_descriptor_type
       type(psb_desc_type), intent(in)    ::  desc_a
       integer,intent(in)                 ::  x(:)  
       integer,intent(out)                ::  y(:)  
       integer, intent(out)               ::  info
       character, intent(in), optional    ::  iact
     end subroutine psb_glob_to_loc2
     subroutine psb_glob_to_loc(x,desc_a,info,iact)
       use psb_descriptor_type
       type(psb_desc_type), intent(in)    ::  desc_a
       integer,intent(inout)              ::  x(:)  
       integer, intent(out)               ::  info
       character, intent(in), optional    ::  iact
     end subroutine psb_glob_to_loc
  end interface

  interface psb_loc_to_glob
     subroutine psb_loc_to_glob2(x,y,desc_a,info,iact)
       use psb_descriptor_type
       type(psb_desc_type), intent(in)    ::  desc_a
       integer,intent(in)                 ::  x(:)  
       integer,intent(out)                ::  y(:)  
       integer, intent(out)               ::  info
       character, intent(in), optional    ::  iact
     end subroutine psb_loc_to_glob2
     subroutine psb_loc_to_glob(x,desc_a,info,iact)
       use psb_descriptor_type
       type(psb_desc_type), intent(in)    ::  desc_a
       integer,intent(inout)              ::  x(:)  
       integer, intent(out)               ::  info
       character, intent(in), optional    ::  iact
     end subroutine psb_loc_to_glob
  end interface


  interface psb_get_boundary
    module procedure psb_get_boundary
  end interface
  
  interface psb_get_overlap
    subroutine psb_get_ovrlap(ovrel,desc,info)
      use psb_descriptor_type
      implicit none 
      integer, allocatable            :: ovrel(:)
      type(psb_desc_type), intent(in) :: desc
      integer, intent(out)            :: info
    end subroutine psb_get_ovrlap
  end interface
  


contains

  subroutine psb_get_boundary(bndel,desc,info)
    use psb_descriptor_type
    use psi_mod
    implicit none 
    integer, allocatable            :: bndel(:)
    type(psb_desc_type), intent(in) :: desc
    integer, intent(out)            :: info
    
    call psi_crea_bnd_elem(bndel,desc,info)

  end subroutine psb_get_boundary
  
  subroutine psb_cdall(ictxt, desc_a, info,mg,ng,parts,vg,vl,flag,nl)
    use psb_descriptor_type
    use psb_serial_mod
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit None
    include 'parts.fh'
    Integer, intent(in)               :: mg,ng,ictxt, vg(:), vl(:),nl
    integer, intent(in)               :: flag
    integer, intent(out)              :: info
    type(psb_desc_type), intent(out)  :: desc_a

    optional :: mg,ng,parts,vg,vl,flag,nl

    interface 
      subroutine psb_cdals(m, n, parts, ictxt, desc_a, info)
        use psb_descriptor_type
        include 'parts.fh'
        Integer, intent(in)                 :: m,n,ictxt
        Type(psb_desc_type), intent(out)    :: desc_a
        integer, intent(out)                :: info
      end subroutine psb_cdals
      subroutine psb_cdalv(v, ictxt, desc_a, info, flag)
        use psb_descriptor_type
        Integer, intent(in)               :: ictxt, v(:)
        integer, intent(in), optional     :: flag
        integer, intent(out)              :: info
        Type(psb_desc_type), intent(out)  :: desc_a
      end subroutine psb_cdalv
      subroutine psb_cd_inloc(v, ictxt, desc_a, info)
        use psb_descriptor_type
        implicit None
        Integer, intent(in)               :: ictxt, v(:)
        integer, intent(out)              :: info
        type(psb_desc_type), intent(out)  :: desc_a
      end subroutine psb_cd_inloc
    end interface
    character(len=20)   :: name, char_err
    integer :: err_act, n_, flag_, i, me, np, nlp
    integer, allocatable :: itmpsz(:) 



    if(psb_get_errstatus() /= 0) return 
    info=0
    name = 'psb_cdall'
    call psb_erractionsave(err_act)
    
    call psb_info(ictxt, me, np)

    if (count((/ present(vg),present(vl),present(parts),present(nl) /)) /= 1) then 
      info=581
      call psb_errpush(info,name,a_err=" vg, vl, parts, nl")
      goto 999 
    endif
    
    if (present(parts)) then 
      if (.not.present(mg)) then 
        info=581
        call psb_errpush(info,name)
        goto 999 
      end if
      if (present(ng)) then 
        n_ = ng
      else
        n_ = mg 
      endif
      call  psb_cdals(mg, n_, parts, ictxt, desc_a, info)

    else if (present(vg)) then 
      if (present(flag)) then 
        flag_=flag
      else
        flag_=0
      endif
      call psb_cdalv(vg, ictxt, desc_a, info, flag_)

    else if (present(vl)) then 
      call psb_cd_inloc(vl,ictxt,desc_a,info)
      
    else if (present(nl)) then 
      allocate(itmpsz(0:np-1),stat=info)
      if (info /= 0) then 
        info = 4000 
        call psb_errpush(info,name)
        goto 999
      endif

      itmpsz = 0
      itmpsz(me) = nl
      call psb_sum(ictxt,itmpsz)
      nlp=0 
      do i=0, me-1
        nlp = nlp + itmpsz(me)
      end do
      call psb_cd_inloc((/(i,i=nlp+1,nlp+nl)/),ictxt,desc_a,info)
      
    endif
    call psb_erractionrestore(err_act)
    return
    
999 continue
    call psb_erractionrestore(err_act)
    if (err_act == act_abort) then
      call psb_error(ictxt)
      return
    end if
    return
    

  end subroutine psb_cdall





end module psb_tools_mod
