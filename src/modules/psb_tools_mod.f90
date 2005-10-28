Module psb_tools_mod
  use psb_const_mod

  interface  psb_alloc
     ! 2-D double precision version
     subroutine psb_dalloc(m, n, x, desc_a, info, js)
       use psb_descriptor_type
       implicit none
       integer, intent(in)                   :: m,n
       real(kind(1.d0)), pointer             :: x(:,:)
       type(psb_desc_type), intent(in)       :: desc_a
       integer                               :: info
       integer, optional, intent(in)         :: js
     end subroutine psb_dalloc
     ! 1-D double precision version
     subroutine psb_dallocv(m, x, desc_a,info)
       use psb_descriptor_type
       integer, intent(in)            :: m
       real(kind(1.d0)), pointer      :: x(:)
       type(psb_desc_type), intent(in):: desc_a
       integer                        :: info
     end subroutine psb_dallocv
     ! 2-D integer version
     subroutine psb_ialloc(m, n, x, desc_a, info,js)
       use psb_descriptor_type
       integer, intent(in)                   :: m,n
       integer, pointer                      :: x(:,:)
       type(psb_desc_type), intent(inout)    :: desc_a
       integer, intent(out)                  :: info
       integer, optional, intent(in)         :: js
     end subroutine psb_ialloc
     subroutine psb_iallocv(m, x, desc_a,info)
       use psb_descriptor_type
       integer, intent(in)            :: m
       integer, pointer               :: x(:)
       type(psb_desc_type), intent(in):: desc_a
       integer                        :: info
     end subroutine psb_iallocv
  end interface


  interface psb_asb
     ! 2-D double precision version
     subroutine psb_dasb(x, desc_a, info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in) ::  desc_a
       real(kind(1.d0)), pointer       ::  x(:,:)
       integer, intent(out)            ::  info
     end subroutine psb_dasb
     ! 1-D double precision version
     subroutine psb_dasbv(x, desc_a, info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in) ::  desc_a
       real(kind(1.d0)), pointer   ::  x(:)
       integer, intent(out)        ::  info
     end subroutine psb_dasbv
     ! 2-D integer version
     subroutine psb_iasb(x, desc_a, info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in) ::  desc_a
       integer, pointer                ::  x(:,:)
       integer, intent(out)            ::  info
     end subroutine psb_iasb
     ! 1-D integer version
     subroutine psb_iasbv(x, desc_a, info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in) ::  desc_a
       integer, pointer   ::  x(:)
       integer, intent(out)        ::  info
     end subroutine psb_iasbv
  end interface

  interface psb_csrovr
     Subroutine psb_dcsrovr(a,desc_a,blk,info,rwcnv,clcnv,outfmt)
       use psb_descriptor_type
       use psb_spmat_type
       Type(psb_dspmat_type),Intent(in)    :: a
       Type(psb_dspmat_type),Intent(inout) :: blk
       Type(psb_desc_type),Intent(in)      :: desc_a
       integer, intent(out)                :: info
       logical, optional, intent(in)       :: rwcnv,clcnv
       character(len=5), optional          :: outfmt 
     end Subroutine psb_dcsrovr
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


  interface psb_descasb
     Subroutine psb_descasb(n_ovr,desc_p,desc_a,a,&
          &       l_tmp_halo,l_tmp_ovr_idx,lworks,lworkr,info)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_dspmat_type),intent(in)  :: a
       type(psb_desc_type),intent(in)    :: desc_a
       type(psb_desc_type),intent(inout) :: desc_p
       integer,intent(in)                :: n_ovr
       Integer, Intent(in)               :: l_tmp_halo,l_tmp_ovr_idx
       Integer, Intent(inout)            :: lworks, lworkr
       integer, intent(out)              :: info
     end Subroutine psb_descasb
  end interface


  interface psb_descprt
     subroutine psb_descprt(iout,desc_p,glob,short)
       use psb_const_mod
       use psb_descriptor_type
       implicit none 
       type(psb_desc_type), intent(in)    :: desc_p
       integer, intent(in)                :: iout
       logical, intent(in), optional      :: glob,short
     end subroutine psb_descprt
  end interface


  interface psb_free
     ! 2-D double precision version
     subroutine psb_dfree(x, desc_a, info)
       use psb_descriptor_type
       real(kind(1.d0)),pointer        :: x(:,:)
       type(psb_desc_type), intent(in) :: desc_a
       integer                         :: info
     end subroutine psb_dfree
     ! 1-D double precision version
     subroutine psb_dfreev(x, desc_a, info)
       use psb_descriptor_type
       real(kind(1.d0)),pointer        :: x(:)
       type(psb_desc_type), intent(in) :: desc_a
       integer                         :: info
     end subroutine psb_dfreev
     ! 2-D integer version
     subroutine psb_ifree(x, desc_a, info)
       use psb_descriptor_type
       integer,pointer                 :: x(:,:)
       type(psb_desc_type), intent(in) :: desc_a
       integer                         :: info
     end subroutine psb_ifree
     ! 1-D integer version
     subroutine psb_ifreev(x, desc_a, info)
       use psb_descriptor_type
       integer, pointer                :: x(:)
       type(psb_desc_type), intent(in) :: desc_a
       integer                         :: info
     end subroutine psb_ifreev
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
  end interface


  interface psb_ins
     ! 2-D double precision version
     subroutine psb_dins(m, n, x, ix, jx, blck, desc_a, info,&
          & iblck, jblck)
       use psb_descriptor_type
       integer, intent(in)                ::  m,n
       type(psb_desc_type), intent(in)    ::  desc_a
       real(kind(1.d0)),pointer           ::  x(:,:)
       integer, intent(in)                ::  ix,jx
       real(kind(1.d0)), intent(in)       ::  blck(:,:)
       integer,intent(out)                ::  info
       integer, optional, intent(in)      ::  iblck,jblck
     end subroutine psb_dins
     ! 2-D double precision square version
     subroutine psb_dinsvm(m, x, ix, jx, blck, desc_a,info,&
          & iblck)
       use psb_descriptor_type
       integer, intent(in)                ::  m
       type(psb_desc_type), intent(in) ::  desc_a
       real(kind(1.d0)),pointer           ::  x(:,:)
       integer, intent(in)                ::  ix,jx
       real(kind(1.d0)), intent(in)       ::  blck(:)
       integer, intent(out)               ::  info
       integer, optional, intent(in)      ::  iblck
     end subroutine psb_dinsvm
     ! 1-D double precision version
     subroutine psb_dinsvv(m, x, ix, blck, desc_a, info,&
          & iblck,insflag)
       use psb_descriptor_type
       integer, intent(in)                ::  m
       type(psb_desc_type), intent(in)    ::  desc_a
       real(kind(1.d0)),pointer           ::  x(:)
       integer, intent(in)                ::  ix
       real(kind(1.d0)), intent(in)       ::  blck(:)
       integer, intent(out)               ::  info
       integer, optional, intent(in)      ::  iblck
       integer, optional, intent(in)      ::  insflag
     end subroutine psb_dinsvv
     ! 2-D integer version
     subroutine psb_iins(m, n, x, ix, jx, blck, desc_a, info,&
          & iblck, jblck)
       use psb_descriptor_type
       integer, intent(in)                ::  m,n
       type(psb_desc_type), intent(in)    ::  desc_a
       integer,pointer                    ::  x(:,:)
       integer, intent(in)                ::  ix,jx
       integer, intent(in)                ::  blck(:,:)
       integer,intent(out)                ::  info
       integer, optional, intent(in)      ::  iblck,jblck
     end subroutine psb_iins
     ! 2-D integer square version
     subroutine psb_iinsvm(m, x, ix, jx, blck, desc_a,info,&
          & iblck)
       use psb_descriptor_type
       integer, intent(in)                ::  m
       type(psb_desc_type), intent(in)    ::  desc_a
       integer, pointer                   ::  x(:,:)
       integer, intent(in)                ::  ix,jx
       integer, intent(in)                ::  blck(:)
       integer, intent(out)               ::  info
       integer, optional, intent(in)      ::  iblck
     end subroutine psb_iinsvm
     ! 1-D integer version
     subroutine psb_iinsvv(m, x, ix, blck, desc_a, info,&
          & iblck,insflag)
       use psb_descriptor_type
       integer, intent(in)                ::  m
       type(psb_desc_type), intent(in)    ::  desc_a
       integer, pointer                   ::  x(:)
       integer, intent(in)                ::  ix
       integer, intent(in)                ::  blck(:)
       integer, intent(out)               ::  info
       integer, optional, intent(in)      ::  iblck
       integer, optional, intent(in)      ::  insflag
     end subroutine psb_iinsvv
  end interface



  interface psb_ptins
     subroutine psb_dptins(ia,ja,blck,desc_a,info)
       use psb_descriptor_type
       use psb_spmat_type
       implicit none
       type(psb_desc_type), intent(inout)    ::  desc_a
       integer, intent(in)                   ::  ia,ja
       type(psb_dspmat_type), intent(in)             ::  blck
       integer,intent(out)                   ::  info
     end subroutine psb_dptins
  end interface

  interface psb_dscall
     subroutine psb_dscall(m, n, parts, icontxt, desc_a, info)
       use psb_descriptor_type
       include 'parts.fh'
       Integer, intent(in)                 :: m,n,icontxt
       Type(psb_desc_type), intent(out)    :: desc_a
       integer, intent(out)                :: info
     end subroutine psb_dscall
     subroutine psb_dscalv(m, v, icontxt, desc_a, info, flag)
       use psb_descriptor_type
       Integer, intent(in)               :: m,icontxt, v(:)
       integer, intent(in), optional     :: flag
       integer, intent(out)              :: info
       Type(psb_desc_type), intent(out)  :: desc_a
     end subroutine psb_dscalv
  end interface
  

  interface psb_dscasb
     subroutine psb_dscasb(desc_a,info)
       use psb_descriptor_type
       Type(psb_desc_type), intent(inout) :: desc_a
       integer, intent(out)               :: info
     end subroutine psb_dscasb
  end interface



  interface psb_dsccpy
     subroutine psb_dsccpy(desc_out, desc_a, info)
       use psb_descriptor_type
       type(psb_desc_type), intent(out) :: desc_out
       type(psb_desc_type), intent(in)  :: desc_a
       integer, intent(out)             :: info
     end subroutine psb_dsccpy
  end interface
  
 
  interface psb_dscfree
     subroutine psb_dscfree(desc_a,info)
       use psb_descriptor_type
       type(psb_desc_type), intent(inout) :: desc_a
       integer, intent(out)               :: info
     end subroutine psb_dscfree
  end interface
  
  interface psb_dscins
     subroutine psb_dscins(nz,ia,ja,desc_a,info,is,js)
       use psb_descriptor_type
       type(psb_desc_type), intent(inout) :: desc_a
       Integer, intent(in)                :: nz,IA(:),JA(:)
       integer, intent(out)               :: info
       integer, intent(in), optional      :: is,js
     end subroutine psb_dscins
  end interface


  interface psb_dscov
     Subroutine psb_dscov(a,desc_a,novr,desc_ov,info)
       use psb_descriptor_type
       Use psb_spmat_type
       integer, intent(in)                :: novr
       Type(psb_dspmat_type), Intent(in)  ::  a
       Type(psb_desc_type), Intent(in)    :: desc_a
       Type(psb_desc_type), Intent(inout) :: desc_ov
       integer, intent(out)               :: info
     end Subroutine psb_dscov
  end interface
       
       
  interface psb_dscren
     subroutine psb_dscren(trans,iperm,desc_a,info)
       use psb_descriptor_type
       type(psb_desc_type), intent(inout)    ::  desc_a
       integer, intent(inout)                ::  iperm(:)
       character, intent(in)                 :: trans
       integer, intent(out)                  :: info
     end subroutine psb_dscren
  end interface
  
  interface psb_spalloc
     subroutine psb_dspalloc(a, desc_a, info, nnz)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_desc_type), intent(inout) :: desc_a
       type(psb_dspmat_type), intent(out) :: a
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: nnz
     end subroutine psb_dspalloc
  end interface

  interface psb_spasb
     subroutine psb_dspasb(a,desc_a, info, afmt, up, dup)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_dspmat_type), intent (inout)   :: a
       type(psb_desc_type), intent(in)         :: desc_a
       integer, intent(out)                    :: info
       integer,optional, intent(in)            :: dup
       character, optional, intent(in)         :: afmt*5, up
     end subroutine psb_dspasb
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
  end interface


  interface psb_spfree
     subroutine psb_dspfree(a, desc_a,info)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_desc_type), intent(in) :: desc_a
       type(psb_dspmat_type), intent(inout)       ::a
       integer, intent(out)        :: info
     end subroutine psb_dspfree
     subroutine psb_dspfrees(a,info)
       use psb_spmat_type
       type(psb_dspmat_type), intent(inout)       ::a
       integer, intent(out)        :: info
     end subroutine psb_dspfrees
  end interface


  interface psb_spins
     subroutine psb_dspins(nz,ia,ja,val,a,desc_a,info,is,js)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_desc_type), intent(inout)   :: desc_a
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(in)                  :: nz,ia(:),ja(:)
       real(kind(1.d0)), intent(in)         :: val(:)
       integer, intent(out)                 :: info
       integer, intent(in), optional        :: is,js
     end subroutine psb_dspins
  end interface


  interface psb_sprn
     subroutine psb_dsprn(a, desc_a,info)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_desc_type), intent(in)      :: desc_a
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(out)                 :: info
     end subroutine psb_dsprn
  end interface


  interface psb_spupdate
     subroutine psb_dspupdate(a, ia, ja, blck, desc_a,info,ix,jx,updflag)
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_desc_type), intent(in)      ::  desc_a
       type(psb_dspmat_type), intent(inout) ::  a
       integer, intent(in)                  ::  ia,ja
       type(psb_dspmat_type), intent(in)    ::  blck
       integer, intent(out)                 ::  info
       integer, optional, intent(in)        ::  ix,jx
       integer, optional, intent(in)        ::  updflag
     end subroutine psb_dspupdate
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


  interface psb_ptasb
     subroutine psb_ptasb(desc_a,info)
       use psb_descriptor_type
       type(psb_desc_type), intent(inout) :: desc_a
       integer,intent(out)                :: info
     end subroutine psb_ptasb
  end interface


  interface psb_dscrep
     subroutine psb_dscrep(m, icontxt, desc_a,info)
       use psb_descriptor_type
       Integer, intent(in)               :: m,icontxt
       Type(psb_desc_type), intent(out)  :: desc_a
       integer, intent(out)              :: info
     end subroutine psb_dscrep
  end interface

  interface psb_dscdec
     subroutine psb_dscdec(nloc, icontxt, desc_a,info)
       use psb_descriptor_type
       Integer, intent(in)               :: nloc,icontxt
       Type(psb_desc_type), intent(out)  :: desc_a
       integer, intent(out)              :: info
     end subroutine psb_dscdec
  end interface


end module psb_tools_mod
