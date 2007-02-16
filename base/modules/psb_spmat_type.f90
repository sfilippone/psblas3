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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	Module to   define D_SPMAT, structure   !!
!!      for sparse matrix.                      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module psb_spmat_type
  use psb_error_mod
  use psb_realloc_mod
  use psb_const_mod
  implicit none 
  ! Typedef: psb_dspmat_type
  !    Contains a sparse matrix

  !
  !     Queries into spmat%info
  !     
  integer, parameter :: psb_nztotreq_=1, psb_nzrowreq_=2
  integer, parameter :: psb_nzsizereq_=3
  !
  !     Entries and values for  spmat%info
  !     

  integer, parameter :: psb_nnz_=1
  integer, parameter :: psb_del_bnd_=7, psb_srtd_=8
  integer, parameter :: psb_state_=9
  integer, parameter :: psb_upd_pnt_=10
  integer, parameter :: psb_dupl_=11,  psb_upd_=12
  integer, parameter :: psb_ifasize_=16
  integer, parameter :: psb_spmat_null_=0, psb_spmat_bld_=1
  integer, parameter :: psb_spmat_asb_=2, psb_spmat_upd_=4
  integer, parameter :: psb_ireg_flgs_=10, psb_ip2_=0
  integer, parameter :: psb_iflag_=2, psb_ichk_=3
  integer, parameter :: psb_nnzt_=4, psb_zero_=5,psb_ipc_=6
  integer, parameter :: psb_dupl_ovwrt_ = 0
  integer, parameter :: psb_dupl_add_   = 1
  integer, parameter :: psb_dupl_err_   = 2
  integer, parameter :: psb_dupl_def_   = psb_dupl_ovwrt_
  integer, parameter :: psb_upd_dflt_   = 0
  integer, parameter :: psb_upd_srch_   = 98764
  integer, parameter :: psb_upd_perm_   = 98765
  integer, parameter :: psb_isrtdcoo_   = 98761
  integer, parameter :: psb_maxjdrows_=8, psb_minjdrows_=4
  integer, parameter :: psb_dbleint_=2
  character(len=5)   :: psb_fidef_='CSR'



  type psb_dspmat_type
    ! Rows & columns 
    integer     :: m, k
    ! Identify the representation method. Es: CSR, JAD, ...
    character(len=5) :: fida
    ! describe some chacteristics of sparse matrix
    character(len=11) :: descra
    ! Contains some additional informations on sparse matrix
    integer     :: infoa(psb_ifasize_)
    ! Contains sparse matrix coefficients
    real(kind(1.d0)), allocatable  :: aspk(:)
    ! Contains indeces that describes sparse matrix structure
    integer, allocatable :: ia1(:), ia2(:)
    ! Permutations matrix
    integer, allocatable :: pl(:), pr(:)
  end type psb_dspmat_type
  type psb_zspmat_type
    ! Rows & columns 
    integer     :: m, k
    ! Identify the representation method. Es: CSR, JAD, ...
    character(len=5) :: fida
    ! describe some chacteristics of sparse matrix
    character(len=11) :: descra
    ! Contains some additional informations on sparse matrix
    integer     :: infoa(psb_ifasize_)
    ! Contains sparse matrix coefficients
    complex(kind(1.d0)), allocatable  :: aspk(:)
    ! Contains indeces that describes sparse matrix structure
    integer, allocatable  :: ia1(:), ia2(:)
    ! Permutations matrix
    integer, allocatable  :: pl(:), pr(:)
  end type psb_zspmat_type

  interface psb_nullify_sp
    module procedure psb_nullify_dsp, psb_nullify_zsp
  end interface

  interface psb_sp_clone
    module procedure psb_dspclone, psb_zspclone
  end interface

  interface psb_sp_setifld
    module procedure psb_dsp_setifld, psb_zsp_setifld
  end interface

  interface psb_sp_getifld
    module procedure psb_dsp_getifld, psb_zsp_getifld
  end interface

  interface psb_sp_transfer
    module procedure psb_dsp_transfer, psb_zsp_transfer
  end interface

  interface psb_sp_trimsize
    module procedure psb_dsp_trimsize, psb_zsp_trimsize
  end interface

  interface psb_sp_reall
    module procedure psb_dspreallocate, psb_dspreall3, &
         & psb_zspreall3, psb_zspreallocate
  end interface

  interface psb_sp_all
    module procedure psb_dspallocate, psb_dspall3, psb_dspallmk, psb_dspallmknz,&
         &  psb_zspallocate, psb_zspall3, psb_zspallmk, psb_zspallmknz
  end interface

  interface psb_sp_free
    module procedure psb_dsp_free, psb_zsp_free
  end interface

  interface psb_sp_reinit
    module procedure psb_dspreinit,  psb_zspreinit
  end interface

  interface psb_sp_sizeof
    module procedure psb_dspsizeof,  psb_zspsizeof
  end interface

  interface psb_sp_get_nrows
    module procedure psb_get_dsp_nrows, psb_get_zsp_nrows
  end interface

  interface psb_sp_get_ncols
    module procedure psb_get_dsp_ncols, psb_get_zsp_ncols
  end interface

  interface psb_sp_get_nnzeros
    module procedure psb_get_dsp_nnzeros, psb_get_zsp_nnzeros
  end interface

  interface psb_sp_get_nzsize
    module procedure psb_get_dsp_nzsize, psb_get_zsp_nzsize
  end interface

  interface psb_sp_get_nnz_row
    module procedure psb_get_dsp_nnz_row, psb_get_zsp_nnz_row
  end interface



  interface psb_sp_info
    module procedure psb_dspinfo, psb_zspinfo
  end interface



contains

  integer function psb_get_dsp_nrows(a)
    type(psb_dspmat_type), intent(in) :: a
    psb_get_dsp_nrows = a%m

    return
  end function psb_get_dsp_nrows

  integer function psb_get_dsp_ncols(a)
    type(psb_dspmat_type), intent(in) :: a
    psb_get_dsp_ncols = a%k

    return
  end function psb_get_dsp_ncols
  integer function psb_get_zsp_nrows(a)
    type(psb_zspmat_type), intent(in) :: a
    psb_get_zsp_nrows = a%m

    return
  end function psb_get_zsp_nrows

  integer function psb_get_zsp_ncols(a)
    type(psb_zspmat_type), intent(in) :: a
    psb_get_zsp_ncols = a%k

    return
  end function psb_get_zsp_ncols
  

  integer function psb_get_dsp_nnzeros(a)
    type(psb_dspmat_type), intent(in) :: a  
    integer :: ires,info
    
    call psb_sp_info(psb_nztotreq_,a,ires,info)
    if (info == 0) then 
      psb_get_dsp_nnzeros = ires
    else
      psb_get_dsp_nnzeros = 0
    end if
  end function psb_get_dsp_nnzeros

  integer function psb_get_zsp_nnzeros(a)
    type(psb_zspmat_type), intent(in) :: a  
    integer :: ires,info
    
    call psb_sp_info(psb_nztotreq_,a,ires,info)
    if (info == 0) then 
      psb_get_zsp_nnzeros = ires
    else
      psb_get_zsp_nnzeros = 0
    end if
  end function psb_get_zsp_nnzeros

  integer function psb_get_dsp_nzsize(a)
    type(psb_dspmat_type), intent(in) :: a  
    integer :: ires,info
    
    call psb_sp_info(psb_nzsizereq_,a,ires,info)
    if (info == 0) then 
      psb_get_dsp_nzsize = ires
    else
      psb_get_dsp_nzsize = 0
    end if
  end function psb_get_dsp_nzsize

  integer function psb_get_zsp_nzsize(a)
    type(psb_zspmat_type), intent(in) :: a  
    integer :: ires,info
    
    call psb_sp_info(psb_nzsizereq_,a,ires,info)
    if (info == 0) then 
      psb_get_zsp_nzsize = ires
    else
      psb_get_zsp_nzsize = 0
    end if
  end function psb_get_zsp_nzsize


  integer function psb_get_dsp_nnz_row(ir,a)
    integer, intent(in)               :: ir
    type(psb_dspmat_type), intent(in) :: a  
    integer :: ires,info
    
    call psb_sp_info(psb_nzrowreq_,a,ires,info,iaux=ir)
    if (info == 0) then 
      psb_get_dsp_nnz_row = ires
    else
      psb_get_dsp_nnz_row = 0
    end if
  end function psb_get_dsp_nnz_row
  integer function psb_get_zsp_nnz_row(ir,a)
    integer, intent(in)               :: ir
    type(psb_zspmat_type), intent(in) :: a  
    integer :: ires,info
    
    call psb_sp_info(psb_nzrowreq_,a,ires,info,iaux=ir)
    if (info == 0) then 
      psb_get_zsp_nnz_row = ires
    else
      psb_get_zsp_nnz_row = 0
    end if
  end function psb_get_zsp_nnz_row


  subroutine psb_nullify_dsp(mat)
    implicit none
    type(psb_dspmat_type), intent(inout) :: mat

!!$    nullify(mat%aspk,mat%ia1,mat%ia2,mat%pl,mat%pr)

    mat%infoa(:)=0
    mat%m=0
    mat%k=0
    mat%fida=''
    mat%descra=''

  end subroutine psb_nullify_dsp

  Subroutine psb_dspreinit(a,info,clear)

    Implicit None

    !....Parameters...
    Type(psb_dspmat_type), intent(inout) :: a
    integer, intent(out)                 :: info
    logical, intent(in), optional        :: clear

    !locals
    logical, parameter  :: debug=.false.
    logical             :: clear_
    character(len=20)   :: name
    
    info = 0
    name = 'psb_sp_reinit'

    if (present(clear)) then 
      clear_ = clear
    else
      clear_ = .true.
    end if
    
    select case(psb_sp_getifld(psb_state_,a,info))
    case(psb_spmat_asb_) 

      if (clear_) a%aspk(:) = dzero

      if (psb_sp_getifld(psb_upd_,a,info)==psb_upd_perm_) then 
        if(a%fida(1:3).eq.'JAD') then
          a%ia1(a%infoa(psb_upd_pnt_)+psb_nnz_) = 0
        else
          a%ia2(a%infoa(psb_upd_pnt_)+psb_nnz_) = 0
        endif
      endif
      a%infoa(psb_state_) = psb_spmat_upd_
    case(psb_spmat_bld_) 
      ! in this case do nothing. this allows sprn to be called 
      ! right after allocate, with spins doing the right thing.
      ! hopefully :-)

    case( psb_spmat_upd_) 

    case default
      info=591     
      call psb_errpush(info,name)
    end select

  end Subroutine psb_dspreinit

  Subroutine psb_dspallocate(a, nnz,info)
    implicit none
    !....Parameters...
    Type(psb_dspmat_type), intent(inout) :: A
    Integer, intent(in)          :: nnz
    integer, intent(out)         :: info

    !locals
    logical, parameter  :: debug=.false.

    info  = 0
    if (nnz.lt.0) then
      info=45
      return
    Endif
    if (debug) write(0,*) 'SPALL : NNZ ',nnz,a%m,a%k
    call psb_sp_reall(a,nnz,info)

    a%pl(1)=0
    a%pr(1)=0
    ! set INFOA fields
    a%fida   = 'COO'
    a%descra = 'GUN'
    a%infoa(:) = 0
    a%m      = 0
    a%k      = 0
    if (debug) write(0,*) 'SPALL : end'
    Return

  End Subroutine psb_dspallocate

  Subroutine psb_dspallmk(m,k,a,info)
    implicit none
    !....Parameters...

    Type(psb_dspmat_type), intent(inout) :: A
    Integer, intent(in)          :: m,k
    Integer, intent(out)         :: info

    !locals
    logical, parameter  :: debug=.false.
    integer  :: nnz

    INFO  = 0
    nnz = 2*max(1,m,k)
    if (debug) write(0,*) 'SPALL : NNZ ',nnz,a%m,a%k
    a%m=max(0,m)
    a%k=max(0,k)
    call psb_sp_reall(a,nnz,info)

    a%pl(1)=0
    a%pr(1)=0
    ! set INFOA fields
    a%fida   = 'COO'
    a%descra = 'GUN'
    a%infoa(:) = 0
    if (debug) write(0,*) 'SPALL : end'
    Return

  end subroutine psb_dspallmk

  Subroutine psb_dspallmknz(m,k,a, nnz,info)
    implicit none
    !....parameters...

    type(psb_dspmat_type), intent(inout) :: a
    integer, intent(in)                  :: m,k,nnz
    integer, intent(out)                 :: info

    !locals
    logical, parameter  :: debug=.false.

    info  = 0
    if (nnz.lt.0) then
      info=45
      return
    endif
    if (debug) write(0,*) 'spall : nnz ',nnz,a%m,a%k
    a%m=max(0,m)
    a%k=max(0,k)
    call psb_sp_reall(a,nnz,info)
    if (debug) write(0,*) 'Check in ALLOCATE ',info,allocated(a%pl),allocated(a%pr)
    a%pl(1)=0
    a%pr(1)=0
    ! set infoa fields
    a%fida   = 'COO'
    a%descra = 'GUN'
    a%infoa(:) = 0
    if (debug) write(0,*) 'spall : end'
    return

  end subroutine psb_dspallmknz


  subroutine psb_dspall3(a, ni1,ni2,nd,info)
    implicit none
    !....Parameters...
    Type(psb_dspmat_type), intent(inout) :: A
    Integer, intent(in)          :: ni1,ni2,nd
    Integer, intent(out)         :: info

    !locals
    logical, parameter  :: debug=.false.

    info  = 0

    call psb_sp_reall(a, ni1,ni2,nd,info)

    a%pl(1)=0
    a%pr(1)=0
    ! set INFOA fields
    a%fida   = 'COO'
    a%descra = 'GUN'
    a%infoa(:) = 0
    a%m      = 0
    a%k      = 0
    if (debug) write(0,*) 'SPALL : end'
    Return

  End Subroutine psb_dspall3


  subroutine psb_dspreallocate(a, nnz,info,ifc)
    implicit none
    !....Parameters...
    Type(psb_dspmat_type), intent(inout)  :: A
    Integer, intent(in)           :: NNZ
    Integer, intent(out)          :: info
    !
    ! ifc is used here to allocate space in IA1 for smart 
    ! regeneration. This probably ought to be changed, 
    ! by adding a new component to d_spmat, or by making
    ! infoa a pointer.    
    !
    Integer, intent(in), optional :: ifc
    integer                       :: ifc_

    !locals
    logical, parameter  :: debug=.false.

    info  = 0
    if (nnz.lt.0) then
      info=45
      return
    endif
    if (present(ifc)) then 
      ifc_ = max(1,ifc)
    else
      ifc_ = 1
    endif

    if (ifc_ == 1) then 
      call psb_realloc(max(nnz,a%m+1,a%k+1),a%ia1,a%ia2,a%aspk,info)
    else
      call psb_realloc(max(nnz,a%m+1,a%k+1),a%aspk,info)
      if (info /= 0) return 
      call psb_realloc(max(nnz,a%m+1,a%k+1),a%ia2,info)
      if (info /= 0) return 
      call psb_realloc(ifc*nnz+200,a%ia1,info)
      if (info /= 0) return 
    end if
    if (info /= 0) return
    call psb_realloc(max(1,a%m),a%pl,info)
    if (info /= 0) return
    call psb_realloc(max(1,a%k),a%pr,info)
    if (debug) write(0,*) allocated(a%ia1),allocated(a%ia2),&
         & allocated(a%aspk),allocated(a%pl),allocated(a%pr),info
    if (info /= 0) return

    Return

  End Subroutine psb_dspreallocate

  subroutine psb_dspreall3(a, ni1,ni2,nd,info)
    implicit none
    !....Parameters...
    Type(psb_dspmat_type), intent(inout)  :: A
    Integer, intent(in)                   :: ni1,ni2,nd
    Integer, intent(inout)                :: info

    !locals
    logical, parameter  :: debug=.false.

    info  = 0
    if (debug) write(0,*) 'Before realloc',nd,size(a%aspk),ni1,ni2
    call psb_realloc(nd,a%aspk,info)
    if (debug) write(0,*) 'After realloc',nd,size(a%aspk),info
    if (info /= 0) return 
    if (debug) write(0,*) 'Before realloc2',ni2,allocated(a%ia2),size(a%ia2)
    call psb_realloc(ni2,a%ia2,info)
    if (info /= 0) return 
    if (debug) write(0,*) 'Before realloc3',ni1,allocated(a%ia1),size(a%ia1)
    call psb_realloc(ni1,a%ia1,info)
    if (info /= 0) return
    if (debug) write(0,*) 'Before realloc4',max(1,a%m),allocated(a%pl),size(a%pl)
    call psb_realloc(max(1,a%m),a%pl,info)
    if (info /= 0) return
    if (debug) write(0,*) 'Before realloc5',max(1,a%k),allocated(a%pr),size(a%pr)
    call psb_realloc(max(1,a%k),a%pr,info)
    if (info /= 0) return

    Return

  End Subroutine psb_dspreall3


  subroutine psb_dspclone(a, b,info)
    implicit none
    !....Parameters...
    Type(psb_dspmat_type), intent(in)   :: A
    Type(psb_dspmat_type), intent(out)  :: B
    Integer, intent(out)                :: info

    !locals

    INFO  = 0
    call psb_nullify_sp(b)
    call psb_safe_cpy(a%aspk,b%aspk,info)
    if (info == 0) call psb_safe_cpy(a%ia1,b%ia1,info)
    if (info == 0) call psb_safe_cpy(a%ia2,b%ia2,info)
    if (info == 0) call psb_safe_cpy(a%pl,b%pl,info)
    if (info == 0) call psb_safe_cpy(a%pr,b%pr,info)
    if (info /= 0) then
      info=2023
      return
    Endif
    b%infoa(:) = a%infoa(:)
    b%fida     = a%fida
    b%descra   = a%descra
    b%m        = a%m
    b%k        = a%k

    Return

  End Subroutine psb_dspclone


  
  ! Will be changed to use MOVE_ALLOC 
  subroutine psb_dsp_transfer(a, b,info)
    implicit none
    !....Parameters...
    Type(psb_dspmat_type), intent(inout)  :: A
    Type(psb_dspmat_type), intent(inout)  :: B
    Integer, intent(out)                  :: info

    !locals
    logical, parameter  :: debug=.false.

    info  = 0


    call psb_transfer( a%aspk,     b%aspk  , info)
    call psb_transfer( a%ia1 ,     b%ia1   , info)
    call psb_transfer( a%ia2 ,     b%ia2   , info)
    call psb_transfer( a%pl  ,     b%pl    , info)
    call psb_transfer( a%pr  ,     b%pr    , info)
    b%infoa(:) = a%infoa(:)
    b%fida     = a%fida
    b%descra   = a%descra
    b%m        = a%m
    b%k        = a%k

    call psb_nullify_sp(a)

    Return

  End Subroutine psb_dsp_transfer


  Subroutine psb_dsp_setifld(val,field,a,info)
    implicit none
    !....Parameters...

    Type(psb_dspmat_type), intent(inout) :: A
    Integer, intent(in)          :: field,val
    Integer, intent(out)         :: info
    
    !locals
    logical, parameter  :: debug=.false.
    info  = 0

!!$    call psb_realloc(psb_ifasize_,a%infoa,info)
    
    if (info == 0) &
         & call psb_setifield(val,field,a%infoa,psb_ifasize_,info)
    

    Return

  end subroutine psb_dsp_setifld


  subroutine psb_dsp_trimsize(a, i1,i2,ia,info)
    use psb_string_mod
    implicit none
    !....Parameters...
    Type(psb_dspmat_type), intent(in) :: A
    Integer, intent(out)              :: i1, i2, ia, info

    !locals
    Integer             :: nza
    logical, parameter  :: debug=.false.

    info  = 0
    if (psb_sp_getifld(psb_upd_,a,info) == psb_upd_perm_) then 
      info = -1 
      i1 = size(a%ia1)
      i2 = size(a%ia2)
      ia = size(a%aspk)
      return
    endif
    select case(toupper(a%fida))
    case('CSR')
      nza = a%ia2(a%m+1)-1
      ia  = nza
      i1  = nza
      i2  = a%m + 1
    case('COO')
      nza = a%infoa(psb_nnz_)
      i1  = nza
      i2  = nza
      ia  = nza
    case('JAD')
      ! Feeling lazy today
      i1 = size(a%ia1)
      i2 = size(a%ia2)
      ia = size(a%aspk)
    case default
      info = -2
      i1 = size(a%ia1)
      i2 = size(a%ia2)
      ia = size(a%aspk)
    end select

    Return

  End Subroutine psb_dsp_trimsize
  
  function psb_dsp_getifld(field,a,info)
    implicit none
    !....Parameters...

    Type(psb_dspmat_type), intent(in) :: A
    Integer, intent(in)          :: field
    Integer                      :: psb_dsp_getifld
    Integer, intent(out)         :: info

    !locals
    logical, parameter  :: debug=.false.
    integer :: val

    info  = 0
    val   = -1

    if ((field < 1).or.(field > psb_ifasize_)) then
      info = -1
      psb_dsp_getifld = val
      return
    endif

    call psb_getifield(val,field,a%infoa,psb_ifasize_,info)
    
    psb_dsp_getifld = val
    Return

  end function psb_dsp_getifld

  function psb_dspsizeof(a)
    implicit none
    !....Parameters...

    Type(psb_dspmat_type), intent(in) :: A
    Integer                      :: psb_dspsizeof

    !locals
    logical, parameter  :: debug=.false.
    integer :: val

    val   = 4*size(a%infoa)
    
    if (allocated(a%aspk)) then 
      val = val + 8 * size(a%aspk)
    endif

    if (allocated(a%ia1)) then 
      val = val + 4 * size(a%ia1)
    endif
    if (allocated(a%ia2)) then 
      val = val + 4 * size(a%ia2)
    endif
    if (allocated(a%pl)) then 
      val = val + 4 * size(a%pl)
    endif
    if (allocated(a%pr)) then 
      val = val + 4 * size(a%pr)
    endif

    
    psb_dspsizeof = val
    Return

  end function psb_dspsizeof


  subroutine psb_dsp_free(a,info)
    implicit none
    !....Parameters...
    Type(psb_dspmat_type), intent(inout)  :: A
    Integer, intent(out)        :: info
    !locals
    logical, parameter  :: debug=.false.
    integer             :: iret
    info  = 0

    if (allocated(a%aspk)) then
!!$      write(0,*) 'Deallocating aspk'
      deallocate(a%aspk,STAT=IRET)
!!$      write(0,*) 'Deallocated  aspk',iret
      if (iret /= 0) info = max(info,1)
    endif
    if (allocated(a%ia1)) then
      deallocate(a%ia1,STAT=IRET)
      if (iret /= 0) info = max(info,2)
    endif
    if (allocated(a%ia2)) then
      deallocate(a%ia2,STAT=IRET)
      if (iret /= 0) info = max(info,3)
    endif
    if (allocated(a%pr)) then
      deallocate(a%pr,STAT=IRET)
      if (iret /= 0) info = max(info,4)
    endif
    if (allocated(a%pl)) then
      deallocate(a%pl,STAT=IRET)
      if (iret /= 0) info = max(info,5)
    endif
    call psb_nullify_sp(a)
!!$    write(0,*) 'End of sp_free ',info
    Return
  End Subroutine psb_dsp_free


  subroutine psb_nullify_zsp(mat)
    implicit none
    type(psb_zspmat_type), intent(inout) :: mat

    mat%infoa(:)=0
    mat%m=0
    mat%k=0
    mat%fida=''
    mat%descra=''

  end subroutine psb_nullify_zsp

  Subroutine psb_zspreinit(a,info,clear)

    Implicit None

    !....Parameters...
    Type(psb_zspmat_type), intent(inout) :: a
    integer, intent(out)                 :: info
    logical, intent(in), optional        :: clear

    !locals
    logical, parameter  :: debug=.false.
    logical             :: clear_
    character(len=20)   :: name
    
    info = 0
    name = 'psb_sp_reinit'

    if (present(clear)) then 
      clear_ = clear
    else
      clear_ = .true.
    end if
    
    select case(psb_sp_getifld(psb_state_,a,info))
    case(psb_spmat_asb_) 

      if (clear_) a%aspk(:) = zzero

      if (psb_sp_getifld(psb_upd_,a,info)==psb_upd_perm_) then 
        if(a%fida(1:3).eq.'JAD') then
          a%ia1(a%infoa(psb_upd_pnt_)+psb_nnz_) = 0
        else
          a%ia2(a%infoa(psb_upd_pnt_)+psb_nnz_) = 0
        endif
      endif
      a%infoa(psb_state_) = psb_spmat_upd_
    case(psb_spmat_bld_) 
      ! in this case do nothing. this allows sprn to be called 
      ! right after allocate, with spins doing the right thing.
      ! hopefully :-)

    case( psb_spmat_upd_) 

    case default
      info=591     
      call psb_errpush(info,name)
    end select

  end Subroutine psb_zspreinit

  Subroutine psb_zspallocate(a, nnz,info)
    implicit none
    !....Parameters...
    Type(psb_zspmat_type), intent(inout) :: A
    Integer, intent(in)          :: nnz
    integer, intent(out)         :: info

    !locals
    logical, parameter  :: debug=.false.

    info  = 0
    if (nnz.lt.0) then
      info=45
      return
    Endif
    if (debug) write(0,*) 'SPALL : NNZ ',nnz,a%m,a%k
    call psb_sp_reall(a,nnz,info)

    a%pl(1)=0
    a%pr(1)=0
    ! set INFOA fields
    a%fida   = 'COO'
    a%descra = 'GUN'
    a%infoa(:) = 0
    a%m      = 0
    a%k      = 0
    if (debug) write(0,*) 'SPALL : end'
    Return

  End Subroutine psb_zspallocate

  Subroutine psb_zspallmk(m,k,a,info)
    implicit none
    !....Parameters...

    Type(psb_zspmat_type), intent(inout) :: A
    Integer, intent(in)          :: m,k
    Integer, intent(out)         :: info

    !locals
    logical, parameter  :: debug=.false.
    integer  :: nnz

    INFO  = 0
    nnz = 2*max(1,m,k)
    if (debug) write(0,*) 'SPALL : NNZ ',nnz,a%m,a%k
    a%m=max(0,m)
    a%k=max(0,k)
    call psb_sp_reall(a,nnz,info)

    a%pl(1)=0
    a%pr(1)=0
    ! set INFOA fields
    a%fida   = 'COO'
    a%descra = 'GUN'
    a%infoa(:) = 0
    if (debug) write(0,*) 'SPALL : end'
    Return

  end subroutine psb_zspallmk

  Subroutine psb_zspallmknz(m,k,a, nnz,info)
    implicit none
    !....parameters...

    type(psb_zspmat_type), intent(inout) :: a
    integer, intent(in)                  :: m,k,nnz
    integer, intent(out)                 :: info

    !locals
    logical, parameter  :: debug=.false.

    info  = 0
    if (nnz.lt.0) then
      info=45
      return
    endif
    if (debug) write(0,*) 'spall : nnz ',nnz,a%m,a%k
    a%m=max(0,m)
    a%k=max(0,k)
    call psb_sp_reall(a,nnz,info)

    a%pl(1)=0
    a%pr(1)=0
    ! set infoa fields
    a%fida   = 'COO'
    a%descra = 'GUN'
    a%infoa(:) = 0
    if (debug) write(0,*) 'spall : end'
    return

  end subroutine psb_zspallmknz


  subroutine psb_zspall3(a, ni1,ni2,nd,info)
    implicit none
    !....Parameters...
    Type(psb_zspmat_type), intent(inout) :: A
    Integer, intent(in)          :: ni1,ni2,nd
    Integer, intent(out)         :: info

    !locals
    logical, parameter  :: debug=.false.

    info  = 0

    call psb_sp_reall(a, ni1,ni2,nd,info)

    a%pl(1)=0
    a%pr(1)=0
    ! set INFOA fields
    a%fida   = 'COO'
    a%descra = 'GUN'
    a%infoa(:) = 0
    a%m      = 0
    a%k      = 0
    if (debug) write(0,*) 'SPALL : end'
    Return

  End Subroutine psb_zspall3

  subroutine psb_zspreall3(a, ni1,ni2,nz,info)
    implicit none
    !....Parameters...
    Type(psb_zspmat_type), intent(inout)  :: A
    Integer, intent(in)                   :: ni1,ni2,nz
    Integer, intent(inout)                :: info

    !locals
    logical, parameter  :: debug=.false.

    info  = 0
    call psb_realloc(nz,a%aspk,info)
    if (info /= 0) return 
    call psb_realloc(ni2,a%ia2,info)
    if (info /= 0) return 
    call psb_realloc(ni1,a%ia1,info)
    if (info /= 0) return
    call psb_realloc(max(1,a%m),a%pl,info)
    if (info /= 0) return
    call psb_realloc(max(1,a%k),a%pr,info)
    if (info /= 0) return

    Return

  End Subroutine psb_zspreall3


  subroutine psb_zspreallocate(a, nnz,info,ifc)
    implicit none
    !....Parameters...
    Type(psb_zspmat_type), intent(inout)  :: A
    Integer, intent(in)           :: NNZ
    Integer, intent(out)          :: info
    !
    ! ifc is used here to allocate space in IA1 for smart 
    ! regeneration. This probably ought to be changed, 
    ! by adding a new component to d_spmat, or by making
    ! infoa a pointer.    
    !
    Integer, intent(in), optional :: ifc
    integer                       :: ifc_

    !locals
    logical, parameter  :: debug=.false.

    info  = 0
    if (nnz.lt.0) then
      info=45
      return
    endif
    if (present(ifc)) then 
      ifc_ = max(1,ifc)
    else
      ifc_ = 1
    endif

    if (ifc_ == 1) then 
      call psb_realloc(max(nnz,a%m+1,a%k+1),a%ia1,a%ia2,a%aspk,info)
    else
      call psb_realloc(max(nnz,a%m+1,a%k+1),a%aspk,info)
      if (info /= 0) return 
      call psb_realloc(max(nnz,a%m+1,a%k+1),a%ia2,info)
      if (info /= 0) return 
      call psb_realloc(ifc*nnz+200,a%ia1,info)
      if (info /= 0) return 
    end if
    if (info /= 0) return
    call psb_realloc(max(1,a%m),a%pl,info)
    if (info /= 0) return
    call psb_realloc(max(1,a%k),a%pr,info)
    if (info /= 0) return

    Return

  End Subroutine psb_zspreallocate

  subroutine psb_zspclone(a, b,info)
    implicit none
    !....Parameters...
    Type(psb_zspmat_type), intent(in)   :: A
    Type(psb_zspmat_type), intent(out)  :: B
    Integer, intent(out)                :: info

    !locals

    logical, parameter  :: debug=.false.


    INFO  = 0
    call psb_nullify_sp(b)
    call psb_safe_cpy(a%aspk,b%aspk,info)
    if (info == 0) call psb_safe_cpy(a%ia1,b%ia1,info)
    if (info == 0) call psb_safe_cpy(a%ia2,b%ia2,info)
    if (info == 0) call psb_safe_cpy(a%pl,b%pl,info)
    if (info == 0) call psb_safe_cpy(a%pr,b%pr,info)
    if (info /= 0) then
      info=2023
      return
    Endif
    b%infoa(:) = a%infoa(:)
    b%fida     = a%fida
    b%descra   = a%descra
    b%m        = a%m
    b%k        = a%k

    Return

  End Subroutine psb_zspclone


  ! This is done with pointer assignments, but it 
  ! will be feasible with MOVE_ALLOC when we move 
  ! to ALLOCATABLE components. 
  subroutine psb_zsp_transfer(a, b,info)
    implicit none
    !....Parameters...
    Type(psb_zspmat_type), intent(inout)  :: A
    Type(psb_zspmat_type), intent(inout)  :: B
    Integer, intent(out)                  :: info

    !locals
    logical, parameter  :: debug=.false.

    info  = 0

    call psb_transfer( a%aspk,     b%aspk  , info)
    call psb_transfer( a%ia1 ,     b%ia1   , info)
    call psb_transfer( a%ia2 ,     b%ia2   , info)
    call psb_transfer( a%pl  ,     b%pl    , info)
    call psb_transfer( a%pr  ,     b%pr    , info)
    b%infoa(:) = a%infoa(:)
    b%fida     = a%fida
    b%descra   = a%descra
    b%m        = a%m
    b%k        = a%k

    call psb_nullify_sp(a)

    Return

  End Subroutine psb_zsp_transfer

  Subroutine psb_zsp_setifld(val,field,a,info)
    implicit none
    !....Parameters...

    Type(psb_zspmat_type), intent(inout) :: A
    Integer, intent(in)          :: field,val
    Integer, intent(out)         :: info

    !locals
    logical, parameter  :: debug=.false.

    info  = 0


!!$    call psb_realloc(psb_ifasize_,a%infoa,info)
    
    if (info == 0) &
         & call psb_setifield(val,field,a%infoa,psb_ifasize_,info)
    

    Return

  end subroutine psb_zsp_setifld


  subroutine psb_zsp_trimsize(a, i1,i2,ia,info)
    use psb_string_mod
    implicit none
    !....Parameters...
    Type(psb_zspmat_type), intent(in) :: A
    Integer, intent(out)              :: i1, i2, ia, info

    !locals
    Integer             :: nza
    logical, parameter  :: debug=.false.

    info  = 0
    if (psb_sp_getifld(psb_upd_,a,info) == psb_upd_perm_) then 
      info = -1 
      i1 = size(a%ia1)
      i2 = size(a%ia2)
      ia = size(a%aspk)
      return
    endif
    select case(toupper(a%fida))
    case('CSR')
      nza = a%ia2(a%m+1)-1
      ia  = nza
      i1  = nza
      i2  = a%m + 1
    case('COO')
      nza = a%infoa(psb_nnz_)
      i1  = nza
      i2  = nza
      ia  = nza
    case('JAD')
      ! Feeling lazy today
      i1 = size(a%ia1)
      i2 = size(a%ia2)
      ia = size(a%aspk)
    case default
      info = -2
      i1 = size(a%ia1)
      i2 = size(a%ia2)
      ia = size(a%aspk)
    end select

    Return

  End Subroutine psb_zsp_trimsize

  function psb_zsp_getifld(field,a,info)
    implicit none
    !....Parameters...

    Type(psb_zspmat_type), intent(in) :: A
    Integer, intent(in)          :: field
    Integer                      :: psb_zsp_getifld
    Integer, intent(out)         :: info

    !locals
    logical, parameter  :: debug=.false.
    integer :: val

    info  = 0
    val   = -1

    if ((field < 1).or.(field > psb_ifasize_)) then
      info = -1
      psb_zsp_getifld = val
      return
    endif

    call psb_getifield(val,field,a%infoa,psb_ifasize_,info)
    
    psb_zsp_getifld = val
    Return

  end function psb_zsp_getifld

  function psb_zspsizeof(a)
    implicit none
    !....Parameters...

    Type(psb_zspmat_type), intent(in) :: A
    Integer                      :: psb_zspsizeof

    !locals
    logical, parameter  :: debug=.false.
    integer :: val

    val   = 4*size(a%infoa)
    
    if (allocated(a%aspk)) then 
      val = val + 16 * size(a%aspk)
    endif

    if (allocated(a%ia1)) then 
      val = val + 4 * size(a%ia1)
    endif
    if (allocated(a%ia2)) then 
      val = val + 4 * size(a%ia2)
    endif
    if (allocated(a%pl)) then 
      val = val + 4 * size(a%pl)
    endif
    if (allocated(a%pr)) then 
      val = val + 4 * size(a%pr)
    endif
    
    
    psb_zspsizeof = val
    Return

  end function psb_zspsizeof




  subroutine psb_zsp_free(a,info)
    implicit none
    !....Parameters...
    Type(psb_zspmat_type), intent(inout)  :: A
    Integer, intent(out)        :: info
    !locals
    logical, parameter  :: debug=.false.

    info  = 0

    if (allocated(a%aspk)) then
      deallocate(a%aspk,STAT=INFO)
    endif
    if (allocated(a%ia1)) then
      deallocate(a%ia1,STAT=INFO)
    endif
    if ( allocated(a%ia2)) then
      deallocate(a%ia2,STAT=INFO)
    endif
    if ( allocated(a%pr)) then
      deallocate(a%pr,STAT=INFO)
    endif
    if ( allocated(a%pl)) then
      deallocate(a%pl,STAT=INFO)
    endif
    call psb_nullify_sp(a)
    Return
  End Subroutine psb_zsp_free

  subroutine psb_dspinfo(ireq,a,ires,info,iaux)
    use psb_const_mod
    use psb_error_mod
    use psb_string_mod
    implicit none

    type(psb_dspmat_type), intent(in), target :: a
    integer, intent(in)               :: ireq
    integer, intent(out)              :: ires, info
    integer, intent(in), optional     :: iaux

    integer :: i,j,k,ip,jp,nr,irw,nz, err_act, row, ipx, pia, pja, rb,idx, nc
    integer, pointer :: ia1(:), ia2(:), ia3(:), ja(:)
    character(len=20)                 :: name, ch_err

    name='psb_dspinfo'
    info  = 0
    call psb_erractionsave(err_act)


    if (ireq == psb_nztotreq_) then 
      ! The number of nonzeroes
      if (toupper(a%fida) == 'CSR') then 
        nr   = a%m
        ires = a%ia2(nr+1)-1
      else if ((toupper(a%fida) == 'COO').or.(toupper(a%fida) == 'COI')) then 
        ires = a%infoa(psb_nnz_)
      else if (toupper(a%fida) == 'JAD') then 
        ires = a%infoa(psb_nnz_)
      else if (toupper(a%fida) == 'CSC') then 
        nc   = a%k
        ires = a%ia2(nc+1)-1
      else
        ires=-1
        info=136
        ch_err=a%fida(1:3)
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

    else if (ireq == psb_nzrowreq_) then 
      ! The number of nonzeroes in row iaux
      if (.not.present(iaux)) then 
        write(0,*) 'Need IAUX when ireq=nzrowreq'
        ires=-1
        return
      endif
      irw = iaux
      if (irw > a%m) then 
        write(0,*) 'SPINFO: Accessing out of bounds? ',irw,a%m
        ires = 0
        return
      endif
      if (toupper(a%fida) == 'CSR') then 
        ires = a%ia2(irw+1)-a%ia2(irw)
      else if ((toupper(a%fida) == 'COO').or.(toupper(a%fida) == 'COI')) then 

        if (a%infoa(psb_srtd_) == psb_isrtdcoo_) then 
          ! In this case we can do a binary search. 
          nz = a%infoa(psb_nnz_)
          call ibsrch(ip,irw,nz,a%ia1)
          jp = ip
          ! expand [ip,jp] to contain all row entries.
          do 
            if (ip < 2) exit
            if (a%ia1(ip-1) == irw) then  
              ip = ip -1 
            else 
              exit
            end if
          end do

          do
            if (jp > nz) exit
            if (a%ia1(jp) == irw) then
              jp =jp + 1
            else
              exit
            endif
          end do
          ires = jp-ip
        else
          ires = count(a%ia1(1:a%infoa(psb_nnz_))==irw)
        endif
!!$      ires = 0
!!$      do i=1, a%infoa(psb_nnz_) 
!!$        if (a%ia1(i) == irw) ires = ires + 1
!!$      enddo
      else if (toupper(a%fida) == 'JAD') then 
        pia = a%ia2(2) ! points to the beginning of ia(3,png)
        pja = a%ia2(3) ! points to the beginning of ja(:)
        ja  => a%ia2(pja:)             ! the array containing the pointers to ka and aspk
        ia1 => a%ia2(pia:pja-1:3)      ! the array containing the first row index of each block
        ia2 => a%ia2(pia+1:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first jad column
        ia3 => a%ia2(pia+2:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first csr column

        idx=a%pl(irw)
        j=0
        nz=0
        blkfnd: do
          j=j+1
          if(ia1(j).eq.idx) then
            nz=nz+ia3(j)-ia2(j)
            ipx = ia1(j)         ! the first row index of the block
            rb  = idx-ipx        ! the row offset within the block
            row = ia3(j)+rb
            nz  = nz+ja(row+1)-ja(row)
            exit blkfnd
          else if(ia1(j).gt.idx) then
            nz=nz+ia3(j-1)-ia2(j-1)
            ipx = ia1(j-1)         ! the first row index of the block
            rb  = idx-ipx          ! the row offset within the block
            row = ia3(j-1)+rb
            nz  = nz+ja(row+1)-ja(row)
            exit blkfnd
          end if
        end do blkfnd
        ires=nz
      else
        ires=-1
        info=136
        ch_err=a%fida(1:3)
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

    else  if (ireq == psb_nzsizereq_) then 
      if (toupper(a%fida) == 'CSR') then 
        ires = size(a%aspk)
      else if ((toupper(a%fida) == 'COO').or.(toupper(a%fida) == 'COI')) then 
        ires = size(a%aspk)
      else if (toupper(a%fida) == 'JAD') then 
        ires = a%infoa(psb_nnz_)
      else
        ires=-1
        info=136
        ch_err=a%fida(1:3)
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

    else 
      write(0,*) 'Unknown request into SPINFO'
      ires=-1
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine psb_dspinfo


  subroutine psb_zspinfo(ireq,a,ires,info,iaux)
    use psb_const_mod
    use psb_error_mod
    use psb_string_mod
    implicit none

    type(psb_zspmat_type), intent(in), target :: a
    integer, intent(in)               :: ireq
    integer, intent(out)              :: ires, info
    integer, intent(in), optional     :: iaux

    integer :: i,j,k,ip,jp,nr,irw,nz, err_act, row, ipx, pia, pja, rb,idx, nc
    integer, pointer :: ia1(:), ia2(:), ia3(:), ja(:)
    character(len=20)                 :: name, ch_err

    name='psb_zspinfo'
    info  = 0
    call psb_erractionsave(err_act)


    if (ireq == psb_nztotreq_) then 
      ! The number of nonzeroes
      if (toupper(a%fida) == 'CSR') then 
        nr   = a%m
        ires = a%ia2(nr+1)-1
      else if ((toupper(a%fida) == 'COO').or.(toupper(a%fida) == 'COI')) then 
        ires = a%infoa(psb_nnz_)
      else if (toupper(a%fida) == 'JAD') then 
        ires = a%infoa(psb_nnz_)
      else if (toupper(a%fida) == 'CSC') then 
        nc   = a%k
        ires = a%ia2(nc+1)-1
      else
        ires=-1
        info=136
        ch_err=a%fida(1:3)
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

    else if (ireq == psb_nzrowreq_) then 
      ! The number of nonzeroes in row iaux
      if (.not.present(iaux)) then 
        write(0,*) 'Need IAUX when ireq=nzrowreq'
        ires=-1
        return
      endif
      irw = iaux
      if (toupper(a%fida) == 'CSR') then 
        ires = a%ia2(irw+1)-a%ia2(irw)
      else if ((toupper(a%fida) == 'COO').or.(toupper(a%fida) == 'COI')) then 

        if (a%infoa(psb_srtd_) == psb_isrtdcoo_) then 
          ! In this case we can do a binary search. 
          nz = a%infoa(psb_nnz_)
          call ibsrch(ip,irw,nz,a%ia1)
          jp = ip
          ! expand [ip,jp] to contain all row entries.
          do 
            if (ip < 2) exit
            if (a%ia1(ip-1) == irw) then  
              ip = ip -1 
            else 
              exit
            end if
          end do

          do
            if (jp > nz) exit
            if (a%ia1(jp) == irw) then
              jp =jp + 1
            else
              exit
            endif
          end do
          ires = jp-ip
        else
          ires = count(a%ia1(1:a%infoa(psb_nnz_))==irw)
        endif
!!$      ires = 0
!!$      do i=1, a%infoa(psb_nnz_) 
!!$        if (a%ia1(i) == irw) ires = ires + 1
!!$      enddo
      else if (toupper(a%fida) == 'JAD') then 
        pia = a%ia2(2) ! points to the beginning of ia(3,png)
        pja = a%ia2(3) ! points to the beginning of ja(:)
        ja  => a%ia2(pja:)             ! the array containing the pointers to ka and aspk
        ia1 => a%ia2(pia:pja-1:3)      ! the array containing the first row index of each block
        ia2 => a%ia2(pia+1:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first jad column
        ia3 => a%ia2(pia+2:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first csr column

        idx=a%pl(irw)
        j=0
        nz=0
        blkfnd: do
          j=j+1
          if(ia1(j).eq.idx) then
            nz=nz+ia3(j)-ia2(j)
            ipx = ia1(j)         ! the first row index of the block
            rb  = idx-ipx        ! the row offset within the block
            row = ia3(j)+rb
            nz  = nz+ja(row+1)-ja(row)
            exit blkfnd
          else if(ia1(j).gt.idx) then
            nz=nz+ia3(j-1)-ia2(j-1)
            ipx = ia1(j-1)         ! the first row index of the block
            rb  = idx-ipx          ! the row offset within the block
            row = ia3(j-1)+rb
            nz  = nz+ja(row+1)-ja(row)
            exit blkfnd
          end if
        end do blkfnd
        ires=nz
      else
        ires=-1
        info=136
        ch_err=a%fida(1:3)
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

    else  if (ireq == psb_nzsizereq_) then 
      if (toupper(a%fida) == 'CSR') then 
        ires = size(a%aspk)
      else if ((toupper(a%fida) == 'COO').or.(toupper(a%fida) == 'COI')) then 
        ires = size(a%aspk)
      else if (toupper(a%fida) == 'JAD') then 
        ires = a%infoa(psb_nnz_)
      else
        ires=-1
        info=136
        ch_err=a%fida(1:3)
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

    else 
      write(0,*) 'Unknown request into SPINFO'
      ires=-1
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine psb_zspinfo

end module psb_spmat_type

