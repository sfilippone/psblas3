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
  ! Typedef: psb_dspmat_type
  !    Contains a sparse matrix
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
    real(kind(1.d0)), pointer :: aspk(:)=>null()
    ! Contains indeces that describes sparse matrix structure
    integer, pointer :: ia1(:)=>null(), ia2(:)=>null()
    ! Permutations matrix
    integer, pointer :: pl(:)=>null(), pr(:)=>null()
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
    complex(kind(1.d0)), pointer :: aspk(:)=>null()
    ! Contains indeces that describes sparse matrix structure
    integer, pointer :: ia1(:)=>null(), ia2(:)=>null()
    ! Permutations matrix
    integer, pointer :: pl(:)=>null(), pr(:)=>null()
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

contains 

  subroutine psb_nullify_dsp(mat)
    implicit none
    type(psb_dspmat_type), intent(inout) :: mat

    nullify(mat%aspk,mat%ia1,mat%ia2,mat%pl,mat%pr)
    mat%m=0
    mat%k=0
    mat%fida=''
    mat%descra=''

  end subroutine psb_nullify_dsp

  Subroutine psb_dspreinit(a)
    implicit none
    !....Parameters...
    Type(psb_dspmat_type), intent(inout) :: A

    !locals
    logical, parameter  :: debug=.false.

    if (debug) write(0,*) 'spreinit init ',a%fida,a%infoa(psb_nnz_)
    if (a%fida=='COO') a%infoa(psb_nnz_) = 0
    if (associated(a%aspk)) a%aspk(:) = 0.d0
    if (debug) write(0,*) 'spreinit end ',a%fida,a%infoa(psb_nnz_)

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
    if (debug) write(0,*) 'Check in ALLOCATE ',info,associated(a%pl),associated(a%pr)
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
      call psb_realloc(nnz,a%ia1,a%ia2,a%aspk,info)
    else
      call psb_realloc(nnz,a%aspk,info)
      if (info /= 0) return 
      call psb_realloc(nnz,a%ia2,info)
      if (info /= 0) return 
      call psb_realloc(ifc*nnz+200,a%ia1,info)
      if (info /= 0) return 
    end if
    if (info /= 0) return
    call psb_realloc(max(1,a%m),a%pl,info)
    if (info /= 0) return
    call psb_realloc(max(1,a%k),a%pr,info)
    if (debug) write(0,*) associated(a%ia1),associated(a%ia2),&
         & associated(a%aspk),associated(a%pl),associated(a%pr),info
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
    if (debug) write(0,*) 'Before realloc',nd,size(a%aspk)
    call psb_realloc(nd,a%aspk,info)
    if (debug) write(0,*) 'After realloc',nd,size(a%aspk),info
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

  End Subroutine psb_dspreall3


  subroutine psb_dspclone(a, b,info)
    implicit none
    !....Parameters...
    Type(psb_dspmat_type), intent(in)   :: A
    Type(psb_dspmat_type), intent(out)  :: B
    Integer, intent(out)                :: info

    !locals
    Integer             :: nza,nz1, nz2, nzl, nzr
    logical, parameter  :: debug=.false.

    INFO  = 0

    nza = size(a%aspk)
    nz1 = size(a%ia1)
    nz2 = size(a%ia2)
    nzl = size(a%pl)
    nzr = size(a%pr)
    allocate(b%aspk(nza),b%ia1(nz1),b%ia2(nz2),&
         & b%pl(nzl),b%pr(nzr),stat=info)
    if (info /= 0) then
      info=2023
      return
    Endif
    b%aspk(:)  = a%aspk(:)
    b%ia1(:)   = a%ia1(:)
    b%ia2(:)   = a%ia2(:)
    b%pl(:)    = a%pl(:)
    b%pr(:)    = a%pr(:)    
    b%infoa(:) = a%infoa(:)
    b%fida     = a%fida
    b%descra   = a%descra
    b%m        = a%m
    b%k        = a%k

    Return

  End Subroutine psb_dspclone


  ! This is done with pointer assignments, but it 
  ! will be feasible with MOVE_ALLOC when we move 
  ! to ALLOCATABLE components. 
  subroutine psb_dsp_transfer(a, b,info)
    implicit none
    !....Parameters...
    Type(psb_dspmat_type), intent(inout)  :: A
    Type(psb_dspmat_type), intent(inout)  :: B
    Integer, intent(out)                  :: info

    !locals
    Integer             :: nza,nz1, nz2, nzl, nzr
    logical, parameter  :: debug=.false.

    INFO  = 0

    if (associated(b%pr)) then  
      deallocate(b%pr,stat=info)
    end if

    if (associated(b%pl)) then 
      deallocate(b%pl,stat=info)
    end if
    if (associated(b%ia2)) then
      deallocate(b%ia2,stat=info)
    end if
    if (associated(b%ia1))  then
      deallocate(b%ia1,stat=info)
    endif
    if (associated(b%aspk))  then
      deallocate(b%aspk,stat=info)
    endif

    b%aspk  => a%aspk
    b%ia1   => a%ia1
    b%ia2   => a%ia2
    b%pl    => a%pl
    b%pr    => a%pr    
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
    Integer             :: nza,nz1, nz2, nzl, nzr
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

!!$    if (.not.associated(a%infoa)) then 
!!$      info = -2 
!!$      return
!!$    endif
    
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
    
    if (associated(a%aspk)) then 
      val = val + 8 * size(a%aspk)
    endif

    if (associated(a%ia1)) then 
      val = val + 4 * size(a%ia1)
    endif
    if (associated(a%ia2)) then 
      val = val + 4 * size(a%ia2)
    endif
    if (associated(a%pl)) then 
      val = val + 4 * size(a%pl)
    endif
    if (associated(a%pr)) then 
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

    info  = 0

    if (associated(a%aspk)) then
      deallocate(a%aspk,STAT=INFO)
    endif
    if ((info == 0) .and. associated(a%ia1)) then
      deallocate(a%ia1,STAT=INFO)
    endif
    if ((info == 0) .and. associated(a%ia2)) then
      deallocate(a%ia2,STAT=INFO)
    endif
    if ((info == 0) .and. associated(a%pr)) then
      deallocate(a%pr,STAT=INFO)
    endif
    if ((info == 0) .and. associated(a%pl)) then
      deallocate(a%pl,STAT=INFO)
    endif
    call psb_nullify_sp(a)
    Return
  End Subroutine psb_dsp_free


  subroutine psb_nullify_zsp(mat)
    implicit none
    type(psb_zspmat_type), intent(inout) :: mat

    nullify(mat%aspk,mat%ia1,mat%ia2,mat%pl,mat%pr)
    mat%m=0
    mat%k=0
    mat%fida=''
    mat%descra=''

  end subroutine psb_nullify_zsp

  Subroutine psb_zspreinit(a)
    implicit none
    !....Parameters...
    Type(psb_zspmat_type), intent(inout) :: A

    !locals
    logical, parameter  :: debug=.false.

    if (debug) write(0,*) 'spreinit init ',a%fida,a%infoa(psb_nnz_)
    if (a%fida=='COO') a%infoa(psb_nnz_) = 0
    if (associated(a%aspk)) a%aspk(:) = 0.d0
    if (debug) write(0,*) 'spreinit end ',a%fida,a%infoa(psb_nnz_)

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
      call psb_realloc(nnz,a%ia1,a%ia2,a%aspk,info)
    else
      call psb_realloc(nnz,a%aspk,info)
      if (info /= 0) return 
      call psb_realloc(nnz,a%ia2,info)
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
    Integer             :: nza,nz1, nz2, nzl, nzr
    logical, parameter  :: debug=.false.

    INFO  = 0

    nza = size(a%aspk)
    nz1 = size(a%ia1)
    nz2 = size(a%ia2)
    nzl = size(a%pl)
    nzr = size(a%pr)
    allocate(b%aspk(nza),b%ia1(nz1),b%ia2(nz2),&
         & b%pl(nzl),b%pr(nzr),stat=info)
    if (info /= 0) then
      info=2023
      return
    Endif
    b%aspk(:)  = a%aspk(:)
    b%ia1(:)   = a%ia1(:)
    b%ia2(:)   = a%ia2(:)
    b%pl(:)    = a%pl(:)
    b%pr(:)    = a%pr(:)    
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
    Integer             :: nza,nz1, nz2, nzl, nzr
    logical, parameter  :: debug=.false.

    INFO  = 0

    if (associated(b%pr)) then  
      deallocate(b%pr,stat=info)
    end if

    if (associated(b%pl)) then 
      deallocate(b%pl,stat=info)
    end if
    if (associated(b%ia2)) then
      deallocate(b%ia2,stat=info)
    end if
    if (associated(b%ia1))  then
      deallocate(b%ia1,stat=info)
    endif
    if (associated(b%aspk))  then
      deallocate(b%aspk,stat=info)
    endif

    b%aspk  => a%aspk
    b%ia1   => a%ia1
    b%ia2   => a%ia2
    b%pl    => a%pl
    b%pr    => a%pr    
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
    Integer             :: nza,nz1, nz2, nzl, nzr
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

!!$    if (.not.associated(a%infoa)) then 
!!$      info = -2 
!!$      return
!!$    endif
    
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
    
    if (associated(a%aspk)) then 
      val = val + 16 * size(a%aspk)
    endif

    if (associated(a%ia1)) then 
      val = val + 4 * size(a%ia1)
    endif
    if (associated(a%ia2)) then 
      val = val + 4 * size(a%ia2)
    endif
    if (associated(a%pl)) then 
      val = val + 4 * size(a%pl)
    endif
    if (associated(a%pr)) then 
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

    if (associated(a%aspk)) then
      deallocate(a%aspk,STAT=INFO)
    endif
    if ((info == 0) .and. associated(a%ia1)) then
      deallocate(a%ia1,STAT=INFO)
    endif
    if ((info == 0) .and. associated(a%ia2)) then
      deallocate(a%ia2,STAT=INFO)
    endif
    if ((info == 0) .and. associated(a%pr)) then
      deallocate(a%pr,STAT=INFO)
    endif
    if ((info == 0) .and. associated(a%pl)) then
      deallocate(a%pl,STAT=INFO)
    endif
    call psb_nullify_sp(a)
    Return
  End Subroutine psb_zsp_free


end module psb_spmat_type

