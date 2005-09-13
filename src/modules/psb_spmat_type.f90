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
    integer     :: infoa(10)
    ! Contains sparse matrix coefficients
    real(kind(1.d0)), pointer :: aspk(:)=>null()
    ! Contains indeces that describes sparse matrix structure
    integer, pointer :: ia1(:)=>null(), ia2(:)=>null()
    ! Permutations matrix
    integer, pointer :: pl(:)=>null(), pr(:)=>null()
 end type psb_dspmat_type
  
  interface psb_nullify_sp
    module procedure psb_nullify_dsp
  end interface

  interface psb_spclone
    module procedure psb_dspclone
  end interface

  interface psb_spreall
    module procedure psb_dspreallocate, psb_dspreall3
  end interface

  interface psb_spall
    module procedure psb_dspallocate, psb_dspall3, psb_dspallmk, psb_dspallmknz
  end interface

  interface psb_spfree
    module procedure psb_dspfree
  end interface

  interface psb_spreinit
    module procedure psb_dspreinit
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
    call psb_spreall(a,nnz,info)

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
    call psb_spreall(a,nnz,info)

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
    call psb_spreall(a,nnz,info)

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
    
    call psb_spreall(a, ni1,ni2,nd,info)
    
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
    call psb_realloc(nd,a%aspk,info)
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


  subroutine psb_dspfree(a,info)
    implicit none
    !....Parameters...
    Type(psb_dspmat_type), intent(inout)  :: A
    Integer, intent(out)        :: info

    !locals
    logical, parameter  :: debug=.false.

    INFO  = 0

    deallocate(a%aspk,a%ia1,a%ia2,a%pr,a%pl,STAT=INFO)
    
    call psb_nullify_sp(a)

    Return

  End Subroutine psb_dspfree


end module psb_spmat_type

