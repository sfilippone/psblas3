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
!
!
! package: psb_d_linmap_mod
!    Defines data types and interfaces for mapping between vectors belonging
!    to different spaces.
!
module psb_d_linmap_mod

  use psb_const_mod
  use psb_d_mat_mod, only : psb_dspmat_type
  use psb_desc_mod, only : psb_desc_type
  use psb_base_linmap_mod


  type, extends(psb_base_linmap_type) ::  psb_dlinmap_type 
    type(psb_dspmat_type) :: map_X2Y, map_Y2X
  contains
    procedure, pass(map) :: sizeof   => d_map_sizeof
    procedure, pass(map) :: is_asb   => d_is_asb
    procedure, pass(map) :: free     => d_free
    procedure, pass(map) :: clone    => d_clone
  end type psb_dlinmap_type


  interface psb_map_X2Y
    subroutine psb_d_map_X2Y(alpha,x,beta,y,map,info,work)
      import :: psb_ipk_, psb_dpk_, psb_dlinmap_type
      implicit none 
      type(psb_dlinmap_type), intent(in) :: map
      real(psb_dpk_), intent(in)     :: alpha,beta
      real(psb_dpk_), intent(inout)  :: x(:)
      real(psb_dpk_), intent(out)    :: y(:)
      integer(psb_ipk_), intent(out)           :: info 
      real(psb_dpk_), optional       :: work(:)
    end subroutine psb_d_map_X2Y
    subroutine psb_d_map_X2Y_vect(alpha,x,beta,y,map,info,work)
      use psb_d_vect_mod, only : psb_d_vect_type
      import :: psb_ipk_, psb_dpk_, psb_dlinmap_type
      implicit none 
      type(psb_dlinmap_type), intent(in) :: map
      real(psb_dpk_), intent(in)     :: alpha,beta
      type(psb_d_vect_type), intent(inout)  :: x,y
      integer(psb_ipk_), intent(out)           :: info 
      real(psb_dpk_), optional       :: work(:)
    end subroutine psb_d_map_X2Y_vect
  end interface

  interface psb_map_Y2X
    subroutine psb_d_map_Y2X(alpha,x,beta,y,map,info,work)
      import :: psb_ipk_, psb_dpk_, psb_dlinmap_type
      implicit none 
      type(psb_dlinmap_type), intent(in) :: map
      real(psb_dpk_), intent(in)     :: alpha,beta
      real(psb_dpk_), intent(inout)  :: x(:)
      real(psb_dpk_), intent(out)    :: y(:)
      integer(psb_ipk_), intent(out)           :: info 
      real(psb_dpk_), optional       :: work(:)
    end subroutine psb_d_map_Y2X
    subroutine psb_d_map_Y2X_vect(alpha,x,beta,y,map,info,work)
      use psb_d_vect_mod, only : psb_d_vect_type
      import :: psb_ipk_, psb_dpk_, psb_dlinmap_type
      implicit none 
      type(psb_dlinmap_type), intent(in) :: map
      real(psb_dpk_), intent(in)     :: alpha,beta
      type(psb_d_vect_type), intent(inout)  :: x,y
      integer(psb_ipk_), intent(out)           :: info 
      real(psb_dpk_), optional       :: work(:)
    end subroutine psb_d_map_Y2X_vect
  end interface


  interface psb_map_cscnv
    module procedure psb_d_map_cscnv
  end interface

  interface psb_linmap_sub
    module procedure psb_d_linmap_sub
  end interface

  interface psb_move_alloc
    module procedure  psb_dlinmap_transfer
  end interface

  interface psb_linmap
    function psb_d_linmap(map_kind,desc_X, desc_Y, map_X2Y, map_Y2X,iaggr,naggr)
      use psb_d_mat_mod, only : psb_dspmat_type
      import :: psb_ipk_, psb_dlinmap_type, psb_desc_type
      implicit none 
      type(psb_dlinmap_type)                  :: psb_d_linmap    
      type(psb_desc_type), target             :: desc_X, desc_Y
      type(psb_dspmat_type), intent(inout)    :: map_X2Y, map_Y2X
      integer(psb_ipk_), intent(in)           :: map_kind
      integer(psb_ipk_), intent(in), optional :: iaggr(:), naggr(:)
    end function psb_d_linmap
  end interface

  private :: d_map_sizeof, d_is_asb, d_free





contains

  function d_map_sizeof(map) result(val)
    use psb_desc_mod
    use psb_d_mat_mod
    implicit none 
    class(psb_dlinmap_type), intent(in) :: map
    integer(psb_long_int_k_) :: val

    val = map%psb_base_linmap_type%sizeof()
    val = val + map%map_X2Y%sizeof()
    val = val + map%map_Y2X%sizeof()

  end function d_map_sizeof


  function d_is_asb(map) result(val)
    use psb_desc_mod
    implicit none 
    class(psb_dlinmap_type), intent(in) :: map
    logical  :: val

    val = map%psb_base_linmap_type%is_asb() .and. &
         & map%map_X2Y%is_asb() .and.map%map_Y2X%is_asb() 
    
  end function d_is_asb


  subroutine psb_d_map_cscnv(map,info,type,mold)    
    use psb_d_mat_mod
    implicit none
    type(psb_dlinmap_type), intent(inout)  :: map
    integer(psb_ipk_), intent(out)                   :: info
    character(len=*), intent(in), optional :: type
    class(psb_d_base_sparse_mat), intent(in), optional :: mold

    call map%map_X2Y%cscnv(info,type=type,mold=mold)
    if (info == psb_success_)&
         & call map%map_Y2X%cscnv(info,type=type,mold=mold)

  end subroutine psb_d_map_cscnv

  subroutine psb_d_linmap_sub(out_map,map_kind,desc_X, desc_Y,&
       & map_X2Y, map_Y2X,iaggr,naggr)
    use psb_d_mat_mod
    implicit none 
    type(psb_dlinmap_type), intent(out)     :: out_map    
    type(psb_desc_type), target             :: desc_X, desc_Y
    type(psb_dspmat_type), intent(inout)    :: map_X2Y, map_Y2X
    integer(psb_ipk_), intent(in)           :: map_kind
    integer(psb_ipk_), intent(in), optional :: iaggr(:), naggr(:)
    out_map = psb_linmap(map_kind,desc_X,desc_Y,map_X2Y,map_Y2X,iaggr,naggr)
  end subroutine psb_d_linmap_sub

  subroutine  psb_dlinmap_transfer(mapin,mapout,info)
    use psb_realloc_mod
    use psb_desc_mod
    use psb_mat_mod, only : psb_move_alloc
    implicit none 
    type(psb_dlinmap_type) :: mapin,mapout
    integer(psb_ipk_), intent(out)      :: info 
    
    call psb_move_alloc(mapin%psb_base_linmap_type, &
         & mapout%psb_base_linmap_type,info)
    call psb_move_alloc(mapin%map_X2Y,mapout%map_X2Y,info)
    call psb_move_alloc(mapin%map_Y2X,mapout%map_Y2X,info)

  end subroutine psb_dlinmap_transfer

  subroutine  d_free(map,info)
    use psb_desc_mod
    implicit none 
    class(psb_dlinmap_type) :: map
    integer(psb_ipk_), intent(out)      :: info 
    
    call map%psb_base_linmap_type%free(info)
    
    call map%map_X2Y%free()
    call map%map_Y2X%free()

  end subroutine d_free
  

  subroutine  d_clone(map,mapout,info)
    use psb_desc_mod
    implicit none 
    class(psb_dlinmap_type), intent(inout) :: map
    class(psb_dlinmap_type), intent(inout) :: mapout
    integer(psb_ipk_)     :: info 
    
    call mapout%free(info)
    ! Base clone!    
    if (info == 0) call &
         & map%psb_base_linmap_type%clone(mapout%psb_base_linmap_type,info)
    if (info == 0) call map%map_X2Y%clone(mapout%map_X2Y,info)
    if (info == 0) call map%map_Y2X%clone(mapout%map_Y2X,info)

  end subroutine d_clone
  

end module psb_d_linmap_mod

