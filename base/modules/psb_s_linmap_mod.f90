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
! package: psb_s_linmap_mod
!    Defines data types and interfaces for mapping between vectors belonging
!    to different spaces.
!
module psb_s_linmap_mod

  use psb_const_mod
  use psb_s_mat_mod, only : psb_sspmat_type
  use psb_desc_mod, only : psb_desc_type
  use psb_base_linmap_mod


  type, extends(psb_base_linmap_type) ::  psb_slinmap_type 
    type(psb_sspmat_type) :: map_X2Y, map_Y2X
  contains
    procedure, pass(map) :: sizeof   => s_map_sizeof
    procedure, pass(map) :: is_asb   => s_is_asb
    procedure, pass(map) :: free     => s_free
    procedure, pass(map) :: clone    => s_clone
  end type psb_slinmap_type


  interface psb_map_X2Y
    subroutine psb_s_map_X2Y(alpha,x,beta,y,map,info,work)
      import :: psb_ipk_, psb_spk_, psb_slinmap_type
      implicit none 
      type(psb_slinmap_type), intent(in) :: map
      real(psb_spk_), intent(in)     :: alpha,beta
      real(psb_spk_), intent(inout)  :: x(:)
      real(psb_spk_), intent(out)    :: y(:)
      integer(psb_ipk_), intent(out)           :: info 
      real(psb_spk_), optional       :: work(:)
    end subroutine psb_s_map_X2Y
    subroutine psb_s_map_X2Y_vect(alpha,x,beta,y,map,info,work)
      use psb_s_vect_mod, only : psb_s_vect_type
      import :: psb_ipk_, psb_spk_, psb_slinmap_type
      implicit none 
      type(psb_slinmap_type), intent(in) :: map
      real(psb_spk_), intent(in)     :: alpha,beta
      type(psb_s_vect_type), intent(inout)  :: x,y
      integer(psb_ipk_), intent(out)           :: info 
      real(psb_spk_), optional       :: work(:)
    end subroutine psb_s_map_X2Y_vect
  end interface

  interface psb_map_Y2X
    subroutine psb_s_map_Y2X(alpha,x,beta,y,map,info,work)
      import :: psb_ipk_, psb_spk_, psb_slinmap_type
      implicit none 
      type(psb_slinmap_type), intent(in) :: map
      real(psb_spk_), intent(in)     :: alpha,beta
      real(psb_spk_), intent(inout)  :: x(:)
      real(psb_spk_), intent(out)    :: y(:)
      integer(psb_ipk_), intent(out)           :: info 
      real(psb_spk_), optional       :: work(:)
    end subroutine psb_s_map_Y2X
    subroutine psb_s_map_Y2X_vect(alpha,x,beta,y,map,info,work)
      use psb_s_vect_mod, only : psb_s_vect_type
      import :: psb_ipk_, psb_spk_, psb_slinmap_type
      implicit none 
      type(psb_slinmap_type), intent(in) :: map
      real(psb_spk_), intent(in)     :: alpha,beta
      type(psb_s_vect_type), intent(inout)  :: x,y
      integer(psb_ipk_), intent(out)           :: info 
      real(psb_spk_), optional       :: work(:)
    end subroutine psb_s_map_Y2X_vect
  end interface


  interface psb_map_cscnv
    module procedure psb_s_map_cscnv
  end interface

  interface psb_linmap_sub
    module procedure psb_s_linmap_sub
  end interface

  interface psb_move_alloc
    module procedure  psb_slinmap_transfer
  end interface

  interface psb_linmap
    function psb_s_linmap(map_kind,desc_X, desc_Y, map_X2Y, map_Y2X,iaggr,naggr)
      use psb_s_mat_mod, only : psb_sspmat_type
      import :: psb_ipk_, psb_slinmap_type, psb_desc_type
      implicit none 
      type(psb_slinmap_type)                  :: psb_s_linmap    
      type(psb_desc_type), target             :: desc_X, desc_Y
      type(psb_sspmat_type), intent(inout)    :: map_X2Y, map_Y2X
      integer(psb_ipk_), intent(in)           :: map_kind
      integer(psb_ipk_), intent(in), optional :: iaggr(:), naggr(:)
    end function psb_s_linmap
  end interface

  private :: s_map_sizeof, s_is_asb, s_free





contains

  function s_map_sizeof(map) result(val)
    use psb_desc_mod
    use psb_s_mat_mod
    implicit none 
    class(psb_slinmap_type), intent(in) :: map
    integer(psb_long_int_k_) :: val

    val = map%psb_base_linmap_type%sizeof()
    val = val + map%map_X2Y%sizeof()
    val = val + map%map_Y2X%sizeof()

  end function s_map_sizeof


  function s_is_asb(map) result(val)
    use psb_desc_mod
    implicit none 
    class(psb_slinmap_type), intent(in) :: map
    logical  :: val

    val = map%psb_base_linmap_type%is_asb() .and. &
         & map%map_X2Y%is_asb() .and.map%map_Y2X%is_asb() 
    
  end function s_is_asb


  subroutine psb_s_map_cscnv(map,info,type,mold)    
    use psb_s_mat_mod
    implicit none
    type(psb_slinmap_type), intent(inout)  :: map
    integer(psb_ipk_), intent(out)                   :: info
    character(len=*), intent(in), optional :: type
    class(psb_s_base_sparse_mat), intent(in), optional :: mold

    call map%map_X2Y%cscnv(info,type=type,mold=mold)
    if (info == psb_success_)&
         & call map%map_Y2X%cscnv(info,type=type,mold=mold)

  end subroutine psb_s_map_cscnv

  subroutine psb_s_linmap_sub(out_map,map_kind,desc_X, desc_Y,&
       & map_X2Y, map_Y2X,iaggr,naggr)
    use psb_s_mat_mod
    implicit none 
    type(psb_slinmap_type), intent(out)     :: out_map    
    type(psb_desc_type), target             :: desc_X, desc_Y
    type(psb_sspmat_type), intent(inout)    :: map_X2Y, map_Y2X
    integer(psb_ipk_), intent(in)           :: map_kind
    integer(psb_ipk_), intent(in), optional :: iaggr(:), naggr(:)
    out_map = psb_linmap(map_kind,desc_X,desc_Y,map_X2Y,map_Y2X,iaggr,naggr)
  end subroutine psb_s_linmap_sub

  subroutine  psb_slinmap_transfer(mapin,mapout,info)
    use psb_realloc_mod
    use psb_desc_mod
    use psb_mat_mod, only : psb_move_alloc
    implicit none 
    type(psb_slinmap_type) :: mapin,mapout
    integer(psb_ipk_), intent(out)      :: info 
    
    call psb_move_alloc(mapin%psb_base_linmap_type, &
         & mapout%psb_base_linmap_type,info)
    call psb_move_alloc(mapin%map_X2Y,mapout%map_X2Y,info)
    call psb_move_alloc(mapin%map_Y2X,mapout%map_Y2X,info)

  end subroutine psb_slinmap_transfer

  subroutine  s_free(map,info)
    use psb_desc_mod
    implicit none 
    class(psb_slinmap_type) :: map
    integer(psb_ipk_), intent(out)      :: info 
    
    call map%psb_base_linmap_type%free(info)
    
    call map%map_X2Y%free()
    call map%map_Y2X%free()

  end subroutine s_free
  

  subroutine  s_clone(map,mapout,info)
    use psb_desc_mod
    implicit none 
    class(psb_slinmap_type), intent(inout) :: map
    class(psb_slinmap_type), intent(inout) :: mapout
    integer(psb_ipk_)     :: info 
    
    call mapout%free(info)
    ! Base clone!    
    if (info == 0) call &
         & map%psb_base_linmap_type%clone(mapout%psb_base_linmap_type,info)
    if (info == 0) call map%map_X2Y%clone(mapout%map_X2Y,info)
    if (info == 0) call map%map_Y2X%clone(mapout%map_Y2X,info)

  end subroutine s_clone
  

end module psb_s_linmap_mod

