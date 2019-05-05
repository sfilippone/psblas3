!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!    
!
!
! package: psb_d_linmap_mod
!    Defines data types and interfaces for mapping between vectors belonging
!    to different spaces.
!
module psb_d_linmap_mod

  use psb_const_mod
  use psb_d_mat_mod
  use psb_desc_mod
  use psb_base_linmap_mod


  type, extends(psb_base_linmap_type) ::  psb_dlinmap_type 
    type(psb_dspmat_type) :: mat_U2V, mat_V2U
  contains
    procedure, pass(map) :: map_U2V_a  => psb_d_map_U2V_a
    procedure, pass(map) :: map_U2V_v  => psb_d_map_U2V_v
    generic, public      :: map_U2V => map_U2V_a, map_U2V_v

    procedure, pass(map) :: map_V2U_a  => psb_d_map_V2U_a
    procedure, pass(map) :: map_V2U_v  => psb_d_map_V2U_v
    generic, public      :: map_V2U => map_V2U_a, map_V2U_v
    
    procedure, pass(map) :: sizeof   => d_map_sizeof
    procedure, pass(map) :: is_asb   => d_is_asb
    procedure, pass(map) :: free     => d_free
    procedure, pass(map) :: clone    => d_clone
    procedure, pass(map) :: cnv      => psb_d_map_cscnv
  end type psb_dlinmap_type


  interface psb_map_U2V
    subroutine psb_d_map_U2V_a(alpha,x,beta,y,map,info,work)
      import :: psb_ipk_, psb_dpk_, psb_dlinmap_type
      implicit none 
      class(psb_dlinmap_type), intent(in) :: map
      real(psb_dpk_), intent(in)     :: alpha,beta
      real(psb_dpk_), intent(inout)  :: x(:)
      real(psb_dpk_), intent(out)    :: y(:)
      integer(psb_ipk_), intent(out)           :: info 
      real(psb_dpk_), optional       :: work(:)
    end subroutine psb_d_map_U2V_a
    subroutine psb_d_map_U2V_v(alpha,x,beta,y,map,info,work,vtx,vty)
      use psb_d_vect_mod, only : psb_d_vect_type
      import :: psb_ipk_, psb_dpk_, psb_dlinmap_type
      implicit none 
      class(psb_dlinmap_type), intent(in) :: map
      real(psb_dpk_), intent(in)     :: alpha,beta
      type(psb_d_vect_type), intent(inout)  :: x,y
      integer(psb_ipk_), intent(out)           :: info 
      real(psb_dpk_), optional       :: work(:)
      type(psb_d_vect_type), optional, target, intent(inout)  :: vtx,vty
    end subroutine psb_d_map_U2V_v
  end interface

  interface psb_map_V2U
    subroutine psb_d_map_V2U_a(alpha,x,beta,y,map,info,work)
      import :: psb_ipk_, psb_dpk_, psb_dlinmap_type
      implicit none 
      class(psb_dlinmap_type), intent(in) :: map
      real(psb_dpk_), intent(in)     :: alpha,beta
      real(psb_dpk_), intent(inout)  :: x(:)
      real(psb_dpk_), intent(out)    :: y(:)
      integer(psb_ipk_), intent(out)           :: info 
      real(psb_dpk_), optional       :: work(:)
    end subroutine psb_d_map_V2U_a
    subroutine psb_d_map_V2U_v(alpha,x,beta,y,map,info,work,vtx,vty)
      use psb_d_vect_mod, only : psb_d_vect_type
      import :: psb_ipk_, psb_dpk_, psb_dlinmap_type
      implicit none 
      class(psb_dlinmap_type), intent(in) :: map
      real(psb_dpk_), intent(in)     :: alpha,beta
      type(psb_d_vect_type), intent(inout)  :: x,y
      integer(psb_ipk_), intent(out)           :: info 
      real(psb_dpk_), optional       :: work(:)
      type(psb_d_vect_type), optional, target, intent(inout)  :: vtx,vty
    end subroutine psb_d_map_V2U_v
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
    function psb_d_linmap(map_kind,desc_U, desc_V, mat_U2V, mat_V2U,iaggr,naggr)
      use psb_d_mat_mod, only : psb_dspmat_type
      import :: psb_ipk_, psb_dlinmap_type, psb_desc_type, psb_lpk_
      implicit none 
      type(psb_dlinmap_type)                  :: psb_d_linmap    
      type(psb_desc_type), target             :: desc_U, desc_V
      type(psb_dspmat_type), intent(inout)    :: mat_U2V, mat_V2U
      integer(psb_ipk_), intent(in)           :: map_kind
      integer(psb_lpk_), intent(in), optional :: iaggr(:), naggr(:)
    end function psb_d_linmap
  end interface

  private :: d_map_sizeof, d_is_asb, d_free





contains

  function d_map_sizeof(map) result(val)
    implicit none 
    class(psb_dlinmap_type), intent(in) :: map
    integer(psb_epk_) :: val

    val = map%psb_base_linmap_type%sizeof()
    val = val + map%mat_U2V%sizeof()
    val = val + map%mat_V2U%sizeof()

  end function d_map_sizeof


  function d_is_asb(map) result(val)
    implicit none 
    class(psb_dlinmap_type), intent(in) :: map
    logical  :: val

    val = map%psb_base_linmap_type%is_asb() .and. &
         & map%mat_U2V%is_asb() .and.map%mat_V2U%is_asb() 
    
  end function d_is_asb


  subroutine psb_d_map_cscnv(map,info,type,mold,imold)    
    implicit none
    class(psb_dlinmap_type), intent(inout)  :: map
    integer(psb_ipk_), intent(out)                   :: info
    character(len=*), intent(in), optional :: type
    class(psb_d_base_sparse_mat), intent(in), optional :: mold
    class(psb_i_base_vect_type), intent(in), optional  :: imold

    if (map%mat_U2V%is_asb())&
         & call map%mat_U2V%cscnv(info,type=type,mold=mold)
    if (info == psb_success_ .and.map%mat_V2U%is_asb())&
         & call map%mat_V2U%cscnv(info,type=type,mold=mold)
    if (present(imold)) then 
      call map%desc_U%cnv(mold=imold)
      call map%desc_V%cnv(mold=imold)
    end if

  end subroutine psb_d_map_cscnv

  subroutine psb_d_linmap_sub(out_map,map_kind,desc_U, desc_V,&
       & mat_U2V, mat_V2U,iaggr,naggr)
    use psb_d_mat_mod
    implicit none 
    type(psb_dlinmap_type), intent(out)     :: out_map    
    type(psb_desc_type), target             :: desc_U, desc_V
    type(psb_dspmat_type), intent(inout)    :: mat_U2V, mat_V2U
    integer(psb_ipk_), intent(in)           :: map_kind
    integer(psb_lpk_), intent(in), optional :: iaggr(:), naggr(:)
    out_map = psb_linmap(map_kind,desc_U,desc_V,mat_U2V,mat_V2U,iaggr,naggr)
  end subroutine psb_d_linmap_sub

  subroutine  psb_dlinmap_transfer(mapin,mapout,info)
    use psb_realloc_mod
    implicit none 
    type(psb_dlinmap_type) :: mapin,mapout
    integer(psb_ipk_), intent(out)      :: info 
    
    call psb_move_alloc(mapin%psb_base_linmap_type, &
         & mapout%psb_base_linmap_type,info)
    call psb_move_alloc(mapin%mat_U2V,mapout%mat_U2V,info)
    call psb_move_alloc(mapin%mat_V2U,mapout%mat_V2U,info)

  end subroutine psb_dlinmap_transfer

  subroutine  d_free(map,info)
    implicit none 
    class(psb_dlinmap_type) :: map
    integer(psb_ipk_), intent(out)      :: info 
    
    call map%psb_base_linmap_type%free(info)
    
    call map%mat_U2V%free()
    call map%mat_V2U%free()

  end subroutine d_free
  

  subroutine  d_clone(map,mapout,info)
    use psb_error_mod
    implicit none 
    class(psb_dlinmap_type), intent(inout) :: map
    class(psb_base_linmap_type), intent(inout) :: mapout
    integer(psb_ipk_)     :: info 
    
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='clone'

    info = 0
    select type(mout => mapout)
    class is (psb_dlinmap_type)
      call mout%free(info)
      ! Base clone!    
      if (info == 0) call &
           & map%psb_base_linmap_type%clone(mout%psb_base_linmap_type,info)
      if (info == 0) call map%mat_U2V%clone(mout%mat_U2V,info)
      if (info == 0) call map%mat_V2U%clone(mout%mat_V2U,info)
    class default
      info = psb_err_invalid_dynamic_type_
      info = psb_err_missing_override_method_
      call psb_errpush(info,name,m_err=(/2/))
      call psb_erractionsave(err_act)

      call psb_error_handler(err_act)
    end select

      
  end subroutine d_clone
  

end module psb_d_linmap_mod

