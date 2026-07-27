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
!         software without specific prior written permission.
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
module psb_i_mmio_mod


  use psb_base_mod, only :  psb_ipk_, psb_lpk_,&
       & psb_i_vect_type, psb_l_vect_type

  !public mm_mat_read, mm_mat_write, mm_array_read, mm_array_write

  interface mm_array_read
    subroutine mm_ivet_read(b, info, iunit, filename)   
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_), allocatable, intent(out)  :: b(:)
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_ivet_read
    subroutine mm_ivet2_read(b, info, iunit, filename)   
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_), allocatable, intent(out)  :: b(:,:)
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_ivet2_read
    subroutine mm_ivect_read(b, info, iunit, filename)   
      import :: psb_ipk_, psb_i_vect_type
      implicit none
      type(psb_i_vect_type), intent(inout)  :: b
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in) :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_ivect_read
    subroutine mm_lvect_read(b, info, iunit, filename)   
      import :: psb_ipk_, psb_l_vect_type
      implicit none
      type(psb_l_vect_type), intent(inout)  :: b
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in) :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_lvect_read
#if defined(PSB_IPK4) && defined(PSB_LPK8) 
    subroutine mm_lvet_read(b, info, iunit, filename)   
      import :: psb_ipk_, psb_lpk_
      implicit none
      integer(psb_lpk_), allocatable, intent(out)  :: b(:)
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_lvet_read
    subroutine mm_lvet2_read(b, info, iunit, filename)   
      import :: psb_ipk_, psb_lpk_
      implicit none
      integer(psb_lpk_), allocatable, intent(out)  :: b(:,:)
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_lvet2_read
#endif
  end interface mm_array_read

  interface mm_array_write
    subroutine mm_ivet2_write(b, header, info, iunit, filename)   
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_), intent(in)  :: b(:,:)
      character(len=*), intent(in) :: header
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_ivet2_write
    subroutine mm_ivet1_write(b, header, info, iunit, filename)   
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_), intent(in)  :: b(:)
      character(len=*), intent(in) :: header
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_ivet1_write
#if defined(PSB_IPK4) && defined(PSB_LPK8) 
    subroutine mm_lvet2_write(b, header, info, iunit, filename)   
      import :: psb_ipk_, psb_lpk_
      implicit none
      integer(psb_lpk_), intent(in)  :: b(:,:)
      character(len=*), intent(in) :: header
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_lvet2_write
    subroutine mm_lvet1_write(b, header, info, iunit, filename)   
      import :: psb_ipk_, psb_lpk_
      implicit none
      integer(psb_lpk_), intent(in)  :: b(:)
      character(len=*), intent(in) :: header
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_lvet1_write
    subroutine mm_ivect_write(b, header, info, iunit, filename)   
      import :: psb_ipk_,psb_i_vect_type
      implicit none
      type(psb_i_vect_type), intent(inout)  :: b
      character(len=*), intent(in) :: header
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_ivect_write
    subroutine mm_lvect_write(b, header, info, iunit, filename)   
      import :: psb_ipk_,psb_l_vect_type
      implicit none
      type(psb_l_vect_type), intent(inout)  :: b
      character(len=*), intent(in) :: header
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_lvect_write
#endif
  end interface mm_array_write

  interface mm_vet_write
    procedure mm_ivet1_write, mm_ivet2_write
  end interface

#if 0
#if ! defined(PSB_HAVE_BUGGY_GENERICS)
  public mm_vet_read, mm_vet_write
#endif
#endif
end module psb_i_mmio_mod
