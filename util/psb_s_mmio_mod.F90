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
module psb_s_mmio_mod


  use psb_base_mod, only :  psb_ipk_, psb_lpk_, psb_spk_,&
       & psb_s_vect_type, psb_sspmat_type, psb_lsspmat_type

  !public mm_mat_read, mm_mat_write, mm_array_read, mm_array_write
  interface mm_array_read
    subroutine mm_svet_read(b, info, iunit, filename)   
      import :: psb_spk_, psb_ipk_
      implicit none
      real(psb_spk_), allocatable, intent(out)  :: b(:)
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_svet_read
    subroutine mm_svet2_read(b, info, iunit, filename)   
      import :: psb_spk_, psb_ipk_
      implicit none
      real(psb_spk_), allocatable, intent(out)  :: b(:,:)
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_svet2_read
    subroutine mm_svect_read(b, info, iunit, filename)   
      import :: psb_ipk_,psb_s_vect_type
      implicit none
      type(psb_s_vect_type), intent(inout)  :: b
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in) :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_svect_read
  end interface mm_array_read


#if 0
#if ! defined(PSB_HAVE_BUGGY_GENERICS)
  interface mm_vet_read
    procedure mm_svet_read, mm_svet2_read
  end interface
#endif
#endif
  
  interface mm_array_write
    subroutine mm_svet2_write(b, header, info, iunit, filename)   
      import :: psb_spk_, psb_ipk_
      implicit none
      real(psb_spk_), intent(in)  :: b(:,:)
      character(len=*), intent(in) :: header
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_svet2_write
    subroutine mm_svet1_write(b, header, info, iunit, filename)   
      import :: psb_spk_, psb_ipk_
      implicit none
      real(psb_spk_), intent(in)  :: b(:)
      character(len=*), intent(in) :: header
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine mm_svet1_write
    subroutine mm_svect_write(b, header, info, iunit, filename)   
      import :: psb_ipk_, psb_s_vect_type
      implicit none
      type(psb_s_vect_type), intent(inout)    :: b
      character(len=*), intent(in)            :: header
      integer(psb_ipk_), intent(out)          :: info
      integer(psb_ipk_), optional, intent(in) :: iunit
      character(len=*), optional, intent(in)  :: filename
    end subroutine mm_svect_write
  end interface

#if ! defined(PSB_HAVE_BUGGY_GENERICS)
  interface mm_vet_write
    procedure mm_svet1_write, mm_svet2_write
  end interface
#endif

  interface mm_mat_read
    subroutine smm_mat_read(a, info, iunit, filename)   
      import :: psb_sspmat_type, psb_ipk_
      implicit none
      type(psb_sspmat_type), intent(out)  :: a
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine smm_mat_read
    subroutine lsmm_mat_read(a, info, iunit, filename)   
      import :: psb_lsspmat_type, psb_ipk_
      implicit none
      type(psb_lsspmat_type), intent(out)  :: a
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine lsmm_mat_read
  end interface

  interface mm_mat_write
    subroutine smm_mat_write(a,mtitle,info,iunit,filename)
      import :: psb_sspmat_type, psb_ipk_ 
      implicit none
      type(psb_sspmat_type), intent(in)  :: a
      integer(psb_ipk_), intent(out)        :: info
      character(len=*), intent(in) :: mtitle
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine smm_mat_write
    subroutine lsmm_mat_write(a,mtitle,info,iunit,filename)
      import :: psb_lsspmat_type, psb_ipk_ 
      implicit none
      type(psb_lsspmat_type), intent(in)  :: a
      integer(psb_ipk_), intent(out)        :: info
      character(len=*), intent(in) :: mtitle
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
    end subroutine lsmm_mat_write
  end interface
#if 0
#if ! defined(PSB_HAVE_BUGGY_GENERICS)
  public mm_vet_read, mm_vet_write
#endif
#endif
end module psb_s_mmio_mod
