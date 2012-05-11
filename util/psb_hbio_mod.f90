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
module psb_hbio_mod

  use psb_base_mod, only :  psb_ipk_, psb_spk_, psb_dpk_,&
       & psb_sspmat_type, psb_cspmat_type, &
       & psb_dspmat_type, psb_zspmat_type


  public hb_read, hb_write
  interface hb_read
    subroutine shb_read(a, iret, iunit, filename,b,g,x,mtitle)   
      import :: psb_sspmat_type, psb_spk_, psb_ipk_
      implicit none
      type(psb_sspmat_type), intent(out)     :: a
      integer(psb_ipk_), intent(out)                   :: iret
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
      real(psb_spk_), optional, allocatable, intent(out)  :: b(:,:), g(:,:), x(:,:) 
      character(len=72), optional, intent(out) :: mtitle
    end subroutine shb_read
    subroutine dhb_read(a, iret, iunit, filename,b,g,x,mtitle)   
      import :: psb_dspmat_type, psb_dpk_, psb_ipk_
      implicit none
      type(psb_dspmat_type), intent(out)     :: a
      integer(psb_ipk_), intent(out)                   :: iret
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
      real(psb_dpk_), optional, allocatable, intent(out)  :: b(:,:), g(:,:), x(:,:) 
      character(len=72), optional, intent(out) :: mtitle
    end subroutine dhb_read
    subroutine chb_read(a, iret, iunit, filename,b,g,x,mtitle)   
      import :: psb_cspmat_type, psb_spk_, psb_ipk_
      implicit none
      type(psb_cspmat_type), intent(out)     :: a
      integer(psb_ipk_), intent(out)                   :: iret
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
      complex(psb_spk_), optional, allocatable, intent(out)  :: b(:,:), g(:,:), x(:,:) 
      character(len=72), optional, intent(out) :: mtitle
    end subroutine chb_read
    subroutine zhb_read(a, iret, iunit, filename,b,g,x,mtitle)   
      import :: psb_zspmat_type, psb_dpk_, psb_ipk_
      implicit none
      type(psb_zspmat_type), intent(out)     :: a
      integer(psb_ipk_), intent(out)                   :: iret
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
      complex(psb_dpk_), optional, allocatable, intent(out)  :: b(:,:), g(:,:), x(:,:) 
      character(len=72), optional, intent(out) :: mtitle
    end subroutine zhb_read
  end interface

  interface hb_write
    subroutine shb_write(a,iret,iunit,filename,key,rhs,g,x,mtitle)
      import :: psb_sspmat_type, psb_spk_, psb_ipk_
      implicit none
      type(psb_sspmat_type), intent(inout)  :: a
      integer(psb_ipk_), intent(out)        :: iret
      character(len=*), optional, intent(in) :: mtitle
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
      character(len=*), optional, intent(in) :: key
      real(psb_spk_), optional             :: rhs(:), g(:), x(:)
    end subroutine shb_write
    subroutine dhb_write(a,iret,iunit,filename,key,rhs,g,x,mtitle)
      import :: psb_dspmat_type, psb_dpk_, psb_ipk_
      implicit none
      type(psb_dspmat_type), intent(inout)  :: a
      integer(psb_ipk_), intent(out)        :: iret
      character(len=*), optional, intent(in) :: mtitle
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
      character(len=*), optional, intent(in) :: key
      real(psb_dpk_), optional             :: rhs(:), g(:), x(:)
    end subroutine dhb_write
    subroutine chb_write(a,iret,iunit,filename,key,rhs,g,x,mtitle)
      import :: psb_cspmat_type, psb_spk_, psb_ipk_
      implicit none
      type(psb_cspmat_type), intent(inout)  :: a
      integer(psb_ipk_), intent(out)        :: iret
      character(len=*), optional, intent(in) :: mtitle
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
      character(len=*), optional, intent(in) :: key
      complex(psb_spk_), optional             :: rhs(:), g(:), x(:)
    end subroutine chb_write
    subroutine zhb_write(a,iret,iunit,filename,key,rhs,g,x,mtitle)
      import :: psb_zspmat_type, psb_dpk_, psb_ipk_
      implicit none
      type(psb_zspmat_type), intent(inout)  :: a
      integer(psb_ipk_), intent(out)        :: iret
      character(len=*), optional, intent(in) :: mtitle
      integer(psb_ipk_), optional, intent(in)          :: iunit
      character(len=*), optional, intent(in) :: filename
      character(len=*), optional, intent(in) :: key
      complex(psb_dpk_), optional             :: rhs(:), g(:), x(:)
    end subroutine zhb_write
  end interface

end module psb_hbio_mod
