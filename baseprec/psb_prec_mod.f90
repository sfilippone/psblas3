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

module psb_prec_mod
  use psb_prec_type

  interface psb_precbld
    subroutine psb_dprecbld(a,desc_a,prec,info,upd)
      use psb_base_mod
      use psb_prec_type
      implicit none
      type(psb_dspmat_type), intent(in), target  :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(psb_dprec_type), intent(inout)        :: prec
      integer, intent(out)                       :: info
      character, intent(in),optional             :: upd
    end subroutine psb_dprecbld
    subroutine psb_zprecbld(a,desc_a,prec,info,upd)
      use psb_base_mod
      use psb_prec_type
      implicit none
      type(psb_zspmat_type), intent(in), target  :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(psb_zprec_type), intent(inout)        :: prec
      integer, intent(out)                       :: info
      character, intent(in),optional             :: upd
    end subroutine psb_zprecbld
  end interface

  interface psb_precset
    subroutine psb_dprecset(prec,ptype,info,iv,rs,rv)
      use psb_base_mod
      use psb_prec_type
      implicit none
      type(psb_dprec_type), intent(inout)    :: prec
      character(len=*), intent(in)           :: ptype
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: iv(:)
      real(kind(1.d0)), optional, intent(in) :: rs
      real(kind(1.d0)), optional, intent(in) :: rv(:)
    end subroutine psb_dprecset
    subroutine psb_zprecset(prec,ptype,info,iv,rs,rv)
      use psb_base_mod
      use psb_prec_type
      implicit none
      type(psb_zprec_type), intent(inout)    :: prec
      character(len=*), intent(in)           :: ptype
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: iv(:)
      real(kind(1.d0)), optional, intent(in) :: rs
      real(kind(1.d0)), optional, intent(in) :: rv(:)
    end subroutine psb_zprecset
  end interface


  interface psb_precfree
    subroutine psb_dprecfree(p,info)
      use psb_base_mod
      use psb_prec_type
      type(psb_dprec_type), intent(inout) :: p
      integer, intent(out)                :: info
    end subroutine psb_dprecfree
    subroutine psb_zprecfree(p,info)
      use psb_base_mod
      use psb_prec_type
      type(psb_zprec_type), intent(inout) :: p
      integer, intent(out)                :: info
    end subroutine psb_zprecfree
  end interface

  interface psb_prc_aply
    subroutine psb_dprc_aply(prec,x,y,desc_data,info,trans,work)
      use psb_base_mod
      use psb_prec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_dprec_type), intent(in)  :: prec
      real(kind(0.d0)),intent(inout)    :: x(:), y(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
      real(kind(0.d0)),intent(inout), optional, target :: work(:)
    end subroutine psb_dprc_aply
    subroutine psb_dprc_aply1(prec,x,desc_data,info,trans)
      use psb_base_mod
      use psb_prec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_dprec_type), intent(in)  :: prec
      real(kind(0.d0)),intent(inout)    :: x(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
    end subroutine psb_dprc_aply1
    subroutine psb_zprc_aply(prec,x,y,desc_data,info,trans,work)
      use psb_base_mod
      use psb_prec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_zprec_type), intent(in)  :: prec
      complex(kind(0.d0)),intent(inout) :: x(:), y(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
      complex(kind(0.d0)),intent(inout), optional, target :: work(:)
    end subroutine psb_zprc_aply
    subroutine psb_zprc_aply1(prec,x,desc_data,info,trans)
      use psb_base_mod
      use psb_prec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_zprec_type), intent(in)  :: prec
      complex(kind(0.d0)),intent(inout) :: x(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
    end subroutine psb_zprc_aply1
  end interface

end module psb_prec_mod
