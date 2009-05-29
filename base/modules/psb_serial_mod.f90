!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
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
module psb_serial_mod
  use psb_const_mod
  use psb_spmat_type
  use psb_string_mod
  use psb_sort_mod

  use psi_serial_mod, &
       & psb_gth => psi_gth,&
       & psb_sct => psi_sct

  interface psb_csrws
    subroutine psb_dcsrws(rw,a,info,trans)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type) :: a
      real(psb_dpk_), allocatable   :: rw(:) 
      integer :: info
      character, optional :: trans
    end subroutine psb_dcsrws
    subroutine psb_zcsrws(rw,a,info,trans)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type) :: a
      complex(psb_dpk_), allocatable :: rw(:) 
      integer :: info
      character, optional :: trans
    end subroutine psb_zcsrws
  end interface

!!$  interface psb_cssm
!!$    subroutine psb_scssm(alpha,t,b,beta,c,info,trans,unitd,d)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_sspmat_type) :: t
!!$      real(psb_spk_) :: alpha, beta, b(:,:), c(:,:)
!!$      integer :: info
!!$      character, optional :: trans, unitd
!!$      real(psb_spk_), optional, target :: d(:)
!!$    end subroutine psb_scssm
!!$    subroutine psb_scssv(alpha,t,b,beta,c,info,trans,unitd,d)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_sspmat_type) :: t
!!$      real(psb_spk_) :: alpha, beta, b(:), c(:)
!!$      integer :: info
!!$      character, optional :: trans, unitd
!!$      real(psb_spk_), optional, target :: d(:)
!!$    end subroutine psb_scssv
!!$    subroutine psb_dcssm(alpha,t,b,beta,c,info,trans,unitd,d)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_dspmat_type) :: t
!!$      real(psb_dpk_) :: alpha, beta, b(:,:), c(:,:)
!!$      integer :: info
!!$      character, optional :: trans, unitd
!!$      real(psb_dpk_), optional, target :: d(:)
!!$    end subroutine psb_dcssm
!!$    subroutine psb_dcssv(alpha,t,b,beta,c,info,trans,unitd,d)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_dspmat_type) :: t
!!$      real(psb_dpk_) :: alpha, beta, b(:), c(:)
!!$      integer :: info
!!$      character, optional :: trans, unitd
!!$      real(psb_dpk_), optional, target :: d(:)
!!$    end subroutine psb_dcssv
!!$    subroutine psb_ccssm(alpha,t,b,beta,c,info,trans,unitd,d)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_cspmat_type) :: t
!!$      complex(psb_spk_) :: alpha, beta, b(:,:), c(:,:)
!!$      integer :: info
!!$      character, optional :: trans, unitd
!!$      complex(psb_spk_), optional, target :: d(:)
!!$    end subroutine psb_ccssm
!!$    subroutine psb_ccssv(alpha,t,b,beta,c,info,trans,unitd,d)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_cspmat_type) :: t
!!$      complex(psb_spk_) :: alpha, beta, b(:), c(:)
!!$      integer :: info
!!$      character, optional :: trans, unitd
!!$      complex(psb_spk_), optional, target :: d(:)
!!$    end subroutine psb_ccssv
!!$    subroutine psb_zcssm(alpha,t,b,beta,c,info,trans,unitd,d)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_zspmat_type) :: t
!!$      complex(psb_dpk_) :: alpha, beta, b(:,:), c(:,:)
!!$      integer :: info
!!$      character, optional :: trans, unitd
!!$      complex(psb_dpk_), optional, target :: d(:)
!!$    end subroutine psb_zcssm
!!$    subroutine psb_zcssv(alpha,t,b,beta,c,info,trans,unitd,d)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_zspmat_type) :: t
!!$      complex(psb_dpk_) :: alpha, beta, b(:), c(:)
!!$      integer :: info
!!$      character, optional :: trans, unitd
!!$      complex(psb_dpk_), optional, target :: d(:)
!!$    end subroutine psb_zcssv
!!$  end interface

!!$  interface psb_csmm
!!$    module procedure psb_scsmm, psb_scsmv, psb_dcsmm, psb_dcsmv,&
!!$         & psb_ccsmm, psb_ccsmv, psb_zcsmm, psb_zcsmv
!!$    subroutine psb_scsmv(alpha,a,b,beta,c,info,trans)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_sspmat_type) :: a
!!$      real(psb_spk_) :: alpha, beta, b(:), c(:)
!!$      integer :: info
!!$      character, optional :: trans
!!$    end subroutine psb_scsmv
!!$    subroutine psb_scsmm(alpha,a,b,beta,c,info,trans)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_sspmat_type) :: a
!!$      real(psb_spk_) :: alpha, beta, b(:,:), c(:,:)
!!$      integer :: info
!!$      character, optional :: trans
!!$    end subroutine psb_scsmm
!!$    subroutine psb_dcsmv(alpha,a,b,beta,c,info,trans)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_dspmat_type) :: a
!!$      real(psb_dpk_) :: alpha, beta, b(:), c(:)
!!$      integer :: info
!!$      character, optional :: trans
!!$    end subroutine psb_dcsmv
!!$    subroutine psb_dcsmm(alpha,a,b,beta,c,info,trans)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_dspmat_type) :: a
!!$      real(psb_dpk_) :: alpha, beta, b(:,:), c(:,:)
!!$      integer :: info
!!$      character, optional :: trans
!!$    end subroutine psb_dcsmm
!!$    subroutine psb_ccsmv(alpha,a,b,beta,c,info,trans)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_cspmat_type) :: a
!!$      complex(psb_spk_) :: alpha, beta, b(:), c(:)
!!$      integer :: info
!!$      character, optional :: trans
!!$    end subroutine psb_ccsmv
!!$    subroutine psb_ccsmm(alpha,a,b,beta,c,info,trans)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_cspmat_type) :: a
!!$      complex(psb_spk_) :: alpha, beta, b(:,:), c(:,:)
!!$      integer :: info
!!$      character, optional :: trans
!!$    end subroutine psb_ccsmm
!!$    subroutine psb_zcsmv(alpha,a,b,beta,c,info,trans)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_zspmat_type) :: a
!!$      complex(psb_dpk_) :: alpha, beta, b(:), c(:)
!!$      integer :: info
!!$      character, optional :: trans
!!$    end subroutine psb_zcsmv
!!$    subroutine psb_zcsmm(alpha,a,b,beta,c,info,trans)
!!$      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
!!$           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
!!$      type(psb_zspmat_type) :: a
!!$      complex(psb_dpk_) :: alpha, beta, b(:,:), c(:,:)
!!$      integer :: info
!!$      character, optional :: trans
!!$    end subroutine psb_zcsmm
!!$  end interface

  interface psb_cest
    subroutine psb_cest(afmt, m,n,nnz, lia1, lia2, lar, iup, info)
      integer, intent(in) ::  m,n,nnz,iup
      integer, intent(out) :: lia1, lia2, lar, info
      character(len=*), intent(inout) :: afmt
    end subroutine psb_cest
  end interface

  interface psb_spcnv
    subroutine psb_sspcnv2(ain, a, info, afmt, upd, dupl)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_sspmat_type), intent (in)     :: ain
      type(psb_sspmat_type), intent (out)    :: a
      integer, intent(out)                   :: info
      integer,optional, intent(in)           :: dupl, upd
      character(len=*), optional, intent(in) :: afmt
    end subroutine psb_sspcnv2
    subroutine psb_sspcnv1(a, info, afmt, upd, dupl)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_sspmat_type), intent (inout)   :: a
      integer, intent(out)                    :: info
      integer,optional, intent(in)            :: dupl, upd
      character(len=*), optional, intent(in)  :: afmt
    end subroutine psb_sspcnv1
    subroutine psb_dspcnv2(ain, a, info, afmt, upd, dupl)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type), intent (in)     :: ain
      type(psb_dspmat_type), intent (out)    :: a
      integer, intent(out)                   :: info
      integer,optional, intent(in)           :: dupl, upd
      character(len=*), optional, intent(in) :: afmt
    end subroutine psb_dspcnv2
    subroutine psb_dspcnv1(a, info, afmt, upd, dupl)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type), intent (inout)   :: a
      integer, intent(out)                    :: info
      integer,optional, intent(in)            :: dupl, upd
      character(len=*), optional, intent(in)  :: afmt
    end subroutine psb_dspcnv1
    subroutine psb_cspcnv2(ain, a, info, afmt, upd, dupl)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_cspmat_type), intent (in)     :: ain
      type(psb_cspmat_type), intent (out)    :: a
      integer, intent(out)                   :: info
      integer,optional, intent(in)           :: dupl, upd
      character(len=*), optional, intent(in) :: afmt
    end subroutine psb_cspcnv2
    subroutine psb_cspcnv1(a, info, afmt, upd, dupl)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_cspmat_type), intent (inout)   :: a
      integer, intent(out)                    :: info
      integer,optional, intent(in)            :: dupl, upd
      character(len=*), optional, intent(in)  :: afmt
    end subroutine psb_cspcnv1
    subroutine psb_zspcnv2(ain, a, info, afmt, upd, dupl)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type), intent (in)     :: ain
      type(psb_zspmat_type), intent (out)    :: a
      integer, intent(out)                   :: info
      integer,optional, intent(in)           :: dupl, upd
      character(len=*), optional, intent(in) :: afmt
    end subroutine psb_zspcnv2
    subroutine psb_zspcnv1(a, info, afmt, upd, dupl)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type), intent (inout)   :: a
      integer, intent(out)                    :: info
      integer,optional, intent(in)            :: dupl, upd
      character(len=*), optional, intent(in)  :: afmt
    end subroutine psb_zspcnv1
  end interface



  interface psb_fixcoo
    subroutine psb_sfixcoo(a,info,idir)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_sspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      integer, intent(in), optional :: idir
    end subroutine psb_sfixcoo
    subroutine psb_dfixcoo(a,info,idir)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      integer, intent(in), optional :: idir
    end subroutine psb_dfixcoo
    subroutine psb_cfixcoo(a,info,idir)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_cspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      integer, intent(in), optional :: idir
    end subroutine psb_cfixcoo
    subroutine psb_zfixcoo(a,info,idir)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      integer, intent(in), optional :: idir
    end subroutine psb_zfixcoo
  end interface

  interface psb_ipcoo2csr
    subroutine psb_sipcoo2csr(a,info,rwshr)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_sspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      logical, optional :: rwshr
    end subroutine psb_sipcoo2csr
    subroutine psb_dipcoo2csr(a,info,rwshr)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      logical, optional :: rwshr
    end subroutine psb_dipcoo2csr
    subroutine psb_cipcoo2csr(a,info,rwshr)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_cspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      logical, optional :: rwshr
    end subroutine psb_cipcoo2csr
    subroutine psb_zipcoo2csr(a,info,rwshr)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      logical, optional :: rwshr
    end subroutine psb_zipcoo2csr
  end interface

  interface psb_ipcoo2csc
    subroutine psb_sipcoo2csc(a,info,clshr)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_sspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      logical, optional :: clshr
    end subroutine psb_sipcoo2csc
    subroutine psb_dipcoo2csc(a,info,clshr)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      logical, optional :: clshr
    end subroutine psb_dipcoo2csc
    subroutine psb_cipcoo2csc(a,info,clshr)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_cspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      logical, optional :: clshr
    end subroutine psb_cipcoo2csc
    subroutine psb_zipcoo2csc(a,info,clshr)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      logical, optional :: clshr
    end subroutine psb_zipcoo2csc
  end interface

  interface psb_ipcsr2coo
    subroutine psb_sipcsr2coo(a,info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_sspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
    end subroutine psb_sipcsr2coo
    subroutine psb_dipcsr2coo(a,info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
    end subroutine psb_dipcsr2coo
    subroutine psb_cipcsr2coo(a,info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_cspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
    end subroutine psb_cipcsr2coo
    subroutine psb_zipcsr2coo(a,info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
    end subroutine psb_zipcsr2coo
  end interface

  interface psb_csprt
    subroutine psb_scsprt(iout,a,iv,irs,ics,head,ivr,ivc)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      integer, intent(in)       :: iout
      type(psb_sspmat_type), intent(in) :: a
      integer, intent(in), optional :: iv(:)
      integer, intent(in), optional :: irs,ics
      character(len=*), optional    :: head
      integer, intent(in), optional :: ivr(:),ivc(:)
    end subroutine psb_scsprt
    subroutine psb_dcsprt(iout,a,iv,irs,ics,head,ivr,ivc)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      integer, intent(in)       :: iout
      type(psb_dspmat_type), intent(in) :: a
      integer, intent(in), optional :: iv(:)
      integer, intent(in), optional :: irs,ics
      character(len=*), optional    :: head
      integer, intent(in), optional :: ivr(:),ivc(:)
    end subroutine psb_dcsprt
    subroutine psb_ccsprt(iout,a,iv,irs,ics,head,ivr,ivc)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      integer, intent(in)       :: iout
      type(psb_cspmat_type), intent(in) :: a
      integer, intent(in), optional :: iv(:)
      integer, intent(in), optional :: irs,ics
      character(len=*), optional    :: head
      integer, intent(in), optional :: ivr(:),ivc(:)
    end subroutine psb_ccsprt
    subroutine psb_zcsprt(iout,a,iv,irs,ics,head,ivr,ivc)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      integer, intent(in)       :: iout
      type(psb_zspmat_type), intent(in) :: a
      integer, intent(in), optional :: iv(:)
      integer, intent(in), optional :: irs,ics
      character(len=*), optional    :: head
      integer, intent(in), optional :: ivr(:),ivc(:)
    end subroutine psb_zcsprt
  end interface

  interface psb_neigh
    subroutine psb_sneigh(a,idx,neigh,n,info,lev)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_sspmat_type), intent(in) :: a   
      integer, intent(in)       :: idx 
      integer, intent(out)      :: n   
      integer, allocatable          :: neigh(:)
      integer, intent(out)  :: info
      integer, optional, intent(in) :: lev 
    end subroutine psb_sneigh
    subroutine psb_dneigh(a,idx,neigh,n,info,lev)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type), intent(in) :: a   
      integer, intent(in)       :: idx 
      integer, intent(out)      :: n   
      integer, allocatable          :: neigh(:)
      integer, intent(out)  :: info
      integer, optional, intent(in) :: lev 
    end subroutine psb_dneigh
    subroutine psb_cneigh(a,idx,neigh,n,info,lev)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_cspmat_type), intent(in) :: a   
      integer, intent(in)       :: idx 
      integer, intent(out)      :: n   
      integer, allocatable      :: neigh(:)
      integer, intent(out)  :: info
      integer, optional, intent(in) :: lev 
    end subroutine psb_cneigh
    subroutine psb_zneigh(a,idx,neigh,n,info,lev)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type), intent(in) :: a   
      integer, intent(in)       :: idx 
      integer, intent(out)      :: n   
      integer, allocatable      :: neigh(:)
      integer, intent(out)  :: info
      integer, optional, intent(in) :: lev 
    end subroutine psb_zneigh
  end interface

  interface psb_coins
    subroutine psb_scoins(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl,rebuild)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      integer, intent(in) :: nz, imin,imax,jmin,jmax
      integer, intent(in) :: ia(:),ja(:)
      real(psb_spk_), intent(in) :: val(:)
      type(psb_sspmat_type), intent(inout) :: a
      integer, intent(out) :: info
      integer, intent(in), optional :: gtl(:)
      logical, optional, intent(in) :: rebuild
    end subroutine psb_scoins
    subroutine psb_dcoins(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl,rebuild)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      integer, intent(in) :: nz, imin,imax,jmin,jmax
      integer, intent(in) :: ia(:),ja(:)
      real(psb_dpk_), intent(in) :: val(:)
      type(psb_dspmat_type), intent(inout) :: a
      integer, intent(out) :: info
      integer, intent(in), optional :: gtl(:)
      logical, optional, intent(in) :: rebuild
    end subroutine psb_dcoins
    subroutine psb_ccoins(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl,rebuild)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      integer, intent(in) :: nz, imin,imax,jmin,jmax
      integer, intent(in) :: ia(:),ja(:)
      complex(psb_spk_), intent(in) :: val(:)
      type(psb_cspmat_type), intent(inout) :: a
      integer, intent(out) :: info
      integer, intent(in), optional :: gtl(:)
      logical, optional, intent(in) :: rebuild
    end subroutine psb_ccoins
    subroutine psb_zcoins(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl,rebuild)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      integer, intent(in) :: nz, imin,imax,jmin,jmax
      integer, intent(in) :: ia(:),ja(:)
      complex(psb_dpk_), intent(in) :: val(:)
      type(psb_zspmat_type), intent(inout) :: a
      integer, intent(out) :: info
      integer, intent(in), optional :: gtl(:)
      logical, optional, intent(in) :: rebuild
    end subroutine psb_zcoins
  end interface


  interface psb_symbmm
    subroutine psb_ssymbmm(a,b,c,info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_sspmat_type) :: a,b,c
      integer               :: info
    end subroutine psb_ssymbmm
    subroutine psb_dsymbmm(a,b,c,info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type) :: a,b,c
      integer               :: info
    end subroutine psb_dsymbmm
    subroutine psb_csymbmm(a,b,c,info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_cspmat_type) :: a,b,c
      integer               :: info
    end subroutine psb_csymbmm
    subroutine psb_zsymbmm(a,b,c,info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type) :: a,b,c
      integer               :: info
    end subroutine psb_zsymbmm
  end interface

  interface psb_numbmm
    subroutine psb_snumbmm(a,b,c)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_sspmat_type) :: a,b,c
    end subroutine psb_snumbmm
    subroutine psb_dnumbmm(a,b,c)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type) :: a,b,c
    end subroutine psb_dnumbmm
    subroutine psb_cnumbmm(a,b,c)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_cspmat_type) :: a,b,c
    end subroutine psb_cnumbmm
    subroutine psb_znumbmm(a,b,c)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type) :: a,b,c
    end subroutine psb_znumbmm
  end interface

  interface psb_transp
    subroutine psb_stransp(a,b,c,fmt)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_sspmat_type), intent(in) :: a
      type(psb_sspmat_type), intent(out)   :: b
      integer, optional :: c
      character(len=*), optional :: fmt
    end subroutine psb_stransp
    subroutine psb_dtransp(a,b,c,fmt)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type), intent(in) :: a
      type(psb_dspmat_type), intent(out)   :: b
      integer, optional :: c
      character(len=*), optional :: fmt
    end subroutine psb_dtransp
    subroutine psb_ctransp(a,b,c,fmt)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_cspmat_type), intent(in) :: a
      type(psb_cspmat_type), intent(out)   :: b
      integer, optional :: c
      character(len=*), optional :: fmt
    end subroutine psb_ctransp
    subroutine psb_ztransp(a,b,c,fmt)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type), intent(in) :: a
      type(psb_zspmat_type), intent(out)   :: b
      integer, optional :: c
      character(len=*), optional :: fmt
    end subroutine psb_ztransp
    subroutine psb_stransp1(a,c,fmt)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_sspmat_type), intent(inout) :: a
      integer, optional :: c
      character(len=*), optional :: fmt
    end subroutine psb_stransp1
    subroutine psb_dtransp1(a,c,fmt)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type), intent(inout) :: a
      integer, optional :: c
      character(len=*), optional :: fmt
    end subroutine psb_dtransp1
    subroutine psb_ctransp1(a,c,fmt)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_cspmat_type), intent(inout) :: a
      integer, optional :: c
      character(len=*), optional :: fmt
    end subroutine psb_ctransp1
    subroutine psb_ztransp1(a,c,fmt)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type), intent(inout) :: a
      integer, optional :: c
      character(len=*), optional :: fmt
    end subroutine psb_ztransp1
  end interface

  interface psb_transc
    subroutine psb_ctransc(a,b,c,fmt)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_cspmat_type), intent(in) :: a
      type(psb_cspmat_type), intent(out)   :: b
      integer, optional :: c
      character(len=*), optional :: fmt
    end subroutine psb_ctransc
    subroutine psb_ztransc(a,b,c,fmt)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type), intent(in) :: a
      type(psb_zspmat_type), intent(out)   :: b
      integer, optional :: c
      character(len=*), optional :: fmt
    end subroutine psb_ztransc
  end interface

  interface psb_rwextd
    subroutine psb_srwextd(nr,a,info,b,rowscale)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      integer, intent(in) :: nr
      type(psb_sspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      type(psb_sspmat_type), intent(in), optional  :: b
      logical, intent(in), optional :: rowscale
    end subroutine psb_srwextd
    subroutine psb_drwextd(nr,a,info,b,rowscale)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      integer, intent(in) :: nr
      type(psb_dspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      type(psb_dspmat_type), intent(in), optional  :: b
      logical, intent(in), optional :: rowscale
    end subroutine psb_drwextd
    subroutine psb_crwextd(nr,a,info,b,rowscale)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      integer, intent(in) :: nr
      type(psb_cspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      type(psb_cspmat_type), intent(in), optional  :: b
      logical, intent(in), optional :: rowscale
    end subroutine psb_crwextd
    subroutine psb_zrwextd(nr,a,info,b,rowscale)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      integer, intent(in) :: nr
      type(psb_zspmat_type), intent(inout) :: a
      integer, intent(out)         :: info
      type(psb_zspmat_type), intent(in), optional  :: b
      logical, intent(in), optional :: rowscale
    end subroutine psb_zrwextd
  end interface

  interface psb_csnmi
    function psb_scsnmi(a,info,trans)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_sspmat_type), intent(in)  :: a
      integer, intent(out)       :: info
      character, optional        :: trans
      real(psb_spk_)             :: psb_scsnmi
    end function psb_scsnmi
    function psb_dcsnmi(a,info,trans)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type), intent(in)  :: a
      integer, intent(out)       :: info
      character, optional        :: trans
      real(psb_dpk_)             :: psb_dcsnmi
    end function psb_dcsnmi
    function psb_ccsnmi(a,info,trans)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_cspmat_type), intent(in)  :: a
      integer, intent(out)       :: info
      character, optional        :: trans
      real(psb_spk_)             :: psb_ccsnmi
    end function psb_ccsnmi
    function psb_zcsnmi(a,info,trans)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type), intent(in)  :: a
      integer, intent(out)       :: info
      character, optional        :: trans
      real(psb_dpk_)             :: psb_zcsnmi
    end function psb_zcsnmi
  end interface

  interface psb_sp_clip
    subroutine psb_sspclip(a,b,info,imin,imax,jmin,jmax,rscale,cscale)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      implicit none 
      type(psb_sspmat_type), intent(in)  :: a
      type(psb_sspmat_type), intent(out) :: b
      integer, intent(out)               :: info
      integer, intent(in), optional      :: imin,imax,jmin,jmax
      logical, intent(in), optional      :: rscale,cscale
    end subroutine psb_sspclip
    subroutine psb_dspclip(a,b,info,imin,imax,jmin,jmax,rscale,cscale)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      implicit none 
      type(psb_dspmat_type), intent(in)  :: a
      type(psb_dspmat_type), intent(out) :: b
      integer, intent(out)               :: info
      integer, intent(in), optional      :: imin,imax,jmin,jmax
      logical, intent(in), optional      :: rscale,cscale
    end subroutine psb_dspclip
    subroutine psb_cspclip(a,b,info,imin,imax,jmin,jmax,rscale,cscale)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      implicit none 
      type(psb_cspmat_type), intent(in)  :: a
      type(psb_cspmat_type), intent(out) :: b
      integer, intent(out)               :: info
      integer, intent(in), optional      :: imin,imax,jmin,jmax
      logical, intent(in), optional      :: rscale,cscale
    end subroutine psb_cspclip
    subroutine psb_zspclip(a,b,info,imin,imax,jmin,jmax,rscale,cscale)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      implicit none 
      type(psb_zspmat_type), intent(in)  :: a
      type(psb_zspmat_type), intent(out) :: b
      integer, intent(out)               :: info
      integer, intent(in), optional      :: imin,imax,jmin,jmax
      logical, intent(in), optional      :: rscale,cscale
    end subroutine psb_zspclip
  end interface

  interface psb_sp_getdiag
     subroutine psb_sspgtdiag(a,d,info)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_sspmat_type), intent(in)     :: a
       real(psb_spk_), intent(inout) :: d(:) 
       integer, intent(out)  :: info
     end subroutine psb_sspgtdiag
     subroutine psb_dspgtdiag(a,d,info)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_dspmat_type), intent(in)     :: a
       real(psb_dpk_), intent(inout) :: d(:) 
       integer, intent(out)  :: info
     end subroutine psb_dspgtdiag
     subroutine psb_cspgtdiag(a,d,info)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_cspmat_type), intent(in)     :: a
       complex(psb_spk_), intent(inout) :: d(:) 
       integer, intent(out)  :: info
     end subroutine psb_cspgtdiag
     subroutine psb_zspgtdiag(a,d,info)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_zspmat_type), intent(in)     :: a
       complex(psb_dpk_), intent(inout) :: d(:) 
       integer, intent(out)  :: info
     end subroutine psb_zspgtdiag
  end interface

  interface psb_sp_scal
     subroutine psb_sspscals(a,d,info)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_sspmat_type), intent(inout) :: a
       real(psb_spk_), intent(in) :: d 
       integer, intent(out)  :: info
     end subroutine psb_sspscals
     subroutine psb_sspscal(a,d,info)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_sspmat_type), intent(inout) :: a
       real(psb_spk_), intent(in) :: d(:) 
       integer, intent(out)  :: info
     end subroutine psb_sspscal
     subroutine psb_dspscals(a,d,info)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_dspmat_type), intent(inout) :: a
       real(psb_dpk_), intent(in) :: d 
       integer, intent(out)  :: info
     end subroutine psb_dspscals
     subroutine psb_dspscal(a,d,info)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_dspmat_type), intent(inout) :: a
       real(psb_dpk_), intent(in) :: d(:) 
       integer, intent(out)  :: info
     end subroutine psb_dspscal
     subroutine psb_cspscals(a,d,info)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_cspmat_type), intent(inout) :: a
       complex(psb_spk_), intent(in) :: d
       integer, intent(out)  :: info
     end subroutine psb_cspscals
     subroutine psb_cspscal(a,d,info)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_cspmat_type), intent(inout) :: a
       complex(psb_spk_), intent(in) :: d(:) 
       integer, intent(out)  :: info
     end subroutine psb_cspscal
     subroutine psb_zspscals(a,d,info)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_zspmat_type), intent(inout) :: a
       complex(psb_dpk_), intent(in) :: d
       integer, intent(out)  :: info
     end subroutine psb_zspscals
     subroutine psb_zspscal(a,d,info)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_zspmat_type), intent(inout) :: a
       complex(psb_dpk_), intent(in) :: d(:) 
       integer, intent(out)  :: info
     end subroutine psb_zspscal
  end interface


  interface psb_sp_setbld
    subroutine psb_dspsetbld1(a,info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type), intent(inout) :: a
      integer, intent(out)               :: info      
    end subroutine psb_dspsetbld1
    subroutine psb_dspsetbld2(a,b,info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type), intent(in)  :: a
      type(psb_dspmat_type), intent(out) :: b
      integer, intent(out)               :: info      
    end subroutine psb_dspsetbld2
    subroutine psb_zspsetbld1(a,info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type), intent(inout) :: a
      integer, intent(out)               :: info      
    end subroutine psb_zspsetbld1
    subroutine psb_zspsetbld2(a,b,info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type), intent(in)  :: a
      type(psb_zspmat_type), intent(out) :: b
      integer, intent(out)               :: info      
    end subroutine psb_zspsetbld2
  end interface

  interface psb_sp_shift
    subroutine psb_dspshift(alpha,a,beta,b,info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type), intent(in)  :: a
      type(psb_dspmat_type), intent(out) :: b
      real(psb_dpk_), intent(in)         :: alpha, beta
      integer, intent(out)               :: info      
    end subroutine psb_dspshift
    subroutine psb_zspshift(alpha,a,beta,b,info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type), intent(in)  :: a
      type(psb_zspmat_type), intent(out) :: b
      complex(psb_dpk_), intent(in)      :: alpha, beta
      integer, intent(out)               :: info      
    end subroutine psb_zspshift
  end interface

  interface psb_sp_getblk
     subroutine psb_sspgtblk(irw,a,b,info,append,iren,lrw,srt)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_sspmat_type), intent(in) :: a
       integer, intent(in)       :: irw
       type(psb_sspmat_type), intent(inout)    :: b
       integer, intent(out)  :: info
       logical, intent(in), optional :: append
       integer, intent(in), target, optional :: iren(:)
       integer, intent(in), optional :: lrw
       logical, intent(in), optional :: srt
     end subroutine psb_sspgtblk
     subroutine psb_dspgtblk(irw,a,b,info,append,iren,lrw,srt)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_dspmat_type), intent(in) :: a
       integer, intent(in)       :: irw
       type(psb_dspmat_type), intent(inout)    :: b
       integer, intent(out)  :: info
       logical, intent(in), optional :: append
       integer, intent(in), target, optional :: iren(:)
       integer, intent(in), optional :: lrw
       logical, intent(in), optional :: srt
     end subroutine psb_dspgtblk
     subroutine psb_cspgtblk(irw,a,b,info,append,iren,lrw,srt)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_cspmat_type), intent(in) :: a
       integer, intent(in)       :: irw
       type(psb_cspmat_type), intent(inout)    :: b
       integer, intent(out)  :: info
       logical, intent(in), optional :: append
       integer, intent(in), target, optional :: iren(:)
       integer, intent(in), optional :: lrw
       logical, intent(in), optional :: srt
     end subroutine psb_cspgtblk
     subroutine psb_zspgtblk(irw,a,b,info,append,iren,lrw,srt)
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       type(psb_zspmat_type), intent(in) :: a
       integer, intent(in)       :: irw
       type(psb_zspmat_type), intent(inout)    :: b
       integer, intent(out)  :: info
       logical, intent(in), optional :: append
       integer, intent(in), target, optional :: iren(:)
       integer, intent(in), optional :: lrw
       logical, intent(in), optional :: srt
     end subroutine psb_zspgtblk
  end interface

  interface psb_sp_getrow
     subroutine  psb_sspgetrow(irw,a,nz,ia,ja,val,info,iren,lrw,append,nzin)
       ! Output is always in  COO format 
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       implicit none

       type(psb_sspmat_type), intent(in)    :: a
       integer, intent(in)                  :: irw
       integer, intent(out)                 :: nz
       integer, allocatable, intent(inout)  :: ia(:), ja(:)
       real(psb_spk_), allocatable,  intent(inout)    :: val(:)
       integer,intent(out)                  :: info
       logical, intent(in), optional        :: append
       integer, intent(in), optional        :: iren(:)
       integer, intent(in), optional        :: lrw, nzin
     end subroutine psb_sspgetrow
     subroutine  psb_dspgetrow(irw,a,nz,ia,ja,val,info,iren,lrw,append,nzin)
       ! Output is always in  COO format 
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       implicit none

       type(psb_dspmat_type), intent(in)    :: a
       integer, intent(in)                  :: irw
       integer, intent(out)                 :: nz
       integer, allocatable, intent(inout)  :: ia(:), ja(:)
       real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
       integer,intent(out)                  :: info
       logical, intent(in), optional        :: append
       integer, intent(in), optional        :: iren(:)
       integer, intent(in), optional        :: lrw, nzin
     end subroutine psb_dspgetrow
     subroutine  psb_cspgetrow(irw,a,nz,ia,ja,val,info,iren,lrw,append,nzin)
       ! Output is always in  COO format 
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       implicit none

       type(psb_cspmat_type), intent(in)    :: a
       integer, intent(in)                  :: irw
       integer, intent(out)                 :: nz
       integer, allocatable, intent(inout)  :: ia(:), ja(:)
       complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
       integer,intent(out)                  :: info
       logical, intent(in), optional        :: append
       integer, intent(in), optional        :: iren(:)
       integer, intent(in), optional        :: lrw, nzin
     end subroutine psb_cspgetrow
     subroutine  psb_zspgetrow(irw,a,nz,ia,ja,val,info,iren,lrw,append,nzin)
       ! Output is always in  COO format 
       use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
            & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
       implicit none

       type(psb_zspmat_type), intent(in)    :: a
       integer, intent(in)                  :: irw
       integer, intent(out)                 :: nz
       integer, allocatable, intent(inout)  :: ia(:), ja(:)
       complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
       integer,intent(out)                  :: info
       logical, intent(in), optional        :: append
       integer, intent(in), optional        :: iren(:)
       integer, intent(in), optional        :: lrw, nzin
     end subroutine psb_zspgetrow
  end interface



  interface psb_csrp
    subroutine psb_dcsrp(trans,iperm,a, info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_dspmat_type), intent(inout)  ::  a
      integer, intent(inout)                :: iperm(:), info
      character, intent(in)                 :: trans
    end subroutine psb_dcsrp
    subroutine psb_zcsrp(trans,iperm,a, info)
      use psb_spmat_type, only : psb_sspmat_type, psb_dspmat_type,&
           & psb_cspmat_type, psb_zspmat_type, psb_spk_, psb_dpk_
      type(psb_zspmat_type), intent(inout)  ::  a
      integer, intent(inout)                :: iperm(:), info
      character, intent(in)                 :: trans
    end subroutine psb_zcsrp
  end interface


  interface psb_gelp
    ! 2-D version
    subroutine psb_sgelp(trans,iperm,x,info)
      use psb_const_mod
      real(psb_spk_), intent(inout)      ::  x(:,:)
      integer, intent(in)                  ::  iperm(:)
      integer, intent(out)                 ::  info
      character, intent(in)                :: trans
    end subroutine psb_sgelp
    ! 1-D version
    subroutine psb_sgelpv(trans,iperm,x,info)
      use psb_const_mod
      real(psb_spk_), intent(inout)    ::  x(:)
      integer, intent(in)                  ::  iperm(:)
      integer, intent(out)                 ::  info
      character, intent(in)              :: trans
    end subroutine psb_sgelpv
    subroutine psb_dgelp(trans,iperm,x,info)
      use psb_const_mod
      real(psb_dpk_), intent(inout)      ::  x(:,:)
      integer, intent(in)                  ::  iperm(:)
      integer, intent(out)                 ::  info
      character, intent(in)                :: trans
    end subroutine psb_dgelp
    ! 1-D version
    subroutine psb_dgelpv(trans,iperm,x,info)
      use psb_const_mod
      real(psb_dpk_), intent(inout)    ::  x(:)
      integer, intent(in)                  ::  iperm(:)
      integer, intent(out)                 ::  info
      character, intent(in)              :: trans
    end subroutine psb_dgelpv
    ! 2-D version
    subroutine psb_cgelp(trans,iperm,x,info)
      use psb_const_mod
      complex(psb_spk_), intent(inout)      ::  x(:,:)
      integer, intent(in)                  ::  iperm(:)
      integer, intent(out)                 ::  info
      character, intent(in)                :: trans
    end subroutine psb_cgelp
    ! 1-D version
    subroutine psb_cgelpv(trans,iperm,x,info)
      use psb_const_mod
      complex(psb_spk_), intent(inout)    ::  x(:)
      integer, intent(in)                  ::  iperm(:)
      integer, intent(out)                 ::  info
      character, intent(in)              :: trans
    end subroutine psb_cgelpv
    ! 2-D version
    subroutine psb_zgelp(trans,iperm,x,info)
      use psb_const_mod
      complex(psb_dpk_), intent(inout)      ::  x(:,:)
      integer, intent(in)                  ::  iperm(:)
      integer, intent(out)                 ::  info
      character, intent(in)                :: trans
    end subroutine psb_zgelp
    ! 1-D version
    subroutine psb_zgelpv(trans,iperm,x,info)
      use psb_const_mod
      complex(psb_dpk_), intent(inout)    ::  x(:)
      integer, intent(in)                  ::  iperm(:)
      integer, intent(out)                 ::  info
      character, intent(in)              :: trans
    end subroutine psb_zgelpv
  end interface


end module psb_serial_mod

