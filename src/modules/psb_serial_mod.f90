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
module psb_serial_mod
  use psb_spmat_type
  use psb_string_mod

  interface psb_csdp
     subroutine psb_dcsdp(a, b,info,ifc,check,trans,unitd)
       use psb_spmat_type
       type(psb_dspmat_type), intent(in)   :: a
       type(psb_dspmat_type), intent(inout)  :: b
       integer, intent(out)        :: info
       integer, intent(in), optional :: ifc
       character, intent(in), optional :: check,trans,unitd
     end subroutine psb_dcsdp
  end interface

  interface psb_csrws
     subroutine psb_dcsrws(rw,a,info,trans)
       use psb_spmat_type
       type(psb_dspmat_type) :: a
       real(kind(1.d0)), pointer  :: rw(:) 
       integer :: info
       character, optional :: trans
     end subroutine psb_dcsrws
  end interface



  interface psb_cssm
     subroutine psb_dcssm(alpha,t,b,beta,c,info,trans,unitd,d)
       use psb_spmat_type
       type(psb_dspmat_type) :: t
       real(kind(1.d0)) :: alpha, beta, b(:,:), c(:,:)
       integer :: info
       character, optional :: trans, unitd
       real(kind(1.d0)), optional, target :: d(:)
     end subroutine psb_dcssm
     subroutine psb_dcssv(alpha,t,b,beta,c,info,trans,unitd,d)
       use psb_spmat_type
       type(psb_dspmat_type) :: t
       real(kind(1.d0)) :: alpha, beta, b(:), c(:)
       integer :: info
       character, optional :: trans, unitd
       real(kind(1.d0)), optional, target :: d(:)
     end subroutine psb_dcssv
  end interface

  interface psb_csmm
     subroutine psb_dcsmv(alpha,a,b,beta,c,info,trans)
       use psb_spmat_type
       type(psb_dspmat_type) :: a
       real(kind(1.d0)) :: alpha, beta, b(:), c(:)
       integer :: info
       character, optional :: trans
     end subroutine psb_dcsmv
     subroutine psb_dcsmm(alpha,a,b,beta,c,info,trans)
       use psb_spmat_type
       type(psb_dspmat_type) :: a
       real(kind(1.d0)) :: alpha, beta, b(:,:), c(:,:)
       integer :: info
       character, optional :: trans
     end subroutine psb_dcsmm
  end interface

  interface psb_fixcoo
     subroutine psb_dfixcoo(a,info,idir)
       use psb_spmat_type
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(out)         :: info
       integer, intent(in), optional :: idir
     end subroutine psb_dfixcoo
  end interface

  interface psb_ipcoo2csr
     subroutine psb_dipcoo2csr(a,info,rwshr)
       use psb_spmat_type
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(out)         :: info
       logical, optional :: rwshr
     end subroutine psb_dipcoo2csr
  end interface

  interface psb_ipcoo2csc
     subroutine psb_dipcoo2csc(a,info,clshr)
       use psb_spmat_type
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(out)         :: info
       logical, optional :: clshr
     end subroutine psb_dipcoo2csc
  end interface

  interface psb_ipcsr2coo
     subroutine psb_dipcsr2coo(a,info)
       use psb_spmat_type
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(out)         :: info
     end subroutine psb_dipcsr2coo
  end interface

  interface psb_csprt
     subroutine psb_dcsprt(iout,a,iv,irs,ics,head,ivr,ivc)
       use psb_spmat_type
       integer, intent(in)       :: iout
       type(psb_dspmat_type), intent(in) :: a
       integer, intent(in), optional :: iv(:)
       integer, intent(in), optional :: irs,ics
       character(len=*), optional    :: head
       integer, intent(in), optional :: ivr(:),ivc(:)
     end subroutine psb_dcsprt
  end interface

  interface psb_spgtdiag
     subroutine psb_dspgtdiag(a,d,info)
       use psb_spmat_type
       type(psb_dspmat_type), intent(in)     :: a
       real(kind(1.d0)), intent(inout) :: d(:) 
       integer, intent(out)  :: info
     end subroutine psb_dspgtdiag
  end interface

  interface psb_spscal
     subroutine psb_dspscal(a,d,info)
       use psb_spmat_type
       type(psb_dspmat_type), intent(inout) :: a
       real(kind(1.d0)), intent(in) :: d(:) 
       integer, intent(out)  :: info
     end subroutine psb_dspscal
  end interface


  interface psb_spinfo
     subroutine psb_dspinfo(ireq,a,ires,info,iaux)
       use psb_spmat_type
       type(psb_dspmat_type), intent(in) :: a
       integer, intent(in)       :: ireq
       integer, intent(out)      :: ires
       integer, intent(out)  :: info
       integer, intent(in), optional :: iaux
     end subroutine psb_dspinfo
  end interface

  interface psb_spgtrow
     subroutine psb_dspgtrow(irw,a,b,info,append,iren,lrw)
       use psb_spmat_type
       type(psb_dspmat_type), intent(in) :: a
       integer, intent(in)       :: irw
       type(psb_dspmat_type), intent(inout)    :: b
       logical, intent(in), optional :: append
       integer, intent(in), target, optional :: iren(:)
       integer, intent(in), optional :: lrw
       integer, intent(out)  :: info
     end subroutine psb_dspgtrow
  end interface

  interface psb_neigh
     subroutine psb_dneigh(a,idx,neigh,n,info,lev)
       use psb_spmat_type
       type(psb_dspmat_type), intent(in) :: a   
       integer, intent(in)       :: idx 
       integer, intent(out)      :: n   
       integer, pointer          :: neigh(:)
       integer, intent(out)  :: info
       integer, optional, intent(in) :: lev 
     end subroutine psb_dneigh
  end interface

  interface psb_coins
     subroutine psb_dcoins(nz,ia,ja,val,a,gtl,imin,imax,jmin,jmax,info)
       use psb_spmat_type
       integer, intent(in) :: nz, imin,imax,jmin,jmax
       integer, intent(in) :: ia(:),ja(:),gtl(:)
       real(kind(1.d0)), intent(in) :: val(:)
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(out) :: info
     end subroutine psb_dcoins
  end interface


  interface psb_symbmm
     subroutine psb_dsymbmm(a,b,c)
       use psb_spmat_type
       type(psb_dspmat_type) :: a,b,c
     end subroutine psb_dsymbmm
  end interface

  interface psb_numbmm
     subroutine psb_dnumbmm(a,b,c)
       use psb_spmat_type
       type(psb_dspmat_type) :: a,b,c
     end subroutine psb_dnumbmm
  end interface

  interface psb_transp
     subroutine psb_dtransp(a,b,c,fmt)
       use psb_spmat_type
       type(psb_dspmat_type) :: a,b
       integer, optional :: c
       character(len=*), optional :: fmt
     end subroutine psb_dtransp
  end interface

  interface psb_rwextd
     subroutine psb_drwextd(nr,a,info,b)
       use psb_spmat_type
       integer, intent(in) :: nr
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(out)         :: info
       type(psb_dspmat_type), intent(in), optional  :: b
     end subroutine psb_drwextd
  end interface

  interface psb_csnmi
     real(kind(1.d0)) function psb_dcsnmi(a,info,trans)
       use psb_spmat_type
       type(psb_dspmat_type), intent(in)  :: a
       integer, intent(out)       :: info
       character, optional        :: trans
     end function psb_dcsnmi
  end interface

end module psb_serial_mod

