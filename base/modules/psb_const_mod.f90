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

module psb_const_mod
  ! This should be integer(8), and normally different from default integer. 
  integer, parameter  :: longndig=12
  integer, parameter  :: psb_long_int_k_ = selected_int_kind(longndig)
  !
  ! These must be the kind parameter corresponding to MPI_DOUBLE_PRECISION
  ! and MPI_REAL
  !
  integer, parameter  :: psb_dpk_ = kind(1.d0)
  integer, parameter  :: psb_spk_ = kind(1.e0)
  integer             :: psb_sizeof_dp, psb_sizeof_sp
  integer             :: psb_sizeof_int, psb_sizeof_long_int

  !
  !     Handy & miscellaneous constants
  !
  integer, parameter             :: izero=0, ione=1
  integer, parameter             :: itwo=2, ithree=3,mone=-1, psb_root_=0
  real(psb_spk_), parameter      :: szero=0.e0, sone=1.e0
  real(psb_dpk_), parameter      :: dzero=0.d0, done=1.d0
  complex(psb_spk_), parameter   :: czero=(0.e0,0.0e0)
  complex(psb_spk_), parameter   :: cone=(1.e0,0.0e0)
  complex(psb_dpk_), parameter   :: zzero=(0.d0,0.0d0)
  complex(psb_dpk_), parameter   :: zone=(1.d0,0.0d0)
  real(psb_dpk_), parameter      :: epstol=1.d-32, psb_percent_=0.7
  real(psb_spk_), parameter      :: s_epstol=1.e-16 ! How to choose this?
  character, parameter           :: psb_all_='A',  psb_topdef_=' '

end module psb_const_mod
