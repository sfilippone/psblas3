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
module psb_comm_mod

  use psb_i2_comm_a_mod
  use psb_m_comm_a_mod
  use psb_e_comm_a_mod
  use psb_s_comm_a_mod
  use psb_d_comm_a_mod
  use psb_c_comm_a_mod
  use psb_z_comm_a_mod

  use psb_i_comm_mod
  use psb_l_comm_mod
  use psb_s_comm_mod
  use psb_d_comm_mod
  use psb_c_comm_mod
  use psb_z_comm_mod

  ! Import scheme symbols and re-export them as public symbols from this module.
  use psb_comm_schemes_mod, only: psb_comm_handle_type, psb_comm_isend_irecv_, &
    & psb_comm_ineighbor_alltoallv_, psb_comm_persistent_ineighbor_alltoallv_, &
    & psb_comm_rma_pull_, psb_comm_rma_push_, psb_comm_unknown_, &
    & psb_comm_status_unknown_, psb_comm_status_start_, psb_comm_status_wait_, &
    & psb_comm_status_sync_

  public :: psb_comm_handle_type, psb_comm_isend_irecv_, psb_comm_ineighbor_alltoallv_, &
    & psb_comm_persistent_ineighbor_alltoallv_, psb_comm_rma_pull_, psb_comm_rma_push_, &
    & psb_comm_unknown_, psb_comm_status_unknown_, psb_comm_status_start_, &
    & psb_comm_status_wait_, psb_comm_status_sync_

end module psb_comm_mod
