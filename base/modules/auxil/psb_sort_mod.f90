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
!  The merge-sort and quicksort routines are implemented in the
!  serial/aux directory
!  References:
!  D. Knuth
!  The Art of Computer Programming, vol. 3
!  Addison-Wesley
!  
!  Aho, Hopcroft, Ullman
!  Data Structures and Algorithms
!  Addison-Wesley
!

module psb_sort_mod
  use psb_const_mod
  use psb_ip_reord_mod
  use psi_serial_mod
  
  use psb_m_hsort_mod
  use psb_m_isort_mod
  use psb_m_msort_mod
  use psb_m_qsort_mod
  
  use psb_e_hsort_mod
  use psb_e_isort_mod
  use psb_e_msort_mod
  use psb_e_qsort_mod
  
  use psb_s_hsort_mod
  use psb_s_isort_mod
  use psb_s_msort_mod
  use psb_s_qsort_mod
  
  use psb_d_hsort_mod
  use psb_d_isort_mod
  use psb_d_msort_mod
  use psb_d_qsort_mod
  
  use psb_c_hsort_mod
  use psb_c_isort_mod
  use psb_c_msort_mod
  use psb_c_qsort_mod
  
  use psb_z_hsort_mod
  use psb_z_isort_mod
  use psb_z_msort_mod
  use psb_z_qsort_mod

  use psb_i_hsort_x_mod
  use psb_l_hsort_x_mod
  use psb_s_hsort_x_mod
  use psb_d_hsort_x_mod
  use psb_c_hsort_x_mod
  use psb_z_hsort_x_mod
  
end module psb_sort_mod
