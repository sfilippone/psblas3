!   
!                Parallel Sparse BLAS  version 3.4
!      (C) Copyright 2006, 2010, 2015
!                         Salvatore Filippone    University of Rome Tor Vergata
!                         Alfredo Buttari        CNRS-IRIT, Toulouse
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
module psb_blockpart_mod
  public part_block, bld_partblock
  
contains
  subroutine part_block(global_indx,n,np,pv,nv)
    use psb_base_mod, only : psb_ipk_, psb_mpik_
    implicit none

    integer(psb_ipk_), intent(in)  ::  global_indx, n
    integer(psb_ipk_), intent(in)  ::  np
    integer(psb_ipk_), intent(out) ::  nv
    integer(psb_ipk_), intent(out) ::  pv(*)
    integer(psb_ipk_) :: dim_block

    dim_block = (n + np - 1)/np
    nv = 1  
    pv(nv) = (global_indx - 1) / dim_block

    return
  end subroutine part_block
      



  subroutine bld_partblock(n,np,ivg)      
    use psb_base_mod, only : psb_ipk_
    integer(psb_ipk_) :: n,np,ivg(*)

    integer(psb_ipk_) :: dim_block,i


    dim_block = (n + np - 1)/np
    do i=1,n
      ivg(i) = (i - 1) / dim_block
    enddo

  end subroutine bld_partblock



end module psb_blockpart_mod

