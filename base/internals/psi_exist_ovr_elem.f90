!
!             Parallel Sparse BLAS  version 3.5
!   (C) Copyright 2006, 2010, 2015, 2017
!       Salvatore Filippone     
!       Alfredo Buttari        CNRS-IRIT, Toulouse
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!   1. Redistributions of source code must retain the above copyright
!      notice, this list of conditions and the following disclaimer.
!   2. Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions, and the following disclaimer in the
!      documentation and/or other materials provided with the distribution.
!   3. The name of the PSBLAS group or the names of its contributors may
!      not be used to endorse or promote products derived from this
!      software without specific written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
! 
integer function psi_exist_ovr_elem(ovr_elem, dim_list,elem_searched)
  use psb_const_mod
  !    PURPOSE:
  !    == = ====
  !
  !    If ELEM_SEARCHED exist in the list OVR_ELEM returns its position in
  !    the list, else returns -1
  !
  !
  !     INPUT
  !     == = ===
  !     OVRLAP_ELEMENT_D.: Contains for all overlap points belonging to 
  !                        the current process:
  !                          1. overlap point index
  !                          2. Number of domains sharing that overlap point
  !                        the end is marked by a -1...............................
  !
  !    DIM_LIST..........: Dimension of list OVRLAP_ELEMENT_D
  !
  !    ELEM_SEARCHED.....:point's  Local index identifier to be searched.

  implicit none

  !     ....Scalars parameters....
  integer(psb_ipk_) :: dim_list,elem_searched
  !     ...array parameters....
  integer(psb_ipk_) :: ovr_elem(dim_list,*)

  !     ...local scalars....
  integer(psb_ipk_) :: i

  i=1
  do while ((i.le.dim_list).and.(ovr_elem(i,1).ne.elem_searched))
    i=i+1
  enddo
  if ((i.le.dim_list).and.(ovr_elem(i,1).eq.elem_searched)) then
    psi_exist_ovr_elem=i
  else
    psi_exist_ovr_elem=-1
  endif
end function psi_exist_ovr_elem

