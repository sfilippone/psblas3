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
integer function psi_list_search(list,lenght_list,elem)
  use psb_const_mod
  implicit none 
  !returns position of elem in a array list
  !of lenght lenght_list, if this element does not exist
  !returns -1
  integer(psb_ipk_) :: list(*)
  integer(psb_ipk_) :: lenght_list
  integer(psb_ipk_) :: elem

  integer(psb_ipk_) :: i

  i=1
  do while ((i.le.lenght_list).and.(list(i).ne.elem))
    i=i+1
  enddo
  if (i.le.lenght_list)  then 
    if (list(i).eq.elem) then
      psi_list_search=i
    else
      psi_list_search=-1
    endif
  else
    psi_list_search=-1
  endif
end function psi_list_search

