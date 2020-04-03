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

module psb_z_util_cbind_mod
   use iso_c_binding
   use psb_base_mod
   use psb_util_mod
   use psb_objhandle_mod
   use psb_base_string_cbind_mod

contains

    function psb_c_zmm_mat_write(ah,matrixtitle,filename) bind(c) result(res)
        implicit none
        real(c_double_complex) :: res

        type(psb_c_zspmat)   :: ah
        character(c_char)      :: matrixtitle(*)
        character(c_char)      :: filename(*)

        type(psb_zspmat_type), pointer :: ap
        character(1024)         :: mtitle
        character(1024)         :: fname
        integer(psb_c_ipk_)     :: info


        res = -1.0
        if (c_associated(ah%item)) then
          call c_f_pointer(ah%item,ap)
        else
          return
        end if

        call stringc2f(matrixtitle,mtitle)
        call stringc2f(filename,fname)

        call mm_mat_write(ap,mtitle,info,filename=fname)

        res = info

    end function psb_c_zmm_mat_write

end module
