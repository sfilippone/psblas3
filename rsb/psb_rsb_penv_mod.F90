!                Parallel Sparse BLAS   GPU plugin 
!      (C) Copyright 2013
!  
!                         Salvatore Filippone
!                         Alessandro Fanfarillo
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
  

module psb_rsb_penv_mod
  use psb_const_mod
  use psb_penv_mod
  !use psi_comm_buffers_mod, only : psb_buffer_queue
  use iso_c_binding

!  interface psb_rsb_init
!    module procedure  psb_rsb_init
!  end interface
#if defined(HAVE_RSB)
  interface 
    function psb_C_rsb_init() &
         & result(res) bind(c,name='rsbInit')
      use iso_c_binding
      integer(c_int)		:: res
    end function psb_C_rsb_init
 end interface
 
  interface 
     function psb_C_rsb_exit() &
         & result(res) bind(c,name='rsbExit')
       use iso_c_binding
       integer(c_int)		:: res
     end function psb_C_rsb_exit
  end interface

#endif

contains
  ! !!!!!!!!!!!!!!!!!!!!!!
  !
  ! Environment handling 
  !
  ! !!!!!!!!!!!!!!!!!!!!!!


  subroutine psb_rsb_init()
    use psb_penv_mod
    use psb_const_mod
    use psb_error_mod
    ! type(psb_ctxt_type), intent(in) :: ctxt
    ! integer, intent(in), optional :: dev

    integer :: info

#if defined (HAVE_RSB)
    info = psb_C_rsb_init()
    if (info/=0) write(*,*) 'error during rsb_init'
#endif
  end subroutine psb_rsb_init

  subroutine psb_rsb_exit()
    use psb_penv_mod
    use psb_const_mod
    use psb_error_mod
    ! type(psb_ctxt_type), intent(in) :: ctxt
    ! integer, intent(in), optional :: dev

    integer :: info

#if defined (HAVE_RSB)
    info = psb_C_rsb_exit()
    if (info/=0) write(*,*) 'error during rsb_exit'
#endif
  end subroutine psb_rsb_exit

end module psb_rsb_penv_mod
