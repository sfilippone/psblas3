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
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine psb_c_apply2_vect(prec,x,y,desc_data,info,trans,work)
  use psb_base_mod
  use psb_c_prec_type, psb_protect_name => psb_c_apply2_vect
  implicit none 
  type(psb_desc_type),intent(in)       :: desc_data
  class(psb_cprec_type), intent(inout) :: prec
  type(psb_c_vect_type),intent(inout)  :: x
  type(psb_c_vect_type),intent(inout)  :: y
  integer(psb_ipk_), intent(out)                 :: info
  character(len=1), optional           :: trans
  complex(psb_spk_),intent(inout), optional, target :: work(:)

  character     :: trans_ 
  complex(psb_spk_), pointer :: work_(:)
  integer(psb_ipk_) :: ictxt,np,me
  integer(psb_ipk_) :: err_act
  character(len=20)   :: name

  name = 'psb_c_apply2v'
  info = psb_success_
  call psb_erractionsave(err_act)

  ictxt = desc_data%get_context()
  call psb_info(ictxt, me, np)

  if (present(trans)) then 
    trans_=psb_toupper(trans)
  else
    trans_='N'
  end if

  if (present(work)) then 
    work_ => work
  else
    allocate(work_(4*desc_data%get_local_cols()),stat=info)
    if (info /= psb_success_) then 
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999      
    end if

  end if

  if (.not.allocated(prec%prec)) then 
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner")
    goto 9999
  end if

  call prec%prec%apply(cone,x,czero,y,desc_data,info,&
       & trans=trans_,work=work_)

  if (present(work)) then 
  else
    deallocate(work_,stat=info)
    if (info /= psb_success_) then 
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='DeAllocate')
      goto 9999      
    end if
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_c_apply2_vect

subroutine psb_c_apply1_vect(prec,x,desc_data,info,trans,work)
  use psb_base_mod
  use psb_c_prec_type, psb_protect_name => psb_c_apply1_vect
  implicit none 
  type(psb_desc_type),intent(in)       :: desc_data
  class(psb_cprec_type), intent(inout) :: prec
  type(psb_c_vect_type),intent(inout)  :: x
  integer(psb_ipk_), intent(out)                 :: info
  character(len=1), optional           :: trans
  complex(psb_spk_),intent(inout), optional, target :: work(:)

  type(psb_c_vect_type)       :: ww
  character     :: trans_ 
  complex(psb_spk_), pointer :: work_(:)
  integer(psb_ipk_) :: ictxt,np,me
  integer(psb_ipk_) :: err_act
  character(len=20)   :: name

  name = 'psb_c_apply1v'
  info = psb_success_
  call psb_erractionsave(err_act)

  ictxt = desc_data%get_context()
  call psb_info(ictxt, me, np)

  if (present(trans)) then 
    trans_=psb_toupper(trans)
  else
    trans_='N'
  end if

  if (present(work)) then 
    work_ => work
  else
    allocate(work_(4*desc_data%get_local_cols()),stat=info)
    if (info /= psb_success_) then 
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999      
    end if

  end if

  if (.not.allocated(prec%prec)) then 
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner")
    goto 9999
  end if

  call psb_geasb(ww,desc_data,info,mold=x%v,scratch=.true.)
  if (info == 0) call prec%prec%apply(cone,x,czero,ww,desc_data,info,&
       & trans=trans_,work=work_)
  if (info == 0) call psb_geaxpby(cone,ww,czero,x,desc_data,info)
  call psb_gefree(ww,desc_data,info)
  if (present(work)) then 
  else
    deallocate(work_,stat=info)
    if (info /= psb_success_) then 
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='DeAllocate')
      goto 9999      
    end if
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_c_apply1_vect

subroutine psb_c_apply2v(prec,x,y,desc_data,info,trans,work)
  use psb_base_mod
  use psb_c_prec_type, psb_protect_name => psb_c_apply2v
  implicit none 
  type(psb_desc_type),intent(in)    :: desc_data
  class(psb_cprec_type), intent(inout) :: prec
  complex(psb_spk_),intent(inout)   :: x(:)
  complex(psb_spk_),intent(inout)   :: y(:)
  integer(psb_ipk_), intent(out)              :: info
  character(len=1), optional        :: trans
  complex(psb_spk_),intent(inout), optional, target :: work(:)

  character     :: trans_ 
  complex(psb_spk_), pointer :: work_(:)
  integer(psb_ipk_) :: ictxt,np,me
  integer(psb_ipk_) :: err_act
  character(len=20)   :: name

  name='psb_c_apply2v'
  info = psb_success_
  call psb_erractionsave(err_act)

  ictxt = desc_data%get_context()
  call psb_info(ictxt, me, np)

  if (present(trans)) then 
    trans_=trans
  else
    trans_='N'
  end if

  if (present(work)) then 
    work_ => work
  else
    allocate(work_(4*desc_data%get_local_cols()),stat=info)
    if (info /= psb_success_) then 
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999      
    end if

  end if

  if (.not.allocated(prec%prec)) then 
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner")
    goto 9999
  end if
  call prec%prec%apply(cone,x,czero,y,desc_data,info,trans_,work=work_)
  if (present(work)) then 
  else
    deallocate(work_,stat=info)
    if (info /= psb_success_) then 
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='DeAllocate')
      goto 9999      
    end if
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_c_apply2v

subroutine psb_c_apply1v(prec,x,desc_data,info,trans)
  use psb_base_mod
  use psb_c_prec_type, psb_protect_name => psb_c_apply1v
  implicit none 
  type(psb_desc_type),intent(in)    :: desc_data
  class(psb_cprec_type), intent(inout) :: prec
  complex(psb_spk_),intent(inout)   :: x(:)
  integer(psb_ipk_), intent(out)              :: info
  character(len=1), optional        :: trans

  character     :: trans_
  integer(psb_ipk_) :: ictxt,np,me
  integer(psb_ipk_) :: err_act
  complex(psb_spk_), pointer :: WW(:), w1(:)
  character(len=20)   :: name
  name='psb_c_apply1v'
  info = psb_success_
  call psb_erractionsave(err_act)


  ictxt=desc_data%get_context()
  call psb_info(ictxt, me, np)
  if (present(trans)) then 
    trans_=psb_toupper(trans)
  else
    trans_='N'
  end if

  if (.not.allocated(prec%prec)) then 
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner")
    goto 9999
  end if
  allocate(ww(size(x)),w1(size(x)),stat=info)
  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='Allocate')
    goto 9999      
  end if
  call prec%prec%apply(cone,x,czero,ww,desc_data,info,&
       & trans_,work=w1)
  if(info /= psb_success_) goto 9999
  x(:) = ww(:)
  deallocate(ww,W1,stat=info)
  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='DeAllocate')
    goto 9999      
  end if


  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_error_handler(err_act)
  return

end subroutine psb_c_apply1v

