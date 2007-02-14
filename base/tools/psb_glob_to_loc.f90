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
! File: psb_glob_to_loc.f90
!
! Subroutine: psb_glob_to_loc2
!    Performs global to local indexes translation
! 
! Parameters: 
!    x        - integer, dimension(:).    Array containing the indices to be translated.
!    y        - integer, dimension(:).    Array containing the translated indices.
!    desc_a   - type(<psb_desc_type>).    The communication descriptor.        
!    info     - integer.                  Eventually returns an error code.
!    iact     - integer(optional).        A character defining the behaviour of this subroutine when is found an index not belonging to the calling process
!
subroutine psb_glob_to_loc2(x,y,desc_a,info,iact,owned)

  use psb_descriptor_type
  use psb_const_mod
  use psb_error_mod
  use psb_string_mod
  use psb_penv_mod
  use psi_mod
  implicit none

  !...parameters....
  type(psb_desc_type), intent(in) ::  desc_a
  integer, intent(in)                ::  x(:)  
  integer, intent(out)               ::  y(:), info
  character, intent(in), optional    ::  iact
  logical, intent(in),  optional     :: owned

  !....locals....
  integer                            ::  n, i, tmp
  character                          ::  act
  integer                            ::  int_err(5), err_act
  real(kind(1.d0))                   ::  real_val
  logical                            :: owned_
  integer, parameter                 ::  zero=0
  character(len=20)   :: name

  if(psb_get_errstatus() /= 0) return 
  info=0
  name = 'glob_to_loc'
  call psb_erractionsave(err_act)

  if (present(iact)) then
    act=iact
  else
    act='I'
  endif
  act = toupper(act)
  if (present(owned)) then 
    owned_=owned
  else
    owned_=.false.
  end if
    
  int_err=0
  real_val = 0.d0

  n = size(x)
  call psi_idx_cnv(n,x,y,desc_a,info,owned=owned_)

  select case(act)
  case('E','I')
    call psb_erractionrestore(err_act)
    return
  case('W')
    if ((info /= 0).or.(count(y(1:n)<0) >0)) then
      write(0,'("Error ",i5," in subroutine glob_to_loc")') info
    end if
  case('A')
    if ((info /= 0).or.(count(y(1:n)<0) >0)) then
      call psb_errpush(info,name)
      goto 9999
    end if
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error()
  end if
  return


end subroutine psb_glob_to_loc2


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
! Subroutine: psb_glob_to_loc
!    Performs global to local indexes translation
! 
! Parameters: 
!    x        - integer, dimension(:).    Array containing the indices to be translated.
!    desc_a   - type(<psb_desc_type>).    The communication descriptor.        
!    info     - integer.                  Eventually returns an error code.
!    iact     - integer(optional).        A character defining the behaviour of this subroutine when is found an index not belonging to the calling process
!
subroutine psb_glob_to_loc(x,desc_a,info,iact,owned)

  use psb_penv_mod
  use psb_descriptor_type
  use psb_const_mod
  use psb_error_mod
  use psb_string_mod
  use psi_mod
  implicit none

  !...parameters....
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(inout)           :: x(:)  
  integer, intent(out)             :: info
  logical, intent(in),  optional     :: owned
  character, intent(in), optional  :: iact

  !....locals....
  integer                          :: n, i, tmp, nk, key, idx, ih, nh, lb, ub, lm
  character                        :: act
  integer                          :: int_err(5), err_act, dectype
  real(kind(1.d0))                 :: real_val, t0, t1,t2
  logical                            :: owned_
  integer, parameter               :: zero=0
  character(len=20)   :: name
  integer             :: ictxt, iam, np

  if(psb_get_errstatus() /= 0) return 
  info=0
  name = 'glob_to_loc'
  ictxt = desc_a%matrix_data(psb_ctxt_)
  call psb_info(ictxt,iam,np)
  call psb_erractionsave(err_act)

  dectype  = desc_a%matrix_data(psb_dec_type_)
  if (present(iact)) then
    act=iact
  else
    act='I'
  endif
  if (present(owned)) then 
    owned_=owned
  else
    owned_=.false.
  end if

  act = toupper(act)

  n = size(x)
  call psi_idx_cnv(n,x,desc_a,info,owned=owned_)

  select case(act)
  case('E','I')
    call psb_erractionrestore(err_act)
    return
  case('W')
    if ((info /= 0).or.(count(x(1:n)<0) >0)) then
      write(0,'("Error ",i5," in subroutine glob_to_loc")') info
    end if
  case('A')
    if ((info /= 0).or.(count(x(1:n)<0) >0)) then
      write(0,*) count(x(1:n)<0)
      call psb_errpush(info,name)
      goto 9999
    end if
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error()
  end if
  return

contains 

  subroutine inlbsrch(ipos,key,n,v)
    implicit none
    integer ipos, key, n
    integer v(n)

    integer lb, ub, m


    lb = 1 
    ub = n
    ipos = -1 

    do 
      if (lb > ub) return
      m = (lb+ub)/2
      if (key.eq.v(m))  then
        ipos = m 
        return
      else if (key.lt.v(m))  then
        ub = m-1
      else 
        lb = m + 1
      end if
    enddo
    return
  end subroutine inlbsrch

  subroutine inner_cnv(n,x,hashsize,hashmask,hashv,glb_lc)
    integer :: n, hashsize,hashmask,x(:), hashv(0:),glb_lc(:,:)
    integer :: i, ih, key, idx,nh,tmp,lb,ub,lm
    do i=1, n
      key = x(i) 
      ih  = iand(key,hashmask)
      idx = hashv(ih)
      nh  = hashv(ih+1) - hashv(ih) 
      if (nh > 0) then 
        tmp = -1 
        lb = idx
        ub = idx+nh-1
        do 
          if (lb>ub) exit
          lm = (lb+ub)/2
          if (key==glb_lc(lm,1)) then 
            tmp = lm
            exit
          else if (key<glb_lc(lm,1)) then 
            ub = lm - 1
          else
            lb = lm + 1
          end if
        end do
      else 
        tmp = -1
      end if
      if (tmp > 0) then 
        x(i) = glb_lc(tmp,2)
      else         
        x(i) = tmp 
      end if
    end do
  end subroutine inner_cnv

end subroutine psb_glob_to_loc

