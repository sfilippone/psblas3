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
subroutine psb_dbaseprc_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
  !
  !  Compute   Y <-  beta*Y + alpha*K^-1 X 
  !  where K is a a basic preconditioner stored in prec
  ! 

  use psb_base_mod
  use psb_prec_type
  implicit none 

  type(psb_desc_type),intent(in)      :: desc_data
  type(psb_dprec_type), intent(in) :: prec
  real(kind(0.d0)),intent(inout)      :: x(:), y(:)
  real(kind(0.d0)),intent(in)         :: alpha,beta
  character(len=1)                    :: trans
  real(kind(0.d0)),target             :: work(:)
  integer, intent(out)                :: info

  ! Local variables
  integer :: n_row,n_col, int_err(5)
  real(kind(1.d0)), pointer :: ww(:), aux(:), tx(:),ty(:)
  character     ::diagl, diagu
  integer :: ictxt,np,me,i, isz, nrg, err_act
  real(kind(1.d0)) :: t1, t2, t3, t4, t5, t6, t7
  logical,parameter                 :: debug=.false., debugprt=.false.
  character(len=20)   :: name, ch_err

  interface psb_bjac_aply
     subroutine psb_dbjac_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
       use psb_base_mod
       use psb_prec_type
       type(psb_desc_type), intent(in)       :: desc_data
       type(psb_dprec_type), intent(in)   :: prec
       real(kind(0.d0)),intent(inout)        :: x(:), y(:)
       real(kind(0.d0)),intent(in)           :: alpha,beta
       character(len=1)                      :: trans
       real(kind(0.d0)),target               :: work(:)
       integer, intent(out)                  :: info
     end subroutine psb_dbjac_aply
  end interface

  name='psb_baseprc_aply'
  info = 0
  call psb_erractionsave(err_act)

  ictxt=desc_data%matrix_data(psb_ctxt_)
  call psb_info(ictxt, me, np)

  diagl='U'
  diagu='U'

  select case(trans)
  case('N','n')
  case('T','t','C','c')
  case default
    info=40
    int_err(1)=6
    ch_err(2:2)=trans
    goto 9999
  end select

  select case(prec%iprcparm(p_type_))

  case(noprec_)

    call psb_geaxpby(alpha,x,beta,y,desc_data,info)

  case(diagsc_)
    
    if (size(work) >= size(x)) then 
      ww => work
    else
      allocate(ww(size(x)),stat=info)
      if (info /= 0) then 
        call psb_errpush(4010,name,a_err='Allocate')
        goto 9999      
      end if
    end if

    n_row=desc_data%matrix_data(psb_n_row_)
    ww(1:n_row) = x(1:n_row)*prec%d(1:n_row)
    call psb_geaxpby(alpha,ww,beta,y,desc_data,info)

    if (size(work) < size(x)) then 
      deallocate(ww,stat=info)
      if (info /= 0) then 
        call psb_errpush(4010,name,a_err='Deallocate')
        goto 9999      
      end if
    end if

  case(bja_)

    call psb_bjac_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
    if(info.ne.0) then
      info=4010
      ch_err='psb_bjac_aply'
      goto 9999
    end if

  case default
    write(0,*) 'Invalid PRE%PREC ',prec%iprcparm(p_type_),':',&
         & min_prec_,noprec_,diagsc_,bja_
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name,i_err=int_err,a_err=ch_err)
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_dbaseprc_aply

