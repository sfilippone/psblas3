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
subroutine psb_dbjac_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
  !
  !  Compute   Y <-  beta*Y + alpha*K^-1 X 
  !  where K is a a Block Jacobi  preconditioner stored in prec
  !  Note that desc_data may or may not be the same as prec%desc_data,
  !  but since both are INTENT(IN) this should be legal. 
  ! 

  use psb_base_mod
  use psb_prec_type
  implicit none 

  type(psb_desc_type), intent(in)       :: desc_data
  type(psb_dprec_type), intent(in)   :: prec
  real(kind(0.d0)),intent(inout)        :: x(:), y(:)
  real(kind(0.d0)),intent(in)           :: alpha,beta
  character(len=1)                      :: trans
  real(kind(0.d0)),target               :: work(:)
  integer, intent(out)                  :: info

  ! Local variables
  integer :: n_row,n_col
  real(kind(1.d0)), pointer :: ww(:), aux(:), tx(:),ty(:),tb(:)
  character     ::diagl, diagu
  integer :: ictxt,np,me,i, nrg, err_act, int_err(5)
  real(kind(1.d0)) :: t1, t2, t3, t4, t5, t6, t7
  logical,parameter                 :: debug=.false., debugprt=.false.
  character(len=20)   :: name, ch_err

  name='psb_bjac_aply'
  info = 0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_data)
  call psb_info(ictxt, me, np)

  diagl='U'
  diagu='U'

  select case(trans)
  case('N','n')
  case('T','t','C','c')
  case default
    call psb_errpush(40,name)
    goto 9999
  end select


  n_row=desc_data%matrix_data(psb_n_row_)
  n_col=desc_data%matrix_data(psb_n_col_)

  if (n_col <= size(work)) then 
    ww => work(1:n_col)
    if ((4*n_col+n_col) <= size(work)) then 
      aux => work(n_col+1:)
    else
      allocate(aux(4*n_col),stat=info)
      if (info /= 0) then 
        call psb_errpush(4010,name,a_err='Allocate')
        goto 9999      
      end if

    endif
  else
    allocate(ww(n_col),aux(4*n_col),stat=info)
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if
  endif
    

  if (prec%iprcparm(jac_sweeps_) == 1) then 


    select case(prec%iprcparm(f_type_))
    case(f_ilu_n_,f_ilu_e_) 

      select case(trans)
      case('N','n')

        call psb_spsm(done,prec%av(l_pr_),x,dzero,ww,desc_data,info,&
             & trans='N',unit=diagl,choice=psb_none_,work=aux)
        if(info /=0) goto 9999
        ww(1:n_row) = ww(1:n_row)*prec%d(1:n_row)
        call psb_spsm(alpha,prec%av(u_pr_),ww,beta,y,desc_data,info,&
             & trans='N',unit=diagu,choice=psb_none_, work=aux)
        if(info /=0) goto 9999

      case('T','t','C','c')
        call psb_spsm(done,prec%av(u_pr_),x,dzero,ww,desc_data,info,&
             & trans=trans,unit=diagu,choice=psb_none_, work=aux)
        if(info /=0) goto 9999
        ww(1:n_row) = ww(1:n_row)*prec%d(1:n_row)
        call psb_spsm(alpha,prec%av(l_pr_),ww,beta,y,desc_data,info,&
             & trans=trans,unit=diagl,choice=psb_none_,work=aux)
        if(info /=0) goto 9999

      end select


    case default
      write(0,*) 'Unknown factorization type in bjac_aply',prec%iprcparm(f_type_)
    end select
    if (debugprt) write(0,*)' Y: ',y(:)

  else if (prec%iprcparm(jac_sweeps_) > 1) then 

    ! Note: we have to add TRANS to this one !!!!!!!!! 

    if (size(prec%av) < ap_nd_) then 
      info = 4011
      goto 9999
    endif

    allocate(tx(n_col),ty(n_col),stat=info)
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    tx = dzero
    ty = dzero
    select case(prec%iprcparm(f_type_)) 
    case(f_ilu_n_,f_ilu_e_) 
      do i=1, prec%iprcparm(jac_sweeps_) 
        !   X(k+1) = M^-1*(b-N*X(k))
        ty(1:n_row) = x(1:n_row)
        call psb_spmm(-done,prec%av(ap_nd_),tx,done,ty,&
             &   prec%desc_data,info,work=aux)
        if(info /=0) goto 9999
        call psb_spsm(done,prec%av(l_pr_),ty,dzero,ww,&
             & prec%desc_data,info,&
             & trans='N',unit='U',choice=psb_none_,work=aux)
        if(info /=0) goto 9999
        ww(1:n_row) = ww(1:n_row)*prec%d(1:n_row)
        call psb_spsm(done,prec%av(u_pr_),ww,dzero,tx,&
             & prec%desc_data,info,&
             & trans='N',unit='U',choice=psb_none_,work=aux)
        if(info /=0) goto 9999
      end do

    end select
    
    call psb_geaxpby(alpha,tx,beta,y,desc_data,info)
    

    deallocate(tx,ty)


  else

    goto 9999

  endif

  if (n_col <= size(work)) then 
    if ((4*n_col+n_col) <= size(work)) then 
    else
      deallocate(aux)
    endif
  else
    deallocate(ww,aux)
  endif


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

end subroutine psb_dbjac_aply

