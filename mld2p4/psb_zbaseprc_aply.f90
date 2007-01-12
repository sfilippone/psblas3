!!$ 
!!$ 
!!$                    MD2P4
!!$    Multilevel Domain Decomposition Parallel Preconditioner Package for PSBLAS
!!$                      for 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$                       Daniela di Serafino    Second University of Naples
!!$                       Pasqua D'Ambra         ICAR-CNR                      
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the MD2P4 group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MD2P4 GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
subroutine psb_zbaseprc_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
  !
  !  Compute   Y <-  beta*Y + alpha*K^-1 X 
  !  where K is a a basic preconditioner stored in prec
  ! 
  use psb_base_mod
  use psb_prec_type
  implicit none 

  type(psb_desc_type),intent(in)      :: desc_data
  type(psb_zbaseprc_type), intent(in) :: prec
  complex(kind(0.d0)),intent(inout)      :: x(:), y(:)
  complex(kind(0.d0)),intent(in)         :: alpha,beta
  character(len=1)                    :: trans
  complex(kind(0.d0)),target             :: work(:)
  integer, intent(out)                :: info

  ! Local variables
  integer :: n_row,n_col, int_err(5)
  complex(kind(1.d0)), pointer :: ww(:), aux(:), tx(:),ty(:)
  character     ::diagl, diagu
  integer :: ictxt,np,me,i, isz, nrg, err_act
  real(kind(1.d0)) :: t1, t2, t3, t4, t5, t6, t7, mpi_wtime
  logical,parameter                 :: debug=.false., debugprt=.false.
  external mpi_wtime
  character(len=20)   :: name, ch_err

  interface psb_bjac_aply
     subroutine psb_zbjac_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
       use psb_base_mod
       use psb_prec_type
       type(psb_desc_type), intent(in)       :: desc_data
       type(psb_zbaseprc_type), intent(in)   :: prec
       complex(kind(0.d0)),intent(inout)        :: x(:), y(:)
       complex(kind(0.d0)),intent(in)           :: alpha,beta
       character(len=1)                      :: trans
       complex(kind(0.d0)),target               :: work(:)
       integer, intent(out)                  :: info
     end subroutine psb_zbjac_aply
  end interface
  
  name='psb_zbaseprc_aply'
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

  case(asm_,ras_,ash_,rash_)

    if (prec%iprcparm(n_ovr_)==0) then 
      ! shortcut: this fixes performance for RAS(0) == BJA
      call psb_bjac_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      if(info.ne.0) then
        info=4010
        ch_err='psb_bjacaply'
        goto 9999
      end if

    else
      ! Note: currently trans is unused.
      n_row=prec%desc_data%matrix_data(psb_n_row_)
      n_col=prec%desc_data%matrix_data(psb_n_col_)

      isz=max(n_row,N_COL)
      if ((6*isz) <= size(work)) then 
        ww => work(1:isz)
        tx => work(isz+1:2*isz)
        ty => work(2*isz+1:3*isz)
        aux => work(3*isz+1:)
      else if ((4*isz) <= size(work)) then 
        aux => work(1:)
        allocate(ww(isz),tx(isz),ty(isz),stat=info)
        if (info /= 0) then 
          call psb_errpush(4010,name,a_err='Allocate')
          goto 9999      
        end if
      else if ((3*isz) <= size(work)) then 
        ww => work(1:isz)
        tx => work(isz+1:2*isz)
        ty => work(2*isz+1:3*isz)
        allocate(aux(4*isz),stat=info)
        if (info /= 0) then 
          call psb_errpush(4010,name,a_err='Allocate')
          goto 9999      
        end if

      else 
        allocate(ww(isz),tx(isz),ty(isz),&
             &aux(4*isz),stat=info)
        if (info /= 0) then 
          call psb_errpush(4010,name,a_err='Allocate')
          goto 9999      
        end if

      endif

      if (debugprt) write(0,*)' vdiag: ',prec%d(:)
      if (debug) write(0,*) 'Bi-CGSTAB with Additive Schwarz prec' 

      tx(1:desc_data%matrix_data(psb_n_row_)) = x(1:desc_data%matrix_data(psb_n_row_)) 
      tx(desc_data%matrix_data(psb_n_row_)+1:isz) = zzero

      if (prec%iprcparm(restr_)==psb_halo_) then 
        call psb_halo(tx,prec%desc_data,info,work=aux)
        if(info /=0) then
          info=4010
          ch_err='psb_halo'
          goto 9999
        end if
      else if (prec%iprcparm(restr_) /= psb_none_) then 
        write(0,*) 'Problem in PRC_APLY: Unknown value for restriction ',&
             &prec%iprcparm(restr_)
      end if

      if (prec%iprcparm(iren_)>0) then 
        call zgelp('N',n_row,1,prec%perm,tx,isz,ww,isz,info)
        if(info /=0) then
          info=4010
          ch_err='psb_zgelp'
          goto 9999
        end if
      endif

      call psb_bjac_aply(zone,prec,tx,zzero,ty,prec%desc_data,trans,aux,info)
      if(info.ne.0) then
        info=4010
        ch_err='psb_bjac_aply'
        goto 9999
      end if

      if (prec%iprcparm(iren_)>0) then 
        call zgelp('N',n_row,1,prec%invperm,ty,isz,ww,isz,info)
        if(info /=0) then
          info=4010
          ch_err='psb_zgelp'
          goto 9999
        end if
      endif

      select case (prec%iprcparm(prol_)) 

      case(psb_none_) 
        ! Would work anyway, but since it's supposed to do nothing...
        ! call f90_psovrl(ty,prec%desc_data,update=prec%a_restrict)

      case(psb_sum_,psb_avg_) 
        call psb_ovrl(ty,prec%desc_data,info,&
             & update=prec%iprcparm(prol_),work=aux)
        if(info /=0) then
          info=4010
          ch_err='psb_ovrl'
          goto 9999
        end if

      case default
        write(0,*) 'Problem in PRC_APLY: Unknown value for prolongation ',&
             & prec%iprcparm(prol_)
      end select

      call psb_geaxpby(alpha,ty,beta,y,desc_data,info) 


      if ((6*isz) <= size(work)) then 
      else if ((4*isz) <= size(work)) then 
        deallocate(ww,tx,ty)
      else if ((3*isz) <= size(work)) then 
        deallocate(aux)
      else 
        deallocate(ww,aux,tx,ty)
      endif
    end if
  case default
    write(0,*) 'Invalid PRE%PREC ',prec%iprcparm(p_type_),':',&
         & min_prec_,noprec_,diagsc_,bja_,&
         & asm_,ras_,ash_,rash_
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name,i_err=int_err,a_err=ch_err)
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error()
    return
  end if
  return

end subroutine psb_zbaseprc_aply

