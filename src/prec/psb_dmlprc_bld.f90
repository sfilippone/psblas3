!!$ 
!!$ 
!!$              MPcube: Multilevel Parallel Preconditioners Package 
!!$                      for 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$                       Daniela Di Serafino    II University of Naples
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
!!$    3. The name of the MPCUBE group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MPCUBE GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
subroutine psb_dmlprc_bld(a,desc_a,p,info)

  use psb_serial_mod
  use psb_tools_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_const_mod
  use psb_error_mod
  implicit none 

  type(psb_dspmat_type), intent(in), target :: a
  type(psb_desc_type), intent(in)           :: desc_a
  type(psb_dbase_prec), intent(inout)       :: p
  integer, intent(out)                      :: info

  integer :: i, nrg, nzg, err_act,k
  character(len=20) :: name, ch_err
  
  interface psb_ilu_fct
     subroutine psb_dilu_fct(a,l,u,d,info,blck)
       use psb_spmat_type
       integer, intent(out)                ::     info
       type(psb_dspmat_type),intent(in)    :: a
       type(psb_dspmat_type),intent(inout) :: l,u
       type(psb_dspmat_type),intent(in), optional, target :: blck
       real(kind(1.d0)), intent(inout)     ::  d(:)
     end subroutine psb_dilu_fct
  end interface

  interface psb_genaggrmap
     subroutine psb_dgenaggrmap(aggr_type,a,desc_a,nlaggr,ilaggr,info)
       use psb_spmat_type
       use psb_descriptor_type
       implicit none
       integer, intent(in)               :: aggr_type
       type(psb_dspmat_type), intent(in) :: a
       type(psb_desc_type), intent(in)   :: desc_a
       integer, pointer                  :: ilaggr(:),nlaggr(:)
       integer, intent(out)              :: info
     end subroutine psb_dgenaggrmap
  end interface

  interface psb_bldaggrmat
     subroutine psb_dbldaggrmat(a,desc_a,p,info)
       use psb_prec_type
       use psb_descriptor_type
       use psb_spmat_type
       type(psb_dspmat_type), intent(in), target :: a
       type(psb_dbase_prec), intent(inout)        :: p
       type(psb_desc_type), intent(in)           :: desc_a
       integer, intent(out)                      :: info
     end subroutine psb_dbldaggrmat
  end interface

  integer :: icontxt, nprow, npcol, me, mycol
  
  name='psb_mlprec_bld'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  p%aorig => a
  allocate(p%av(smth_avsz),stat=info)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if
  
  do i=1, smth_avsz
     call psb_nullify_sp(p%av(i))
     call psb_spall(0,0,p%av(i),1,info)
     if(info /= 0) then
        info=4010
        ch_err='psb_spall'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if
  end do
  
  ! Currently this is ignored by gen_aggrmap, but it could be 
  ! changed in the future. Need to package nlaggr & mlia in a 
  ! private data structure? 

  call psb_genaggrmap(p%iprcparm(aggr_alg_),a,desc_a,p%nlaggr,p%mlia,info)
  if(info /= 0) then
     info=4010
     ch_err='psb_gen_aggrmap'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  call psb_bldaggrmat(a,desc_a,p,info)
  if(info /= 0) then
     info=4010
     ch_err='psb_bld_aggrmat'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  nrg = p%av(ac_)%m
  call psb_spinfo(psb_nztotreq_,p%av(ac_),nzg,info)
  call psb_ipcoo2csr(p%av(ac_),info)
  if(info /= 0) then
     info=4011
     ch_err='psb_ipcoo2csr'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  allocate(p%d(nrg),stat=info) 
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if

  select case(p%iprcparm(f_type_)) 
  case(f_ilu_n_,f_ilu_e_) 
     call psb_spreall(p%av(l_pr_),nzg,info)
     call psb_spreall(p%av(u_pr_),nzg,info)
     call psb_ilu_fct(p%av(ac_),p%av(l_pr_),p%av(u_pr_),p%d,info)
     if(info /= 0) then
        info=4011
        ch_err='psb_ilu_fct'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

  case(f_slu_) 
!!$    call psb_spall(0,0,p%av(l_pr_),1,info)
!!$    call psb_spall(0,0,p%av(u_pr_),1,info)
    call psb_ipcsr2coo(p%av(ac_),info)
     if(info /= 0) then
        info=4011
        ch_err='psb_ipcsr2coo'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if
    k=0
    do i=1,p%av(ac_)%infoa(psb_nnz_)
      if (p%av(ac_)%ia2(i) <= p%av(ac_)%m) then 
        k = k + 1
        p%av(ac_)%aspk(k) = p%av(ac_)%aspk(i)
        p%av(ac_)%ia1(k) = p%av(ac_)%ia1(i)
        p%av(ac_)%ia2(k) = p%av(ac_)%ia2(i)
      end if
    end do
    p%av(ac_)%infoa(psb_nnz_) = k
    call psb_ipcoo2csr(p%av(ac_),info)
    call psb_spinfo(psb_nztotreq_,p%av(ac_),nzg,info)
    call psb_slu_factor(nrg,nzg,&
         & p%av(ac_)%aspk,p%av(ac_)%ia2,p%av(ac_)%ia1,p%iprcparm(slu_ptr_),info)
     if(info /= 0) then
        info=4011
        ch_err='psb_slu_factor'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

  case(f_umf_) 
!!$    call psb_spall(0,0,p%av(l_pr_),1,info)
!!$    call psb_spall(0,0,p%av(u_pr_),1,info)
    call psb_ipcsr2coo(p%av(ac_),info)
     if(info /= 0) then
        info=4011
        ch_err='psb_ipcsr2coo'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if
    k=0
    do i=1,p%av(ac_)%infoa(psb_nnz_)
      if (p%av(ac_)%ia2(i) <= p%av(ac_)%m) then 
        k = k + 1
        p%av(ac_)%aspk(k) = p%av(ac_)%aspk(i)
        p%av(ac_)%ia1(k) = p%av(ac_)%ia1(i)
        p%av(ac_)%ia2(k) = p%av(ac_)%ia2(i)
      end if
    end do
    p%av(ac_)%infoa(psb_nnz_) = k
    call psb_ipcoo2csc(p%av(ac_),info)
    call psb_spinfo(psb_nztotreq_,p%av(ac_),nzg,info)
    call psb_umf_factor(nrg,nzg,&
         & p%av(ac_)%aspk,p%av(ac_)%ia1,p%av(ac_)%ia2,&
         & p%iprcparm(umf_symptr_),p%iprcparm(umf_numptr_),info)
     if(info /= 0) then
        info=4011
        ch_err='psb_umf_factor'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

  case default
     write(0,*) 'Invalid fact type for multi level',(p%iprcparm(f_type_)) 
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error()
     return
  end if
  Return

end subroutine psb_dmlprc_bld
