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
subroutine psb_dslu_bld(a,desc_a,p,info)
  use psb_base_mod
  use psb_prec_type

  implicit none 

  type(psb_dspmat_type), intent(inout)      :: a
  type(psb_desc_type), intent(in)        :: desc_a
  type(psb_dbaseprc_type), intent(inout) :: p
  integer, intent(out)                   :: info


  type(psb_dspmat_type)    :: blck, atmp
  character(len=5)         :: fmt
  character                :: upd='F'
  integer                  :: i,j,nza,nzb,nzt,ictxt,me,np,err_act
  logical, parameter :: debug=.false.
  character(len=20)   :: name, ch_err

  interface psb_asmatbld
    Subroutine psb_dasmatbld(ptype,novr,a,blk,desc_data,upd,desc_p,info,outfmt)
      use psb_base_mod
      use psb_prec_type
      integer, intent(in)                  :: ptype,novr
      Type(psb_dspmat_type), Intent(in)    ::  a
      Type(psb_dspmat_type), Intent(inout) ::  blk
      Type(psb_desc_type), Intent(inout)   :: desc_p
      Type(psb_desc_type), Intent(in)      :: desc_data 
      Character, Intent(in)                :: upd
      integer, intent(out)                 :: info
      character(len=5), optional           :: outfmt
    end Subroutine psb_dasmatbld
  end interface

  if(psb_get_errstatus().ne.0) return 
  info=0
  name='psb_slu_bld'
  call psb_erractionsave(err_act)

  ictxt = desc_a%matrix_data(psb_ctxt_)

  call psb_info(ictxt, me, np)

  fmt = 'COO'
  call psb_nullify_sp(blck)    
  call psb_nullify_sp(atmp)    

  atmp%fida='COO'
  if (Debug) then 
    write(0,*) me, 'SPLUBLD: Calling  csdp'
    call psb_barrier(ictxt)
  endif

  call psb_csdp(a,atmp,info)
  if(info /= 0) then
    info=4010
    ch_err='psb_csdp'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  nza = atmp%infoa(psb_nnz_)
  if (Debug) then 
    write(0,*) me, 'SPLUBLD: Done csdp',info,nza,atmp%m,atmp%k
    call psb_barrier(ictxt)
  endif
  call psb_asmatbld(p%iprcparm(p_type_),p%iprcparm(n_ovr_),a,&
       & blck,desc_a,upd,p%desc_data,info,outfmt=fmt)  
  if(info /= 0) then
    info=4010
    ch_err='psb_asmatbld'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  nzb = blck%infoa(psb_nnz_)
  if (Debug) then 
    write(0,*) me, 'SPLUBLD: Done asmatbld',info,nzb,blck%fida
    call psb_barrier(ictxt)
  endif
  if (nzb > 0 ) then 
    if (size(atmp%aspk)<nza+nzb) then 
      call psb_sp_reall(atmp,nza+nzb,info)
      if(info /= 0) then
        info=4010
        ch_err='psb_sp_reall'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    endif
    if (Debug) then 
      write(0,*) me, 'SPLUBLD: Done realloc',info,nza+nzb,atmp%fida
      call psb_barrier(ictxt)
    endif

    do j=1,nzb
      atmp%aspk(nza+j) = blck%aspk(j)
      atmp%ia1(nza+j)  = blck%ia1(j)
      atmp%ia2(nza+j)  = blck%ia2(j)
    end do
    atmp%infoa(psb_nnz_) = nza+nzb
    atmp%m = atmp%m + blck%m
    atmp%k = max(a%k,blck%k)
  else
    atmp%infoa(psb_nnz_) = nza
    atmp%m = a%m 
    atmp%k = a%k
  endif

  i=0
  do j=1, atmp%infoa(psb_nnz_) 
    if (atmp%ia2(j) <= atmp%m) then 
      i = i + 1
      atmp%aspk(i) = atmp%aspk(j)
      atmp%ia1(i) = atmp%ia1(j)
      atmp%ia2(i) = atmp%ia2(j)
    endif
  enddo
  atmp%infoa(psb_nnz_) = i


  call psb_ipcoo2csr(atmp,info)
  if(info /= 0) then
    info=4010
    ch_err='psb_ipcoo2csr'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  nzt = psb_sp_get_nnzeros(atmp)

  if (Debug) then 
    write(0,*) me,'Calling psb_slu_factor ',nzt,atmp%m,&
         & atmp%k,p%desc_data%matrix_data(psb_n_row_)
    call psb_barrier(ictxt)
  endif

  call psb_dslu_factor(atmp%m,nzt,&
       & atmp%aspk,atmp%ia2,atmp%ia1,p%iprcparm(slu_ptr_),info)
  if(info /= 0) then
    ch_err='psb_slu_fact'
    call psb_errpush(4110,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
    goto 9999
  end if

  if (Debug) then 
    write(0,*) me, 'SPLUBLD: Done slu_Factor',info,p%iprcparm(slu_ptr_)
    call psb_barrier(ictxt)
  endif

  call psb_sp_free(blck,info)
  call psb_sp_free(atmp,info)
  if(info /= 0) then
    info=4010
    ch_err='psb_sp_free'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error()
    return
  end if
  return

end subroutine psb_dslu_bld

