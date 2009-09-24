!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
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
subroutine psb_sbjac_bld(a,desc_a,p,upd,info)
  use psb_base_mod
  use psb_prec_mod, psb_protect_name => psb_sbjac_bld
  implicit none
  !                                                                               
  !     .. Scalar Arguments ..                                                    
  integer, intent(out)                      :: info
  !     .. array Arguments ..                                                     
  type(psb_s_sparse_mat), intent(in), target :: a
  type(psb_sprec_type), intent(inout)    :: p
  type(psb_desc_type), intent(in)           :: desc_a
  character, intent(in)                     :: upd

  !     .. Local Scalars ..                                                       
  integer  ::    i, m
  integer  ::    int_err(5)
  character ::        trans, unitd
  type(psb_s_csr_sparse_mat), allocatable  :: lf, uf
  real(psb_dpk_) :: t1,t2,t3,t4,t5,t6, t7, t8
  integer   nztota,  err_act, n_row, nrow_a,n_col, nhalo
  integer :: ictxt,np,me
  character(len=20)      :: name, ch_err


  if(psb_get_errstatus() /= 0) return 
  info=0
  name='psb_sbjac_bld'
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)

  m = a%get_nrows()
  if (m < 0) then
    info = 10
    int_err(1) = 1
    int_err(2) = m
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif
  trans = 'N'
  unitd = 'U'

  call psb_cdcpy(desc_a,p%desc_data,info)
  if(info /= 0) then
    info=4010
    ch_err='psb_cdcpy'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if


  select case(p%iprcparm(psb_f_type_))

  case(psb_f_ilu_n_) 

    if (allocated(p%av)) then 
      if (size(p%av) < psb_bp_ilu_avsz) then 
        do i=1,size(p%av) 
          call p%av(i)%free()
        enddo
        deallocate(p%av,stat=info)
      endif
    end if
    if (.not.allocated(p%av)) then 
      allocate(p%av(psb_max_avsz),stat=info)
      if (info /= 0) then
        call psb_errpush(4000,name)
        goto 9999
      end if
    endif

    nrow_a = psb_cd_get_local_rows(desc_a)
    nztota = a%get_nzeros()

    n_col  = psb_cd_get_local_cols(desc_a)
    nhalo  = n_col-nrow_a
    n_row  = p%desc_data%matrix_data(psb_n_row_)
    
    allocate(lf,uf,stat=info)
    if (info == 0) call lf%allocate(n_row,n_row,nztota)
    if (info == 0) call uf%allocate(n_row,n_row,nztota)
    
    if(info/=0) then
      info=4010
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if (allocated(p%d)) then 
      if (size(p%d) < n_row) then 
        deallocate(p%d)
      endif
    endif
    if (.not.allocated(p%d)) then 
      allocate(p%d(n_row),stat=info)
      if (info /= 0) then 
        call psb_errpush(4010,name,a_err='Allocate')
        goto 9999      
      end if

    endif
    t3 = psb_wtime()
    ! This is where we have no renumbering, thus no need 
    call psb_ilu_fct(a,lf,uf,p%d,info)

    if(info==0) then
      call p%av(psb_l_pr_)%mv_from(lf)
      call p%av(psb_u_pr_)%mv_from(uf)
      call p%av(psb_l_pr_)%set_asb()
      call p%av(psb_u_pr_)%set_asb()
      call p%av(psb_l_pr_)%trim()
      call p%av(psb_u_pr_)%trim()
    else
      info=4010
      ch_err='psb_ilu_fct'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case(psb_f_none_) 
    info=4010
    ch_err='Inconsistent prec  psb_f_none_'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999

  case default
    info=4010
    ch_err='Unknown psb_f_type_'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end select


  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return


end subroutine psb_sbjac_bld


