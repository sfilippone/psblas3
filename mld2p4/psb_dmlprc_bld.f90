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
subroutine psb_dmlprc_bld(a,desc_a,p,info)

  use psb_base_mod
  use psb_prec_type
  implicit none 

  type(psb_dspmat_type), intent(in), target :: a
  type(psb_desc_type), intent(in), target   :: desc_a
  type(psb_dbaseprc_type), intent(inout),target    :: p
  integer, intent(out)                      :: info

  type(psb_desc_type)                       :: desc_ac

  integer :: i, nrg, nzg, err_act,k
  character(len=20) :: name, ch_err
  logical, parameter :: debug=.false.
  type(psb_dspmat_type)                     :: ac

  interface psb_baseprc_bld
    subroutine psb_dbaseprc_bld(a,desc_a,p,info,upd)
      use psb_base_mod
      use psb_prec_type
      type(psb_dspmat_type), target              :: a
      type(psb_desc_type), intent(in)            :: desc_a
      type(psb_dbaseprc_type),intent(inout)      :: p
      integer, intent(out)                       :: info
      character, intent(in), optional            :: upd
    end subroutine psb_dbaseprc_bld
  end interface

  interface psb_genaggrmap
    subroutine psb_dgenaggrmap(aggr_type,a,desc_a,nlaggr,ilaggr,info)
      use psb_base_mod
      use psb_prec_type
      implicit none
      integer, intent(in)               :: aggr_type
      type(psb_dspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer, allocatable              :: ilaggr(:),nlaggr(:)
      integer, intent(out)              :: info
    end subroutine psb_dgenaggrmap
  end interface

  interface psb_bldaggrmat
    subroutine psb_dbldaggrmat(a,desc_a,ac,desc_ac,p,info)
      use psb_base_mod
      use psb_prec_type
      type(psb_dspmat_type), intent(in), target :: a
      type(psb_desc_type), intent(in)           :: desc_a
      type(psb_dspmat_type), intent(out),target :: ac
      type(psb_desc_type), intent(inout)        :: desc_ac
      type(psb_dbaseprc_type), intent(inout), target :: p
      integer, intent(out)                      :: info
    end subroutine psb_dbldaggrmat
  end interface

  integer :: ictxt, np, me

  name='psb_mlprec_bld'
  if (psb_get_errstatus().ne.0) return 
  info = 0
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt,me,np)
  call psb_erractionsave(err_act)
  call psb_nullify_sp(ac)


  if (.not.allocated(p%iprcparm)) then 
    info = 2222
    call psb_errpush(info,name)
    goto 9999
  endif
  call psb_check_def(p%iprcparm(ml_type_),'Multilevel type',&
       &   mult_ml_prec_,is_legal_ml_type)
  call psb_check_def(p%iprcparm(aggr_alg_),'aggregation',&
       &   loc_aggr_,is_legal_ml_aggr_kind)
  call psb_check_def(p%iprcparm(smth_kind_),'Smoother kind',&
       &   smth_omg_,is_legal_ml_smth_kind)
  call psb_check_def(p%iprcparm(coarse_mat_),'Coarse matrix',&
       &   mat_distr_,is_legal_ml_coarse_mat)
  call psb_check_def(p%iprcparm(smth_pos_),'smooth_pos',&
       &   pre_smooth_,is_legal_ml_smooth_pos)


!!$  nullify(p%desc_data)
  select case(p%iprcparm(f_type_))
  case(f_ilu_n_)      
    call psb_check_def(p%iprcparm(ilu_fill_in_),'Level',0,is_legal_ml_lev)
  case(f_ilu_e_)                 
    call psb_check_def(p%dprcparm(fact_eps_),'Eps',dzero,is_legal_ml_eps)
  end select
  call psb_check_def(p%dprcparm(smooth_omega_),'omega',dzero,is_legal_omega)
  call psb_check_def(p%iprcparm(jac_sweeps_),'Jacobi sweeps',&
       & 1,is_legal_jac_sweeps)


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

  if (debug) write(0,*) 'Out from genaggrmap',p%nlaggr

  call psb_nullify_desc(desc_ac)
  call psb_bldaggrmat(a,desc_a,ac,desc_ac,p,info)
  if(info /= 0) then
    info=4010
    ch_err='psb_bld_aggrmat'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  if (debug) write(0,*) 'Out from bldaggrmat',desc_ac%matrix_data(:)



  call psb_baseprc_bld(ac,desc_ac,p,info)
  if (debug) write(0,*) 'Out from baseprcbld',info
  if(info /= 0) then
    info=4010
    ch_err='psb_baseprc_bld'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  

  !
  ! We have used a separate ac because:
  ! 1. We want to reuse the same routines psb_ilu_bld etc.
  ! 2. We do NOT want to pass an argument twice to them 
  !    p%av(ac_) and p, as this would violate the Fortran standard
  ! Hence a separate AC and a TRANSFER function at the end. 
  !
  call psb_sp_transfer(ac,p%av(ac_),info)
  p%base_a => p%av(ac_)
  call psb_cdtransfer(desc_ac,p%desc_ac,info)

  if (info /= 0) then 
    info=4010
    ch_err='psb_cdtransfer'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  p%base_desc => p%desc_ac

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
