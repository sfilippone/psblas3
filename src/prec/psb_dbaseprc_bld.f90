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
subroutine psb_dbaseprc_bld(a,desc_a,p,info,upd)

  use psb_serial_mod
  Use psb_spmat_type
  use psb_descriptor_type
  use psb_prec_type
  use psb_tools_mod
  use psb_comm_mod
  use psb_const_mod
  use psb_psblas_mod
  use psb_error_mod
  Implicit None

  type(psb_dspmat_type), target           :: a
  type(psb_desc_type), intent(in)         :: desc_a
  type(psb_dbaseprc_type),intent(inout)   :: p
  integer, intent(out)                    :: info
  character, intent(in), optional         :: upd

  interface psb_ilu_bld
    subroutine psb_dilu_bld(a,desc_data,p,upd,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      integer, intent(out) :: info
      type(psb_dspmat_type), intent(in), target :: a
      type(psb_desc_type),intent(in)            :: desc_data
      type(psb_dbaseprc_type), intent(inout)    :: p
      character, intent(in)                     :: upd
    end subroutine psb_dilu_bld
  end interface

  interface psb_slu_bld
    subroutine psb_dslu_bld(a,desc_a,p,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      use psb_const_mod
      implicit none 

      type(psb_dspmat_type), intent(in)      :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(psb_dbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine psb_dslu_bld
  end interface

  interface psb_umf_bld
    subroutine psb_dumf_bld(a,desc_a,p,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      use psb_const_mod
      implicit none 

      type(psb_dspmat_type), intent(in)      :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(psb_dbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine psb_dumf_bld
  end interface

  ! Local scalars
  Integer      :: err, nnzero, n_row, n_col,I,j,k,icontxt,&
       & me,mycol,nprow,npcol,mglob,lw, mtype, nrg, nzg, err_act
  real(kind(1.d0))         :: temp, real_err(5)
  real(kind(1.d0)),pointer :: gd(:), work(:)
  integer      :: int_err(5)
  character    :: iupd

  logical, parameter :: debug=.false.   
  integer,parameter  :: iroot=0,iout=60,ilout=40
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  err=0
  call psb_erractionsave(err_act)
  name = 'psb_baseprc_bld'

  if (debug) write(0,*) 'Entering baseprc_bld'
  info = 0
  int_err(1) = 0
  icontxt = desc_a%matrix_data(psb_ctxt_)
  n_row   = desc_a%matrix_data(psb_n_row_)
  n_col   = desc_a%matrix_data(psb_n_col_)
  mglob   = desc_a%matrix_data(psb_m_)
  if (debug) write(0,*) 'Preconditioner Blacs_gridinfo'
  call blacs_gridinfo(icontxt, nprow, npcol, me, mycol)

  if (present(upd)) then 
    if (debug) write(0,*) 'UPD ', upd
    if ((UPD.eq.'F').or.(UPD.eq.'T')) then
      IUPD=UPD
    else
      IUPD='F'
    endif
  else
    IUPD='F'
  endif

  !
  ! Should add check to ensure all procs have the same... 
  !
  ! ALso should define symbolic names for the preconditioners. 
  !

  call psb_check_def(p%iprcparm(p_type_),'base_prec',&
       &  diagsc_,is_legal_base_prec)

  allocate(p%desc_data,stat=info)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if

  call psb_nullify_desc(p%desc_data)

  select case(p%iprcparm(p_type_)) 
  case (noprec_)
    ! Do nothing. 


  case (diagsc_)

    if (debug) write(0,*) 'Precond: Diagonal scaling'
    ! diagonal scaling

    call psb_realloc(n_col,p%d,info)
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psb_realloc')
      goto 9999
    end if

    call psb_csrws(p%d,a,info,trans='N')
    if(info /= 0) then
      info=4010
      ch_err='psb_csrws'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if (debug) write(ilout+me,*) 'VDIAG ',n_row
    do i=1,n_row
      if (p%d(i).eq.0.0d0) then
        p%d(i)=1.d0
      else
        p%d(i) =  1.d0/p%d(i)
      endif

      if (debug) write(ilout+me,*) i,desc_a%loc_to_glob(i), p%d(i)
      if (p%d(i).lt.0.d0) then
        write(0,*) me,'Negative RWS? ',i,p%d(i)
      endif
    end do
    if (a%pl(1) /= 0) then
      allocate(work(n_row),stat=info)
      if (info /= 0) then
        info=4000
        call psb_errpush(info,name)
        goto 9999
      end if
      call  psb_gelp('n',a%pl,p%d,desc_a,info)
      if(info /= 0) then
        info=4010
        ch_err='psb_dgelp'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

      deallocate(work)
    endif

    if (debug) then
      allocate(gd(mglob),stat=info)       
      if (info /= 0) then 
        call psb_errpush(4010,name,a_err='Allocate')
        goto 9999      
      end if

      call  psb_gather(gd, p%d, desc_a, info, iroot=iroot)
      if(info /= 0) then
        info=4010
        ch_err='psb_dgatherm'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

      if (me.eq.iroot) then
        write(iout+nprow,*) 'VDIAG CHECK ',mglob
        do i=1,mglob
          write(iout+nprow,*) i,gd(i)
        enddo
      endif
      deallocate(gd)
    endif
    if (debug) write(*,*) 'Preconditioner DIAG computed OK'


  case (bja_,asm_)

    call psb_check_def(p%iprcparm(n_ovr_),'overlap',&
         &  0,is_legal_n_ovr)
    call psb_check_def(p%iprcparm(restr_),'restriction',&
         &  psb_halo_,is_legal_restrict)
    call psb_check_def(p%iprcparm(prol_),'prolongator',&
         &  psb_none_,is_legal_prolong)
    call psb_check_def(p%iprcparm(iren_),'renumbering',&
         &  renum_none_,is_legal_renum)
    call psb_check_def(p%iprcparm(f_type_),'fact',&
         &  f_ilu_n_,is_legal_ml_fact)

    if (debug) write(0,*)me, ': Calling PSB_ILU_BLD'


    select case(p%iprcparm(f_type_))

    case(f_ilu_n_,f_ilu_e_) 
      call psb_ilu_bld(a,desc_a,p,iupd,info)
      if(debug) write(0,*)me,': out of psb_ilu_bld'
      if(info /= 0) then
        info=4010
        ch_err='psb_ilu_bld'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

    case(f_slu_)

      if(debug) write(0,*)me,': calling slu_bld'
      call psb_slu_bld(a,desc_a,p,info)
      if(info /= 0) then
        info=4010
        ch_err='slu_bld'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

    case(f_umf_)
      if(debug) write(0,*)me,': calling umf_bld'
      call psb_umf_bld(a,desc_a,p,info)
      if(info /= 0) then
        info=4010
        ch_err='umf_bld'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

    case(f_none_) 
      write(0,*) 'Fact=None in BASEPRC_BLD Bja/ASM??'

    case default
      write(0,*) 'Unknown factor type in baseprc_bld bja/asm: ',&
           &p%iprcparm(f_type_)
    end select

  end select



  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error()
    return
  end if
  return

end subroutine psb_dbaseprc_bld

