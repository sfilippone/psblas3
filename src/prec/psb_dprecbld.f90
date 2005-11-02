subroutine psb_dprecbld(a,p,desc_a,info,upd)

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

  integer, intent(out)                       :: info
  type(psb_dspmat_type), target              :: a
  type(psb_dprec_type),intent(inout)         :: p
  type(psb_desc_type), intent(in)            :: desc_a
  character, intent(in), optional            :: upd

  interface psb_cslu
     subroutine psb_dcslu(a,desc_data,p,upd,info)
       use psb_serial_mod
       use psb_descriptor_type
       use psb_prec_type
       integer, intent(out) :: info
       type(psb_dspmat_type), intent(in), target :: a
       type(psb_desc_type),intent(in)            :: desc_data
       type(psb_dbase_prec), intent(inout)       :: p
       character, intent(in)                     :: upd
     end subroutine psb_dcslu
  end interface

  interface 
    subroutine psb_splu_bld(a,desc_a,p,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      use psb_const_mod
      implicit none 
      
      type(psb_dspmat_type), intent(in)   :: a
      type(psb_desc_type), intent(in)     :: desc_a
      type(psb_dbase_prec), intent(inout) :: p
      integer, intent(out)                :: info
    end subroutine psb_splu_bld
  end interface
  interface 
    subroutine psb_umf_bld(a,desc_a,p,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      use psb_const_mod
      implicit none 
      
      type(psb_dspmat_type), intent(in)   :: a
      type(psb_desc_type), intent(in)     :: desc_a
      type(psb_dbase_prec), intent(inout) :: p
      integer, intent(out)                :: info
    end subroutine psb_umf_bld
  end interface

  ! Local scalars
  Integer      :: err, nnzero, n_row, n_col,I,j,icontxt,&
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
  name = 'psb_precbld'

  if (debug) write(0,*) 'Entering precbld',P%prec,desc_a%matrix_data(:)
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

  if (.not.associated(p%baseprecv)) then 
    !! Error 1: should call precset
  end if
  !
  ! Should add check to ensure all procs have the same... 
  !
  ! ALso should define symbolic names for the preconditioners. 
  !

  call psb_check_def(p%baseprecv(1)%iprcparm(p_type_),'base_prec',diagsc_,is_legal_base_prec)
  allocate(p%baseprecv(1)%desc_data)
  call psb_nullify_desc(p%baseprecv(1)%desc_data)

  select case(p%baseprecv(1)%iprcparm(p_type_)) 
  case (noprec_)
     ! Do nothing. 


  case (diagsc_)

     if (debug) write(0,*) 'Precond: Diagonal scaling'
     ! diagonal scaling

     call psb_realloc(n_col,p%baseprecv(1)%d,info)
     if (info /= 0) then
        call psb_errpush(4010,name,a_err='psb_realloc')
        goto 9999
     end if

     call psb_csrws(p%baseprecv(1)%d,a,info,trans='N')
     if(info /= 0) then
        info=4010
        ch_err='psb_csrws'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

     if (debug) write(ilout+me,*) 'VDIAG ',n_row
     do i=1,n_row
        if (p%baseprecv(1)%d(i).eq.0.0d0) then
           p%baseprecv(1)%d(i)=1.d0
        else
           p%baseprecv(1)%d(i) =  1.d0/p%baseprecv(1)%d(i)
        endif

        if (debug) write(ilout+me,*) i,desc_a%loc_to_glob(i),&
             & p%baseprecv(1)%d(i)
        if (p%baseprecv(1)%d(i).lt.0.d0) then
           write(0,*) me,'Negative RWS? ',i,p%baseprecv(1)%d(i)
        endif
     end do
     if (a%pl(1) /= 0) then
        allocate(work(n_row),stat=info)
        if (info /= 0) then
           info=4000
           call psb_errpush(info,name)
           goto 9999
        end if
        call  psb_gelp('n',a%pl,p%baseprecv(1)%d,desc_a,info)
        if(info /= 0) then
           info=4010
           ch_err='psb_dgelp'
           call psb_errpush(info,name,a_err=ch_err)
           goto 9999
        end if

        deallocate(work)
     endif

     if (debug) then
        allocate(gd(mglob))       
        call  psb_gather(gd, p%baseprecv(1)%d, desc_a, info, iroot=iroot)
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

     call psb_check_def(p%baseprecv(1)%iprcparm(n_ovr_),'overlap',&
          &  0,is_legal_n_ovr)
     call psb_check_def(p%baseprecv(1)%iprcparm(restr_),'restriction',&
          &  psb_halo_,is_legal_restrict)
     call psb_check_def(p%baseprecv(1)%iprcparm(prol_),'prolongator',&
          &  psb_none_,is_legal_prolong)

     if ((p%baseprecv(1)%iprcparm(iren_)<0).or.(p%baseprecv(1)%iprcparm(iren_)>2)) then 
        write(0,*) 'Bad PREC%IRENUM value, defaulting to 0', &
             &  p%baseprecv(1)%iprcparm(iren_)
        p%baseprecv(1)%iprcparm(iren_) = 0
     endif

     if (debug) write(0,*)me, ': Calling PSB_DCSLU'


     select case(p%baseprecv(1)%iprcparm(f_type_))

     case(f_ilu_n_,f_ilu_e_) 
        call psb_cslu(a,desc_a,p%baseprecv(1),iupd,info)
        if(debug) write(0,*)me,': out of psb_dcslu'
        if(info /= 0) then
           info=4010
           ch_err='psb_dcslu'
           call psb_errpush(info,name,a_err=ch_err)
           goto 9999
        end if

     case(f_slu_)
        p%baseprecv(1)%av => null()
        if(debug) write(0,*)me,': calling splu_bld'
        call psb_splu_bld(a,desc_a,p%baseprecv(1),info)
        if(info /= 0) then
           info=4010
           ch_err='splu_bld'
           call psb_errpush(info,name,a_err=ch_err)
           goto 9999
        end if

     case(f_umf_)
        p%baseprecv(1)%av => null()
        if(debug) write(0,*)me,': calling umf_bld'
        call psb_umf_bld(a,desc_a,p%baseprecv(1),info)
        if(info /= 0) then
           info=4010
           ch_err='umf_bld'
           call psb_errpush(info,name,a_err=ch_err)
           goto 9999
        end if

     case(f_none_) 
        write(0,*) 'Fact=None in PRECBLD Bja/ASM??'

     case default
        write(0,*) 'Unknown factor type in precbld bja/asm: ',&
             &p%baseprecv(1)%iprcparm(f_type_)
     end select

  end select


  if (size(p%baseprecv) >1) then 
    
    if (.not.associated(p%baseprecv(2)%iprcparm)) then 
      info = 2222
      call psb_errpush(info,name)
      goto 9999
    endif
    call psb_check_def(p%baseprecv(2)%iprcparm(aggr_alg_),'aggregation',&
         &   loc_aggr_,is_legal_ml_aggr_kind)
    call psb_check_def(p%baseprecv(2)%iprcparm(smth_kind_),'Smoother kind',&
         &   smth_omg_,is_legal_ml_smth_kind)
    call psb_check_def(p%baseprecv(2)%iprcparm(coarse_mat_),'Coarse matrix',&
         &   mat_distr_,is_legal_ml_coarse_mat)
    call psb_check_def(p%baseprecv(2)%iprcparm(smth_pos_),'smooth_pos',&
         &   pre_smooth_,is_legal_ml_smooth_pos)
    call psb_check_def(p%baseprecv(2)%iprcparm(f_type_),'fact',f_ilu_n_,is_legal_ml_fact)

    allocate(p%baseprecv(2)%desc_data)
    call psb_nullify_desc(p%baseprecv(2)%desc_data)

    select case(p%baseprecv(2)%iprcparm(f_type_))
    case(f_ilu_n_)      
      call psb_check_def(p%baseprecv(2)%iprcparm(ilu_fill_in_),'Level',0,is_legal_ml_lev)
    case(f_ilu_e_)                 
      call psb_check_def(p%baseprecv(2)%dprcparm(fact_eps_),'Eps',0.0d0,is_legal_ml_eps)
    end select
    call psb_check_def(p%baseprecv(2)%dprcparm(smooth_omega_),'omega',0.0d0,is_legal_omega)
    call psb_check_def(p%baseprecv(2)%iprcparm(jac_sweeps_),'Jacobi sweeps',&
         & 1,is_legal_jac_sweeps)

    call psb_mlprec_bld(a,desc_a,p%baseprecv(2),info)
    if(info /= 0) then
       info=4010
       ch_err='psb_mlprec_bld'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    end if

  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error()
     return
  end if
  return

end subroutine psb_dprecbld


subroutine psb_splu_bld(a,desc_a,p,info)
  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_tools_mod
  use psb_const_mod
  implicit none 

  type(psb_dspmat_type), intent(in)   :: a
  type(psb_desc_type), intent(in)     :: desc_a
  type(psb_dbase_prec), intent(inout) :: p
  integer, intent(out)                :: info


  type(psb_dspmat_type)    :: blck, atmp
  character(len=5)         :: fmt
  character                :: upd='F'
  integer                  :: i,j,nza,nzb,nzt,icontxt, me,mycol,nprow,npcol,err_act
  logical, parameter :: debug=.false.
  character(len=20)   :: name, ch_err

  interface psb_csrsetup
    Subroutine psb_dcsrsetup(ptype,novr,a,blk,desc_data,upd,desc_p,info,outfmt)
      use psb_serial_mod
      Use psb_descriptor_type
      Use psb_prec_type
      integer, intent(in)                  :: ptype,novr
      Type(psb_dspmat_type), Intent(in)    ::  a
      Type(psb_dspmat_type), Intent(inout) ::  blk
      Type(psb_desc_type), Intent(inout)   :: desc_p
      Type(psb_desc_type), Intent(in)      :: desc_data 
      Character, Intent(in)                :: upd
      integer, intent(out)                 :: info
      character(len=5), optional           :: outfmt
    end Subroutine psb_dcsrsetup
 end interface

  if(psb_get_errstatus().ne.0) return 
  info=0
  name='psb_splu_bld'
  call psb_erractionsave(err_act)

  icontxt = desc_A%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt, nprow, npcol, me, mycol)


  fmt = 'COO'
  call psb_nullify_sp(blck)    
  call psb_nullify_sp(atmp)    
  atmp%fida='COO'
  if (Debug) then 
     write(0,*) me, 'SPLUBLD: Calling  csdp'
     call blacs_barrier(icontxt,'All')
  endif

  call psb_dcsdp(a,atmp,info)
  if(info /= 0) then
     info=4010
     ch_err='psb_dcsdp'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if
  nza = atmp%infoa(psb_nnz_)
  if (Debug) then 
     write(0,*) me, 'SPLUBLD: Done csdp',info,nza,atmp%m,atmp%k
     call blacs_barrier(icontxt,'All')
  endif
  call psb_csrsetup(p%iprcparm(p_type_),p%iprcparm(n_ovr_),a,&
       & blck,desc_a,upd,p%desc_data,info,outfmt=fmt)  
  if(info /= 0) then
     info=4010
     ch_err='psb_csrsetup'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  nzb = blck%infoa(psb_nnz_)
  if (Debug) then 
     write(0,*) me, 'SPLUBLD: Done csrsetup',info,nzb,blck%fida
     call blacs_barrier(icontxt,'All')
  endif
  if (nzb > 0 ) then 
     if (size(atmp%aspk)<nza+nzb) then 
        call psb_spreall(atmp,nza+nzb,info)
        if(info /= 0) then
           info=4010
           ch_err='psb_spreall'
           call psb_errpush(info,name,a_err=ch_err)
           goto 9999
        end if
     endif
     if (Debug) then 
        write(0,*) me, 'SPLUBLD: Done realloc',info,nza+nzb,atmp%fida
        call blacs_barrier(icontxt,'All')
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
  call psb_spinfo(psb_nztotreq_,atmp,nzt,info)
  if(info /= 0) then
     info=4010
     ch_err='psb_spinfo'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if
  if (Debug) then 
     write(0,*) me,'Calling fort_slu_factor ',nzt,atmp%m,&
          & atmp%k,p%desc_data%matrix_data(psb_n_row_)
     call blacs_barrier(icontxt,'All')
  endif

  call fort_slu_factor(atmp%m,nzt,&
       & atmp%aspk,atmp%ia2,atmp%ia1,p%iprcparm(slu_ptr_),info)
  if(info /= 0) then
     info=4010
     ch_err='fort_slu_fact'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if (Debug) then 
     write(0,*) me, 'SPLUBLD: Done slu_Factor',info,p%iprcparm(slu_ptr_)
     call blacs_barrier(icontxt,'All')
  endif

  call psb_spfree(blck,info)
  call psb_spfree(atmp,info)
  if(info /= 0) then
     info=4010
     ch_err='psb_spfree'
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

end subroutine psb_splu_bld



subroutine psb_umf_bld(a,desc_a,p,info)
  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_tools_mod
  use psb_const_mod
  implicit none 

  type(psb_dspmat_type), intent(in)   :: a
  type(psb_desc_type), intent(in)     :: desc_a
  type(psb_dbase_prec), intent(inout) :: p
  integer, intent(out)                :: info


  type(psb_dspmat_type)    :: blck, atmp
  character(len=5)         :: fmt
  character                :: upd='F'
  integer                  :: i,j,nza,nzb,nzt,icontxt, me,mycol,nprow,npcol,err_act
  logical, parameter :: debug=.false.
  character(len=20)   :: name, ch_err

  interface psb_csrsetup
    Subroutine psb_dcsrsetup(ptype,novr,a,blk,desc_data,upd,desc_p,info,outfmt)
      use psb_serial_mod
      Use psb_descriptor_type
      Use psb_prec_type
      integer, intent(in)                  :: ptype,novr
      Type(psb_dspmat_type), Intent(in)    ::  a
      Type(psb_dspmat_type), Intent(inout) ::  blk
      Type(psb_desc_type), Intent(inout)   :: desc_p
      Type(psb_desc_type), Intent(in)      :: desc_data 
      Character, Intent(in)                :: upd
      integer, intent(out)                 :: info
      character(len=5), optional           :: outfmt
    end Subroutine psb_dcsrsetup
 end interface

  info=0
  name='psb_umf_bld'
  call psb_erractionsave(err_act)

  icontxt = desc_A%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt, nprow, npcol, me, mycol)


  fmt = 'COO'
  call psb_nullify_sp(blck)    
  call psb_nullify_sp(atmp)    
  atmp%fida='COO'
  if (Debug) then 
     write(0,*) me, 'UMFBLD: Calling  csdp'
     call blacs_barrier(icontxt,'All')
  endif

  call psb_dcsdp(a,atmp,info)
  if(info /= 0) then
     info=4010
     ch_err='psb_dcsdp'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if
  nza = atmp%infoa(psb_nnz_)
  if (Debug) then 
     write(0,*) me, 'UMFBLD: Done csdp',info,nza,atmp%m,atmp%k
     call blacs_barrier(icontxt,'All')
  endif
  call psb_csrsetup(p%iprcparm(p_type_),p%iprcparm(n_ovr_),a,&
       & blck,desc_a,upd,p%desc_data,info,outfmt=fmt)  
  if(info /= 0) then
     info=4010
     ch_err='psb_csrsetup'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  nzb = blck%infoa(psb_nnz_)
  if (Debug) then 
     write(0,*) me, 'UMFBLD: Done csrsetup',info,nzb,blck%fida
     call blacs_barrier(icontxt,'All')
  endif
  if (nzb > 0 ) then 
     if (size(atmp%aspk)<nza+nzb) then 
        call psb_spreall(atmp,nza+nzb,info)
        if(info /= 0) then
           info=4010
           ch_err='psb_spreall'
           call psb_errpush(info,name,a_err=ch_err)
           goto 9999
        end if
     endif
     if (Debug) then 
        write(0,*) me, 'UMFBLD: Done realloc',info,nza+nzb,atmp%fida
        call blacs_barrier(icontxt,'All')
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


  call psb_ipcoo2csc(atmp,info)
  if(info /= 0) then
     info=4010
     ch_err='psb_ipcoo2csc'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if
  call psb_spinfo(psb_nztotreq_,atmp,nzt,info)
  if(info /= 0) then
     info=4010
     ch_err='psb_spinfo'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if
  if (Debug) then 
     write(0,*) me,'Calling fort_slu_factor ',nzt,atmp%m,&
          & atmp%k,p%desc_data%matrix_data(psb_n_row_)
     call blacs_barrier(icontxt,'All')
  endif

  call fort_umf_factor(atmp%m,nzt,&
       & atmp%aspk,atmp%ia1,atmp%ia2,&
       & p%iprcparm(umf_symptr_),p%iprcparm(umf_numptr_),info)
  if(info /= 0) then
     info=4010
     ch_err='fort_umf_fact'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if (Debug) then 
     write(0,*) me, 'UMFBLD: Done umf_Factor',info,p%iprcparm(umf_numptr_)
     call blacs_barrier(icontxt,'All')
  endif
  call psb_spfree(blck,info)
  call psb_spfree(atmp,info)
  if(info /= 0) then
     info=4010
     ch_err='psb_spfree'
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

end subroutine psb_umf_bld




subroutine psb_mlprec_bld(a,desc_a,p,info)

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
  
  interface psb_splu
     subroutine psb_dsplu(a,l,u,d,info,blck)
       use psb_spmat_type
       integer, intent(out)                ::     info
       type(psb_dspmat_type),intent(in)    :: a
       type(psb_dspmat_type),intent(inout) :: l,u
       type(psb_dspmat_type),intent(in), optional, target :: blck
       real(kind(1.d0)), intent(inout)     ::  d(:)
     end subroutine psb_dsplu
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
  allocate(p%av(smth_avsz))

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
     call psb_errpush(info,name)
     goto 9999
  end if

  allocate(p%d(nrg)) 

  select case(p%iprcparm(f_type_)) 
  case(f_ilu_n_,f_ilu_e_) 
     call psb_spreall(p%av(l_pr_),nzg,info)
     call psb_spreall(p%av(u_pr_),nzg,info)
     call psb_splu(p%av(ac_),p%av(l_pr_),p%av(u_pr_),p%d,info)
     if(info /= 0) then
        info=4011
        call psb_errpush(info,name)
        goto 9999
     end if

  case(f_slu_) 
    call psb_spall(0,0,p%av(l_pr_),1,info)
    call psb_spall(0,0,p%av(u_pr_),1,info)
    call psb_ipcsr2coo(p%av(ac_),info)
     if(info /= 0) then
        info=4011
        call psb_errpush(info,name)
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
    call fort_slu_factor(nrg,nzg,&
         & p%av(ac_)%aspk,p%av(ac_)%ia2,p%av(ac_)%ia1,p%iprcparm(slu_ptr_),info)
     if(info /= 0) then
        info=4011
        call psb_errpush(info,name)
        goto 9999
     end if

  case(f_umf_) 
    call psb_spall(0,0,p%av(l_pr_),1,info)
    call psb_spall(0,0,p%av(u_pr_),1,info)
    call psb_ipcsr2coo(p%av(ac_),info)
     if(info /= 0) then
        info=4011
        call psb_errpush(info,name)
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
    call fort_umf_factor(nrg,nzg,&
         & p%av(ac_)%aspk,p%av(ac_)%ia1,p%av(ac_)%ia2,&
         & p%iprcparm(umf_symptr_),p%iprcparm(umf_numptr_),info)
     if(info /= 0) then
        info=4011
        call psb_errpush(info,name)
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

end subroutine psb_mlprec_bld
