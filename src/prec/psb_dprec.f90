<<<<<<< variant A
  /home/buttari/NUMERICAL/psblas-2.0/src/prec:
  total used in directory 1000 available 4183416
  drwxr-xr-x  2 buttari users   4096 Oct 14 11:26 CVS
  -rw-r--r--  1 buttari users   9723 Sep  7 11:16 fort_slu_impl.c
  -rw-r--r--  1 buttari users   3616 Oct 13 17:19 fort_slu_impl.o
  -rw-r--r--  1 buttari users  20777 Sep  7 11:16 gps.f
  -rw-r--r--  1 buttari users  47396 Oct 13 17:19 gps.o
  -rw-r--r--  1 buttari users    623 Sep 12 17:38 Makefile
  -rw-r--r--  1 buttari users  23481 Oct 13 20:40 psb_dbldaggrmat.f90
  -rw-r--r--  1 buttari users 158084 Oct 13 20:40 psb_dbldaggrmat.o
  -rw-r--r--  1 buttari users  19821 Sep 27 11:35 psb_dcslu.f90
  -rw-r--r--  1 buttari users 109004 Oct 13 17:19 psb_dcslu.o
  -rw-r--r--  1 buttari users   6503 Sep 12 14:28 psb_dcsrsetup.f90
  -rw-r--r--  1 buttari users  15096 Oct 13 17:19 psb_dcsrsetup.o
  -rw-r--r--  1 buttari users   6868 Oct 14 10:38 psb_dgenaggrmap.f90
  -rw-r--r--  1 buttari users  37432 Oct 14 10:39 psb_dgenaggrmap.o
  -rw-r--r--  1 buttari users  17221 Oct 14 10:38 psb_dprecbld.f90
  -rw-r--r--  1 buttari users 111376 Oct 14 10:39 psb_dprecbld.o
  -rw-r--r--  1 buttari users  27689 Oct 13 12:40 psb_dprec.f90
  -rw-r--r--  1 buttari users   1669 Oct  7 16:25 psb_dprecfree.f90
  -rw-r--r--  1 buttari users   5632 Oct 13 17:19 psb_dprecfree.o
  -rw-r--r--  1 buttari users 200168 Oct 13 17:19 psb_dprec.o
  -rw-r--r--  1 buttari users   5353 Sep 30 14:29 psb_dprecset.f90
  -rw-r--r--  1 buttari users  72964 Oct 13 17:19 psb_dprecset.o
  -rw-r--r--  1 buttari users  11328 Oct  5 17:28 psb_dsplu.f90
  -rw-r--r--  1 buttari users  44108 Oct 13 17:19 psb_dsplu.o
>>>>>>> variant B
subroutine psb_dprecaply(prec,x,y,desc_data,info,trans, work)

  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_psblas_mod
  use psb_const_mod
  use psb_error_mod
  implicit none

  type(psb_desc_type),intent(in)      :: desc_data
  type(psb_dprec_type), intent(in)    :: prec
  real(kind(0.d0)),intent(inout)      :: x(:), y(:)
  integer, intent(out)                :: info
  character(len=1), optional          :: trans
  real(kind(0.d0)), optional, target  :: work(:)

  ! Local variables
  character     ::trans_ 
  real(kind(1.d0)), pointer :: work_(:)
  integer :: icontxt,nprow,npcol,me,mycol,err_act, int_err(5)
  logical,parameter                 :: debug=.false., debugprt=.false.
  real(kind(1.d0)), parameter       :: one=1.d0, zero=0.d0
  external mpi_wtime
  character(len=20)   :: name, ch_err

  interface
     subroutine psb_dbaseprcaply(prec,x,beta,y,desc_data,trans,work,info)
       use psb_descriptor_type
       use psb_prec_type
       type(psb_desc_type),intent(in)    :: desc_data
       type(psb_dbase_prec), intent(in)  :: prec
       real(kind(0.d0)),intent(inout)    :: x(:), y(:)
       real(kind(0.d0)),intent(in)       :: beta
       character(len=1)                  :: trans
       real(kind(0.d0)),target           :: work(:)
       integer, intent(out)              :: info
     end subroutine psb_dbaseprcaply
  end interface

  interface
     subroutine psb_dmlprcaply(baseprecv,x,beta,y,desc_data,trans,work,info)
       use psb_descriptor_type
       use psb_prec_type
       type(psb_desc_type),intent(in)    :: desc_data
       type(psb_dbase_prec), intent(in)  :: baseprecv(:)
       real(kind(0.d0)),intent(in)       :: beta
       real(kind(0.d0)),intent(inout)    :: x(:), y(:)
       character                         :: trans
       real(kind(0.d0)),target           :: work(:)
       integer, intent(out)              :: info
     end subroutine psb_dmlprcaply
  end interface
  
  name='psb_dprecaply'
  info = 0
  call psb_erractionsave(err_act)

  icontxt=desc_data%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,nprow,npcol,me,mycol)

  if (present(trans)) then 
    trans_=trans
  else
    trans_='N'
  end if

  if (present(work)) then 
    work_ => work
  else
    allocate(work_(4*desc_data%matrix_data(psb_n_col_)))
  end if

  if (.not.(associated(prec%baseprecv))) then 
    write(0,*) 'Inconsistent preconditioner: neither SMTH nor BASE?'      
  end if
  if (size(prec%baseprecv) >1) then 
    call psb_dmlprcaply(prec%baseprecv,x,zero,y,desc_data,trans_,work_,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psb_dmlprcaply')
      goto 9999
    end if

  else  if (size(prec%baseprecv) == 1) then 
    call psb_dbaseprcaply(prec%baseprecv(1),x,zero,y,desc_data,trans_, work_,info)
  else 
    write(0,*) 'Inconsistent preconditioner: size of baseprecv???' 
  endif

  if (present(work)) then 
  else
    deallocate(work_)
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

end subroutine psb_dprecaply



subroutine psb_dbaseprcaply(prec,x,beta,y,desc_data,trans,work,info)
  !
  !  Compute   Y <-  beta*Y + K^-1 X 
  !  where K is a a basic preconditioner stored in prec
  ! 

  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_psblas_mod
  use psb_const_mod
  use psb_error_mod
  implicit none 

  type(psb_desc_type),intent(in)    :: desc_data
  type(psb_dbase_prec), intent(in)  :: prec
  real(kind(0.d0)),intent(inout)    :: x(:), y(:)
  real(kind(0.d0)),intent(in)       :: beta
  character(len=1)                  :: trans
  real(kind(0.d0)),target           :: work(:)
  integer, intent(out)              :: info

  ! Local variables
  integer :: n_row,n_col, int_err(5)
  real(kind(1.d0)), pointer :: ww(:), aux(:), tx(:),ty(:)
  character     ::diagl, diagu
  integer :: icontxt,nprow,npcol,me,mycol,i, isz, nrg, err_act
  real(kind(1.d0)) :: t1, t2, t3, t4, t5, t6, t7, mpi_wtime
  logical,parameter                 :: debug=.false., debugprt=.false.
  real(kind(1.d0)), parameter       :: one=1.d0, zero=0.d0
  external mpi_wtime
  character(len=20)   :: name, ch_err

  interface
     subroutine psb_dbjacaply(prec,x,beta,y,desc_data,trans,work,info)
       use psb_descriptor_type
       use psb_prec_type
       type(psb_desc_type), intent(in)       :: desc_data
       type(psb_dbase_prec), intent(in)      :: prec
       real(kind(0.d0)),intent(inout)        :: x(:), y(:)
       real(kind(0.d0)),intent(in)           :: beta
       character(len=1)                      :: trans
       real(kind(0.d0)),target               :: work(:)
       integer, intent(out)                  :: info
     end subroutine psb_dbjacaply
  end interface
  
  name='psb_dbaseprcaply'
  info = 0
  call psb_erractionsave(err_act)

  icontxt=desc_data%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,nprow,npcol,me,mycol)

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

    n_row=desc_data%matrix_data(psb_n_row_)
    if (beta==zero) then 
      y(1:n_row) = x(1:n_row)
    else if (beta==one) then 
      y(1:n_row) = x(1:n_row) + y(1:n_row)
    else if (beta==-one) then 
      y(1:n_row) = x(1:n_row) - y(1:n_row)
    else 
      y(1:n_row) = x(1:n_row) + beta*y(1:n_row)
    end if

  case(diagsc_)

    n_row=desc_data%matrix_data(psb_n_row_)
    if (beta==zero) then 
      y(1:n_row) = x(1:n_row)*prec%d(1:n_row)
    else if (beta==one) then 
      y(1:n_row) = x(1:n_row)*prec%d(1:n_row) + y(1:n_row)
    else if (beta==-one) then 
      y(1:n_row) = x(1:n_row)*prec%d(1:n_row) - y(1:n_row)
    else 
      y(1:n_row) = x(1:n_row)*prec%d(1:n_row) + beta*y(1:n_row)
    end if

  case(bja_)

    call psb_dbjacaply(prec,x,beta,y,desc_data,trans,work,info)
    if(info.ne.0) then
       info=4010
       ch_err='psb_bjacaply'
       goto 9999
    end if

  case(asm_,ras_,ash_,rash_)

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
      allocate(ww(isz),tx(isz),ty(isz))
    else if ((3*isz) <= size(work)) then 
      ww => work(1:isz)
      tx => work(isz+1:2*isz)
      ty => work(2*isz+1:3*isz)
      allocate(aux(4*isz))
    else 
      allocate(ww(isz),tx(isz),ty(isz),&
           &aux(4*isz))
    endif

    if (debug) write(0,*)' vdiag: ',prec%d(:)
    if (debug) write(0,*) 'Bi-CGSTAB with Additive Schwarz prec' 

    tx(1:desc_data%matrix_data(psb_n_row_)) = x(1:desc_data%matrix_data(psb_n_row_)) 
    tx(desc_data%matrix_data(psb_n_row_)+1:isz) = zero

    if (prec%iprcparm(restr_)==psb_halo_) then 
       call psb_halo(tx,prec%desc_data,info,work=aux)
      if(info /=0) then
         info=4010
         ch_err='psb_halo'
         goto 9999
      end if
    else if (prec%iprcparm(restr_) /= psb_none_) then 
      write(0,*) 'Problem in PRCAPLY: Unknown value for restriction ',&
           &prec%iprcparm(restr_)
    end if

    if (prec%iprcparm(iren_)>0) then 
      call psb_dgelp('N',n_row,1,prec%perm,tx,isz,ww,isz,info)
      if(info /=0) then
         info=4010
         ch_err='psb_dgelp'
         goto 9999
      end if
    endif

    call psb_dbjacaply(prec,tx,zero,ty,prec%desc_data,trans,aux,info)

    if (prec%iprcparm(iren_)>0) then 
      call psb_dgelp('N',n_row,1,prec%invperm,ty,isz,ww,isz,info)
      if(info /=0) then
         info=4010
         ch_err='psb_dgelp'
         goto 9999
      end if
    endif

    select case (prec%iprcparm(prol_)) 

    case(psb_none_) 
      ! Would work anyway, but since it's supposed to do nothing...
      ! call f90_psovrl(ty,prec%desc_data,update_type=prec%a_restrict)

    case(psb_sum_,psb_avg_) 
       call psb_ovrl(ty,prec%desc_data,info,&
            & update_type=prec%iprcparm(prol_),work=aux)
       if(info /=0) then
         info=4010
         ch_err='psb_ovrl'
         goto 9999
      end if

    case default
      write(0,*) 'Problem in PRCAPLY: Unknown value for prolongation ',&
           & prec%iprcparm(prol_)
    end select

    if (beta == zero)  then 
      y(1:desc_data%matrix_data(psb_n_row_)) = ty(1:desc_data%matrix_data(psb_n_row_)) 
    else if (beta == one) then 
      y(1:desc_data%matrix_data(psb_n_row_)) = y(1:desc_data%matrix_data(psb_n_row_)) + &
           & ty(1:desc_data%matrix_data(psb_n_row_)) 
    else if (beta == -one) then 
      y(1:desc_data%matrix_data(psb_n_row_)) = -y(1:desc_data%matrix_data(psb_n_row_)) + &
           & ty(1:desc_data%matrix_data(psb_n_row_)) 
    else 
      y(1:desc_data%matrix_data(psb_n_row_)) = beta*y(1:desc_data%matrix_data(psb_n_row_)) + &
           & ty(1:desc_data%matrix_data(psb_n_row_)) 
    end if


    if ((6*isz) <= size(work)) then 
    else if ((4*isz) <= size(work)) then 
      deallocate(ww,tx,ty)
    else if ((3*isz) <= size(work)) then 
      deallocate(aux)
    else 
      deallocate(ww,aux,tx,ty)
    endif

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

end subroutine psb_dbaseprcaply



subroutine psb_dbjacaply(prec,x,beta,y,desc_data,trans,work,info)
  !
  !  Compute   Y <-  beta*Y + K^-1 X 
  !  where K is a a Block Jacobi  preconditioner stored in prec
  !  Note that desc_data may or may not be the same as prec%desc_data,
  !  but since both are INTENT(IN) this should be legal. 
  ! 

  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_psblas_mod
  use psb_const_mod
  use psb_error_mod
  implicit none 

  type(psb_desc_type), intent(in)       :: desc_data
  type(psb_dbase_prec), intent(in)      :: prec
  real(kind(0.d0)),intent(inout)        :: x(:), y(:)
  real(kind(0.d0)),intent(in)           :: beta
  character(len=1)                      :: trans
  real(kind(0.d0)),target               :: work(:)
  integer, intent(out)                  :: info

  ! Local variables
  integer :: n_row,n_col
  real(kind(1.d0)), pointer :: ww(:), aux(:), tx(:),ty(:),tb(:)
  character     ::diagl, diagu
  integer :: icontxt,nprow,npcol,me,mycol,i, isz, nrg, err_act, int_err(5)
  real(kind(1.d0)) :: t1, t2, t3, t4, t5, t6, t7, mpi_wtime
  logical,parameter                 :: debug=.false., debugprt=.false.
  real(kind(1.d0)), parameter       :: one=1.d0, zero=0.d0
  external mpi_wtime
  character(len=20)   :: name, ch_err

  name='psb_dbjacaply'
  info = 0
  call psb_erractionsave(err_act)

  icontxt=desc_data%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,nprow,npcol,me,mycol)

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
      allocate(aux(4*n_col))
    endif
  else
    allocate(ww(n_col),aux(4*n_col))
  endif


  if (prec%iprcparm(jac_sweeps_) == 1) then 


    select case(prec%iprcparm(f_type_))
    case(f_ilu_n_,f_ilu_e_) 

      select case(trans)
      case('N','n')

        call psb_spsm(one,prec%av(l_pr_),x,zero,ww,desc_data,info,&
             & trans='N',unit=diagl,choice=psb_none_,work=aux)
        if(info /=0) goto 9999
        ww(1:n_row) = ww(1:n_row)*prec%d(1:n_row)
        call psb_spsm(one,prec%av(u_pr_),ww,beta,y,desc_data,info,&
             & trans='N',unit=diagu,choice=psb_none_, work=aux)
        if(info /=0) goto 9999

      case('T','t','C','c')
        call psb_spsm(one,prec%av(u_pr_),x,zero,ww,desc_data,info,&
             & trans=trans,unit=diagu,choice=psb_none_, work=aux)
        if(info /=0) goto 9999
        ww(1:n_row) = ww(1:n_row)*prec%d(1:n_row)
        call psb_spsm(one,prec%av(l_pr_),ww,beta,y,desc_data,info,&
             & trans=trans,unit=diagl,choice=psb_none_,work=aux)
        if(info /=0) goto 9999

      end select

    case(f_slu_)

      ww(1:n_row) = x(1:n_row)

      select case(trans)
      case('N','n')
        call fort_slu_solve(0,n_row,1,ww,n_row,prec%iprcparm(slu_ptr_),info)
      case('T','t','C','c')
        call fort_slu_solve(1,n_row,1,ww,n_row,prec%iprcparm(slu_ptr_),info)
      end select

      if(info /=0) goto 9999

      if (beta == 0.d0) then 
        y(1:n_row) = ww(1:n_row)
      else if (beta==1.d0) then 
        y(1:n_row) = ww(1:n_row) + y(1:n_row) 
      else if (beta==-1.d0) then 
        y(1:n_row) = ww(1:n_row) - y(1:n_row) 
      else 
        y(1:n_row) = ww(1:n_row) + beta*y(1:n_row) 
      endif

    case default
      write(0,*) 'Unknown factorization type in bjac_prcaply',prec%iprcparm(f_type_)
    end select
    if (debug) write(0,*)' Y: ',y(:)

  else if (prec%iprcparm(jac_sweeps_) > 1) then 

    ! Note: we have to add TRANS to this one !!!!!!!!! 

    if (size(prec%av) < ap_nd_) then 
      info = 4011
      goto 9999
    endif

    allocate(tx(n_col),ty(n_col))
    tx = zero
    ty = zero
    select case(prec%iprcparm(f_type_)) 
    case(f_ilu_n_,f_ilu_e_) 
      do i=1, prec%iprcparm(jac_sweeps_) 
        !   X(k+1) = M^-1*(b-N*X(k))
        ty(1:n_row) = x(1:n_row)
        call psb_spmm(-one,prec%av(ap_nd_),tx,one,ty,&
             &   prec%desc_data,info,work=aux)
        if(info /=0) goto 9999
        call psb_spsm(one,prec%av(l_pr_),ty,zero,ww,&
             & prec%desc_data,info,&
             & trans='N',unit='U',choice=psb_none_,work=aux)
        if(info /=0) goto 9999
        ww(1:n_row) = ww(1:n_row)*prec%d(1:n_row)
        call psb_spsm(one,prec%av(u_pr_),ww,zero,tx,&
             & prec%desc_data,info,&
             & trans='N',unit='U',choice=psb_none_,work=aux)
        if(info /=0) goto 9999
      end do

    case(f_slu_) 
      do i=1, prec%iprcparm(jac_sweeps_) 
        !   X(k+1) = M^-1*(b-N*X(k))
        ty(1:n_row) = x(1:n_row)
        call psb_spmm(-one,prec%av(ap_nd_),tx,one,ty,&
             &   prec%desc_data,info,work=aux)
        if(info /=0) goto 9999

        call fort_slu_solve(0,n_row,1,ty,n_row,prec%iprcparm(slu_ptr_),info)
        if(info /=0) goto 9999
        tx(1:n_row) = ty(1:n_row)        
      end do
    end select
    
    if (beta == 0.d0) then 
      y(1:n_row) = tx(1:n_row)
    else if (beta==1.d0) then 
      y(1:n_row) = tx(1:n_row) + y(1:n_row) 
    else if (beta==-1.d0) then 
      y(1:n_row) = tx(1:n_row) - y(1:n_row) 
    else 
      y(1:n_row) = tx(1:n_row) + beta*y(1:n_row) 
    endif

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
  if (err_act.eq.act_abort) then
    call psb_error()
    return
  end if
  return

end subroutine psb_dbjacaply


subroutine psb_dmlprcaply(baseprecv,x,beta,y,desc_data,trans,work,info)
  !
  !  Compute   Y <-  beta*Y + K^-1 X 
  !  where K is a multilevel (actually 2-level) preconditioner stored in prec
  ! 

  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_psblas_mod
  use psb_blacs_mod
  use psb_const_mod
  use psb_error_mod
  implicit none

  type(psb_desc_type),intent(in)    :: desc_data
  type(psb_dbase_prec), intent(in)  :: baseprecv(:)
  real(kind(0.d0)),intent(in)       :: beta
  real(kind(0.d0)),intent(inout)    :: x(:), y(:)
  character                         :: trans
  real(kind(0.d0)),target           :: work(:)
  integer, intent(out)              :: info


  ! Local variables
  integer :: n_row,n_col
  real(kind(1.d0)), allocatable :: tx(:),ty(:),t2l(:),w2l(:),&
       &   x2l(:),b2l(:),tz(:),tty(:)
  character     ::diagl, diagu
  integer :: icontxt,nprow,npcol,me,mycol,i, isz, nrg,nr2l,err_act, iptype, int_err(5)
  real(kind(1.d0)) :: omega
  real(kind(1.d0)) :: t1, t2, t3, t4, t5, t6, t7, mpi_wtime
  logical, parameter          :: debug=.false., debugprt=.false.
  real(kind(1.d0)), parameter :: one=1.d0, zero=0.d0
  integer      :: ismth
  external mpi_wtime
  character(len=20)   :: name, ch_err

  interface
     subroutine psb_dbaseprcaply(prec,x,beta,y,desc_data,trans,work,info)
       use psb_descriptor_type
       use psb_prec_type
       type(psb_desc_type),intent(in)    :: desc_data
       type(psb_dbase_prec), intent(in)  :: prec
       real(kind(0.d0)),intent(inout)    :: x(:), y(:)
       real(kind(0.d0)),intent(in)       :: beta
       character(len=1)                  :: trans
       real(kind(0.d0)),target           :: work(:)
       integer, intent(out)              :: info
     end subroutine psb_dbaseprcaply
  end interface

  name='psb_dmlprcaply'
  info = 0
  call psb_erractionsave(err_act)


  icontxt=desc_data%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,nprow,npcol,me,mycol)

  omega=baseprecv(2)%dprcparm(smooth_omega_)
  ismth=baseprecv(2)%iprcparm(smth_kind_)

  select case(baseprecv(2)%iprcparm(ml_type_)) 
  case(no_ml_) 
    ! Should not really get here.
    write(0,*) 'Smooth preconditioner with no multilevel in MLPRCAPLY????' 

  case(add_ml_prec_)

    ! 
    !  Additive multilevel
    !
    t1 = mpi_wtime()
    n_row = desc_data%matrix_data(psb_n_row_)
    n_col = baseprecv(1)%desc_data%matrix_data(psb_n_col_)
    call psb_dbaseprcaply(baseprecv(1),x,beta,y,desc_data,trans,work,info)
    if(info /=0) goto 9999

    nr2l  = baseprecv(2)%desc_data%matrix_data(psb_n_col_)
    nrg   = baseprecv(2)%desc_data%matrix_data(psb_n_row_)
    allocate(t2l(nr2l),w2l(nr2l))
    t2l(:) = zero
    w2l(:) = zero

    if (ismth  /= no_smth_) then 
      !
      ! Smoothed aggregation
      !
      allocate(tx(max(n_row,n_col)),ty(max(n_row,n_col)),&
           & tz(max(n_row,n_col)))
      tx(1:desc_data%matrix_data(psb_n_row_)) = x(1:desc_data%matrix_data(psb_n_row_)) 
      tx(desc_data%matrix_data(psb_n_row_)+1:max(n_row,n_col)) = zero
      ty(desc_data%matrix_data(psb_n_row_)+1:max(n_row,n_col)) = zero
      tz(desc_data%matrix_data(psb_n_row_)+1:max(n_row,n_col)) = zero


      if (baseprecv(2)%iprcparm(glb_smth_) >0) then 
        call psb_halo(tx,desc_data,info,work=work) 
        if(info /=0) goto 9999
      else
        tx(desc_data%matrix_data(psb_n_row_)+1:max(n_row,n_col)) = zero
      end if

      call psb_csmm(one,baseprecv(2)%av(sm_pr_t_),tx,zero,t2l,info)
      if(info /=0) goto 9999

    else
      !
      ! Raw  aggregation, may take shortcut
      !
      do i=1,desc_data%matrix_data(psb_n_row_)
        t2l(baseprecv(2)%mlia(i)) = t2l(baseprecv(2)%mlia(i)) + x(i)
      end do

    end if

    if (baseprecv(2)%iprcparm(coarse_mat_)==mat_repl_) Then 
      call gsum2d(icontxt,'All',t2l(1:nrg))
    else if (baseprecv(2)%iprcparm(coarse_mat_) /= mat_distr_) Then 
        write(0,*) 'Unknown value for baseprecv(2)%iprcparm(coarse_mat_) ',&
             & baseprecv(2)%iprcparm(coarse_mat_)
    endif

    w2l=t2l
    call psb_dbaseprcaply(baseprecv(2),w2l,zero,t2l,baseprecv(2)%desc_data,'N',work,info)


    if (ismth  /= no_smth_) then 

      call psb_csmm(one,baseprecv(2)%av(sm_pr_),t2l,zero,ty,info)
      if(info /=0) goto 9999
      ! 
      ! Finally add back into Y. 
      ! 
      call psb_axpby(one,ty,one,y,desc_data,info)
      if(info /=0) goto 9999
      deallocate(tx,ty,tz)

    else

      do i=1, desc_data%matrix_data(psb_n_row_)
        y(i) = y(i) + t2l(baseprecv(2)%mlia(i))
      enddo

    end if

    if (debug) write(0,*)' Y2: ',Y(:)

    deallocate(t2l,w2l)

  case(mult_ml_prec_)

    ! 
    !  Multiplicative multilevel
    !  Pre/post smoothing versions. 

    select case(baseprecv(2)%iprcparm(smth_pos_))

    case(post_smooth_)


      t1    = mpi_wtime()
      n_row = desc_data%matrix_data(psb_n_row_)
      n_col = baseprecv(1)%desc_data%matrix_data(psb_n_col_)
      nr2l  = baseprecv(2)%desc_data%matrix_data(psb_n_col_)
      nrg   = baseprecv(2)%desc_data%matrix_data(psb_n_row_)
      allocate(t2l(nr2l),w2l(nr2l),tx(n_col),ty(n_col))
      t2l(:) = zero
      w2l(:) = zero

      !
      ! Need temp copies to handle Y<- betaY + K^-1 X
      ! One of the temp copies is not strictly needed when beta==zero
      !

      if (debug) write(0,*)' mult_ml_apply  omega ',omega
      if (debug) write(0,*)' mult_ml_apply  X: ',X(:)
      call psb_axpby(one,x,zero,tx,desc_data,info)
      if(info /=0) goto 9999

      if (ismth  /= no_smth_) then 
        !
        ! Smoothed aggregation
        !
        allocate(tz(max(n_row,n_col)))

        if (baseprecv(2)%iprcparm(glb_smth_) >0) then 
          call psb_halo(tx,desc_data,info,work=work) 
          if(info /=0) goto 9999
        else
          tx(desc_data%matrix_data(psb_n_row_)+1:max(n_row,n_col)) = zero
        end if

        call psb_csmm(one,baseprecv(2)%av(sm_pr_t_),tx,zero,t2l,info)
        if(info /=0) goto 9999

      else
        !
        ! Raw  aggregation, may take shortcut
        !
        do i=1,desc_data%matrix_data(psb_n_row_)
          t2l(baseprecv(2)%mlia(i)) = t2l(baseprecv(2)%mlia(i)) + tx(i)
        end do
      end if

      if (baseprecv(2)%iprcparm(coarse_mat_)==mat_repl_) Then 
        call gsum2d(icontxt,'All',t2l(1:nrg))
      else if (baseprecv(2)%iprcparm(coarse_mat_) /= mat_distr_) Then 
        write(0,*) 'Unknown value for baseprecv(2)%iprcparm(coarse_mat_) ',&
             & baseprecv(2)%iprcparm(coarse_mat_)
      endif

      t6 = mpi_wtime()
      w2l=t2l
      call psb_dbaseprcaply(baseprecv(2),w2l,zero,t2l,baseprecv(2)%desc_data,'N',work,info)
      if(info /=0) goto 9999

      if (ismth  /= no_smth_) then 
        if (ismth == smth_omg_) &
             & call psb_halo(t2l,baseprecv(2)%desc_data,info,work=work) 
        call psb_csmm(one,baseprecv(2)%av(sm_pr_),t2l,zero,ty,info)
        if(info /=0) goto 9999
        ! 
        ! Finally add back into Y. 
        ! 
        deallocate(tz)
      else
        ty(:) = zero
        do i=1, desc_data%matrix_data(psb_n_row_)
          ty(i) = ty(i) + t2l(baseprecv(2)%mlia(i))
        enddo

      end if
      deallocate(t2l,w2l)

      call psb_spmm(-one,baseprecv(2)%aorig,ty,one,tx,desc_data,info,work=work)
      if(info /=0) goto 9999

      call psb_dbaseprcaply(baseprecv(1),tx,one,ty,desc_data,trans,work,info)
      if(info /=0) goto 9999

      call psb_axpby(one,ty,beta,y,desc_data,info)
      if(info /=0) goto 9999

      deallocate(tx,ty)



    case(pre_smooth_)

      t1 = mpi_wtime()
      n_row = desc_data%matrix_data(psb_n_row_)
      n_col = baseprecv(1)%desc_data%matrix_data(psb_n_col_)
      nr2l  = baseprecv(2)%desc_data%matrix_data(psb_n_col_)
      nrg   = baseprecv(2)%desc_data%matrix_data(psb_n_row_)
      allocate(t2l(nr2l),w2l(nr2l),tx(n_col),ty(n_col),tty(n_col))
      t2l(:) = zero
      w2l(:) = zero

      !
      ! Need temp copies to handle Y<- betaY + K^-1 X
      ! One of the temp copies is not strictly needed when beta==zero
      !
      call psb_axpby(one,x,zero,tx,desc_data,info)
      call psb_axpby(one,y,zero,ty,desc_data,info)
      if(info /=0) goto 9999

      call psb_dbaseprcaply(baseprecv(1),x,zero,tty,desc_data,trans,work,info)
      if(info /=0) goto 9999

      call psb_spmm(-one,baseprecv(2)%aorig,tty,one,tx,desc_data,info,work=work)
      if(info /=0) goto 9999

      if (ismth  /= no_smth_) then 
        allocate(tz(max(n_row,n_col)))

        if (baseprecv(2)%iprcparm(glb_smth_) >0) then 
          call psb_halo(tx,desc_data,info,work=work) 
          if(info /=0) goto 9999
        else
          tx(desc_data%matrix_data(psb_n_row_)+1:max(n_row,n_col)) = zero
        end if

        call psb_csmm(one,baseprecv(2)%av(sm_pr_t_),tx,zero,t2l,info)
        if(info /=0) goto 9999

      else
        !
        ! Raw  aggregation, may take shortcuts
        !
        do i=1,desc_data%matrix_data(psb_n_row_)
          t2l(baseprecv(2)%mlia(i)) = t2l(baseprecv(2)%mlia(i)) + tx(i)
        end do
      end if

      if (baseprecv(2)%iprcparm(coarse_mat_)==mat_repl_) Then 
        call gsum2d(icontxt,'All',t2l(1:nrg))
      else if (baseprecv(2)%iprcparm(coarse_mat_) /= mat_distr_) Then 
        write(0,*) 'Unknown value for baseprecv(2)%iprcparm(coarse_mat_) ',&
             & baseprecv(2)%iprcparm(coarse_mat_)
      endif

      t6 = mpi_wtime()
      w2l=t2l
      call psb_dbaseprcaply(baseprecv(2),w2l,zero,t2l,baseprecv(2)%desc_data,'N',work,info)
      if(info /=0) goto 9999

      if (ismth  /= no_smth_) then 

        if (ismth == smth_omg_) &
             & call psb_halo(t2l,baseprecv(2)%desc_data,info,work=work) 
        call psb_csmm(one,baseprecv(2)%av(sm_pr_),t2l,zero,ty,info)
        if(info /=0) goto 9999

        call psb_axpby(one,ty,one,tty,desc_data,info)
        if(info /=0) goto 9999

        deallocate(tz)
      else

        do i=1, desc_data%matrix_data(psb_n_row_)
          tty(i) = tty(i) + t2l(baseprecv(2)%mlia(i))
        enddo

      end if

      call psb_axpby(one,tty,beta,y,desc_data,info)
      if(info /=0) goto 9999

      deallocate(t2l,w2l,tx,ty,tty)




    case default

      write(0,*) 'Unknown value for ml_smooth_pos',baseprecv(2)%iprcparm(smth_pos_)

    end select

  case default
    write(0,*) me, 'Wrong mltype into PRCAPLY ',&
         & baseprecv(2)%iprcparm(ml_type_)
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error()
    return
  end if
  return

end subroutine psb_dmlprcaply


subroutine psb_dprecaply1(prec,x,desc_data,info,trans)

  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_psblas_mod
  use psb_const_mod
  use psb_error_mod
  implicit none

  type(psb_desc_type),intent(in)    :: desc_data
  type(psb_dprec_type), intent(in)  :: prec
  real(kind(0.d0)),intent(inout)    :: x(:)
  integer, intent(out)              :: info
  character(len=1), optional        :: trans
  logical,parameter                 :: debug=.false., debugprt=.false.
  real(kind(1.d0)), parameter       :: one=1.d0, zero=0.d0


  ! Local variables
  character     :: trans_
  integer :: icontxt,nprow,npcol,me,mycol,i, isz, err_act, int_err(5)
  real(kind(1.d0)), pointer :: WW(:), w1(:)
  character(len=20)   :: name, ch_err
  name='psb_dprec1'
  info = 0
  call psb_erractionsave(err_act)
  

  icontxt=desc_data%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,nprow,npcol,me,mycol)
  if (present(trans)) then 
    trans_=trans
  else
    trans_='N'
  end if

  allocate(ww(size(x)),w1(size(x)))
  call psb_dprecaply(prec,x,ww,desc_data,info,trans_,w1)
  if(info /=0) goto 9999
  x(:) = ww(:)
  deallocate(ww,W1)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error()
     return
  end if
  return
end subroutine psb_dprecaply1

####### Ancestor
subroutine psb_dprecaply(prec,x,y,desc_data,info,trans, work)

  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_psblas_mod
  use psb_const_mod
  use psb_error_mod
  implicit none

  type(psb_desc_type),intent(in)      :: desc_data
  type(psb_dprec_type), intent(in)    :: prec
  real(kind(0.d0)),intent(inout)      :: x(:), y(:)
  integer, intent(out)                :: info
  character(len=1), optional          :: trans
  real(kind(0.d0)), optional, target  :: work(:)

  ! Local variables
  character     ::trans_ 
  real(kind(1.d0)), pointer :: work_(:)
  integer :: icontxt,nprow,npcol,me,mycol,err_act, int_err(5)
  logical,parameter                 :: debug=.false., debugprt=.false.
  real(kind(1.d0)), parameter       :: one=1.d0, zero=0.d0
  external mpi_wtime
  character(len=20)   :: name, ch_err

  interface
     subroutine psb_dbaseprcaply(prec,x,beta,y,desc_data,trans,work,info)
       use psb_descriptor_type
       use psb_prec_type
       type(psb_desc_type),intent(in)    :: desc_data
       type(psb_dbase_prec), intent(in)  :: prec
       real(kind(0.d0)),intent(inout)    :: x(:), y(:)
       real(kind(0.d0)),intent(in)       :: beta
       character(len=1)                  :: trans
       real(kind(0.d0)),target           :: work(:)
       integer, intent(out)              :: info
     end subroutine psb_dbaseprcaply
  end interface

  interface
     subroutine psb_dmlprcaply(baseprecv,x,beta,y,desc_data,trans,work,info)
       use psb_descriptor_type
       use psb_prec_type
       type(psb_desc_type),intent(in)    :: desc_data
       type(psb_dbase_prec), intent(in)  :: baseprecv(:)
       real(kind(0.d0)),intent(in)       :: beta
       real(kind(0.d0)),intent(inout)    :: x(:), y(:)
       character                         :: trans
       real(kind(0.d0)),target           :: work(:)
       integer, intent(out)              :: info
     end subroutine psb_dmlprcaply
  end interface
  
  name='psb_dprecaply'
  info = 0
  call psb_erractionsave(err_act)

  icontxt=desc_data%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,nprow,npcol,me,mycol)

  if (present(trans)) then 
    trans_=trans
  else
    trans_='N'
  end if

  if (present(work)) then 
    work_ => work
  else
    allocate(work_(4*desc_data%matrix_data(psb_n_col_)))
  end if

  if (.not.(associated(prec%baseprecv))) then 
    write(0,*) 'Inconsistent preconditioner: neither SMTH nor BASE?'      
  end if
  if (size(prec%baseprecv) >1) then 
    call psb_dmlprcaply(prec%baseprecv,x,zero,y,desc_data,trans_,work_,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psb_dmlprcaply')
      goto 9999
    end if

  else  if (size(prec%baseprecv) == 1) then 
    call psb_dbaseprcaply(prec%baseprecv(1),x,zero,y,desc_data,trans_, work_,info)
  else 
    write(0,*) 'Inconsistent preconditioner: size of baseprecv???' 
  endif

  if (present(work)) then 
  else
    deallocate(work_)
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

end subroutine psb_dprecaply



subroutine psb_dbaseprcaply(prec,x,beta,y,desc_data,trans,work,info)
  !
  !  Compute   Y <-  beta*Y + K^-1 X 
  !  where K is a a basic preconditioner stored in prec
  ! 

  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_psblas_mod
  use psb_const_mod
  use psb_error_mod
  implicit none 

  type(psb_desc_type),intent(in)    :: desc_data
  type(psb_dbase_prec), intent(in)  :: prec
  real(kind(0.d0)),intent(inout)    :: x(:), y(:)
  real(kind(0.d0)),intent(in)       :: beta
  character(len=1)                  :: trans
  real(kind(0.d0)),target           :: work(:)
  integer, intent(out)              :: info

  ! Local variables
  integer :: n_row,n_col, int_err(5)
  real(kind(1.d0)), pointer :: ww(:), aux(:), tx(:),ty(:)
  character     ::diagl, diagu
  integer :: icontxt,nprow,npcol,me,mycol,i, isz, nrg, err_act
  real(kind(1.d0)) :: t1, t2, t3, t4, t5, t6, t7, mpi_wtime
  logical,parameter                 :: debug=.false., debugprt=.false.
  real(kind(1.d0)), parameter       :: one=1.d0, zero=0.d0
  external mpi_wtime
  character(len=20)   :: name, ch_err

  interface
     subroutine psb_dbjacaply(prec,x,beta,y,desc_data,trans,work,info)
       use psb_descriptor_type
       use psb_prec_type
       type(psb_desc_type), intent(in)       :: desc_data
       type(psb_dbase_prec), intent(in)      :: prec
       real(kind(0.d0)),intent(inout)        :: x(:), y(:)
       real(kind(0.d0)),intent(in)           :: beta
       character(len=1)                      :: trans
       real(kind(0.d0)),target               :: work(:)
       integer, intent(out)                  :: info
     end subroutine psb_dbjacaply
  end interface
  
  name='psb_dbaseprcaply'
  info = 0
  call psb_erractionsave(err_act)

  icontxt=desc_data%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,nprow,npcol,me,mycol)

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

    n_row=desc_data%matrix_data(psb_n_row_)
    if (beta==zero) then 
      y(1:n_row) = x(1:n_row)
    else if (beta==one) then 
      y(1:n_row) = x(1:n_row) + y(1:n_row)
    else if (beta==-one) then 
      y(1:n_row) = x(1:n_row) - y(1:n_row)
    else 
      y(1:n_row) = x(1:n_row) + beta*y(1:n_row)
    end if

  case(diagsc_)

    n_row=desc_data%matrix_data(psb_n_row_)
    if (beta==zero) then 
      y(1:n_row) = x(1:n_row)*prec%d(1:n_row)
    else if (beta==one) then 
      y(1:n_row) = x(1:n_row)*prec%d(1:n_row) + y(1:n_row)
    else if (beta==-one) then 
      y(1:n_row) = x(1:n_row)*prec%d(1:n_row) - y(1:n_row)
    else 
      y(1:n_row) = x(1:n_row)*prec%d(1:n_row) + beta*y(1:n_row)
    end if

  case(bja_)

    call psb_dbjacaply(prec,x,beta,y,desc_data,trans,work,info)
    if(info.ne.0) then
       info=4010
       ch_err='psb_bjacaply'
       goto 9999
    end if

  case(asm_,ras_,ash_,rash_)

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
      allocate(ww(isz),tx(isz),ty(isz))
    else if ((3*isz) <= size(work)) then 
      ww => work(1:isz)
      tx => work(isz+1:2*isz)
      ty => work(2*isz+1:3*isz)
      allocate(aux(4*isz))
    else 
      allocate(ww(isz),tx(isz),ty(isz),&
           &aux(4*isz))
    endif

    if (debug) write(0,*)' vdiag: ',prec%d(:)
    if (debug) write(0,*) 'Bi-CGSTAB with Additive Schwarz prec' 

    tx(1:desc_data%matrix_data(psb_n_row_)) = x(1:desc_data%matrix_data(psb_n_row_)) 
    tx(desc_data%matrix_data(psb_n_row_)+1:isz) = zero

    if (prec%iprcparm(restr_)==psb_halo_) then 
       call psb_halo(tx,prec%desc_data,info,work=aux)
      if(info /=0) then
         info=4010
         ch_err='psb_halo'
         goto 9999
      end if
    else if (prec%iprcparm(restr_) /= psb_none_) then 
      write(0,*) 'Problem in PRCAPLY: Unknown value for restriction ',&
           &prec%iprcparm(restr_)
    end if

    if (prec%iprcparm(iren_)>0) then 
      call psb_dgelp('N',n_row,1,prec%perm,tx,isz,ww,isz,info)
      if(info /=0) then
         info=4010
         ch_err='psb_dgelp'
         goto 9999
      end if
    endif

    call psb_dbjacaply(prec,tx,zero,ty,prec%desc_data,trans,aux,info)

    if (prec%iprcparm(iren_)>0) then 
      call psb_dgelp('N',n_row,1,prec%invperm,ty,isz,ww,isz,info)
      if(info /=0) then
         info=4010
         ch_err='psb_dgelp'
         goto 9999
      end if
    endif

    select case (prec%iprcparm(prol_)) 

    case(psb_none_) 
      ! Would work anyway, but since it's supposed to do nothing...
      ! call f90_psovrl(ty,prec%desc_data,update_type=prec%a_restrict)

    case(psb_sum_,psb_avg_) 
       call psb_ovrl(ty,prec%desc_data,info,&
            & update_type=prec%iprcparm(prol_),work=aux)
       if(info /=0) then
         info=4010
         ch_err='psb_ovrl'
         goto 9999
      end if

    case default
      write(0,*) 'Problem in PRCAPLY: Unknown value for prolongation ',&
           & prec%iprcparm(prol_)
    end select

    if (beta == zero)  then 
      y(1:desc_data%matrix_data(psb_n_row_)) = ty(1:desc_data%matrix_data(psb_n_row_)) 
    else if (beta == one) then 
      y(1:desc_data%matrix_data(psb_n_row_)) = y(1:desc_data%matrix_data(psb_n_row_)) + &
           & ty(1:desc_data%matrix_data(psb_n_row_)) 
    else if (beta == -one) then 
      y(1:desc_data%matrix_data(psb_n_row_)) = -y(1:desc_data%matrix_data(psb_n_row_)) + &
           & ty(1:desc_data%matrix_data(psb_n_row_)) 
    else 
      y(1:desc_data%matrix_data(psb_n_row_)) = beta*y(1:desc_data%matrix_data(psb_n_row_)) + &
           & ty(1:desc_data%matrix_data(psb_n_row_)) 
    end if


    if ((6*isz) <= size(work)) then 
    else if ((4*isz) <= size(work)) then 
      deallocate(ww,tx,ty)
    else if ((3*isz) <= size(work)) then 
      deallocate(aux)
    else 
      deallocate(ww,aux,tx,ty)
    endif

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

end subroutine psb_dbaseprcaply



subroutine psb_dbjacaply(prec,x,beta,y,desc_data,trans,work,info)
  !
  !  Compute   Y <-  beta*Y + K^-1 X 
  !  where K is a a Block Jacobi  preconditioner stored in prec
  !  Note that desc_data may or may not be the same as prec%desc_data,
  !  but since both are INTENT(IN) this should be legal. 
  ! 

  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_psblas_mod
  use psb_const_mod
  use psb_error_mod
  implicit none 

  type(psb_desc_type), intent(in)       :: desc_data
  type(psb_dbase_prec), intent(in)      :: prec
  real(kind(0.d0)),intent(inout)        :: x(:), y(:)
  real(kind(0.d0)),intent(in)           :: beta
  character(len=1)                      :: trans
  real(kind(0.d0)),target               :: work(:)
  integer, intent(out)                  :: info

  ! Local variables
  integer :: n_row,n_col
  real(kind(1.d0)), pointer :: ww(:), aux(:), tx(:),ty(:),tb(:)
  character     ::diagl, diagu
  integer :: icontxt,nprow,npcol,me,mycol,i, isz, nrg, err_act, int_err(5)
  real(kind(1.d0)) :: t1, t2, t3, t4, t5, t6, t7, mpi_wtime
  logical,parameter                 :: debug=.false., debugprt=.false.
  real(kind(1.d0)), parameter       :: one=1.d0, zero=0.d0
  external mpi_wtime
  character(len=20)   :: name, ch_err

  name='psb_dbjacaply'
  info = 0
  call psb_erractionsave(err_act)

  icontxt=desc_data%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,nprow,npcol,me,mycol)

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
      allocate(aux(4*n_col))
    endif
  else
    allocate(ww(n_col),aux(4*n_col))
  endif


  if (prec%iprcparm(jac_sweeps_) == 1) then 


    select case(prec%iprcparm(f_type_))
    case(f_ilu_n_,f_ilu_e_) 

      select case(trans)
      case('N','n')

        call psb_spsm(one,prec%av(l_pr_),x,zero,ww,desc_data,info,&
             & trans='N',unit=diagl,choice=psb_none_,work=aux)
        if(info /=0) goto 9999
        ww(1:n_row) = ww(1:n_row)*prec%d(1:n_row)
        call psb_spsm(one,prec%av(u_pr_),ww,beta,y,desc_data,info,&
             & trans='N',unit=diagu,choice=psb_none_, work=aux)
        if(info /=0) goto 9999

      case('T','t','C','c')
        call psb_spsm(one,prec%av(u_pr_),x,zero,ww,desc_data,info,&
             & trans=trans,unit=diagu,choice=psb_none_, work=aux)
        if(info /=0) goto 9999
        ww(1:n_row) = ww(1:n_row)*prec%d(1:n_row)
        call psb_spsm(one,prec%av(l_pr_),ww,beta,y,desc_data,info,&
             & trans=trans,unit=diagl,choice=psb_none_,work=aux)
        if(info /=0) goto 9999

      end select

    case(f_slu_)

      ww(1:n_row) = x(1:n_row)

      select case(trans)
      case('N','n')
        call fort_slu_solve(0,n_row,1,ww,n_row,prec%iprcparm(slu_ptr_),info)
      case('T','t','C','c')
        call fort_slu_solve(1,n_row,1,ww,n_row,prec%iprcparm(slu_ptr_),info)
      end select

      if(info /=0) goto 9999

      if (beta == 0.d0) then 
        y(1:n_row) = ww(1:n_row)
      else if (beta==1.d0) then 
        y(1:n_row) = ww(1:n_row) + y(1:n_row) 
      else if (beta==-1.d0) then 
        y(1:n_row) = ww(1:n_row) - y(1:n_row) 
      else 
        y(1:n_row) = ww(1:n_row) + beta*y(1:n_row) 
      endif

    case default
      write(0,*) 'Unknown factorization type in bjac_prcaply',prec%iprcparm(f_type_)
    end select
    if (debug) write(0,*)' Y: ',y(:)

  else if (prec%iprcparm(jac_sweeps_) > 1) then 

    ! Note: we have to add TRANS to this one !!!!!!!!! 

    if (size(prec%av) < ap_nd_) then 
      info = 4011
      goto 9999
    endif

    allocate(tx(n_col),ty(n_col))
    tx = zero
    ty = zero
    select case(prec%iprcparm(f_type_)) 
    case(f_ilu_n_,f_ilu_e_) 
      do i=1, prec%iprcparm(jac_sweeps_) 
        !   X(k+1) = M^-1*(b-N*X(k))
        ty(1:n_row) = x(1:n_row)
        call psb_spmm(-one,prec%av(ap_nd_),tx,one,ty,&
             &   prec%desc_data,info,work=aux)
        if(info /=0) goto 9999
        call psb_spsm(one,prec%av(l_pr_),ty,zero,ww,&
             & prec%desc_data,info,&
             & trans='N',unit='U',choice=psb_none_,work=aux)
        if(info /=0) goto 9999
        ww(1:n_row) = ww(1:n_row)*prec%d(1:n_row)
        call psb_spsm(one,prec%av(u_pr_),ww,zero,tx,&
             & prec%desc_data,info,&
             & trans='N',unit='U',choice=psb_none_,work=aux)
        if(info /=0) goto 9999
      end do

    case(f_slu_) 
      do i=1, prec%iprcparm(jac_sweeps_) 
        !   X(k+1) = M^-1*(b-N*X(k))
        ty(1:n_row) = x(1:n_row)
        call psb_spmm(-one,prec%av(ap_nd_),tx,one,ty,&
             &   prec%desc_data,info,work=aux)
        if(info /=0) goto 9999

        call fort_slu_solve(0,n_row,1,ty,n_row,prec%iprcparm(slu_ptr_),info)
        if(info /=0) goto 9999
        tx(1:n_row) = ty(1:n_row)        
      end do
    end select
    
    if (beta == 0.d0) then 
      y(1:n_row) = tx(1:n_row)
    else if (beta==1.d0) then 
      y(1:n_row) = tx(1:n_row) + y(1:n_row) 
    else if (beta==-1.d0) then 
      y(1:n_row) = tx(1:n_row) - y(1:n_row) 
    else 
      y(1:n_row) = tx(1:n_row) + beta*y(1:n_row) 
    endif

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
  if (err_act.eq.act_abort) then
    call psb_error()
    return
  end if
  return

end subroutine psb_dbjacaply


subroutine psb_dmlprcaply(baseprecv,x,beta,y,desc_data,trans,work,info)
  !
  !  Compute   Y <-  beta*Y + K^-1 X 
  !  where K is a multilevel (actually 2-level) preconditioner stored in prec
  ! 

  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_psblas_mod
  use psb_blacs_mod
  use psb_const_mod
  use psb_error_mod
  implicit none

  type(psb_desc_type),intent(in)    :: desc_data
  type(psb_dbase_prec), intent(in)  :: baseprecv(:)
  real(kind(0.d0)),intent(in)       :: beta
  real(kind(0.d0)),intent(inout)    :: x(:), y(:)
  character                         :: trans
  real(kind(0.d0)),target           :: work(:)
  integer, intent(out)              :: info


  ! Local variables
  integer :: n_row,n_col
  real(kind(1.d0)), allocatable :: tx(:),ty(:),t2l(:),w2l(:),&
       &   x2l(:),b2l(:),tz(:),tty(:)
  character     ::diagl, diagu
  integer :: icontxt,nprow,npcol,me,mycol,i, isz, nrg,nr2l,err_act, iptype, int_err(5)
  real(kind(1.d0)) :: omega
  real(kind(1.d0)) :: t1, t2, t3, t4, t5, t6, t7, mpi_wtime
  logical, parameter          :: debug=.false., debugprt=.false.
  real(kind(1.d0)), parameter :: one=1.d0, zero=0.d0
  integer      :: ismth
  external mpi_wtime
  character(len=20)   :: name, ch_err

  interface
     subroutine psb_dbaseprcaply(prec,x,beta,y,desc_data,trans,work,info)
       use psb_descriptor_type
       use psb_prec_type
       type(psb_desc_type),intent(in)    :: desc_data
       type(psb_dbase_prec), intent(in)  :: prec
       real(kind(0.d0)),intent(inout)    :: x(:), y(:)
       real(kind(0.d0)),intent(in)       :: beta
       character(len=1)                  :: trans
       real(kind(0.d0)),target           :: work(:)
       integer, intent(out)              :: info
     end subroutine psb_dbaseprcaply
  end interface

  name='psb_dmlprcaply'
  info = 0
  call psb_erractionsave(err_act)


  icontxt=desc_data%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,nprow,npcol,me,mycol)

  omega=baseprecv(2)%dprcparm(smooth_omega_)
  ismth=baseprecv(2)%iprcparm(smth_kind_)

  select case(baseprecv(2)%iprcparm(ml_type_)) 
  case(no_ml_) 
    ! Should not really get here.
    write(0,*) 'Smooth preconditioner with no multilevel in MLPRCAPLY????' 

  case(add_ml_prec_)

    ! 
    !  Additive multilevel
    !
    t1 = mpi_wtime()
    n_row = desc_data%matrix_data(psb_n_row_)
    n_col = baseprecv(1)%desc_data%matrix_data(psb_n_col_)
    call psb_dbaseprcaply(baseprecv(1),x,beta,y,desc_data,trans,work,info)
    if(info /=0) goto 9999

    nr2l  = baseprecv(2)%desc_data%matrix_data(psb_n_col_)
    nrg   = baseprecv(2)%desc_data%matrix_data(psb_n_row_)
    allocate(t2l(nr2l),w2l(nr2l))
    t2l(:) = zero
    w2l(:) = zero

    if (ismth  /= no_smth_) then 
      !
      ! Smoothed aggregation
      !
      allocate(tx(max(n_row,n_col)),ty(max(n_row,n_col)),&
           & tz(max(n_row,n_col)))
      tx(1:desc_data%matrix_data(psb_n_row_)) = x(1:desc_data%matrix_data(psb_n_row_)) 
      tx(desc_data%matrix_data(psb_n_row_)+1:max(n_row,n_col)) = zero
      ty(desc_data%matrix_data(psb_n_row_)+1:max(n_row,n_col)) = zero
      tz(desc_data%matrix_data(psb_n_row_)+1:max(n_row,n_col)) = zero


      if (baseprecv(2)%iprcparm(glb_smth_) >0) then 
        call psb_halo(tx,desc_data,info,work=work) 
        if(info /=0) goto 9999
      else
        tx(desc_data%matrix_data(psb_n_row_)+1:max(n_row,n_col)) = zero
      end if

      call psb_csmm(one,baseprecv(2)%av(sm_pr_t_),tx,zero,t2l,info)
      if(info /=0) goto 9999

    else
      !
      ! Raw  aggregation, may take shortcut
      !
      do i=1,desc_data%matrix_data(psb_n_row_)
        t2l(baseprecv(2)%mlia(i)) = t2l(baseprecv(2)%mlia(i)) + x(i)
      end do

    end if

    if (baseprecv(2)%iprcparm(coarse_mat_)==mat_repl_) Then 
      call gsum2d(icontxt,'All',t2l(1:nrg))
    else if (baseprecv(2)%iprcparm(coarse_mat_) /= mat_distr_) Then 
        write(0,*) 'Unknown value for baseprecv(2)%iprcparm(coarse_mat_) ',&
             & baseprecv(2)%iprcparm(coarse_mat_)
    endif

    w2l=t2l
    call psb_dbaseprcaply(baseprecv(2),w2l,zero,t2l,baseprecv(2)%desc_data,'N',work,info)


    if (ismth  /= no_smth_) then 

      call psb_csmm(one,baseprecv(2)%av(sm_pr_),t2l,zero,ty,info)
      if(info /=0) goto 9999
      ! 
      ! Finally add back into Y. 
      ! 
      call psb_axpby(one,ty,one,y,desc_data,info)
      if(info /=0) goto 9999
      deallocate(tx,ty,tz)

    else

      do i=1, desc_data%matrix_data(psb_n_row_)
        y(i) = y(i) + t2l(baseprecv(2)%mlia(i))
      enddo

    end if

    if (debug) write(0,*)' Y2: ',Y(:)

    deallocate(t2l,w2l)

  case(mult_ml_prec_)

    ! 
    !  Multiplicative multilevel
    !  Pre/post smoothing versions. 

    select case(baseprecv(2)%iprcparm(smth_pos_))

    case(post_smooth_)


      t1    = mpi_wtime()
      n_row = desc_data%matrix_data(psb_n_row_)
      n_col = baseprecv(1)%desc_data%matrix_data(psb_n_col_)
      nr2l  = baseprecv(2)%desc_data%matrix_data(psb_n_col_)
      nrg   = baseprecv(2)%desc_data%matrix_data(psb_n_row_)
      allocate(t2l(nr2l),w2l(nr2l),tx(n_col),ty(n_col))
      t2l(:) = zero
      w2l(:) = zero

      !
      ! Need temp copies to handle Y<- betaY + K^-1 X
      ! One of the temp copies is not strictly needed when beta==zero
      !

      if (debug) write(0,*)' mult_ml_apply  omega ',omega
      if (debug) write(0,*)' mult_ml_apply  X: ',X(:)
      call psb_axpby(one,x,zero,tx,desc_data,info)
      if(info /=0) goto 9999

      if (ismth  /= no_smth_) then 
        !
        ! Smoothed aggregation
        !
        allocate(tz(max(n_row,n_col)))

        if (baseprecv(2)%iprcparm(glb_smth_) >0) then 
          call psb_halo(tx,desc_data,info,work=work) 
          if(info /=0) goto 9999
        else
          tx(desc_data%matrix_data(psb_n_row_)+1:max(n_row,n_col)) = zero
        end if

        call psb_csmm(one,baseprecv(2)%av(sm_pr_t_),tx,zero,t2l,info)
        if(info /=0) goto 9999

      else
        !
        ! Raw  aggregation, may take shortcut
        !
        do i=1,desc_data%matrix_data(psb_n_row_)
          t2l(baseprecv(2)%mlia(i)) = t2l(baseprecv(2)%mlia(i)) + tx(i)
        end do
      end if

      if (baseprecv(2)%iprcparm(coarse_mat_)==mat_repl_) Then 
        call gsum2d(icontxt,'All',t2l(1:nrg))
      else if (baseprecv(2)%iprcparm(coarse_mat_) /= mat_distr_) Then 
        write(0,*) 'Unknown value for baseprecv(2)%iprcparm(coarse_mat_) ',&
             & baseprecv(2)%iprcparm(coarse_mat_)
      endif

      t6 = mpi_wtime()
      w2l=t2l
      call psb_dbaseprcaply(baseprecv(2),w2l,zero,t2l,baseprecv(2)%desc_data,'N',work,info)
      if(info /=0) goto 9999

      if (ismth  /= no_smth_) then 
        if (ismth == smth_omg_) &
             & call psb_halo(t2l,baseprecv(2)%desc_data,info,work=work) 
        call psb_csmm(one,baseprecv(2)%av(sm_pr_),t2l,zero,ty,info)
        if(info /=0) goto 9999
        ! 
        ! Finally add back into Y. 
        ! 
        deallocate(tz)
      else
        ty(:) = zero
        do i=1, desc_data%matrix_data(psb_n_row_)
          ty(i) = ty(i) + t2l(baseprecv(2)%mlia(i))
        enddo

      end if
      deallocate(t2l,w2l)

      call psb_spmm(-one,baseprecv(2)%aorig,ty,one,tx,desc_data,info,work=work)
      if(info /=0) goto 9999

      call psb_dbaseprcaply(baseprecv(1),tx,one,ty,desc_data,trans,work,info)
      if(info /=0) goto 9999

      call psb_axpby(one,ty,beta,y,desc_data,info)
      if(info /=0) goto 9999

      deallocate(tx,ty)



    case(pre_smooth_)

      t1 = mpi_wtime()
      n_row = desc_data%matrix_data(psb_n_row_)
      n_col = baseprecv(1)%desc_data%matrix_data(psb_n_col_)
      nr2l  = baseprecv(2)%desc_data%matrix_data(psb_n_col_)
      nrg   = baseprecv(2)%desc_data%matrix_data(psb_n_row_)
      allocate(t2l(nr2l),w2l(nr2l),tx(n_col),ty(n_col),tty(n_col))
      t2l(:) = zero
      w2l(:) = zero

      !
      ! Need temp copies to handle Y<- betaY + K^-1 X
      ! One of the temp copies is not strictly needed when beta==zero
      !
      call psb_axpby(one,x,zero,tx,desc_data,info)
      call psb_axpby(one,y,zero,ty,desc_data,info)
      if(info /=0) goto 9999

      call psb_dbaseprcaply(baseprecv(1),x,zero,tty,desc_data,trans,work,info)
      if(info /=0) goto 9999

      call psb_spmm(-one,baseprecv(2)%aorig,tty,one,tx,desc_data,info,work=work)
      if(info /=0) goto 9999

      if (ismth  /= no_smth_) then 
        allocate(tz(max(n_row,n_col)))

        if (baseprecv(2)%iprcparm(glb_smth_) >0) then 
          call psb_halo(tx,desc_data,info,work=work) 
          if(info /=0) goto 9999
        else
          tx(desc_data%matrix_data(psb_n_row_)+1:max(n_row,n_col)) = zero
        end if

        call psb_csmm(one,baseprecv(2)%av(sm_pr_t_),tx,zero,t2l,info)
        if(info /=0) goto 9999

      else
        !
        ! Raw  aggregation, may take shortcuts
        !
        do i=1,desc_data%matrix_data(psb_n_row_)
          t2l(baseprecv(2)%mlia(i)) = t2l(baseprecv(2)%mlia(i)) + tx(i)
        end do
      end if

      if (baseprecv(2)%iprcparm(coarse_mat_)==mat_repl_) Then 
        call gsum2d(icontxt,'All',t2l(1:nrg))
      else if (baseprecv(2)%iprcparm(coarse_mat_) /= mat_distr_) Then 
        write(0,*) 'Unknown value for baseprecv(2)%iprcparm(coarse_mat_) ',&
             & baseprecv(2)%iprcparm(coarse_mat_)
      endif

      t6 = mpi_wtime()
      w2l=t2l
      call psb_dbaseprcaply(baseprecv(2),w2l,zero,t2l,baseprecv(2)%desc_data,'N',work,info)
      if(info /=0) goto 9999

      if (ismth  /= no_smth_) then 

        if (ismth == smth_omg_) &
             & call psb_halo(t2l,baseprecv(2)%desc_data,info,work=work) 
        call psb_csmm(one,baseprecv(2)%av(sm_pr_),t2l,zero,ty,info)
        if(info /=0) goto 9999

        call psb_axpby(one,ty,one,tty,desc_data,info)
        if(info /=0) goto 9999

        deallocate(tz)
      else

        do i=1, desc_data%matrix_data(psb_n_row_)
          tty(i) = tty(i) + t2l(baseprecv(2)%mlia(i))
        enddo

      end if

      call psb_axpby(one,tty,beta,y,desc_data,info)
      if(info /=0) goto 9999

      deallocate(t2l,w2l,tx,ty,tty)




    case default

      write(0,*) 'Unknown value for ml_smooth_pos',baseprecv(2)%iprcparm(smth_pos_)

    end select

  case default
    write(0,*) me, 'Wrong mltype into PRCAPLY ',&
         & baseprecv(2)%iprcparm(ml_type_)
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error()
    return
  end if
  return

end subroutine psb_dmlprcaply


subroutine psb_dprecaply1(prec,x,desc_data,info,trans)

  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_psblas_mod
  use psb_const_mod
  use psb_error_mod
  implicit none

  type(psb_desc_type),intent(in)    :: desc_data
  type(psb_dprec_type), intent(in)  :: prec
  real(kind(0.d0)),intent(inout)    :: x(:)
  integer, intent(out)              :: info
  character(len=1), optional        :: trans
  logical,parameter                 :: debug=.false., debugprt=.false.
  real(kind(1.d0)), parameter       :: one=1.d0, zero=0.d0


  ! Local variables
  character     :: trans_
  integer :: icontxt,nprow,npcol,me,mycol,i, isz, err_act, int_err(5)
  real(kind(1.d0)), pointer :: WW(:), w1(:)
  character(len=20)   :: name, ch_err
  name='psb_dprec1'
  info = 0
  call psb_erractionsave(err_act)
  

  icontxt=desc_data%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,nprow,npcol,me,mycol)
  if (present(trans)) then 
    trans_=trans
  else
    trans_='N'
  end if

  allocate(ww(size(x)),w1(size(x)))
  call psb_dprecaply(prec,x,ww,desc_data,info,trans_,w1)
  if(info /=0) goto 9999
  x(:) = ww(:)
  deallocate(ww,W1)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error()
     return
  end if
  return
end subroutine psb_dprecaply1

======= end
