!*****************************************************************************
!*                                                                           *
!* This routine does two things:                                             *
!*    1. Builds the auxiliary descriptor. This is always done even for       *
!*       Block Jacobi.                                                       *
!*    2. Retrieves the remote matrix pieces.                                 *
!*                                                                           *
!*    All of 1. is done under f90_dscov, which is independent of CSR, and    *
!*    has been placed in the TOOLS directory because it might be used for    *
!*    building a descriptor for an extended stencil in a PDE solver without  *
!*    necessarily applying AS precond.                                       *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*****************************************************************************
Subroutine psb_dcsrsetup(ptype,novr,a,blk,desc_data,upd,desc_p,info,outfmt)

  use psb_serial_mod
  use psb_descriptor_type
  Use psb_prec_type
  use psb_tools_mod
  use psb_const_mod
  use psb_error_mod
  Implicit None

  !     .. Array Arguments ..
  integer, intent(in)                  :: ptype,novr
  Type(psb_dspmat_type), Intent(in)    ::  a
  Type(psb_dspmat_type), Intent(inout) ::  blk
  integer, intent(out)                 :: info
  Type(psb_desc_type), Intent(inout)   :: desc_p
  Type(psb_desc_type), Intent(in)      :: desc_data 
  Character, Intent(in)                :: upd
  character(len=5), optional           :: outfmt


  real(kind(1.d0)) :: t1,t2,t3,mpi_wtime
  external  mpi_wtime
  integer   idscb,idsce,iovrb,iovre, err, irank, icomm

  !     .. Local Scalars ..
  Integer ::  k, tot_elem,proc,&
       &  point,nprow,npcol, me, mycol, start,m,nnzero,&
       &  icontxt, lovr, n_col, linp,ier,n,int_err(5),&
       &  tot_recv, ircode, n_row, nztot,nhalo, nrow_a,err_act
  Logical,Parameter :: debug=.false., debugprt=.false.
  character(len=20) :: name, ch_err
  name='psb_dcsrsetup'
  info=0
  call psb_erractionsave(err_act)

  If(debug) Write(0,*)'IN DCSRSETUP  ', upd
  icontxt=desc_data%matrix_data(psb_ctxt_)
  tot_recv=0

  nrow_a = desc_data%matrix_data(psb_n_row_)
  nnzero = Size(a%aspk)
  n_col  = desc_data%matrix_data(psb_n_col_)
  nhalo  = n_col-nrow_a


  If (ptype == bja_) Then
    !
    ! Block Jacobi. Copy the descriptor, just in case we want to
    ! do the renumbering. 
    !
    call psb_spall(0,0,blk,1,info)
    if(info /= 0) then
      info=4010
      ch_err='psb_spall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    blk%fida        = 'COO'
    blk%infoa(psb_nnz_) = 0

    If (upd == 'F') Then
      call psb_dsccpy(desc_p,desc_data,info)
      if(info /= 0) then
        info=4010
        ch_err='psb_ddsccpy'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    endif

  Else If (ptype == asm_) Then


    !
    ! Additive Schwarz variant. 
    !
    !

    icontxt=desc_data%matrix_data(psb_ctxt_)

    if (novr < 0) then
      info=3
      int_err(1)=novr
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    endif

    if (novr == 0) then 
      !
      ! This is really just Block Jacobi.....
      !
      call psb_spall(0,0,blk,1,info)
      if(info /= 0) then
        info=4010
        ch_err='psb_spall'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
      blk%fida='COO'
      blk%infoa(psb_nnz_)=0
      if (debug) write(0,*) 'Calling desccpy'
      if (upd == 'F') then 
        call psb_dsccpy(desc_p,desc_data,info)
        if(info /= 0) then
          info=4010
          ch_err='psb_dsccpy'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if
        if (debug) write(0,*) 'Early return from dcsrsetup: P>=3 N_OVR=0'
      endif
      return
    endif

    call blacs_get(icontxt,10,icomm )
!!$    call MPI_Comm_rank(icomm,irank,ierr)
!!$    idscb  = mpe_log_get_event_number()
!!$    idsce  = mpe_log_get_event_number()
!!$    iovrb  = mpe_log_get_event_number()
!!$    iovre  = mpe_log_get_event_number()
!!$    if (irank==0) then 
!!$      info = mpe_describe_state(idscb,idsce,"DSCASB ","NavyBlue")
!!$      info = mpe_describe_state(iovrb,iovre,"DSCOVR ","DeepPink")
!!$    endif
!!$

    Call blacs_gridinfo(icontxt,nprow,npcol,me,mycol)
    If(debug)Write(0,*)'BEGIN dcsrsetup',me,upd
    t1 = mpi_wtime()


    If (upd == 'F') Then
      !
      !  Build the  auiliary descriptor',desc_p%matrix_data(psb_n_row_)
      ! 
      call psb_dscov(a,desc_data,novr,desc_p,info)
      if(info /= 0) then
        info=4010
        ch_err='psb_dscov'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    Endif

    if(debug) write(0,*) me,' From dscov _:',desc_p%matrix_data(psb_n_row_),desc_p%matrix_data(psb_n_col_)


    n_row = desc_p%matrix_data(psb_n_row_)
    t2 = mpi_wtime()

    if (debug) write(0,*) 'Before dcsrovr ',blk%fida,blk%m,psb_nnz_,blk%infoa(psb_nnz_)
!!$    ierr = MPE_Log_event( iovrb, 0, "st OVR" )
!!$    blk%m = n_row-nrow_a
!!$    blk%k = n_row

    if (present(outfmt)) then 
      if(debug) write(0,*) me,': Calling CSROVR with ',size(blk%ia2)
      Call psb_csrovr(a,desc_p,blk,info,outfmt=outfmt)
    else
      if(debug) write(0,*) me,': Calling CSROVR with ',size(blk%ia2)
      Call psb_csrovr(a,desc_p,blk,info)
    end if


    if(info /= 0) then
      info=4010
      ch_err='psb_csrovr'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if (debug) write(0,*) 'After psb_dcsrovr ',blk%fida,blk%m,psb_nnz_,blk%infoa(psb_nnz_)
!!$    ierr = MPE_Log_event( iovre, 0, "ed OVR" )

    t3 = mpi_wtime()
    if (debugprt) then 
      open(40+me) 
      call psb_csprt(40+me,blk,head='% Ovrlap rows')
      close(40+me)
    endif


  End If

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error()
    return
  end if
  Return

End Subroutine psb_dcsrsetup

