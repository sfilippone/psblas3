! File: psb_dspmm.f90
!
! Subroutine: psb_dspmm
!     Performs one of the distributed matrix-vector operations
!
!     sub( Y ) := alpha * Pr * A * Pc * sub( X )  + beta * sub( Y ),  or
!
!     sub( Y ) := alpha * Pr * A' * Pr * sub( X )  + beta * sub( Y ),
!
!     where:
!
!        sub( X ) denotes *if* TRANS = 'N',
!
!                       X(1:N,JX:JX+K-1),
!
!                     *else*
!
!                       X(1:M,JX:JX+K-1).
!
!                     *end if*
!
!        sub( Y ) denotes *if* trans = 'N',
!
!                       Y(1:M,JY:JY+K-1),
!
!                     *else*
!
!                       Y(1:N,JY:JY+K-1)
!
!                     *end* *if*
!
!  alpha and beta are scalars, and sub( X ) and sub( Y ) are distributed
!  vectors and A is a M-by-N distributed matrix.
!
! Parameters:   
!    alpha  -  real.                        The scalar alpha.
!    a      -  type(<psb_dspmat_type>).     The sparse matrix containing A.
!    x      -  real,dimension(:,:).         The input vector containing the entries of sub( X ).
!    beta   -  real.                        The scalar beta.
!    y      -  real,dimension(:,:).         The input vector containing the entries of sub( Y ).
!    desc_a -  type(<psb_desc_type>).       The communication descriptor.
!    info   -  integer.                     Eventually returns an error code.
!    trans  -  character(optional).         Whether A or A'. If not present 'N' is assumed.
!    k      -  integer(optional).           The number of right-hand sides.
!    jx     -  integer(optional).           The column offset for sub( X ). If not present 1 is assumed.
!    jy     -  integer(optional).           The column offset for sub( Y ). If not present 1 is assumed.
!    work   -  real,dimension(:)(optional). Working area.
!    doswap -  integer(optional).           Whether to performe halo updates.
! 
subroutine  psb_dspmm(alpha,a,x,beta,y,desc_a,info,&
     & trans, k, jx, jy, work, doswap)   

  use psb_spmat_type
  use psb_serial_mod
  use psb_descriptor_type
  use psb_comm_mod
  use psi_mod
  use psb_error_mod
  implicit none

  real(kind(1.D0)), intent(in)             :: alpha, beta
  real(kind(1.d0)), intent(inout), target  :: x(:,:)
  real(kind(1.d0)), intent(inout), target  :: y(:,:)
  type(psb_dspmat_type), intent(in)        :: a
  type(psb_desc_type), intent(in)          :: desc_a
  integer, intent(out)                     :: info
  real(kind(1.d0)), optional, pointer      :: work(:)
  character, intent(in), optional          :: trans
  integer, intent(in), optional            :: k, jx, jy,doswap

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, n, iix, jjx, ia, ja, iia, jja, temp(2), ix, iy, ik, ijx, ijy,&
       & idoswap, m, nrow, ncol, lldx, lldy, liwork, llwork, iiy, jjy,&
       & i, ib, ib1
  integer, parameter       :: nb=4
  real(kind(1.d0)),pointer :: tmpx(:), xp(:,:), yp(:,:), iwork(:)
  character                :: itrans
  character(len=20)        :: name, ch_err

  name='psb_dspmm'
  info=0
  call psb_erractionsave(err_act)

  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
  if (nprow == -1) then
     info = 2010
     call psb_errpush(info,name)
     goto 9999
  else if (npcol /= 1) then
     info = 2030
     int_err(1) = npcol
     call psb_errpush(info,name)
     goto 9999
  endif

  ia = 1
  ja = 1

  ix = 1
  if (present(jx)) then
     ijx = jx
  else
     ijx = 1
  endif

  iy = 1
  if (present(jy)) then
     ijy = jy
  else
     ijy = 1
  endif

  if (present(doswap)) then
     idoswap = doswap
  else
     idoswap = 1
  endif

  if (present(k)) then     
     ik = min(k,size(x,2)-ijx+1)
     ik = min(ik,size(y,2)-ijy+1)
  else
     ik = min(size(x,2)-ijx+1,size(y,2)-ijy+1)
  endif

  if (present(trans)) then     
     if((trans.eq.'N').or.(trans.eq.'T')) then
        itrans = trans
     else if (trans.eq.'C') then
        info = 3020
        call psb_errpush(info,name)
        goto 9999
     else
        info = 70
        call psb_errpush(info,name)
        goto 9999
     end if
  else
     itrans = 'N'
  endif

  m    = desc_a%matrix_data(psb_m_)
  n    = desc_a%matrix_data(psb_n_)
  nrow = desc_a%matrix_data(psb_n_row_)
  ncol = desc_a%matrix_data(psb_n_col_)
  lldx = size(x,1)
  lldy = size(y,1)

  ! check for presence/size of a work area
  liwork= 2*ncol
  if (a%pr(1) /= 0) llwork = liwork + n * ik
  if (a%pl(1) /= 0) llwork = llwork + m * ik
  if (present(work)) then     
     if(size(work).lt.liwork) then
        call psb_realloc(liwork,work,info)
        if(info.ne.0) then
           info=4010
           ch_err='psb_realloc'
           call psb_errpush(info,name,a_err=ch_err)
           goto 9999
        end if
     end if
     iwork => work
  else
     call psb_realloc(liwork,iwork,info)
     if(info.ne.0) then
        info=4010
        ch_err='psb_realloc'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if
  end if
  iwork(1)=0.d0

  ! checking for matrix correctness
  call psb_chkmat(m,n,ia,ja,desc_a%matrix_data,info,iia,jja)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkmat'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if


  if (itrans.eq.'N') then
     !  Matrix is not transposed
     if((ja.ne.ix).or.(ia.ne.iy)) then
        ! this case is not yet implemented
        info = 3030
        call psb_errpush(info,name)
        goto 9999
     end if

     ! checking for vectors correctness
     call psb_chkvect(n,ik,size(x,1),ix,ijx,desc_a%matrix_data,info,iix,jjx)
     call psb_chkvect(m,ik,size(y,1),iy,ijy,desc_a%matrix_data,info,iiy,jjy)
     if(info.ne.0) then
        info=4010
        ch_err='psb_chkvect'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

     if((iix.ne.1).or.(iiy.ne.1)) then
        ! this case is not yet implemented
        info = 3040
        call psb_errpush(info,name)
        goto 9999
     end if

     if(idoswap.lt.0) x(nrow:ncol,1:ik)=0.d0

     ib1=min(nb,ik)
     xp => x(iix:lldx,jjx:jjx+ib1-1)
     if(idoswap.gt.0)&
          & call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
          & ib1,dzero,xp,desc_a,iwork,info)
!!$          & call PSI_dSwapData(ior(SWAP_SEND,SWAP_RECV),ib1,&
!!$          & dzero,x(iix,jjx),lldx,desc_a%matrix_data,&
!!$          & desc_a%halo_index,iwork,liwork,info)
     

     blk: do i=1, ik, nb
        ib=ib1
        ib1 = max(0,min(nb,(ik)-(i-1+ib)))
        xp => x(iix:lldx,jjx+i+ib-1:jjx+i+ib+ib1-2)
        if((ib1.gt.0).and.(idoswap.gt.0))&
             & call psi_swapdata(psb_swap_send_,ib1,&
             & dzero,xp,desc_a,iwork,info)
!!$             & call PSI_dSwapData(SWAP_SEND,ib1,&
!!$             & dzero,x(iix,jjx+i+ib-1),lldx,desc_a%matrix_data,&
!!$             & desc_a%halo_index,iwork,liwork,info)
        if(info.ne.0) exit blk
        
        !  local Matrix-vector product
        call dcsmm(itrans,nrow,ib,ncol,alpha,a%pr,a%fida,&
             & a%descra,a%aspk,a%ia1,a%ia2,a%infoa,a%pl,&
             & x(iix,jjx+i-1),lldx,beta,y(iiy,jjy+i-1),lldy,&
             & iwork,liwork,info)
        if(info.ne.0) exit blk
        
        if((ib1.gt.0).and.(idoswap.gt.0))&
             & call psi_swapdata(psb_swap_send_,ib1,&
             & dzero,xp,desc_a,iwork,info)
!!$             & call PSI_dSwapData(SWAP_RECV,ib1,&
!!$             & dzero,x(iix,jjx+i+ib-1),lldx,desc_a%matrix_data,&
!!$             & desc_a%halo_index,iwork,liwork,info)
        if(info.ne.0) exit blk
     end do blk

     if(info.ne.0) then
        info = 4011
        call psb_errpush(info,name)
        goto 9999
     end if

  else
     !  Matrix is transposed
     if((ja.ne.iy).or.(ia.ne.ix)) then
        ! this case is not yet implemented
        info = 3030
        call psb_errpush(info,name)
        goto 9999
     end if

     if(desc_a%ovrlap_elem(1).ne.-1) then
        info = 3070
        call psb_errpush(info,name)
        goto 9999
     end if

     ! checking for vectors correctness
     call psb_chkvect(m,ik,size(x,1),ix,ijx,desc_a%matrix_data,info,iix,jjx)
     call psb_chkvect(n,ik,size(y,1),iy,ijy,desc_a%matrix_data,info,iiy,jjy)
     if(info.ne.0) then
        info=4010
        ch_err='psb_chkvect'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

     if((iix.ne.1).or.(iiy.ne.1)) then
        ! this case is not yet implemented
        info = 3040
        call psb_errpush(info,name)
        goto 9999
     end if

     if(idoswap.lt.0) y(nrow:ncol,1:ik)=0.d0

     !  local Matrix-vector product
     call dcsmm(itrans,ncol,ik,nrow,alpha,a%pr,a%fida,&
          & a%descra,a%aspk,a%ia1,a%ia2,a%infoa,a%pl,&
          & x(iix,jjx),lldx,beta,y(iiy,jjy),lldy,&
          & iwork,liwork,info)
     if(info.ne.0) then
        info = 4010
        ch_err='dcsmm'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

     yp => y(iiy:lldy,jjy:jjy+ik-1)
     if(idoswap.gt.0)&
          & call psi_swaptran(ior(psb_swap_send_,psb_swap_recv_),&
          & ik,done,yp,desc_a,iwork,info)
!!$          & call PSI_dSwapTran(ior(SWAP_SEND,SWAP_RECV),&
!!$          & ik,done,y(iiy,jjy),lldy,desc_a%matrix_data,&
!!$          & desc_a%halo_index,iwork,liwork,info)
     if(info.ne.0) then
        info = 4010
        ch_err='PSI_dSwapTran'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

  end if

  if(.not.present(work)) deallocate(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end subroutine psb_dspmm




! Subroutine: psb_dspmmv
!     Performs one of the distributed matrix-vector operations
!
!     Y := alpha * Pr * A * Pc * X  + beta * Y,  or
!
!     Y := alpha * Pr * A' * Pr * X  + beta * Y,
!
!  alpha and beta are scalars, and X and Y are distributed
!  vectors and A is a M-by-N distributed matrix.
!
! Parameters:   
!    alpha  -  real.                        The scalar alpha.
!    a      -  type(<psb_dspmat_type>).     The sparse matrix containing A.
!    x      -  real,dimension(:).           The input vector containing the entries of X.
!    beta   -  real.                        The scalar beta.
!    y      -  real,dimension(:.         The input vector containing the entries of Y.
!    desc_a -  type(<psb_desc_type>).       The communication descriptor.
!    info   -  integer.                     Eventually returns an error code.
!    trans  -  character(optional).         Whether A or A'. If not present 'N' is assumed.
!    work   -  real,dimension(:)(optional). Working area.
!    doswap -  integer(optional).           Whether to performe halo updates.
! 
subroutine  psb_dspmv(alpha,a,x,beta,y,desc_a,info,&
     & trans, work, doswap)   

  use psb_spmat_type
  use psb_serial_mod
  use psb_descriptor_type
  use psb_comm_mod
  use psi_mod
  use psb_error_mod
  implicit none

  real(kind(1.D0)), intent(in)             :: alpha, beta
  real(kind(1.d0)), intent(inout), target  :: x(:)
  real(kind(1.d0)), intent(inout), target  :: y(:)
  type(psb_dspmat_type), intent(in)        :: a
  type(psb_desc_type), intent(in)          :: desc_a
  integer, intent(out)                     :: info
  real(kind(1.d0)), optional, pointer      :: work(:)
  character, intent(in), optional          :: trans
  integer, intent(in), optional            :: doswap

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, n, iix, jjx, ia, ja, iia, jja, temp(2), ix, iy, ik, ijx, ijy,&
       & idoswap, m, nrow, ncol, lldx, lldy, liwork, llwork, jx, jy, iiy, jjy,&
       & i, ib, ib1
  integer, parameter       :: nb=4
  real(kind(1.d0)),pointer :: tmpx(:), iwork(:), xp(:), yp(:)
  character                :: itrans
  character(len=20)        :: name, ch_err

  name='psb_dspmv'
  info=0
  call psb_erractionsave(err_act)

  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
  if (nprow == -1) then
     info = 2010
     call psb_errpush(info,name)
     goto 9999
  else if (npcol /= 1) then
     info = 2030
     int_err(1) = npcol
     call psb_errpush(info,name)
     goto 9999
  endif

  ia = 1
  ja = 1
  ix = 1
  jx = 1
  iy = 1
  jy = 1
  ik = 1

  if (present(doswap)) then
     idoswap = doswap
  else
     idoswap = 1
  endif

  if (present(trans)) then     
     if((trans.eq.'N').or.(trans.eq.'T')) then
        itrans = trans
     else if (trans.eq.'C') then
        info = 3020
        call psb_errpush(info,name)
        goto 9999
     else
        info = 70
        call psb_errpush(info,name)
        goto 9999
     end if
  else
     itrans = 'N'
  endif

  m    = desc_a%matrix_data(psb_m_)
  n    = desc_a%matrix_data(psb_n_)
  nrow = desc_a%matrix_data(psb_n_row_)
  ncol = desc_a%matrix_data(psb_n_col_)
  lldx = size(x,1)
  lldy = size(y,1)

  ! check for presence/size of a work area
  liwork= 2*ncol
  if (a%pr(1) /= 0) llwork = liwork + n * ik
  if (a%pl(1) /= 0) llwork = liwork + m * ik
  if (present(work)) then     
     if(size(work).lt.liwork) then
        call psb_realloc(liwork,work,info)
        if(info.ne.0) then
           info=4010
           ch_err='psb_realloc'
           call psb_errpush(info,name,a_err=ch_err)
           goto 9999
        end if
     end if
     iwork => work
  else
     call psb_realloc(liwork,iwork,info)
     if(info.ne.0) then
        info=4010
        ch_err='psb_realloc'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if
  end if

  ! checking for matrix correctness
  call psb_chkmat(m,n,ia,ja,desc_a%matrix_data,info,iia,jja)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkmat'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if


  if (itrans.eq.'N') then
     !  Matrix is not transposed
     if((ja.ne.ix).or.(ia.ne.iy)) then
        ! this case is not yet implemented
        info = 3030
        call psb_errpush(info,name)
        goto 9999
     end if

     ! checking for vectors correctness
     call psb_chkvect(n,ik,size(x),ix,jx,desc_a%matrix_data,info,iix,jjx)
     call psb_chkvect(m,ik,size(y),iy,jy,desc_a%matrix_data,info,iiy,jjy)
     if(info.ne.0) then
        info=4010
        ch_err='psb_chkvect'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

     if((iix.ne.1).or.(iiy.ne.1)) then
        ! this case is not yet implemented
        info = 3040
        call psb_errpush(info,name)
        goto 9999
     end if

     if(idoswap.lt.0) then
        x(nrow:ncol)=0.d0
     else
        call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
             & dzero,x,desc_a,iwork,info)
!!$        call PSI_dSwapData(ior(SWAP_SEND,SWAP_RECV),1,&
!!$             & dzero,x(iix,jjx),lldx,desc_a%matrix_data,&
!!$             & desc_a%halo_index,iwork,liwork,info)
     end if

     !  local Matrix-vector product
     call dcsmm(itrans,nrow,ib,ncol,alpha,a%pr,a%fida,&
          & a%descra,a%aspk,a%ia1,a%ia2,a%infoa,a%pl,&
          & x(iix),lldx,beta,y(iiy),lldy,&
          & iwork,liwork,info)

     if(info.ne.0) then
        info = 4011
        call psb_errpush(info,name)
        goto 9999
     end if

  else
     !  Matrix is transposed
     if((ja.ne.iy).or.(ia.ne.ix)) then
        ! this case is not yet implemented
        info = 3030
        call psb_errpush(info,name)
        goto 9999
     end if

     if(desc_a%ovrlap_elem(1).ne.-1) then
        info = 3070
        call psb_errpush(info,name)
        goto 9999
     end if

     ! checking for vectors correctness
     call psb_chkvect(m,ik,size(x),ix,jx,desc_a%matrix_data,info,iix,jjx)
     call psb_chkvect(n,ik,size(y),iy,jy,desc_a%matrix_data,info,iiy,jjy)
     if(info.ne.0) then
        info=4010
        ch_err='psb_chkvect'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

     if((iix.ne.1).or.(iiy.ne.1)) then
        ! this case is not yet implemented
        info = 3040
        call psb_errpush(info,name)
        goto 9999
     end if

     xp => x(iix:lldx)
     yp => x(iiy:lldy)

     if(idoswap.lt.0) y(nrow:ncol)=0.d0

     !  local Matrix-vector product
     call dcsmm(itrans,ncol,ik,nrow,alpha,a%pr,a%fida,&
          & a%descra,a%aspk,a%ia1,a%ia2,a%infoa,a%pl,&
          & x(iix),lldx,beta,y(iiy),lldy,&
          & iwork,liwork,info)
     if(info.ne.0) then
        info = 4010
        ch_err='dcsmm'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

     if(idoswap.gt.0)&
          & call psi_swaptran(ior(psb_swap_send_,psb_swap_recv_),&
          & done,yp,desc_a,iwork,info)
!!$          & call PSI_dSwapTran(ior(SWAP_SEND,SWAP_RECV),&
!!$          & ik,done,y(iiy,jjy),lldy,desc_a%matrix_data,&
!!$          & desc_a%halo_index),iwork,liwork,info
     if(info.ne.0) then
        info = 4010
        ch_err='PSI_dSwapTran'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

  end if

  if(.not.present(work)) deallocate(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end subroutine psb_dspmv
