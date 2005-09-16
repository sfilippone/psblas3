! File: psb_dspsm.f90
!
! Subroutine: psb_dspsm
!  Performs one of the distributed matrix-vector operations
!
!     sub( Y ) := alpha * Pr * A-1 * Pc *sub( X ) + beta * sub (Y ),   or
!
!     sub( Y ) := alpha * D * Pr * A-1 * Pc * sub( X ) + beta * sub (Y ), or
!
!     sub( Y ) := alpha * Pr * A-1 * Pc * D * sub( X ) + beta * sub (Y ),  or
!
!     sub( Y ) := alpha * Pr * A-T * Pc * sub( X ) + beta * sub (Y ),   or
!
!     sub( Y ) := alpha * D * Pr * A-T * Pc * sub( X ) + beta * sub (Y ), or
!
!     sub( Y ) := alpha * Pr * A-T * Pc * D * sub( X ) + beta * sub (Y ),  or
!
!  where :
!
!        sub( X ) denotes X(1:M,JX:JX+K-1),
!      
!        sub( Y ) denotes Y(1:M,JY:JY+K-1).
!
!  sub( X ) is a distributed
!  vector and T is a M-by-M distributed triangular matrix.
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
!    unitd  -  character(optional).         Specify some type of operation with the diagonal matrix D.
!    choice -  integer(optional).           The kind of update to perform on overlap elements.
!    d      -  real,dimension(:)(optional). Matrix for diagonal scaling.
!    k      -  integer(optional).           The number of right-hand sides.
!    jx     -  integer(optional).           The column offset for sub( X ). If not present 1 is assumed.
!    jy     -  integer(optional).           The column offset for sub( Y ). If not present 1 is assumed.
!    work   -  real,dimension(:)(optional). Working area.
! 
subroutine  psb_dspsm(alpha,a,x,beta,y,desc_a,info,&
     & trans, unitd, choice, d, k, jx, jy, work)   

  use psb_spmat_type
  use psb_serial_mod
  use psb_descriptor_type
  use psb_comm_mod
  use psi_mod
  use psb_check_mod
  use psb_error_mod
  implicit none

  real(kind(1.D0)), intent(in)              :: alpha, beta
  real(kind(1.d0)), intent(in), target      :: x(:,:)
  real(kind(1.d0)), intent(inout), target   :: y(:,:)
  type (psb_dspmat_type), intent(in)        :: a
  type(psb_desc_type), intent(in)           :: desc_a
  integer, intent(out)                      :: info
  real(kind(1.d0)), intent(in), optional, target      :: d(:)
  real(kind(1.d0)), optional, pointer       :: work(:)
  character, intent(in), optional           :: trans, unitd
  integer, intent(in), optional             :: choice
  integer, intent(in), optional             :: k, jx, jy

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, n, iix, jjx, ia, ja, iia, jja, temp(2), lldx,lldy, lchoice,&
       & ix, iy, ik, ijx, ijy, i, lld,&
       & idoswap, m, nrow, ncol, liwork, llwork, iiy, jjy

  character                :: lunitd
  integer, parameter       :: nb=4
  real(kind(1.d0)),pointer :: iwork(:), xp(:,:), yp(:,:), id(:)
  character                :: itrans
  character(len=20)        :: name, ch_err

  name='psb_dspsm'
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

  ! just this case right now
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

  if (present(k)) then
     ik = min(k,size(x,2)-ijx+1)
     ik = min(ik,size(y,2)-ijy+1)
  else
     ik = min(size(x,2)-ijx+1,size(y,2)-ijy+1)
  endif

  if (present(choice)) then     
    lchoice = choice
  else
    lchoice = psb_avg_
  endif

  if (present(unitd)) then     
    lunitd = unitd
  else
    lunitd = 'U'
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
  nrow = desc_a%matrix_data(psb_n_row_)
  ncol = desc_a%matrix_data(psb_n_col_)
  lldx = size(x,1)
  lldy = size(y,1)

  if((lldx.lt.ncol).or.(lldy.lt.ncol)) then
     info=3010
     call psb_errpush(info,name)
     goto 9999
  end if

  ! check for presence/size of a work area
  liwork= 2*ncol
  if (a%pr(1) /= 0) llwork = liwork + m * ik
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

  if(present(d)) then
     lld = size(d)
     id => d
  else
     lld=1
     allocate(id(1))
     id=1.d0
  end if

  ! checking for matrix correctness
  call psb_chkmat(m,m,ia,ja,desc_a%matrix_data,info,iia,jja)
  ! checking for vectors correctness
  call psb_chkvect(m,ik,size(x,1),ix,ijx,desc_a%matrix_data,info,iix,jjx)
  call psb_chkvect(m,ik,size(y,1),iy,ijy,desc_a%matrix_data,info,iiy,jjy)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect/mat'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if(ja.ne.ix) then
     ! this case is not yet implemented
     info = 3030
  end if

  if((iix.ne.1).or.(iiy.ne.1)) then
     ! this case is not yet implemented
     info = 3040
  end if

  if(info.ne.0) then
     call psb_errpush(info,name)
     goto 9999
  end if

  ! Perform local triangular system solve
  call dcssm(itrans,nrow,ik,alpha,lunitd,id,a%pr,&
       & a%fida,a%descra,a%aspk,a%ia1,a%ia2,a%infoa,&
       & a%pl,x(iix,jjx),lldx,beta,y(iiy,jjy),lldy,&
       & iwork,liwork,info)
  if(info.ne.0) then
     info = 4010
     ch_err='dcssm'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  ! update overlap elements
  if(lchoice.gt.0) then
     yp => y(iiy:lldy,jjy:jjy+ik-1)
     call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),ik,&
          & done,yp,desc_a,iwork,info)
!!$     call PSI_dSwapData(ior(SWAP_SEND,SWAP_RECV),ik,&
!!$          & done,y,lldy,desc_a%matrix_data,desc_a%ovrlap_index,&
!!$          & iwork,liwork,info)

     i=0
     ! switch on update type
     select case (lchoice)
     case(psb_square_root_)
        do while(desc_a%ovrlap_elem(i).ne.-ione)
           y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:) =&
                & y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:)/&
                & sqrt(real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_)))
           i = i+2
        end do
     case(psb_avg_)
        do while(desc_a%ovrlap_elem(i).ne.-ione)
           y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:) =&
                & y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:)/&
                & real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_))
           i = i+2
        end do
     case(psb_sum_)
        ! do nothing
     case default 
        ! wrong value for choice argument
        info = 70
        int_err=(/10,lchoice,0,0,0/)
        call psb_errpush(info,name,i_err=int_err)
        goto 9999
     end select
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
end subroutine psb_dspsm
     
! Subroutine: psb_dspsmv
!  Performs one of the distributed matrix-vector operations
!
!     Y := alpha * Pr * A-1 * Pc * X + beta * Y,   or
!
!     Y := alpha * D * Pr * A-1 * Pc * X + beta * Y, or
!
!     Y := alpha * Pr * A-1 * Pc * D * X + beta * Y,  or
!
!     Y := alpha * Pr * A-T * Pc * X + beta * Y,   or
!
!     Y := alpha * D * Pr * A-T * Pc * X + beta * Y, or
!
!     Y := alpha * Pr * A-T * Pc * D * X + beta * Y,  or
!
!  X is a distributed
!  vector and T is a M-by-M distributed triangular matrix.
!
! Parameters:   
!    alpha  -  real.                        The scalar alpha.
!    a      -  type(<psb_dspmat_type>).     The sparse matrix containing A.
!    x      -  real,dimension(:).           The input vector containing the entries of X.
!    beta   -  real.                        The scalar beta.
!    y      -  real,dimension(:).           The input vector containing the entries of Y.
!    desc_a -  type(<psb_desc_type>).       The communication descriptor.
!    info   -  integer.                     Eventually returns an error code.
!    trans  -  character(optional).         Whether A or A'. If not present 'N' is assumed.
!    unitd  -  character(optional).         Specify some type of operation with the diagonal matrix D.
!    choice -  integer(optional).           The kind of update to perform on overlap elements.
!    d      -  real,dimension(:)(optional). Matrix for diagonal scaling.
!    work   -  real,dimension(:)(optional). Working area.
! 
subroutine  psb_dspsv(alpha,a,x,beta,y,desc_a,info,&
     & trans, unitd, choice, d, work)   
  use psb_spmat_type
  use psb_serial_mod
  use psb_descriptor_type
  use psb_comm_mod
  use psi_mod
  use psb_check_mod
  use psb_error_mod

  real(kind(1.D0)), intent(in)              :: alpha, beta
  real(kind(1.d0)), intent(in), target      :: x(:)
  real(kind(1.d0)), intent(inout), target   :: y(:)
  type(psb_dspmat_type), intent(in)         :: a
  type(psb_desc_type), intent(in)           :: desc_a
  integer, intent(out)                      :: info
  real(kind(1.d0)), intent(in), optional, target    :: d(:)
  real(kind(1.d0)), optional, pointer       :: work(:)
  character, intent(in), optional           :: trans, unitd
  integer, intent(in), optional             :: choice

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, n, iix, jjx, ia, ja, iia, jja, temp(2), lldx,lldy, lchoice,&
       & ix, iy, ik, jx, jy, i, lld,&
       & idoswap, m, nrow, ncol, liwork, llwork, iiy, jjy

  character                :: lunitd
  integer, parameter       :: nb=4
  real(kind(1.d0)),pointer :: iwork(:), xp(:), yp(:), id(:)
  character                :: itrans
  character(len=20)        :: name, ch_err

  name='psb_dspsv'
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

  ! just this case right now
  ia = 1
  ja = 1
  ix = 1
  iy = 1
  ik = 1
  jx= 1
  jy= 1

  if (present(choice)) then     
    lchoice = choice
  else
    lchoice = psb_avg_
  endif

  if (present(unitd)) then     
    lunitd = unitd
  else
    lunitd = 'U'
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
  nrow = desc_a%matrix_data(psb_n_row_)
  ncol = desc_a%matrix_data(psb_n_col_)
  lldx = size(x)
  lldy = size(y)

  if((lldx.lt.ncol).or.(lldy.lt.ncol)) then
     info=3010
     call psb_errpush(info,name)
     goto 9999
  end if

  ! check for presence/size of a work area
  liwork= 2*ncol
  if (a%pr(1) /= 0) llwork = liwork + m * ik
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

  if(present(d)) then
     lld = size(d)
     id => d
  else
     lld=1
     allocate(id(1))
     id=1.d0
  end if

  ! checking for matrix correctness
  call psb_chkmat(m,m,ia,ja,desc_a%matrix_data,info,iia,jja)
  ! checking for vectors correctness
  call psb_chkvect(m,ik,size(x),ix,jx,desc_a%matrix_data,info,iix,jjx)
  call psb_chkvect(m,ik,size(y),iy,jy,desc_a%matrix_data,info,iiy,jjy)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect/mat'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if(ja.ne.ix) then
     ! this case is not yet implemented
     info = 3030
  end if

  if((iix.ne.1).or.(iiy.ne.1)) then
     ! this case is not yet implemented
     info = 3040
  end if

  if(info.ne.0) then
     call psb_errpush(info,name)
     goto 9999
  end if

  ! Perform local triangular system solve
  call dcssm(itrans,nrow,ik,alpha,lunitd,id,a%pr,&
       & a%fida,a%descra,a%aspk,a%ia1,a%ia2,a%infoa,&
       & a%pl,x,lldx,beta,y,lldy,&
       & iwork,liwork,info)
  if(info.ne.0) then
     info = 4010
     ch_err='dcssm'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  ! update overlap elements
  if(lchoice.gt.0) then
     yp => y(iiy:lldy)
     call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
          & done,yp,desc_a,iwork,info)
!!$     call PSI_dSwapData(ior(SWAP_SEND,SWAP_RECV),ik,&
!!$          & done,y,lldy,desc_a%matrix_data,desc_a%ovrlap_index,&
!!$          & iwork,liwork,info)

     i=0
     ! switch on update type
     select case (lchoice)
     case(psb_square_root_)
        do while(desc_a%ovrlap_elem(i).ne.-ione)
           y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_)) =&
                & y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_))/&
                & sqrt(real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_)))
           i = i+2
        end do
     case(psb_avg_)
        do while(desc_a%ovrlap_elem(i).ne.-ione)
           y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_)) =&
                & y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_))/&
                & real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_))
           i = i+2
        end do
     case(psb_sum_)
        ! do nothing
     case default 
        ! wrong value for choice argument
        info = 70
        int_err=(/10,lchoice,0,0,0/)
        call psb_errpush(info,name,i_err=int_err)
        goto 9999
     end select
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
end subroutine psb_dspsv
     
