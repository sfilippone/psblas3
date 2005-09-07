! File: psb_dnrmi.f90
!
! Function: psb_dnrmi
!    Forms the approximated norm of a sparse matrix,                                                                  
!
!    normi := max(abs(sum(A(i,j))))                                                                                   
!
! Parameters:
!    a      -  type(<psb_dspmat_type>).   The sparse matrix containing A.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!
function psb_dnrmi(a,desc_a,info)  
  use psb_descriptor_type
  use psb_serial_mod
  use psb_error_mod
  implicit none

  type(psb_dspmat_type), intent(in)   :: a
  integer, intent(out)                :: info
  type(psb_desc_type), intent(in)     :: desc_a

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, me, mycol,&
       & err_act, n, iia, jja, ia, ja, temp(2)
  real(kind(1.d0))         :: nrmi
  character(len=20)        :: name, ch_err

  name='psb_dnrmi'
  info=0
  call psb_erractionsave(err_act)

  icontxt=desc_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
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
  m = desc_a%matrix_data(m_)
  n = desc_a%matrix_data(n_)

  call psb_chkmat(m,n,ia,ja,desc_a%matrix_data,info,iia,jja)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkmat'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if ((iia.ne.1).or.(jja.ne.1)) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if

  if ((m.ne.0).and.(n.ne.0)) then
     mdim = desc_a%matrix_data(psb_n_row_)
     ndim = desc_a%matrix_data(psb_n_col_)
     nrmi = dcsnmi('N',mdim,ndim,a%fida,&
          & a%descra,a%aspk,a%ia1,a%ia2,&
          & a%infoa,info)

     if(info.ne.0) then
        info=4010
        ch_err='dcsnmi'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

     ! compute global max
     call dgamx2d(icontxt, 'A', ' ', ione, ione, nrmi, ione,&
          &temp ,temp,-ione ,-ione,-ione)
  else
     nrmi = 0.d0
  end if

  psb_nrmi = nrmi

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end function psb_dnrmi
