! File: psb_dsprn.f90
!
! Subroutine: psb_dsprn
!    Reinit sparse matrix structure for psblas routines.
! 
! Parameters: 
!    a        - type(<psb_dspmat_type>).          The sparse matrix to be reinitiated.      
!    desc_a   - type(<psb_desc_type>).            The communication descriptor.
!    info     - integer.                          Eventually returns an error code.
!
Subroutine psb_dsprn(a, desc_a,info)

  use psb_descriptor_type
  use psb_spmat_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  Implicit None

  !....Parameters...
  Type(psb_desc_type), intent(in)      :: desc_a
  Type(psb_dspmat_type), intent(inout) :: a
  integer, intent(out)                 :: info

  !locals
  Integer             :: icontxt
  Integer             :: nprow,npcol,me,mypcol,err,err_act
  integer, parameter  :: ione=1, itwo=2,root=0
  logical, parameter  :: debug=.false.
  integer :: int_err(5)
  real(kind(1.d0)) :: real_err(5)
  character(len=20)  :: name, ch_err

  info = 0
  err  = 0
  int_err(1)=0
  name = 'psb_dsprn'
  call psb_erractionsave(err_act)

  icontxt = desc_a%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
  if (debug) &
       &write(*,*) 'starting spalloc ',icontxt,nprow,npcol,me

  !     ....verify blacs grid correctness..
  if (npcol.ne.1) then
     info = 2030
     call psb_errpush(info,name)
     goto 9999
  endif

  if (debug) &
       &write(*,*) 'got through igamx2d '


  if (.not.is_asb_dec(desc_a%matrix_data(psb_dec_type_))) then
     info=590     
     call psb_errpush(info,name)
     goto 9999
  endif

  if (a%infoa(state_) == spmat_asb) then

     a%aspk(:) = 0.0
     if (ibits(a%infoa(upd_),2,1)==1) then 
        if(a%fida(1:3).eq.'JAD') then
           a%ia1(a%infoa(upd_pnt_)+nnz_) = 0
        else
           a%ia2(a%infoa(upd_pnt_)+nnz_) = 0
        endif
     endif
     a%infoa(state_) = spmat_upd
  else if (a%infoa(state_) == spmat_bld) then
     ! in this case do nothing. this allows sprn to be called 
     ! right after allocate, with spins doing the right thing.
     ! hopefully :-)
  else if (a%infoa(state_) == spmat_upd) then

  else
     info=591     
     call psb_errpush(info,name)
  endif

  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_dsprn
