! File: psb_dspalloc.f90
!
! Subroutine: psb_dspalloc
!    Allocate sparse matrix structure for psblas routines.
! 
! Parameters: 
!    a        - type(<psb_dspmat_type>).       The sparse matrix to be allocated.      
!    desc_a   - type(<psb_desc_type>).         The communication descriptor to be updated.
!    info     - integer.                       Eventually returns an error code.
!    nnz      - integer(optional).             The number of nonzeroes in the matrix.
!
subroutine psb_dspalloc(a, desc_a, info, nnz)

  use psb_descriptor_type
  use psb_spmat_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  implicit none

  !....parameters...
  type(psb_desc_type), intent(inout) :: desc_a
  type(psb_dspmat_type), intent(out) :: a
  integer, intent(out)               :: info
  integer, optional, intent(in)      :: nnz

  !locals
  integer             :: icontxt, dectype
  integer             :: nprow,npcol,myrow,mycol,loc_row,&
       &  length_ia1,length_ia2,err,nprocs, err_act,m,n
  integer             :: int_err(5),temp(1)
  real(kind(1.d0))    :: real_err(5)
  integer, parameter  :: ione=1, itwo=2,root=0
  logical, parameter  :: debug=.false.
  character(len=20)   :: name, ch_err

  info=0
  call psb_erractionsave(err_act)
  name = 'psb_dspalloc'

  icontxt = desc_a%matrix_data(psb_ctxt_)
  dectype=desc_a%matrix_data(psb_dec_type_)
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
!     ....verify blacs grid correctness..
  if (nprow.eq.-1) then
     info = 2010
     call psb_errpush(info,name)
     goto 9999
  else if (npcol.ne.1) then
     info = 2030
     int_err(1) = npcol
     call psb_errpush(info,name,int_err)
     goto 9999
  endif

  !
  ! hmm, not a good idea, not all compilers can rely on any given
  ! value for non initialized pointers. let's avoid this,
  ! and just rely on documentation. 
  ! check if psdalloc is already called for this matrix

  ! set fields in desc_a%matrix_data....
  loc_row = desc_a%matrix_data(psb_n_row_)
  m       = desc_a%matrix_data(psb_m_)
  n       =  desc_a%matrix_data(psb_n_)

  !...allocate matrix data...
  if (present(nnz))then 
    if (nnz.lt.0) then
	info=45
	int_err(1)=7
	int_err(2)=nnz
        call psb_errpush(info,name,int_err)
        goto 9999
     endif
     length_ia1=nnz
     length_ia2=nnz
  else 
     length_ia1=max(1,4*loc_row)
     length_ia2=max(1,4*loc_row)
  endif

  if (debug) write(*,*) 'allocating size:',length_ia1

  !....allocate aspk, ia1, ia2.....
  call psb_spall(loc_row,loc_row,a,length_ia1,info)
  if(info.ne.0) then
     info=4010
     ch_err='spreall'
     call psb_errpush(info,name,int_err)
     goto 9999
  end if

  ! set permutation matrices
  a%pl(1)=0
  a%pr(1)=0
  ! set infoa fields
  a%fida   = 'COO'
  a%descra = 'GUN'
  a%infoa(psb_nnz_)  = 0
  a%infoa(psb_srtd_) = 0
  a%infoa(psb_state_) = psb_spmat_bld_

  if (debug) write(0,*) 'spall: ',  &
       &desc_a%matrix_data(psb_dec_type_),psb_desc_bld_
  desc_a%matrix_data(psb_dec_type_) = psb_desc_bld_
  return
  
  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  
  if (err_act.eq.act_ret) then
     return
  else
     call psb_error(icontxt)
  end if
  return

end subroutine psb_dspalloc
