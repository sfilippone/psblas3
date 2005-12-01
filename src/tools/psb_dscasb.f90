! File: psb_dscasb.f90
!
! Subroutine: psb_dscasb
!   Assembly the psblas communications descriptor.
! 
! Parameters: 
!    desc_a  - type(<psb_desc_type>).         The communication descriptor.
!    info    - integer.                       Eventually returns an error code.
subroutine psb_dscasb(desc_a,info)
  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psi_mod
  use psb_error_mod
  implicit none

  !...Parameters....
  type(psb_desc_type), intent(inout) :: desc_a
  integer, intent(out)               :: info


  !....Locals....
  integer          ::  int_err(5), itemp(2)
  integer,pointer  ::  ovrlap_index(:),halo_index(:)
  integer          ::  i,err,nprow,npcol,me,mypcol,&
       & lovrlap,lhalo,nhalo,novrlap,max_size,max_halo,n_col,ldesc_halo,&
       & ldesc_ovrlap, dectype, err_act
  integer                       :: icontxt,temp(1),n_row
  logical, parameter            :: debug=.false., debugwrt=.false.
  character(len=20)             :: name,ch_err

  info = 0
  int_err(1) = 0
  name = 'psb_dscasb'

  call psb_erractionsave(err_act)

  icontxt = desc_a%matrix_data(psb_ctxt_)
  dectype = desc_a%matrix_data(psb_dec_type_)
  n_row   = desc_a%matrix_data(psb_n_row_)
  n_col   = desc_a%matrix_data(psb_n_col_)

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

  if (.not.psb_is_ok_dec(dectype)) then 
    info = 600
    int_err(1) = dectype
    call psb_errpush(info,name)
    goto 9999
  endif

  if (debug) write (0, *) '   Begin matrix assembly...'

  if (psb_is_bld_dec(dectype)) then 
    if (debug) write(0,*) 'psb_dscasb: Checking rows insertion'
    ! check if all local row are inserted
    do i=1,desc_a%matrix_data(psb_n_col_)
      if (desc_a%loc_to_glob(i) < 0) then
         info=3100
        exit
      endif
    enddo

    if (info /= no_err) then    
       call psb_errpush(info,name,i_err=int_err)
       goto 9999
    endif


    ! comm desc_size is size requested for temporary comm descriptors
    ! (expressed in No of dble element)
    ldesc_halo   = (((3*(n_col-n_row)+1)+1))
    ovrlap_index => desc_a%ovrlap_index
    nullify(desc_a%ovrlap_index)
    halo_index => desc_a%halo_index
    nullify(desc_a%halo_index)

    lhalo = 1
    do while (halo_index(lhalo) /= -1)
      lhalo = lhalo + 1 
    enddo
    nhalo = (lhalo-1)/3
    lovrlap=1
    do while (ovrlap_index(lovrlap) /= -1) 
      lovrlap=lovrlap+1
    enddo
    novrlap = (lovrlap-1)/3

    ! Allocate final comm PSBLAS descriptors
    ! compute necessary dimension of halo index
    max_halo  = max(nhalo,1)
    max_size  = max(1,min(3*desc_a%matrix_data(psb_n_row_),novrlap*3))

    itemp(1) = max_size
    itemp(2) = max_halo
    call igamx2d(icontxt, psb_all_, psb_topdef_, itwo, ione, itemp,&
         & itwo,temp ,temp,-ione ,-ione,-ione)
    max_size = itemp(1) 
    max_halo = itemp(2) 

    ldesc_halo = 3*max_halo+3*nhalo+1

    ! allocate HALO_INDEX field
    call psb_realloc(ldesc_halo, desc_a%halo_index, info)
    ! check on allocate
    if (info /= no_err) then
      info=4010
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    endif

    ! compute necessary dimension of ovrlap index
    ldesc_ovrlap = 2*lovrlap+1
    
    ! allocate OVRLAP_INDEX field
    call psb_realloc(ldesc_ovrlap, desc_a%ovrlap_index, info)
    ! check on allocate
    if (info /= no_err) then
       info=4010
       ch_err='psb_realloc'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    endif
    
    
      
    if (debug) write(0,*) 'psb_dscasb: converting indexes',&
         & nhalo,lhalo,halo_index(lhalo)
    !.... convert comunication stuctures....
    ! first the halo index
    
    call psi_crea_index(desc_a,halo_index,&
         & desc_a%halo_index,.false.,info)
    if(info.ne.0) then
       call psb_errpush(4010,name,a_err='psi_crea_index')
       goto 9999
    end if

    ! then the overlap index
    call psi_crea_index(desc_a,ovrlap_index,&
         & desc_a%ovrlap_index,.true.,info)
    if(info.ne.0) then
       call psb_errpush(4010,name,a_err='psi_crea_index')
       goto 9999
    end if

    ! next is the ovrlap_elem index
    call psi_crea_ovr_elem(desc_a%ovrlap_index,desc_a%ovrlap_elem)

    ! finally bnd_elem
    call psi_crea_bnd_elem(desc_a,info)
    if(info.ne.0) then
       call psb_errpush(4010,name,a_err='psi_crea_bnd_elem')
       goto 9999
    end if

    ! Ok, register into MATRIX_DATA &  free temporary work areas
    desc_a%matrix_data(psb_dec_type_) = psb_desc_asb_

    deallocate(ovrlap_index, stat=info)
    deallocate(halo_index, stat=info)
    if (info /= 0) then
       info =4000
       call psb_errpush(info,name)
       goto 9999
    end if

  else
     info = 600
     call psb_errpush(info,name)
     goto 9999
     if (debug) write(0,*) 'dectype 2 :',dectype,psb_desc_bld_,&
          &psb_desc_asb_,psb_desc_upd_
  endif

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
  
end subroutine psb_dscasb
