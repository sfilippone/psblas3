! File: psb_ptasb.f90
!
! Subroutine: psb_ptasb
!   ???
! 
! Parameters: 
!    desc_a  - type(<psb_desc_type>).         The communication descriptor.
!    info    - integer.                       Eventually returns an error code.
!
subroutine psb_ptasb(desc_a,info)

  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psi_mod
  use psb_error_mod
  implicit none


  !...parameters....
  type(psb_desc_type), intent(inout) :: desc_a
  integer,intent(out)                :: info
  !....locals....
  integer                       ::  int_err(5)
  integer,pointer               ::  ovrlap_index(:),halo_index(:)
  real(kind(1.d0))              ::  real_err(5)
  integer,pointer               ::  work5(:)
  integer                       ::  err_act,&
       & i,nprow,npcol,me,mypcol ,size_req,&
       & lovrlap,lhalo,nhalo,novrlap,max_size,max_size1,&
       & max_halo,size_req1,n_col,lwork5,ldesc_halo,&
       & ldesc_ovrlap, dectype
  integer                       :: icontxt,temp(1),n_row
  integer, parameter            :: ione=1
  real(kind(1.d0))              :: time(10), mpi_wtime
  external mpi_wtime
  logical, parameter :: debug=.false., debugwrt=.false.
  character(len=20)   :: name, ch_err

  info=0
  name = 'psb_ptasb'
  call psb_erractionsave(err_act)

  time(1) = mpi_wtime()


  icontxt = desc_a%matrix_data(psb_ctxt_)
  dectype = desc_a%matrix_data(psb_dec_type_)
  n_row   = desc_a%matrix_data(psb_n_row_)
  n_col   = desc_a%matrix_data(psb_n_col_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
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

  if (.not.is_ok_dec(dectype)) then 
     info = 600
     int_err(1) = dectype
     call psb_errpush(info, name, int_err)
     goto 9999
  endif

  if (debug) write (*, *) '   begin matrix assembly...'

  if (is_bld_dec(dectype)) then 
     if (debug) write(0,*) 'ptasb: checking rows insertion'
     ! check if all local row are inserted
     do i=1,desc_a%matrix_data(psb_n_col_)
        if (desc_a%loc_to_glob(i).lt.0) then
           write(0,*) 'error on index: ',i,desc_a%loc_to_glob(i)
           info=3100
           exit
        endif
     enddo

     if(info.ne.no_err) then
        call psb_errpush(info, name)
        goto 9999
     end if

     ! comm desc_size is size requested for temporary comm descriptors
     ! (expressed in no of dble element)
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

     if (debug) write(0,*) 'ptasb: from asbx',&
          & nhalo,lhalo,halo_index(lhalo)

     ! allocate final comm psblas descriptors

     ! compute necessary dimension of halo index
     max_halo=max(nhalo,1)
     max_size= max(1,min(3*desc_a%matrix_data(psb_n_row_),novrlap*3))
     max_size1=max_size

     call igamx2d(icontxt, all, topdef, ione, ione, max_size,&
          & ione,temp ,temp,-ione ,-ione,-ione)
     call igamx2d(icontxt, all, topdef, ione, ione,max_halo,&
          & ione,temp ,temp,-ione ,-ione,-ione)

     ldesc_halo=3*max_halo+3*nhalo+1

     ! allocate halo_index field
     allocate(desc_a%halo_index(ldesc_halo),stat=info)
     ! check on allocate
     if (info.ne.0) then
        info=2023
        int_err(1)=ldesc_halo
        call psb_errpush(info, name, int_err)
        goto 9999
     endif

     ! compute necessary dimension of ovrlap index
     ldesc_ovrlap=2*lovrlap+1

     ! allocate ovrlap_index field
     allocate(desc_a%ovrlap_index(ldesc_ovrlap),stat=info)
     ! check on allocate
     if (info.ne.0) then
        info=2023
        int_err(1)=ldesc_ovrlap
        call psb_errpush(info, name, int_err)
        goto 9999
     endif

     size_req1=max((max_size+max_size1),(nprow*4)*(nprow+1)*5)

     size_req=((nprow*4)*(nprow+1)*5)+(nprow+1)+(max&
          & (nhalo,size_req1)+1)
     if (info.ne.0)  info =2040
     allocate(work5(size_req),stat=info)
     if (info.ne.0) then
        info=2025
        int_err(1)=size_req
        call psb_errpush(info, name, int_err)
        goto 9999
     endif
     lwork5=size(work5)

     if (debug) write(0,*) 'ptasb: calling convert_comm',&
          & nhalo,lhalo,halo_index(lhalo)
     !.... convert comunication stuctures....
     call psi_convert_comm(desc_a%matrix_data,&
          & halo_index, ovrlap_index,&
          & desc_a%halo_index,size(desc_a%halo_index),&
          & desc_a%ovrlap_index,size(desc_a%ovrlap_index),&
          & desc_a%ovrlap_elem,size(desc_a%ovrlap_elem),&
          & desc_a%bnd_elem,&
          & desc_a%loc_to_glob,desc_a%glob_to_loc,info)
     if(info /= 0) then
        info=4010
        ch_err='psi_convert_comm'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if
     
     ! ok, register into matrix_data &  free temporary work areas
     desc_a%matrix_data(psb_dec_type_) = desc_asb
     deallocate(halo_index,ovrlap_index,&
          & work5, stat=info)
     if (info.ne.0)  then
        info =2040
        call psb_errpush(info, name, int_err)
        goto 9999
     end if

  else
     info = 600
     if (debug) write(0,*) 'dectype 2 :',dectype,desc_bld,&
          &desc_asb,desc_upd
     call psb_errpush(info, name, int_err)
     goto 9999
  endif

  time(4) = mpi_wtime()
  time(4) = time(4) - time(3)
  if (debug) then 
     call dgamx2d(icontxt, all, topdef, ione, ione, time(4),&
          & ione,temp ,temp,-ione ,-ione,-ione)

     write (*, *) '         comm structs assembly: ', time(4)*1.d-3
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_ptasb
