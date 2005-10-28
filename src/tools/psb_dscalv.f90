! File: psb_dscalv.f90
!
! Subroutine: psb_dscalv
!    Allocate descriptor
!    and checks correctness of PARTS subroutine
! 
! Parameters: 
!    m       - integer.                       The number of rows.
!    v       - integer, dimension(:).         The array containg the partitioning scheme.
!    icontxt - integer.                       The communication context.
!    desc_a  - type(<psb_desc_type>).         The communication descriptor.
!    info    - integer.                       Eventually returns an error code
!    flag    - integer.                       ???
subroutine psb_dscalv(m, v, icontxt, desc_a, info, flag)
  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  implicit None
  !....Parameters...
  Integer, intent(in)               :: m,icontxt, v(:)
  integer, intent(in), optional     :: flag
  integer, intent(out)              :: info
  type(psb_desc_type), intent(out)  :: desc_a

  !locals
  Integer             :: counter,i,j,nprow,npcol,myrow,mycol,&
       & loc_row,err,loc_col,nprocs,n,itmpov, k,&
       & l_ov_ix,l_ov_el,idx, flag_, err_act
  Integer             :: INT_ERR(5),TEMP(1),EXCH(2)
  Real(Kind(1.d0))    :: REAL_ERR(5)
  Integer, Pointer    :: temp_ovrlap(:), ov_idx(:),ov_el(:)
  logical, parameter  :: debug=.false.
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  err=0
  name = 'psb_dscalv'

  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
  if (debug) write(*,*) 'psb_dscall: ',nprow,npcol,myrow,mycol
  !     ....verify blacs grid correctness..
  if (npcol /= 1) then
     info = 2030
     int_err(1) = npcol
     call psb_errpush(info,name,i_err=int_err)
     goto 9999
  endif

  n = m
  !... check m and n parameters....
  if (m < 1) then
    info = 10
    int_err(1) = 1
    int_err(2) = m
  else if (n < 1) then
    info = 10
    int_err(1) = 2
    int_err(2) = n
  else if (size(v)<m) then 
    info = 10
    int_err(1) = 2
    int_err(2) = size(v)
  endif

  if (info /= 0) then 
     call psb_errpush(info,name,i_err=int_err)
     goto 9999
  end if

  if (debug) write(*,*) 'psb_dscall:  doing global checks'  
  !global check on m and n parameters
  if (myrow.eq.psb_root_) then
    exch(1)=m
    exch(2)=n
    call igebs2d(icontxt,psb_all_,psb_topdef_, itwo,ione, exch, itwo)
  else
    call igebr2d(icontxt,psb_all_,psb_topdef_, itwo,ione, exch, itwo, psb_root_,&
         & 0)
    if (exch(1) /= m) then
      info=550
      int_err(1)=1
    else if (exch(2) /= n) then
      info=550
      int_err(1)=2
    endif
  endif

  if (info /= 0) then 
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if


  call psb_nullify_desc(desc_a)

  if (present(flag)) then
    flag_=flag
  else
    flag_=0
  endif
  
  if ((flag_<0).or.(flag_>1)) then
    info = 6
    err=info
    call psb_errpush(info,name)
    goto 9999
  end if

  !count local rows number
  ! allocate work vector
  allocate(desc_a%glob_to_loc(m),desc_a%matrix_data(psb_mdata_size_),&
       &temp_ovrlap(m),stat=info)
  if (info /= 0) then     
    info=2025
    int_err(1)=m
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif


  if (debug) write(*,*) 'PSB_DSCALL:  starting main loop' ,info
  counter = 0
  itmpov  = 0
  temp_ovrlap(:) = -1
  do i=1,m
    
    if (((v(i)-flag_) > nprow-1).or.((v(i)-flag_) < 0)) then
      info=580
      int_err(1)=3
      int_err(2)=v(i) - flag_
      int_err(3)=i
      exit
    end if

    if ((v(i)-flag_) == myrow) then
      ! this point belongs to me
      counter=counter+1
      desc_a%glob_to_loc(i) = counter
    else
      desc_a%glob_to_loc(i) = -(nprow+(v(i)-flag_)+1)
    end if
  enddo

  loc_row=counter
  ! check on parts function
  if (debug) write(*,*) 'PSB_DSCALL:  End main loop:' ,loc_row,itmpov,info

  if (info /= 0) then 
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if
  if (debug) write(*,*) 'PSB_DSCALL:  error check:' ,err

  l_ov_ix=0
  l_ov_el=0
  i = 1
  do while (temp_ovrlap(i) /= -1) 
    idx = temp_ovrlap(i)
    i=i+1
    nprocs = temp_ovrlap(i)
    i = i + 1
    l_ov_ix = l_ov_ix+3*(nprocs-1)
    l_ov_el = l_ov_el + 2
    i = i + nprocs     
  enddo

  l_ov_ix = l_ov_ix+3  
  l_ov_el = l_ov_el+3

  if (debug) write(*,*) 'PSB_DSCALL: Ov len',l_ov_ix,l_ov_el
  allocate(ov_idx(l_ov_ix),ov_el(l_ov_el), stat=info)
  if (info /= 0) then
    info=2025
    int_err(1)=loc_col
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  l_ov_ix=0
  l_ov_el=0
  i = 1
  do while (temp_ovrlap(i) /= -1) 
    idx = temp_ovrlap(i)
    i   = i+1
    nprocs = temp_ovrlap(i)
    ov_el(l_ov_el+1)  = idx
    ov_el(l_ov_el+2)  = nprocs
    l_ov_el           = l_ov_el+2
    do j=1, nprocs
      if (temp_ovrlap(i+j) /= myrow) then
        ov_idx(l_ov_ix+1) = temp_ovrlap(i+j)
        ov_idx(l_ov_ix+2) = 1
        ov_idx(l_ov_ix+3) = idx
        l_ov_ix = l_ov_ix+3
      endif
    enddo
    i = i + nprocs +1
  enddo
  l_ov_el         = l_ov_el + 1
  ov_el(l_ov_el)  = -1
  l_ov_ix         = l_ov_ix + 1
  ov_idx(l_ov_ix) = -1

  desc_a%ovrlap_index => ov_idx
  desc_a%ovrlap_elem  => ov_el
  deallocate(temp_ovrlap,stat=info)
  if (info /= 0) then 
    info=4000
    call psb_errpush(info,name)
    goto 9999
  endif

  ! estimate local cols number 
  loc_col=int((psb_colrow_+1.d0)*loc_row)+1  
  allocate(desc_a%loc_to_glob(loc_col),&
       &desc_a%lprm(1),stat=info)  
  if (info /= 0) then
    info=2025
    int_err(1)=loc_col
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  ! set LOC_TO_GLOB array to all "-1" values
  desc_a%lprm(1) = 0
  desc_a%loc_to_glob(:) = -1
  do i=1,m
    k = desc_a%glob_to_loc(i) 
    if (k.gt.0) then 
      desc_a%loc_to_glob(k) = i
    endif
  enddo
  nullify(desc_a%bnd_elem,desc_a%halo_index)

!!$  if (debug) write(*,*) 'PSB_DSCALL:  Last bits in desc_a', loc_row,k
  ! set fields in desc_a%MATRIX_DATA....
  desc_a%matrix_data(psb_n_row_)  = loc_row
  desc_a%matrix_data(psb_n_col_)  = loc_row

  allocate(desc_a%halo_index(1),stat=info)
  if (info /= 0) then 
    info=4000
    call psb_errpush(info,name)
    goto 9999
  endif

  desc_a%halo_index(:) = -1


  desc_a%matrix_data(psb_m_)        = m
  desc_a%matrix_data(psb_n_)        = n
  desc_a%matrix_data(psb_dec_type_) = psb_desc_bld_
  desc_a%matrix_data(psb_ctxt_)     = icontxt
  call blacs_get(icontxt,10,desc_a%matrix_data(psb_mpi_c_))

  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_dscalv
