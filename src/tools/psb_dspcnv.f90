! File: psb_dspcnv.f90
!
! Subroutine: psb_dspcnv
!    converts sparse matrix a into b
! 
! Parameters: 
!    a        - type(<psb_dspmat_type>).          The sparse input matrix.      
!    b        - type(<psb_dspmat_type>).          The sparse output matrix.
!    desc_a   - type(<psb_desc_type>).            The communication descriptor.
!    info     - integer.                          Eventually returns an error code.
!
subroutine psb_dspcnv(a,b,desc_a,info)

  use psb_descriptor_type
  use psb_spmat_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  implicit none
  interface dcsdp

     subroutine dcsdp(check,trans,m,n,unitd,d,&
          & fida,descra,a,ia1,ia2,infoa,&
          & pl,fidh,descrh,h,ih1,ih2,infoh,pr,lh,lh1,lh2,&
          & work,lwork,ierror)
       integer, intent(in)   :: lh, lwork, lh1, lh2, m, n                 
       integer, intent(out)  :: ierror                 
       character, intent(in) :: check, trans, unitd                               
       real(kind(1.d0)), intent(in)  :: d(*), a(*)
       real(kind(1.d0)), intent(out) :: h(*)
       real(kind(1.d0)), intent(inout) :: work(*)
       integer, intent(in)  :: ia1(*), ia2(*), infoa(*)
       integer, intent(out) :: ih1(*), ih2(*), pl(*),pr(*), infoh(*) 
       character, intent(in) ::  fida*5, descra*11
       character, intent(out) :: fidh*5, descrh*11
     end subroutine dcsdp
  end interface


  interface dcsrp

     subroutine dcsrp(trans,m,n,fida,descra,ia1,ia2,&
          & infoa,p,work,lwork,ierror)
       integer, intent(in)  :: m, n, lwork
       integer, intent(out) :: ierror
       character, intent(in) ::       trans
       real(kind(1.d0)), intent(inout) :: work(*)                     
       integer, intent(in)    :: p(*)
       integer, intent(inout) :: ia1(*), ia2(*), infoa(*) 
       character, intent(in)  :: fida*5, descra*11
     end subroutine dcsrp
  end interface

  interface dcsprt
     subroutine dcsprt(m,n,fida,descra,a,ia1,ia2,infoa ,iout,ierror)
       integer, intent(in)  ::  iout,m, n                 
       integer, intent(out) ::  ierror                 
       real(kind(1.d0)), intent(in) :: a(*)
       integer, intent(in)   :: ia1(*), ia2(*), infoa(*)
       character, intent(in) :: fida*5, descra*11
     end subroutine dcsprt
  end interface

  !...parameters....
  type(psb_dspmat_type), intent(in)   :: a
  type(psb_dspmat_type), intent(out)  :: b
  type(psb_desc_type), intent(in)     :: desc_a
  integer, intent(out)                :: info
  !....locals....
  integer                       ::  int_err(5)
  integer,pointer               ::  ovrlap_elem(:),ovrlap_index(:)
  real(kind(1.d0))              ::  d(1)
  integer,pointer               ::  i_temp(:)
  real(kind(1.d0)),pointer      ::  work_dcsdp(:)
  integer                       ::  ia1_size,ia2_size,aspk_size,err_act&
       & ,i,err,nprow,npcol,myrow,mycol,n_col,l_dcsdp, iout, nrow
  integer                       ::  lwork_dcsdp,dectype
  integer                       ::  icontxt,temp(1),n_row
  character                     ::  check*1, trans*1, unitd*1
  integer, parameter            ::  ione=1

  real(kind(1.d0))              :: time(10), mpi_wtime
  external mpi_wtime
  logical, parameter :: debug=.false.
  character(len=20)   :: name, ch_err

  info=0
  name = 'psb_dspcnv'
  call psb_erractionsave(err_act)

  time(1) = mpi_wtime()


  icontxt = desc_a%matrix_data(psb_ctxt_)
  dectype = desc_a%matrix_data(psb_dec_type_)
  n_row   = desc_a%matrix_data(psb_n_row_)
  n_col   = desc_a%matrix_data(psb_n_col_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
  if (nprow.eq.-1) then
     info = 2010
     call psb_errpush(info,name)
     goto 9999
  else if (npcol.ne.1) then
     info = 2030
     int_err(1) = npcol
     call psb_errpush(info,name,i_err=int_err)
     goto 9999
  endif

  if (.not.psb_is_ok_dec((dectype))) then
     info = 600
     int_err(1) = dectype
     call psb_errpush(info,name,i_err=int_err)
     goto 9999
  endif

  if (debug) write (0, *) name,'   begin matrix assembly...'

  ia1_size  = size(a%ia1)
  ia2_size  = size(a%ia2)
  aspk_size = size(a%aspk)

  if (debug) write (0, *) name,'  sizes',ia1_size,ia2_size,aspk_size

  ! convert only without check
  check='N'
  trans='N'
  unitd='U'

  ! l_dcsdp is the size requested for dcsdp procedure
  l_dcsdp=(ia1_size+100)

  b%m=nrow
  b%k=n_col
  call psb_spall(b,ia1_size,ia2_size,aspk_size,info)
  allocate(work_dcsdp(l_dcsdp),stat=info)
  if (info.ne.0) then
     info=2025
     int_err(1)=l_dcsdp
     call psb_errpush(info, name, i_err=int_err)
     goto 9999
  endif

  lwork_dcsdp=size(work_dcsdp)
  ! set infoa(1) to nnzero
  b%pl(:)  = 0
  b%pr(:)  = 0

  if (debug) write (0, *) name,'  calling dcsdp',lwork_dcsdp,&
       &size(work_dcsdp)
  ! convert aspk,ia1,ia2 in requested representation mode
  if (debug) then

  endif
  ! result is put in b
  call dcsdp(check,trans,n_row,n_col,unitd,d,a%fida,a%descra,&
       & a%aspk,a%ia1,a%ia2,a%infoa,&
       & b%pl,b%fida,b%descra,b%aspk,b%ia1,b%ia2,b%infoa,b%pr,&
       & size(b%aspk),size(b%ia1),size(b%ia2),&
       & work_dcsdp,size(work_dcsdp),info)

  if(info.ne.no_err) then
     info=4010
     ch_err='spclone'
     call psb_errpush(info, name, a_err=ch_err)
     goto 9999
  end if

  !
  !  hmmm, have to fix b%pl and b%pr according to a%pl and a%pr!!! 
  !  should work (crossed fingers :-)
  if (a%pr(1).ne.0) then 
     if (b%pr(1).ne.0) then 
        allocate(i_temp(n_col))
        do i=1,  n_col
           i_temp(i) = b%pr(a%pr(i))
        enddo
        deallocate(b%pr)
        b%pr => i_temp
     else
        allocate(i_temp(n_col))
        do i=1,  n_col
           i_temp(i) = a%pr(i)
        enddo
        deallocate(b%pr)
        b%pr => i_temp
     endif
  endif
  if (a%pl(1).ne.0) then 
     if (b%pr(1).ne.0) then 
        allocate(i_temp(n_row))
        do i=1,  n_row
           i_temp(i) = a%pl(b%pl(i))
        enddo
        deallocate(b%pl)
        b%pl => i_temp
     else
        allocate(i_temp(n_row))
        do i=1,  n_row
           i_temp(i) = a%pl(i)
        enddo
        deallocate(b%pl)
        b%pl => i_temp
     endif
  endif


  if (debug) write (0, *) myrow,name,'  from dcsdp ',&
       &b%fida,' pl ', b%pl(:),'pr',b%pr(:)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_dspcnv
