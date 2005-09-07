! File: psb_dcsrp.f90
!
! Subroutine: psb_dcsrp
!    Apply a right permutation to a sparse matrix, i.e. permute the column 
!    indices. 
! 
! Parameters: 
!    trans   - character.                       Whether iperm or its transpose should be applied
!    iperm   - integer, pointer, dimension(:).  A permutation vector; its size must be either N_ROW or N_COL
!    a       - type(<psb_dspmat_type).          The matrix to be permuted
!    desc_a  - type(<psb_desc_type>).           The communication descriptor.
!    info    - integer.                         Eventually returns an error code
subroutine psb_dcsrp(trans,iperm,a, desc_a, info)
  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  !  implicit none

  interface dcsrp

     subroutine dcsrp(trans,m,n,fida,descra,ia1,ia2,&
          & infoa,p,work,lwork,ierror)
       integer, intent(in)  :: m, n, lwork
       integer, intent(out) :: ierror
       character, intent(in) ::       trans
       double precision, intent(inout) :: work(*)                     
       integer, intent(in)    :: p(*)
       integer, intent(inout) :: ia1(*), ia2(*), infoa(*) 
       character, intent(in)  :: fida*5, descra*11
     end subroutine dcsrp
  end interface


  interface isaperm

    logical function isaperm(n,ip)
      integer, intent(in)    :: n   
      integer, intent(inout) :: ip(*)
    end function isaperm
  end interface

  !...parameters....
  type(psb_dspmat_type), intent(inout)  ::  a
  type(psb_desc_type), intent(in)       ::  desc_a
  integer, intent(inout)                :: iperm(:), info
  character, intent(in)                 :: trans
  !....locals....
  integer                               ::  int_err(5),p(1),infoa(10)
  real(kind(1.d0))                      ::  real_err(5)
  integer,pointer                       ::  ipt(:)
  integer                               ::  i,err,nprow,npcol,me,&
       & mypcol ,ierror ,n_col,l_dcsdp, iout, ipsize
  integer                               ::  dectype
  real(kind(1.d0)), pointer             ::  work_dcsdp(:)
  integer                               ::  icontxt,temp(1),n_row,err_act
  character(len=20)                     ::  name, char_err

  real(kind(1.d0))                      ::  time(10), mpi_wtime
  external mpi_wtime
  logical, parameter :: debug=.false.

  time(1) = mpi_wtime()

  icontxt=desc_a%matrix_data(psb_ctxt_)
  dectype=desc_a%matrix_data(psb_dec_type_)
  n_row = desc_a%matrix_data(psb_n_row_)
  n_col = desc_a%matrix_data(psb_n_col_)
     
  info=0
  call psb_erractionsave(err_act)
  name = 'psd_csrp'

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


  if (.not.is_asb_dec(dectype)) then 
    info = 600
    int_err(1) = dectype
    call psb_errpush(info,name,int_err)
    goto 9999
 endif

  ipsize = size(iperm)
  if (.not.((ipsize.eq.n_col).or.(ipsize.eq.n_row) )) then 
    info = 35
    int_err(1) = 1
    int_err(2) = ipsize
    call psb_errpush(info,name,int_err)
    goto 9999
  else
    if (.not.isaperm(ipsize,iperm)) then
      info = 70
      int_err(1) = 1      
      call psb_errpush(info,name,int_err)
      goto 9999
    endif
  endif

  l_dcsdp = (n_col)
  
  call psb_realloc(l_dcsdp,work_dcsdp,info)
  call psb_realloc(n_col,ipt,info)
  if(info /= no_err) then
     info=4010
     char_err='psrealloc'
     call psb_errpush(info,name,a_err=char_err)
     goto 9999
  end if

  if (ipsize.eq.n_col) then 
    do i=1, n_col
      ipt(i) = iperm(i)
    enddo
  else    
    do i=1, n_row
      ipt(i) = iperm(i)
    enddo
    do i=n_row+1,n_col
      ipt(i) = i
    enddo
  endif
  ! crossed fingers.....
  ! fix glob_to_loc/loc_to_glob  mappings, then indices lists
  ! hmm, maybe we should just move all of this onto a different level,
  ! have a specialized subroutine, and do it in the solver context???? 
  if (debug) write(0,*) 'spasb: calling dcsrp',size(work_dcsdp)
  call dcsrp(trans,n_row,n_col,a%fida,a%descra,a%ia1,a%ia2,a%infoa,&
       & ipt,work_dcsdp,size(work_dcsdp),info)
  if(info /= no_err) then
     info=4010
     char_err='dcsrp'
     call psb_errpush(info,name,a_err=char_err)
     goto 9999
  end if
  
  deallocate(ipt,work_dcsdp)
  
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
     call psb_error()
     return
  end if
  return

end subroutine psb_dcsrp
