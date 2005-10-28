! File: psb_dscren.f90
!
! Subroutine: psb_dscren
!    Updates a communication descriptor according to a renumbering scheme.
! 
! Parameters: 
!    trans    - character.                     Whether iperm or its transpose should be applied.
!    iperm    - integer,dimension(:).          The renumbering scheme.
!    desc_a   - type(<psb_desc_type>).         The communication descriptor to be updated.
!    info     - integer.                       Eventually returns an error code.
!
subroutine psb_dscren(trans,iperm,desc_a,info)
  use psb_descriptor_type
  use psb_const_mod
  use psb_error_mod
  implicit none

  interface isaperm
    logical function isaperm(n,ip)
      integer, intent(in)    :: n   
      integer, intent(inout) :: ip(*)
    end function isaperm
  end interface

  !...parameters....
  type(psb_desc_type), intent(inout)    ::  desc_a
  integer, intent(inout)                ::  iperm(:)
  character, intent(in)                 :: trans
  integer, intent(out)                  :: info
  !....locals....
  integer                       :: i,j,err,nprow,npcol,myrow,mycol, n_col, kh, nh
  integer                       :: dectype
  integer                       :: icontxt,temp(1),n_row, int_err(5), err_act
  real(kind(1.d0))              :: time(10), mpi_wtime, real_err(6)
  external mpi_wtime
  logical, parameter            :: debug=.false.
  character(len=20)             :: name, char_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)
  name = 'psb_dcren'

  time(1) = mpi_wtime()

  icontxt=desc_a%matrix_data(psb_ctxt_)
  dectype=desc_a%matrix_data(psb_dec_type_)
  n_row = desc_a%matrix_data(psb_n_row_)
  n_col = desc_a%matrix_data(psb_n_col_)
     
  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
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

  if (.not.psb_is_asb_dec(dectype)) then 
    info = 600
    int_err(1) = dectype
     call psb_errpush(info,name,int_err)
     goto 9999
  endif

  if (iperm(1) /= 0) then 
    if (.not.isaperm(n_row,iperm)) then
      info = 610
      int_err(1) = iperm(1)
      call psb_errpush(info,name,int_err)
      goto 9999
    endif
 endif

  if (debug) write (*, *) '   begin matrix assembly...'

  !check on errors encountered in psdspins
  
  if ((iperm(1) /= 0))   then 

    if (debug) write(0,*) 'spasb: here we go with ',iperm(1) 
    deallocate(desc_a%lprm)
    allocate(desc_a%lprm(n_col))
    if (trans.eq.'n') then 
      do i=1, n_row
        desc_a%lprm(iperm(i)) = i
      enddo
      do i=n_row+1,n_col
        desc_a%lprm(i) = i
      enddo
    else if (trans.eq.'t') then 
      do i=1, n_row
        desc_a%lprm(i) = iperm(i)
      enddo
      do i=n_row+1,n_col
        desc_a%lprm(i) = i
      enddo
    endif
    ! crossed fingers.....
    ! fix glob_to_loc/loc_to_glob  mappings, then indices lists
    ! hmm, maybe we should just moe all of this onto a different level,
    ! have a specialized subroutine, and do it in the solver context???? 
    if (debug) write(0,*) 'spasb: renumbering glob_to_loc'
    do i=1, n_col
      desc_a%glob_to_loc(desc_a%loc_to_glob(desc_a%lprm(i))) = i  
    enddo
    if (debug) write(0,*) 'spasb: renumbering loc_to_glob'
    do i=1,desc_a%matrix_data(psb_m_) 
      j = desc_a%glob_to_loc(i)
      if (j>0) then 
        desc_a%loc_to_glob(j) = i
      endif
    enddo
    if (debug) write(0,*) 'spasb: renumbering halo_index'
    i=1
    kh=desc_a%halo_index(i)
    do while (kh /= -1) 
      i = i+1
      nh = desc_a%halo_index(i)
      do j = i+1, i+nh
        desc_a%halo_index(j) = &
             &desc_a%lprm(desc_a%halo_index(j))
      enddo
      i = i + nh + 1
      nh = desc_a%halo_index(i)
      do j= i+1, i+nh
        desc_a%halo_index(j) = &
             &desc_a%lprm(desc_a%halo_index(j))
      enddo
      i = i + nh + 1
      kh=desc_a%halo_index(i)
    enddo
    if (debug) write(0,*) 'spasb: renumbering ovrlap_index'
    i=1
    kh=desc_a%ovrlap_index(i)
    do while (kh /= -1) 
      i = i + 1
      nh = desc_a%ovrlap_index(i)
      do j= i+1, i+nh
        desc_a%ovrlap_index(j) = &
             &desc_a%lprm(desc_a%ovrlap_index(j))
      enddo
      i = i + nh + 1
      kh=desc_a%ovrlap_index(i)
    enddo
    if (debug) write(0,*) 'spasb: renumbering ovrlap_elem'
    i = 1
    kh=desc_a%ovrlap_elem(i)
    do while (kh /= -1)          
      desc_a%ovrlap_elem(i) = &
           &desc_a%lprm(desc_a%ovrlap_elem(i))
      i = i+2
      kh=desc_a%ovrlap_elem(i)         
    enddo
    if (debug) write(0,*) 'spasb: done renumbering'
    if (debug) then
      write(60+myrow,*) 'n_row ',n_row,' n_col',n_col, ' trans: ',trans
      do i=1,n_col
        write(60+myrow,*)i, ' lprm ', desc_a%lprm(i), ' iperm',iperm(i)
      enddo
      i=1
      kh = desc_a%halo_index(i)
      do while (kh /= -1) 
        write(60+myrow,*) i, kh 
        i = i+1
        kh = desc_a%halo_index(i)
      enddo
      close(60+myrow)
    end if
    
!!$    iperm(1) = 0
  else 
!!$    allocate(desc_a%lprm(1))
!!$    desc_a%lprm(1) = 0       
  endif    


  time(4) = mpi_wtime()
  time(4) = time(4) - time(3)
  if (debug) then 
    call dgamx2d(icontxt, psb_all_, psb_topdef_, ione, ione, time(4),&
         & ione,temp ,temp,-ione ,-ione,-ione)

    write (*, *) '         comm structs assembly: ', time(4)*1.d-3
  end if

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

end subroutine psb_dscren
