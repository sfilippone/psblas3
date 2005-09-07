! File: psb_dptins.f90
!
! Subroutine: psb_dptins
!    insert sparse submatrix to sparse matrix structure for psblas
!    routines
! 
! Parameters: 
!    ia      - integer.                       a global-row corresponding to position at which blck submatrix must be inserted.
!    ja      - integer.                       a global-col corresponding to position at which blck submatrix must be inserted.
!    blck    - type(<psb_dspmat_type>).       The source sparse submatrix.  
!    desc_a  - type(<psb_desc_type>).         The communication descriptor.
!    info    - integer.                       Eventually returns an error code
subroutine psb_dptins(ia,ja,blck,desc_a,info)
  use psb_descriptor_type
  use psb_spmat_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  implicit none

  !....parameters...
  type(psb_desc_type), intent(inout)    ::  desc_a
  integer, intent(in)                   ::  ia,ja
  type(psb_dspmat_type), intent(in)     ::  blck
  integer,intent(out)                   ::  info

  !locals.....

  interface
    subroutine dcsins(m,n,fida,descra,a,ia1,ia2,infoa,&
         & ia,ja,latot,lia1tot,lia2tot,&
         &fidh,descrh,h,ih1,ih2,infoh,ih,jh,work,lwork,ierror)
      implicit none
      !      .. scalar arguments ..
      integer, intent(in) :: m, n, lwork, latot,lia1tot,lia2tot,ia,ja,ih,jh
      integer, intent(out) :: ierror
      !      .. array arguments ..
      double precision, intent(in)    :: h(*)
      double precision, intent(inout) :: a(*), work(*)
      integer, intent(in) :: ih1(*), ih2(*), infoh(10)                         
      integer, intent(inout) :: ia1(*), ia2(*), infoa(10)
      character, intent(in)  :: fida*5, fidh*5,descra*11, descrh*11
    end subroutine dcsins
  end interface

  integer :: i,icontxt,nprocs ,glob_row,row,&
       & k ,start_row,end_row,int_err(5),&
       & first_loc_row,n_row,j, ierror,locix,locjx,&
       & allocated_prcv,dectype,mglob, nnza, err_act
  integer,pointer        :: prcv(:), tia1(:),tia2(:), temp(:)
  integer                :: nprow,npcol, me ,mypcol, iflag, isize, irlc
  integer                ::  m,n, pnt_halo,ncol, nh, ip
  type(psb_dspmat_type)          ::  a
  real(kind(1.d0)),pointer       :: workarea(:),taspk(:)
  logical, parameter :: debug=.false.
  integer, parameter :: nrlcthr=3
  integer, save :: irlcmin,nrlc
  data irlcmin/500/,nrlc/0/
  character(len=20)   :: name, ch_err

  info=0
  call psb_erractionsave(err_act)
  name = 'psb_dptins'

  locix=1
  locjx=1
  icontxt = desc_a%matrix_data(psb_ctxt_)
  dectype = desc_a%matrix_data(psb_dec_type_)
  mglob   = desc_a%matrix_data(psb_m_)
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

  if (.not.psb_is_bld_dec(dectype)) then 
    info = 3110
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(workarea(1),prcv(nprow),stat=info)
  if (info.ne.0) then
    info = 2023
    call psb_errpush(info,name)
    goto 9999
  end if
  call psb_spall(a,size(blck%aspk),info)
  if (info.ne.0) then
    info = 4010
    ch_err='spall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  allocated_prcv = 1
  a%infoa(:)     = 0
  a%fida         = 'coo'
  a%descra       = 'gun'
  n_row          = desc_a%matrix_data(psb_n_row_)
  nnza           = a%infoa(psb_nnz_) 
  m              = blck%m
  n              = blck%k
  row            = ia
  i              = 1
  do while (i.le.m)
    !loop over all blck's rows
    ! row actual block row 

    row      = locix+i-1
    glob_row = ia+i-1
    if (glob_row > mglob) exit 
    if (debug) then 
      write(0,*) 'ptins: inserting ',glob_row
    endif
    k = desc_a%glob_to_loc(glob_row)
    if (k.gt.0) then
      start_row     = row
      first_loc_row = k
!!$      do while ((i.lt.m).and.&
!!$           & (desc_a%glob_to_loc(ia+i).gt.0))
!!$         i=i+1
!!$      enddo
      do 
        if (i>=m) exit
        if ((ia+i)>mglob) exit
        if (desc_a%glob_to_loc(ia+i) <=0 ) exit
        i=i+1
      enddo

      end_row=locix+i-1      
      ! insert blck submatrix in 'coo' format
      call dcsins(end_row-start_row+1,n,a%fida,a%descra,a%aspk,&
           & a%ia1,a%ia2,a%infoa,first_loc_row, ja,&
           & size(a%aspk),size(a%ia1),size(a%ia2),&
           & blck%fida,blck%descra,blck%aspk,blck%ia1,blck%ia2,&
           & blck%infoa,start_row,locjx,workarea,size(workarea),&
           & info)

      if (info.ne.0) then

        if (info.eq.60) then 
          ! try reallocating
          irlc = irlcmin
          do while (info.eq.60)
            if (debug) write(*,*) "attempting reallocation with",irlc

            isize = size(a%ia1)                    
            allocate(tia1(isize+irlc),stat=info)
            if (info.ne.0) goto 9998
            tia1(1:isize) = a%ia1(1:isize)
            deallocate(a%ia1,stat=info)
            if (info.ne.0) goto 9998
            a%ia1 => tia1
            nullify(tia1)

            isize = size(a%ia2)                    
            allocate(tia2(isize+irlc),stat=info)
            if (info.ne.0) goto 9998
            tia2(1:isize) = a%ia2(1:isize)
            deallocate(a%ia2,stat=info)
            if (info.ne.0) goto 9998
            a%ia2 => tia2
            nullify(tia2)

            isize = size(a%aspk)                    
            allocate(taspk(isize+irlc),stat=info)
            if (info.ne.0) goto 9998
            taspk(1:isize) = a%aspk(1:isize)
            deallocate(a%aspk,stat=info)
            if (info.ne.0) goto 9998
            a%aspk => taspk
            nullify(taspk)

9998        if (info.ne.0) then
               info=4000
               call psb_errpush(info,name)
               goto 9999
            end if

            ! insert blck submatrix in 'coo' format
            call dcsins(end_row-start_row+1,n,a%fida,a%descra,a%aspk,&
                 & a%ia1,a%ia2,a%infoa,first_loc_row, ja,&
                 & size(a%aspk), size(a%ia1),size(a%ia2),&
                 & blck%fida,blck%descra,blck%aspk,blck%ia1,blck%ia2,&
                 & blck%infoa,start_row, locjx,workarea,size(workarea),&
                 & info)

            if (info.eq.60) irlc = irlc*2                    
          enddo
          ! if we get here, it means we succesfully reallocated. 
          nrlc = nrlc+1
          if (nrlc .ge. nrlcthr) then
            nrlc = 0
            irlcmin = irlcmin * 2
          endif

       else
          info = 4010
          ch_err='spall'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
       endif
    endif
   endif
   ! next blck's row
   i=i+1
  enddo


  if (.not.associated(desc_a%halo_index)) then
    allocate(desc_a%halo_index(irlcmin))
    desc_a%halo_index(:) = -1
  endif
  pnt_halo=1
  do while (desc_a%halo_index(pnt_halo) .ne.  -1 )
    pnt_halo = pnt_halo + 1
  end do 
  ncol = desc_a%matrix_data(psb_n_col_)

  isize = size(desc_a%halo_index)
  do i = nnza+1,a%infoa(psb_nnz_)
    ip = a%ia2(i)
    k  = desc_a%glob_to_loc(ip)
    if (k.lt.-nprow) then
      k = k + nprow
      k = - k - 1
      ncol = ncol + 1      
      desc_a%glob_to_loc(ip)   = ncol
      isize = size(desc_a%loc_to_glob)
      if (ncol > isize) then 
        nh = ncol + irlcmin
        call psb_realloc(nh,desc_a%loc_to_glob,info,pad=-1)
        if (me==0) then 
          if (debug) write(0,*) 'done realloc ',nh
        end if
        if (info /= 0) then
           info=4000
           call psb_errpush(info,name)
           goto 9999
        end if
        isize = nh
      endif
      desc_a%loc_to_glob(ncol) = ip
      isize = size(desc_a%halo_index)
      if ((pnt_halo+3).gt.isize) then
        nh = isize + irlcmin
        call psb_realloc(nh,desc_a%halo_index,info,pad=-1)
        if (info /= 0) then
           info=4000
           call psb_errpush(info,name)
           goto 9999
        end if
        isize = nh 
      endif
      desc_a%halo_index(pnt_halo)   = k
      desc_a%halo_index(pnt_halo+1) = 1
      desc_a%halo_index(pnt_halo+2) = ncol
      pnt_halo                      = pnt_halo + 3
    endif
  enddo

  desc_a%matrix_data(psb_n_col_) = ncol



  if (allocated_prcv.eq.1) then
    call psb_spfree(a,info)
        if (info /= 0) then
           info=4000
           call psb_errpush(info,name)
           goto 9999
        end if
    deallocate(prcv,workarea,stat=info)
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
end subroutine psb_dptins

