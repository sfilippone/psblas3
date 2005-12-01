subroutine psi_compute_size(desc_data,&
     & index_in, dl_lda, info)

  use psb_const_mod
  use psb_error_mod
  implicit none

  !     ....scalars parameters....
  integer  :: info, dl_lda
  !     .....array parameters....
  integer  :: desc_data(:), index_in(:)
  !     ....local scalars....      
  integer  :: i,npcol,nprow,mycol,myrow,proc,counter, max_index
  integer  :: icontxt, err, err_act, np
  !     ...local array...
  integer  :: exch(2)
  integer  :: int_err(5)
  integer, pointer :: counter_recv(:), counter_dl(:)

  !     ...parameters
  logical, parameter :: debug=.false.
  character(len=20)  :: name

  name='psi_compute_size'
  call psb_get_erraction(err_act)

  info = 0
  icontxt = desc_data(psb_ctxt_)

  call blacs_gridinfo(icontxt,nprow,npcol,myrow,mycol)
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

  np=nprow
  allocate(counter_dl(0:np-1),counter_recv(0:np-1))
  !     ..initialize counters...
  do i=0,np-1
     counter_recv(i)=0
     counter_dl(i)=0
  enddo

  !     ....verify local correctness of halo_in....
  i=1
  do while (index_in(i).ne.-1)
     proc=index_in(i)
     if ((proc.gt.np-1).or.(proc.lt.0)) then
        info = 115
        int_err(1) = 11
        int_err(2) = proc
        call psb_errpush(info,name,i_err=int_err)
        goto 9999
     endif
     counter_dl(proc)=1

     !        ..update no of elements to receive from proc proc..         
     counter_recv(proc)=counter_recv(proc)+&
          & index_in(i+1)

     i=i+index_in(i+1)+2
  enddo

  !     ...computing max_halo: max halo points to be received from
  !                            same processor
  max_index=0
  dl_lda=0

  do i=0,np-1
     if (counter_recv(i).gt.max_index) max_index = counter_recv(i)
     if (counter_dl(i).eq.1) dl_lda = dl_lda+1
  enddo

  !     computing max global value of dl_lda
  call igamx2d(icontxt, psb_all_, psb_topdef_, 1, ione, dl_lda, &
       &1, counter, counter, -ione ,-ione,-ione)

  if (debug) then 
     write(0,*) 'psi_compute_size: ',dl_lda
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psi_compute_size

         

