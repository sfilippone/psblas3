subroutine psb_cest(afmt, nnz, lia1, lia2, lar, up, info)

  use psb_error_mod
  implicit none

  !     .. scalar arguments ..
  integer           ::  nnz, lia1, lia2, lar, info
  character         ::  up
  !     .. array arguments..
  character(len=5)  ::  afmt
  integer           ::  int_val(5), err_act
  character(len=20) ::  name

  name = 'cest'      
  call psb_erractionsave(err_act)

  if (afmt.eq.'???') then 
     afmt = fidef
  endif

  if (up.eq.'y') then
     if (afmt.eq.'JAD') then 
        lia1 = 2*(nnz + nnz/5) +1000
        lia2 = 2*(nnz + nnz/5) +1000
        lar = nnz + nnz/5
     else if (afmt.eq.'COO') then 
        lia1 = nnz
        lia2 = 2*nnz + 1000
        lar = nnz
     else if(afmt.eq.'CSR') then
        lia1 = nnz
        lia2 = 2*nnz + 1000
        lar = nnz
     else
        info = 3012
        call psb_errpush(info,name)
        goto 9999
     endif

  else if (up.eq.'n') then

     if (afmt.eq.'jad') then 
        lia1 = nnz + nnz/5
        lia2 = nnz + nnz/5
        lar = nnz + nnz/5
     else if (afmt.eq.'coo') then 
        lia1 = nnz
        lia2 = nnz
        lar = nnz
     else if(afmt.eq.'csr') then
        lia1 = nnz
        lia2 = nnz
        lar = nnz
     else
        info = 3012
        call psb_errpush(info,name)
        goto 9999
     endif

  else

     info = 3012
     call psb_errpush(info,name,int_val)
     goto 9999

  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if ( err_act .ne. 0 ) then 
     call psb_error()
     return
  endif

  return

end subroutine psb_cest

