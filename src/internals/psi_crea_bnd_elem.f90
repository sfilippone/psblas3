subroutine psi_crea_bnd_elem(desc_a,info)

  use psb_descriptor_type
  use psb_error_mod
  implicit none

  type(psb_desc_type)  :: desc_a
  integer, intent(out) :: info

  integer, pointer :: work(:)
  integer :: i, j, nr, ns, k, irv, err_act
  character(len=20)    :: name, ch_err

  info = 0
  name='psi_crea_bnd_elem'
  call psb_erractionsave(err_act)
  
  allocate(work(size(desc_a%halo_index)),stat=info)
  if (info /= 0 ) then 
     info = 4000
     call psb_errpush(info,name)
     goto 9999
  end if

  i=0
  j=1
  do while(desc_a%halo_index(j) /= -1) 

     nr = desc_a%halo_index(j+1)
     ns = desc_a%halo_index(j+1+nr+1)
     do k=1, ns
        i = i + 1
        work(i) = desc_a%halo_index(j+1+nr+1+k)
     enddo
     j  = j + 1 + ns + 1 + nr + 1
  enddo

  if (i>0) then 
     call isr(i,work)
     j=1
     irv = work(1)
     do k=2, i
        if (work(k) /= irv) then
           irv = work(k)
           j = j + 1 
           work(j) = work(k)
        endif
     enddo
  else 
     j = 0
  endif

  allocate(desc_a%bnd_elem(j+1))
  if (.false.) then 
     desc_a%bnd_elem(1) = j 
     desc_a%bnd_elem(2:j+1) = work(1:j)
  else
     desc_a%bnd_elem(1:j) = work(1:j)
     desc_a%bnd_elem(j+1) = -1
  endif

  deallocate(work)
  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error()
     return
  end if
  return
end subroutine psi_crea_bnd_elem
