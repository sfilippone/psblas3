! File:  psb_dneigh.f90 
! Subroutine: 
! Parameters:

subroutine psb_dneigh(a,idx,neigh,n,info,lev)

  use psb_realloc_mod
  use psb_const_mod
  use psb_spmat_type
  implicit none


  type(psb_dspmat_type), intent(in) :: a   ! the sparse matrix
  integer, intent(in)               :: idx ! the index whose neighbours we want to find
  integer, intent(out)              :: n, info   ! the number of neighbours and the info
  integer, pointer                  :: neigh(:) ! the neighbours
  integer, optional                 :: lev ! level of neighbours to find

  interface psb_spgtrow
     subroutine psb_dspgtrow(irw,a,b,info,append,iren,lrw)
       use psb_spmat_type
       type(psb_dspmat_type), intent(in) :: a
       integer, intent(in)       :: irw
       type(psb_dspmat_type), intent(inout)    :: b
       logical, intent(in), optional :: append
       integer, intent(in), target, optional :: iren(:)
       integer, intent(in), optional :: lrw
       integer, intent(out)  :: info
     end subroutine psb_dspgtrow
  end interface

  type(psb_dspmat_type) :: atmp
  integer, pointer :: tmpn(:)
  integer :: level, dim, i, j, k, r, c, brow,&
       & elem_pt, ii, n1, col_idx, ne, err_act
  character(len=20)                 :: name, ch_err

  name='psb_dneigh'
  info  = 0
  call psb_erractionsave(err_act)

  n    = 0
  info = 0
  if(present(lev)) then 
    if(lev.le.2) then
      level=lev
    else
      write(0,'("Too many levels!!!")')
      return
    endif
  else
    level=1
  end if

  if ((a%fida /= 'CSR')) then 

     call psb_spall(atmp,1,info)
     if(info.ne.0) then
        call psb_errpush(4010,name)
        goto 9999
     end if
     call psb_spgtrow(idx,a,atmp,info)
     if(info.ne.0) then
        call psb_errpush(4010,name)
        goto 9999
     end if
     dim=atmp%infoa(psb_nnz_)
     allocate(tmpn(dim))
     tmpn(1:dim)=atmp%ia2(1:dim)
     
     if(level.eq.2) then
        do i=1,dim
           if ((1<=tmpn(i)).and.(tmpn(i)<=a%m).and.(tmpn(i).ne.idx)) then
              call psb_spgtrow(tmpn(i),a,atmp,info,append=.true.)
              if(info.ne.0) then
                 call psb_errpush(4010,name)
                 goto 9999
              end if
           end if
        end do
     end if

     
     dim=atmp%infoa(psb_nnz_)
     if(dim.gt.size(neigh)) &
          & call psb_realloc(dim,neigh,info)
     if(info.ne.0) then
        call psb_errpush(4010,name)
        goto 9999
     end if
     call psb_spfree(atmp,info)
     if(info.ne.0) then
        call psb_errpush(4010,name)
        goto 9999
     end if
     call psb_nullify_sp(atmp)
     deallocate(tmpn)
     
  else if(a%fida.eq.'CSR') then

    dim=0
    if(level.eq.1) dim=(a%ia2(idx+1)-a%ia2(idx))
    if(dim >size(neigh)) call psb_realloc(dim,neigh,info)
    if(info.ne.izero) then
       info=4010
       ch_err='psrealloc'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    endif
    
    n=0
    if(level.eq.1) then
      do i=a%ia2(idx), a%ia2(idx+1)-1
        n=n+1
        neigh(n)=a%ia1(i)
      end do

    else

      do i=a%ia2(idx), a%ia2(idx+1)-1

        j=a%ia1(i)
        if ((1<=j).and.(j<=a%m).and.(j.ne.idx)) then
           
           dim=dim+ a%ia2(j+1)-a%ia2(j)
           if(dim.gt.size(neigh)) call psb_realloc(dim,neigh,info)
           if(info.ne.izero) then
              info=4010
              ch_err='psrealloc'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
           endif
           
          do k=a%ia2(j), a%ia2(j+1)-1
            n=n+1
            neigh(n)=a%ia1(k)
          end do
        end if
      end do

    end if

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

end subroutine psb_dneigh
