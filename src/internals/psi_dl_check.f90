subroutine psi_dl_check(dep_list,dl_lda,np,length_dl)

  use psb_const_mod
  implicit none

  integer  :: np,dl_lda,length_dl(0:np)
  integer  :: dep_list(dl_lda,0:np)
  ! locals
  integer  :: proc, proc2, i, j

! ...i must order communication in in halo

! ...if in dep_list of process i there is j
!     and in dep_list of process j there isn't i,
!     add to it process i...

  do proc=0,np-1
     i=1
     do while (i.le.length_dl(proc))
        proc2=dep_list(i,proc)
        if (proc2.ne.psb_no_comm_) then
           ! ...search proc in proc2's dep_list....
           j=1
           do while ((j.le.length_dl(proc2).and.&
                & dep_list(j,proc2).ne.proc))
              j=j+1
           enddo
           if ((dep_list(j,proc2).ne.proc).or.&
                & (j.gt.length_dl(proc2))) then

              ! ...proc not found...
              ! ...add proc to proc2's dep_list.....
              length_dl(proc2)=length_dl(proc2)+1
              if (length_dl(proc2).gt.size(dep_list,1)) then
                 write(0,*)'error in crea_halo', proc2,&
                      & length_dl(proc2),'>',size(dep_list,1)
              endif
              dep_list(length_dl(proc2),proc2)=proc
           endif
        endif
        i=i+1
     enddo
  enddo

end subroutine psi_dl_check
