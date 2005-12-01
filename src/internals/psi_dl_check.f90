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
    outer: do 
      if (i >length_dl(proc)) exit outer
      proc2=dep_list(i,proc)
      if (proc2.ne.psb_no_comm_) then
        ! ...search proc in proc2's dep_list....
        j=1
        p2loop:do 
          if (j > length_dl(proc2)) exit p2loop
          if (dep_list(j,proc2) == proc) exit p2loop
          j=j+1
        enddo p2loop

        if (j > length_dl(proc2)) then
          ! ...add proc to proc2 s dep_list.....',proc,proc2
          length_dl(proc2)     = length_dl(proc2)+1
          if (length_dl(proc2) > size(dep_list,1)) then
            write(0,*)'error in crea_halo', proc2,proc,&
                 & length_dl(proc2),'>',size(dep_list,1)
          endif
          dep_list(length_dl(proc2),proc2) = proc
        else if (dep_list(j,proc2) /= proc) then 
          write(0,*) 'PSI_DL_CHECK This should not happen!!! ',&
               & j,proc2,dep_list(j,proc2),proc
        endif
      endif
      i=i+1
    enddo outer
  enddo

end subroutine psi_dl_check
