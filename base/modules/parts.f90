module psb_parts_mod
  interface 
     subroutine psb_parts(glob_index,nrow,np,pv,nv)
       integer, intent (in)  :: glob_index,np,nrow
       integer, intent (out) :: nv, pv(*)
     end subroutine psb_parts
  end interface
end module psb_parts_mod
