module psb_comm_mod

  interface psb_ovrl
     subroutine  psb_dovrlm(x,desc_a,info,jx,ik,work,choice,update_type)
       use psb_descriptor_type
       real(kind(1.d0)), intent(inout)           :: x(:,:)
       type(psb_desc_type), intent(in)           :: desc_a
       integer, intent(out)                      :: info
       real(kind(1.d0)), intent(inout), optional :: work(:)
       logical, intent(in), optional             :: choice
       integer, intent(in), optional             :: update_type,jx,ik
     end subroutine psb_dovrlm
     subroutine  psb_dovrlv(x,desc_a,info,work,choice,update_type)
       use psb_descriptor_type
       real(kind(1.d0)), intent(inout)           :: x(:)
       type(psb_desc_type), intent(in)           :: desc_a
       integer, intent(out)                      :: info
       real(kind(1.d0)), intent(inout), optional :: work(:)
       logical, intent(in), optional             :: choice
       integer, intent(in), optional             :: update_type
     end subroutine psb_dovrlv
  end interface

  interface psb_halo
     subroutine  psb_dhalom(x,desc_a,info,alpha,jx,ik,work,tran,mode)
       use psb_descriptor_type
       real(kind(1.d0)), intent(inout)           :: x(:,:)
       type(psb_desc_type), intent(in)           :: desc_a
       integer, intent(out)                      :: info
       real(kind(1.d0)), intent(in), optional    :: alpha
       real(kind(1.d0)), intent(inout), optional :: work(:)
       integer, intent(in), optional             :: mode,jx,ik
       character, intent(in), optional           :: tran
     end subroutine psb_dhalom
     subroutine  psb_dhalov(x,desc_a,info,alpha,work,tran,mode)
       use psb_descriptor_type
       real(kind(1.d0)), intent(inout)           :: x(:)
       type(psb_desc_type), intent(in)           :: desc_a
       integer, intent(out)                      :: info
       real(kind(1.d0)), intent(in), optional    :: alpha
       real(kind(1.d0)), intent(inout), optional :: work(:)
       integer, intent(in), optional             :: mode
       character, intent(in), optional           :: tran
     end subroutine psb_dhalov
     subroutine  psb_ihalom(x,desc_a,info,alpha,jx,ik,work,tran,mode)
       use psb_descriptor_type
       integer, intent(inout) :: x(:,:)
       type(psb_desc_type), intent(in)        :: desc_a
       integer, intent(out)                   :: info
       real(kind(1.d0)), intent(in), optional :: alpha
       integer, intent(inout), optional       :: work(:)
       integer, intent(in), optional          :: mode,jx,ik
       character, intent(in), optional        :: tran
     end subroutine psb_ihalom
     subroutine  psb_ihalov(x,desc_a,info,alpha,work,tran,mode)
       use psb_descriptor_type
       integer, intent(inout)                 :: x(:)
       type(psb_desc_type), intent(in)        :: desc_a
       integer, intent(out)                   :: info
       real(kind(1.d0)), intent(in), optional :: alpha
       integer, intent(inout), optional       :: work(:)
       integer, intent(in), optional          :: mode
       character, intent(in), optional        :: tran
     end subroutine psb_ihalov
  end interface


  interface psb_dscatter
     subroutine  psb_dscatterm(globx, locx, desc_a, info, iroot,&
          & iiglobx, ijglobx, iilocx,ijlocx,ik)
       use psb_descriptor_type
       real(kind(1.d0)), intent(out)    :: locx(:,:)
       real(kind(1.d0)), intent(in)     :: globx(:,:)
       type(psb_desc_type), intent(in)  :: desc_a
       integer, intent(out)             :: info
       integer, intent(in), optional    :: iroot,iiglobx,&
            & ijglobx,iilocx,ijlocx,ik
     end subroutine psb_dscatterm
     subroutine  psb_dscatterv(globx, locx, desc_a, info, iroot)
       use psb_descriptor_type
       real(kind(1.d0)), intent(out)    :: locx(:)
       real(kind(1.d0)), intent(in)     :: globx(:)
       type(psb_desc_type), intent(in)  :: desc_a
       integer, intent(out)             :: info
       integer, intent(in), optional    :: iroot
     end subroutine psb_dscatterv
  end interface

  interface psb_dgather
     subroutine  psb_dgatherm(globx, locx, desc_a, info, iroot,&
          & iiglobx, ijglobx, iilocx,ijlocx,ik)
       use psb_descriptor_type
       real(kind(1.d0)), intent(in)    :: locx(:,:)
       real(kind(1.d0)), intent(out)   :: globx(:,:)
       type(psb_desc_type), intent(in) :: desc_a
       integer, intent(out)            :: info
       integer, intent(in), optional   :: iroot, iiglobx, ijglobx, iilocx, ijlocx, ik
     end subroutine psb_dgatherm
     subroutine  psb_dgatherv(globx, locx, desc_a, info, iroot,&
          & iiglobx, iilocx)
       use psb_descriptor_type
       real(kind(1.d0)), intent(in)    :: locx(:,:)
       real(kind(1.d0)), intent(out)   :: globx(:,:)
       type(psb_desc_type), intent(in) :: desc_a
       integer, intent(out)            :: info
       integer, intent(in), optional   :: iroot, iiglobx, iilocx
     end subroutine psb_dgatherv
  end interface
  
end module psb_comm_mod
