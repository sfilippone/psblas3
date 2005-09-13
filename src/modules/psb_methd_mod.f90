Module psb_methd_mod

  interface psb_cg
     subroutine psb_dcg(a,prec,b,x,eps,&
	  & desc_a,info,itmax,iter,err,itrace,istop)
       use psb_serial_mod
       use psb_descriptor_type
       use psb_prec_type
       type(psb_dspmat_type), intent(in)  :: a
       type(psb_desc_type), intent(in)    :: desc_a
       real(kind(1.d0)), intent(in)       :: b(:)
       real(kind(1.d0)), intent(inout)    :: x(:)
       real(kind(1.d0)), intent(in)       :: eps
       type(psb_dprec_type), intent(in)   :: prec
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: itmax, itrace,istop
       integer, optional, intent(out)     :: iter
       real(kind(1.d0)), optional, intent(out) :: err
     end subroutine psb_dcg
  end interface

  interface spb_bicg
     subroutine psb_dbicg(a,prec,b,x,eps,&
	  & desc_a,info,itmax,iter,err,itrace,istop)
       use psb_serial_mod
       use psb_descriptor_type
       use psb_prec_type
       type(psb_dspmat_type), intent(in)  :: a
       type(psb_desc_type), intent(in)    :: desc_a
       real(kind(1.d0)), intent(in)       :: b(:)
       real(kind(1.d0)), intent(inout)    :: x(:)
       real(kind(1.d0)), intent(in)       :: eps
       type(psb_dprec_type), intent(in)   :: prec
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: itmax, itrace,istop
       integer, optional, intent(out)     :: iter
       real(kind(1.d0)), optional, intent(out) :: err
     end subroutine psb_dbicg
  end interface

  interface psb_bicgstab
     subroutine psb_dcgstab(a,prec,b,x,eps,&
	  & desc_a,info,itmax,iter,err,itrace,istop)
       use psb_serial_mod
       use psb_descriptor_type
       use psb_prec_type
       type(psb_dspmat_type), intent(in)  :: a
       type(psb_desc_type), intent(in)    :: desc_a
       real(kind(1.d0)), intent(in)       :: b(:)
       real(kind(1.d0)), intent(inout)    :: x(:)
       real(kind(1.d0)), intent(in)       :: eps
       type(psb_dprec_type), intent(in)   :: prec
       integer, intent(out)               :: info
       integer, optional, intent(in)      :: itmax, itrace,istop
       integer, optional, intent(out)     :: iter
       real(kind(1.d0)), optional, intent(out) :: err
     end subroutine psb_dcgstab
  end interface

  interface psb_bicgstabl
    Subroutine psb_dcgstabl(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err, itrace,irst,istop)
      use psb_serial_mod
      use psb_descriptor_type
      Use psb_prec_type
!!$  parameters 
      Type(psb_dspmat_type), Intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      type(psb_dprec_type), intent(in)   :: prec
      Real(Kind(1.d0)), Intent(in)       :: b(:)
      Real(Kind(1.d0)), Intent(inout)    :: x(:)
      Real(Kind(1.d0)), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(Kind(1.d0)), Optional, Intent(out) :: err
    end subroutine psb_dcgstabl
  end interface

  interface psb_rgmres
    Subroutine psb_dgmresr(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
      use psb_serial_mod
      use psb_descriptor_type
      Use psb_prec_type
!!$  parameters 
      Type(psb_dspmat_type), Intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      type(psb_dprec_type), intent(in)   :: prec 
      Real(Kind(1.d0)), Intent(in)       :: b(:)
      Real(Kind(1.d0)), Intent(inout)    :: x(:)
      Real(Kind(1.d0)), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(Kind(1.d0)), Optional, Intent(out) :: err
    end subroutine psb_dgmresr
  end interface

  interface psb_cgs
    subroutine psb_dcgs(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
!!$  parameters 
      type(psb_dspmat_type), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a 
      type(psb_dprec_type), intent(in)   :: prec 
      real(kind(1.d0)), intent(in)       :: b(:)
      real(kind(1.d0)), intent(inout)    :: x(:)
      real(kind(1.d0)), intent(in)       :: eps
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(kind(1.d0)), optional, intent(out) :: err
    end subroutine psb_dcgs
  end interface

end module psb_methd_mod


  
