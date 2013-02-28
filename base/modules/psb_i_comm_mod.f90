!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
module psb_i_comm_mod

  interface psb_ovrl
    subroutine  psb_iovrlm(x,desc_a,info,jx,ik,work,update,mode)
      use psb_desc_mod
      integer(psb_ipk_), intent(inout), target   :: x(:,:)
      type(psb_desc_type), intent(in)            :: desc_a
      integer(psb_ipk_), intent(out)                       :: info
      integer(psb_ipk_), intent(inout), optional, target :: work(:)
      integer(psb_ipk_), intent(in), optional              :: update,jx,ik,mode
    end subroutine psb_iovrlm
    subroutine  psb_iovrlv(x,desc_a,info,work,update,mode)
      use psb_desc_mod
      integer(psb_ipk_), intent(inout), target   :: x(:)
      type(psb_desc_type), intent(in)            :: desc_a
      integer(psb_ipk_), intent(out)                       :: info
      integer(psb_ipk_), intent(inout), optional, target :: work(:)
      integer(psb_ipk_), intent(in), optional              :: update,mode
    end subroutine psb_iovrlv
    subroutine  psb_iovrl_vect(x,desc_a,info,work,update,mode)
      use psb_desc_mod
      use psb_i_vect_mod
      type(psb_i_vect_type), intent(inout)    :: x
      type(psb_desc_type), intent(in)         :: desc_a
      integer(psb_ipk_), intent(out)                    :: info
      integer(psb_ipk_), intent(inout), optional, target :: work(:)
      integer(psb_ipk_), intent(in), optional           :: update,mode
    end subroutine psb_iovrl_vect
  end interface psb_ovrl

  interface psb_halo
    subroutine  psb_ihalom(x,desc_a,info,alpha,jx,ik,work,tran,mode,data)
      use psb_desc_mod
      integer(psb_ipk_), intent(inout), target :: x(:,:)
      type(psb_desc_type), intent(in)          :: desc_a
      integer(psb_ipk_), intent(out)                     :: info
      integer(psb_ipk_), intent(in), optional  :: alpha
      integer(psb_ipk_), target, optional, intent(inout) :: work(:)
      integer(psb_ipk_), intent(in), optional           :: mode,jx,ik,data
      character, intent(in), optional         :: tran
    end subroutine psb_ihalom
    subroutine  psb_ihalov(x,desc_a,info,alpha,work,tran,mode,data)
      use psb_desc_mod
      integer(psb_ipk_), intent(inout)        :: x(:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer(psb_ipk_), intent(out)                    :: info
      integer(psb_ipk_), intent(in), optional :: alpha
      integer(psb_ipk_), target, optional, intent(inout) :: work(:)
      integer(psb_ipk_), intent(in), optional           :: mode,data
      character, intent(in), optional         :: tran
    end subroutine psb_ihalov
    subroutine  psb_ihalo_vect(x,desc_a,info,alpha,work,tran,mode,data)
      use psb_desc_mod
      use psb_i_vect_mod
      type(psb_i_vect_type), intent(inout)   :: x
      type(psb_desc_type), intent(in)         :: desc_a
      integer(psb_ipk_), intent(out)                    :: info
      integer(psb_ipk_), intent(in), optional    :: alpha
      integer(psb_ipk_), target, optional, intent(inout) :: work(:)
      integer(psb_ipk_), intent(in), optional           :: mode,data
      character, intent(in), optional         :: tran
    end subroutine psb_ihalo_vect
  end interface psb_halo


  interface psb_scatter
    subroutine  psb_iscatterm(globx, locx, desc_a, info, root)
      use psb_desc_mod
      integer(psb_ipk_), intent(out) :: locx(:,:)
      integer(psb_ipk_), intent(in)  :: globx(:,:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), intent(in), optional    :: root
    end subroutine psb_iscatterm
    subroutine  psb_iscatterv(globx, locx, desc_a, info, root)
      use psb_desc_mod
      integer(psb_ipk_), intent(out) :: locx(:)
      integer(psb_ipk_), intent(in)  :: globx(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), intent(in), optional    :: root
    end subroutine psb_iscatterv
  end interface psb_scatter

  interface psb_gather
!!$    subroutine  psb_isp_allgather(globa, loca, desc_a, info, root, dupl,keepnum,keeploc)
! !$      use psb_desc_mod
! !$      use psb_mat_mod
! !$      implicit none
! !$      type(psb_ispmat_type), intent(inout) :: loca
! !$      type(psb_ispmat_type), intent(out)   :: globa
! !$      type(psb_desc_type), intent(in) :: desc_a
! !$      integer(psb_ipk_), intent(out)            :: info
! !$      integer(psb_ipk_), intent(in), optional   :: root,dupl
! !$      logical, intent(in), optional   :: keepnum,keeploc
! !$    end subroutine psb_isp_allgather
    subroutine psb_igatherm(globx, locx, desc_a, info, root)
      use psb_desc_mod
      integer(psb_ipk_), intent(in)  :: locx(:,:)
      integer(psb_ipk_), intent(out), allocatable  :: globx(:,:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), intent(in), optional    :: root
    end subroutine psb_igatherm
    subroutine  psb_igatherv(globx, locx, desc_a, info, root)
      use psb_desc_mod
      integer(psb_ipk_), intent(in)  :: locx(:)
      integer(psb_ipk_), intent(out), allocatable  :: globx(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), intent(in), optional    :: root
    end subroutine psb_igatherv
    subroutine  psb_igather_vect(globx, locx, desc_a, info, root)
      use psb_desc_mod
      use psb_i_vect_mod
      type(psb_i_vect_type), intent(inout) :: locx
      integer(psb_ipk_), intent(out), allocatable :: globx(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: root
    end subroutine psb_igather_vect
  end interface psb_gather

end module psb_i_comm_mod
