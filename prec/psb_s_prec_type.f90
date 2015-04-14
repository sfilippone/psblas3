!!$ 
!!$              Parallel Sparse BLAS  version 3.4
!!$    (C) Copyright 2006, 2010, 2015
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	Module to   define PREC_DATA,           !!
!!      structure for preconditioning.          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module psb_s_prec_type

  use psb_prec_const_mod
  use psb_s_base_prec_mod

  type psb_sprec_type
    class(psb_s_base_prec_type), allocatable :: prec
  contains
    procedure, pass(prec)               :: psb_s_apply1_vect
    procedure, pass(prec)               :: psb_s_apply2_vect
    procedure, pass(prec)               :: psb_s_apply2v
    procedure, pass(prec)               :: psb_s_apply1v
    generic, public                     :: apply => psb_s_apply2v, psb_s_apply1v,&
         & psb_s_apply1_vect, psb_s_apply2_vect
    procedure, pass(prec)               :: sizeof => psb_sprec_sizeof
    procedure, pass(prec)               :: clone  => psb_s_prec_clone
    procedure, pass(prec)               :: free   => psb_s_prec_free
  end type psb_sprec_type

  interface psb_precfree
    module procedure psb_s_precfree
  end interface


  interface psb_precdescr
    module procedure psb_sfile_prec_descr
  end interface

  interface psb_precdump
    module procedure psb_s_prec_dump
  end interface

  interface psb_sizeof
    module procedure psb_sprec_sizeof
  end interface

  interface 
    subroutine psb_s_apply2_vect(prec,x,y,desc_data,info,trans,work)
      import :: psb_ipk_, psb_desc_type, psb_sprec_type, psb_s_vect_type, psb_spk_
      type(psb_desc_type),intent(in)       :: desc_data
      class(psb_sprec_type), intent(inout) :: prec
      type(psb_s_vect_type),intent(inout)  :: x
      type(psb_s_vect_type),intent(inout)  :: y
      integer(psb_ipk_), intent(out)                 :: info
      character(len=1), optional           :: trans
      real(psb_spk_),intent(inout), optional, target :: work(:)
    end subroutine psb_s_apply2_vect
  end interface
  
  interface 
    subroutine psb_s_apply1_vect(prec,x,desc_data,info,trans,work)
      import :: psb_ipk_, psb_desc_type, psb_sprec_type, psb_s_vect_type, psb_spk_
        type(psb_desc_type),intent(in)       :: desc_data
      class(psb_sprec_type), intent(inout) :: prec
      type(psb_s_vect_type),intent(inout)  :: x
      integer(psb_ipk_), intent(out)                 :: info
      character(len=1), optional           :: trans
      real(psb_spk_),intent(inout), optional, target :: work(:)
    end subroutine psb_s_apply1_vect
  end interface
  
  interface
    subroutine psb_s_apply2v(prec,x,y,desc_data,info,trans,work)
      import :: psb_ipk_, psb_desc_type, psb_sprec_type, psb_s_vect_type, psb_spk_
      type(psb_desc_type),intent(in)    :: desc_data
      class(psb_sprec_type), intent(in) :: prec
      real(psb_spk_),intent(inout)   :: x(:)
      real(psb_spk_),intent(inout)   :: y(:)
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), optional        :: trans
      real(psb_spk_),intent(inout), optional, target :: work(:)
    end subroutine psb_s_apply2v
  end interface
  
  interface 
    subroutine psb_s_apply1v(prec,x,desc_data,info,trans)
      import :: psb_ipk_, psb_desc_type, psb_sprec_type, psb_s_vect_type, psb_spk_
      type(psb_desc_type),intent(in)    :: desc_data
      class(psb_sprec_type), intent(in) :: prec
      real(psb_spk_),intent(inout)   :: x(:)
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), optional        :: trans
    end subroutine psb_s_apply1v
  end interface
  
contains

  subroutine psb_sfile_prec_descr(p,iout)
    use psb_base_mod
    implicit none 
    type(psb_sprec_type), intent(in) :: p
    integer(psb_ipk_), intent(in), optional    :: iout
    integer(psb_ipk_) :: iout_,info
    character(len=20) :: name='prec_descr' 
    
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if

    if (.not.allocated(p%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
    end if
    call p%prec%precdescr(iout)
    
  end subroutine psb_sfile_prec_descr

  subroutine psb_s_prec_dump(prec,info,prefix,head)
    implicit none 
    type(psb_sprec_type), intent(in) :: prec
    integer(psb_ipk_), intent(out)             :: info
    character(len=*), intent(in), optional :: prefix,head
    !  len of prefix_ 

    info = 0

    if (.not.allocated(prec%prec)) then 
      info = -1
      write(psb_err_unit,*) 'Trying to dump a non-built preconditioner'
      return
    end if
    
    call prec%prec%dump(info,prefix,head)
    
    
  end subroutine psb_s_prec_dump


  subroutine psb_s_precfree(p,info)
    use psb_base_mod
    implicit none 
    type(psb_sprec_type), intent(inout) :: p
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_) :: me, err_act,i
    character(len=20)   :: name
    if(psb_get_errstatus() /= 0) return 
    info=psb_success_
    name = 'psb_precfree'
    call psb_erractionsave(err_act)

    me=-1
    call p%free(info)

    if (info /= 0) goto 9999
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    
    return

  end subroutine psb_s_precfree

  subroutine psb_s_prec_free(prec,info)
    use psb_base_mod
    implicit none 
    class(psb_sprec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)         :: info
    integer(psb_ipk_) :: me, err_act,i
    character(len=20)   :: name
    if(psb_get_errstatus() /= 0) return 
    info=psb_success_
    name = 'psb_precfree'
    call psb_erractionsave(err_act)

    me=-1

    if (allocated(prec%prec)) then 
      call prec%prec%free(info)
      if (info /= psb_success_) goto 9999
      deallocate(prec%prec,stat=info)
      if (info /= psb_success_) goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_s_prec_free

  function psb_sprec_sizeof(prec) result(val)
    class(psb_sprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer(psb_ipk_) :: i
    
    val = 0
    if (allocated(prec%prec)) then 
      val = val + prec%prec%sizeof()
    end if
    
  end function psb_sprec_sizeof

  subroutine psb_s_prec_clone(prec,precout,info)
    implicit none 
    class(psb_sprec_type), intent(inout) :: prec
    class(psb_sprec_type), intent(inout) :: precout
    integer(psb_ipk_), intent(out)             :: info

    info = psb_success_
    call prec%free(info)
    if (allocated(prec%prec)) then 
      call prec%prec%clone(precout%prec,info)
    end if
    
  end subroutine psb_s_prec_clone

end module psb_s_prec_type
