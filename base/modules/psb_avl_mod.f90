!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
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

module psb_item_mod
  type psb_item_int2
    integer :: key, val 
  end type psb_item_int2
  interface psb_sizeof
    module procedure psb_item_int2_size
  end interface
contains
  function  psb_item_int2_size(node)
    use psb_const_mod 
    type(psb_item_int2) :: node
    psb_item_int2_size  = psb_sizeof_int * 2
  end function psb_item_int2_size

  subroutine psb_print_item_key_val(iout,item)
    integer, intent(in)             :: iout
    type(psb_item_int2), intent(in) :: item
    write(iout,*) 'Value: ',item%key,item%val
    call flush(iout)
  end subroutine psb_print_item_key_val
end module psb_item_mod

module psb_avl_mod

  use psb_item_mod
  
  integer, parameter  :: LeftHigh = -1, EqualHeight=0, RightHigh=1
  integer, parameter  :: AVLTreeDuplicate = -123, AVLTreeOK=0, &
       & AVLTreeOutOfMemory=-512, AVLTreeFatalError=-1024
  integer :: level,outlev
  integer, parameter  :: poolsize = 1024

  type psb_treenode_int2
    type(psb_item_int2)              :: item
    type(psb_treenode_int2), pointer :: left=>null(), right=>null()
    integer                          :: balance
  end type psb_treenode_int2

  type psb_treevect_int2
    type(psb_treenode_int2) :: pool(poolsize)
    integer                 :: avail
    type(psb_treevect_int2), pointer :: next=>null(), prev=>null()
  end type psb_treevect_int2
  
  type psb_tree_int2
    type(psb_treevect_int2), pointer :: head=>null(), current=>null()
    type(psb_treenode_int2), pointer :: root=>null()
    integer :: nnodes
  end type psb_tree_int2

  interface psb_sizeof
    module procedure psb_Sizeof_Tree_int2, psb_sizeof_node_int2
  end interface

  interface InitSearchTree
    module procedure InitSearchTree_int2
  end interface

  interface FreeSearchTree
    module procedure FreeSearchTree_int2
  end interface

  interface SearchKey
    module procedure SearchKey_int2
  end interface

  interface SearchInsKey
    module procedure SearchInsKey_int2
  end interface

  interface GetAVLTree
    module procedure GetAVLTree_int2
  end interface

  interface CloneSearchTree
    module procedure CloneSearchTree_int2
  end interface

  interface CloneAVLTree
    module procedure CloneAVLTree_int2
  end interface

  interface GetAVLNode
    module procedure GetAVLNode_int2
  end interface
  interface UnGetAVLNode
    module procedure UnGetAVLNode_int2
  end interface

  interface VisitAVLTree
    module procedure VisitAVLTree_int2, VisitAVLTreeNode_int2
  end interface

  interface VisitAVLTreeLev
    module procedure VisitAVLTreeLev_int2, VisitAVLTreeNodeLev_int2
  end interface
  
  interface AVLTreeLeftBalance
    module procedure AVLTreeLeftBalance_int2
  end interface

  interface AVLTreeRightBalance
    module procedure AVLTreeRightBalance_int2
  end interface

  interface AVLTreeRotateLeft
    module procedure AVLTreeRotateLeft_int2
  end interface

  interface AVLTreeRotateRight
    module procedure AVLTreeRotateRight_int2
  end interface

  interface AVLSearchKey
    module procedure AVLSearchKey_int2
  end interface

  interface AVLSearchInsKey
    module procedure AVLSearchInsKey_int2
  end interface

  interface AVLSearchInsNode
    module procedure AVLSearchInsNode_int2
  end interface

contains

  subroutine InitSearchTree_int2(tree, info)
    type(psb_tree_int2), pointer  :: tree
    integer                       :: info
    
    if (associated(tree)) then 
      call FreeSearchTree(tree,info)
    end if
    call GetAVLTree(tree,info)
    
  end subroutine InitSearchTree_int2

  subroutine CloneSearchTree_int2(treein, treeout)
    type(psb_tree_int2), pointer  :: treein,treeout
    integer             :: info
    if (.not.associated(treein)) then 
      treeout => null()
      return
    endif
    call GetAVLTree(treeout,info)
    call CloneAVLTree(treein%root,treeout)

  end subroutine CloneSearchTree_int2

  recursive subroutine CloneAVLTree_int2(root, tree)
    type(psb_treenode_int2), pointer :: root
    type(psb_tree_int2), pointer     :: tree
    integer                          :: info, key,val,next
    if (.not.associated(root)) return
    key  = root%item%key
    next = root%item%val
    call SearchInsKey(tree,key,val,next,info)
    call CloneAVLTree(root%left,tree)
    call CloneAVLTree(root%right,tree)
  end subroutine CloneAVLTree_int2

  subroutine FreeSearchTree_int2(tree, info)
    type(psb_tree_int2), pointer  :: tree
    integer             :: info
    type(psb_treevect_int2), pointer :: current,next

    if (.not.associated(tree)) return
    current => tree%head
    do 
      if (.not.associated(current)) exit
      next => current%next
      deallocate(current,stat=info) 
      if (info /= 0) then 
        info = AVLTreeFatalError
        return
      end if
      current => next
    end do
    deallocate(tree,stat=info) 
    if (info /= 0) then 
      info = AVLTreeFatalError
      return
    end if
    
  end subroutine FreeSearchTree_int2
  
  function psb_Sizeof_Tree_int2(tree)
    use psb_const_mod 
    type(psb_tree_int2), pointer  :: tree
    integer                       :: psb_Sizeof_Tree_int2
    integer             :: val
    type(psb_treevect_int2), pointer :: current,next

    val = 0 
    if (associated(tree)) then 
      current => tree%head
      do 
        if (.not.associated(current)) exit
        val  = val + 3*psb_sizeof_int + poolsize*psb_sizeof(current%pool(1))
        current => current%next
      end do
    end if
    psb_Sizeof_Tree_int2 = val 
  end function psb_Sizeof_Tree_int2

  function psb_sizeof_node_int2(node)

    use psb_const_mod 
    type(psb_treenode_int2) :: node
    integer                 :: psb_sizeof_node_int2
    integer                 :: val

    
    psb_sizeof_node_int2 = 3*psb_sizeof_int + psb_sizeof(node%item)

  end function psb_sizeof_node_int2

  subroutine SearchKey_int2(tree,key,val,info)
    type(psb_tree_int2), target :: tree
    integer                     :: key,val,info
    type(psb_item_int2), pointer :: retval
    info = 0
    call AVLSearchKey(tree,key,retval,info)
    if (associated(retval)) then 
      val = retval%val
    else
      val = -1
    end if
  end subroutine SearchKey_int2

  subroutine SearchInsKey_int2(tree,key,val, nextval,info)
    type(psb_tree_int2), target :: tree
    integer                     :: key,val,nextval,info

    call AVLSearchInsKey(tree,key,val,nextval,info)
    
  end subroutine SearchInsKey_int2

  subroutine GetAVLTree_int2(tree, info)
    type(psb_tree_int2), pointer  :: tree
    integer                       :: info

    allocate(tree, stat=info) 
    if (info == 0) allocate(tree%head,stat=info)
    if (info == 0) then 
      tree%current => tree%head
      tree%head%avail = 0
      tree%nnodes=0
    end if
    
    if (info /= 0) then 
      write(0,*) 'Failed allocation 1 GetAVLTree '
      info = AVLTreeOutOfMemory

      return
    end if
    
  end subroutine GetAVLTree_int2

  subroutine VisitAVLTree_int2(tree, info,iout)
    type(psb_tree_int2), pointer  :: tree
    integer                       :: info
    integer, optional             :: iout

    info = 0
    if (.not.associated(tree)) return
    call VisitAVLTree(tree%root,iout)

  end subroutine VisitAVLTree_int2

  recursive subroutine VisitAVLTreeNode_int2(root,iout)
    type(psb_treenode_int2), pointer  :: root
    integer, optional             :: iout
    integer                       :: info

    if (.not.associated(root)) return
    call VisitAVLTree(root%left,iout)
    if (present(iout)) then 
      call psb_print_item_key_val(iout,root%item)
    else
      call psb_print_item_key_val(6,root%item)
    end if
    call VisitAVLTree(root%right,iout)
  end subroutine VisitAVLTreeNode_int2  

  subroutine VisitAVLTreeLev_int2(tree, info)
    type(psb_tree_int2), pointer  :: tree
    integer                       :: info

    if (.not.associated(tree)) return
    do outlev = 0, 3
      write(6,*) 'Tree level : ',outlev
      call VisitAVLTreeLev(tree%root,0)
    end do

  end subroutine VisitAVLTreeLev_int2

  recursive subroutine VisitAVLTreeNodeLev_int2(root,level)
    type(psb_treenode_int2), pointer  :: root
    integer                       :: info,level

    if (.not.associated(root)) return
    call VisitAVLTreeLev(root%left,level+1)
    if (level == outlev) call psb_print_item_key_val(6,root%item)
    call VisitAVLTreeLev(root%right,level+1)
  end subroutine VisitAVLTreeNodeLev_int2
  

  function GetAVLNode_int2(tree, info)
    type(psb_tree_int2), target      :: tree
    type(psb_treenode_int2), pointer :: GetAVLNode_int2
    integer                          :: info
    type(psb_treevect_int2), pointer :: current, temp

    GetAVLNode_int2 => null() 
    
    if (.not.associated(tree%current)) then 
      allocate(tree%head,stat=info)
      if (info /= 0) then
        info = AVLTreeOutOfMemory
        return
      end if
      tree%current => tree%head
      tree%current%avail = 0
    end if
    current => tree%current
    do 
      if (current%avail < poolsize) exit
      if (.not.(associated(current%next))) then 
        allocate(temp,stat=info)
        if (info /= 0) then
          info = AVLTreeOutOfMemory
          return
        end if
        temp%avail = 0
        temp%prev => current
        current%next => temp
      end if
      current => current%next
    end do
    tree%current  => current
    current%avail = current%avail + 1
    GetAVLNode_int2 => current%pool(current%avail)
    
  end function GetAVLNode_int2

  subroutine  UnGetAVLNode_int2(tree, info)
    type(psb_tree_int2), target      :: tree
    integer                          :: info
    
    
    if (.not.associated(tree%current)) then 
      return
    end if
    if (tree%current%avail > 0) &
         & tree%current%avail = tree%current%avail - 1
    return
  end subroutine UnGetAVLNode_int2
  
  subroutine AVLSearchKey_int2(tree,key,retval,info)
    type(psb_tree_int2), target      :: tree
    integer                          :: key,info
    type(psb_item_int2), pointer     :: retval
    type(psb_treenode_int2), pointer :: root
    
    retval => null()
    root   => tree%root
    do 
      if (.not.associated(root)) exit
      if (key < root%item%key) then 
        root   => root%left
      else if (key == root%item%key) then 
        retval => root%item
        exit
      else if (key > root%item%key) then 
        root   => root%right
      end if
    end do

  end subroutine AVLSearchKey_int2

  subroutine AVLSearchInsKey_int2(tree,key,val,nextval,info)
    type(psb_tree_int2),  target :: tree
    integer                      :: key,val,nextval,info
    type(psb_treenode_int2), pointer :: itemp 
    logical                      :: taller

    itemp => GetAVLNode(tree,info)  
    if (info /=0) then 
      return
    end if
    if (.not.associated(itemp)) then 
      info = -5
      return
    endif
    itemp%item%key = key
    itemp%item%val = nextval
    itemp%left => null()
    itemp%right => null()

    call AVLSearchInsNode(tree%root,itemp,taller,info)
    val = itemp%item%val
    if (info == AVLTreeDuplicate) then 
      call UnGetAVLNode(tree,info)
!!$      write(0,*) 'From searchInsNode ',key,val,nextval
      info = 0
      return
    else if (info == AVLTreeOK) then 
      tree%nnodes = tree%nnodes + 1 
      info = 0
      return
    else 
      write(0,*) 'Error from inner SearchInsNode '
    endif    
    
  end subroutine AVLSearchInsKey_int2


  recursive subroutine AVLSearchInsNode_int2(root,node,taller,info)
    type(psb_treenode_int2), pointer :: root, node    
    integer                          :: info
    logical                          :: taller    
    
    info = AVLTreeOK
    taller = .false.
    if (.not.associated(root)) then 
      root         => node
      node%balance = EqualHeight
      node%left    => null()
      node%right   => null()
      taller       = .true.
    else if (node%item%key == root%item%key) then 
!!$      write(0,*) 'SearchInsNode : found key',node%item%key,node%item%val,&
!!$           &root%item%key,root%item%val
      info          = AVLTreeDuplicate
      node%item%val = root%item%val
      return

    else if (node%item%key < root%item%key) then 

      call AVLSearchInsNode(root%left,node,taller,info)
      if (info == AVLTreeDuplicate) return
      if (info == AVLTreeFatalError) return
      if (taller) then 
        select case(root%balance)
        case(LeftHigh)
          call AVLTreeLeftBalance(root,taller)
        case(EqualHeight)
          root%balance = LeftHigh
        case(RightHigh)
          root%balance = EqualHeight
          taller       = .false.
        case default
          info = AVLTreeFatalError
        end select
      end if
    else if (node%item%key > root%item%key) then 
      call AVLSearchInsNode(root%right,node,taller,info)
      if (info == AVLTreeDuplicate) return
      if (info == AVLTreeFatalError) return
      if (taller) then 
        select case(root%balance)
        case(LeftHigh)
          root%balance = EqualHeight
          taller       = .false.
        case(EqualHeight)
          root%balance = RightHigh
        case(RightHigh)
          call AVLTreeRightBalance(root,taller)
        case default
          info = AVLTreeFatalError
        end select
      end if
    end if

  end subroutine AVLSearchInsNode_int2



  recursive subroutine AVLTreeLeftBalance_int2(root,taller)
    type(psb_treenode_int2), pointer :: root
    logical                          :: taller    

    type(psb_treenode_int2), pointer :: rs, ls
    
    ls => root%left
    select case (ls%balance)
    case(LeftHigh)
      root%balance = EqualHeight
      ls%balance   = EqualHeight
      call AVLTreeRotateRight(root)
      taller       = .false.
    case(EqualHeight)
      write(0,*) 'Warning: balancing and already balanced left tree? '
    case(RightHigh) 
      rs => ls%right
      select case(rs%balance) 
      case(LeftHigh)
        root%balance = RightHigh
        ls%balance   = EqualHeight
      case(EqualHeight)
        root%balance = EqualHeight
        ls%balance   = EqualHeight
      case(RightHigh)
        root%balance = EqualHeight
        ls%balance   = LeftHigh
      end select
      rs%balance = EqualHeight
      call AVLTreeRotateLeft(root%left)
      call AVLTreeRotateRight(root)
      taller = .false.
    end select

  end subroutine AVLTreeLeftBalance_int2


  recursive subroutine AVLTreeRightBalance_int2(root,taller)
    type(psb_treenode_int2), pointer :: root
    logical                          :: taller    
    type(psb_treenode_int2), pointer :: rs, ls
    
    rs => root%right
    select case (rs%balance)
    case(RightHigh)
      root%balance = EqualHeight
      rs%balance   = EqualHeight
      call AVLTreeRotateLeft(root)
      taller       = .false.
    case(EqualHeight)
      write(0,*) 'Warning: balancing and already balanced right tree? '
    case(LeftHigh) 
      ls => rs%left
      select case(ls%balance) 
      case(RightHigh)
        root%balance = LeftHigh
        rs%balance   = EqualHeight
      case(EqualHeight)
        root%balance = EqualHeight
        rs%balance   = EqualHeight
      case(LeftHigh)
        root%balance = EqualHeight
        rs%balance   = RightHigh
      end select
      ls%balance = EqualHeight
      call AVLTreeRotateRight(root%right)
      call AVLTreeRotateLeft(root)
      taller = .false.
    end select
  end subroutine AVLTreeRightBalance_int2
  


  subroutine AVLTreeRotateLeft_int2(root)
    type(psb_treenode_int2), pointer :: root
    type(psb_treenode_int2), pointer :: temp
    if (.not.associated(root)) then 
      return
    endif
    if (.not.associated(root%right)) then 
      return
    endif
    temp       => root%right
    root%right => temp%left
    temp%left  => root
    root       => temp
    
  end subroutine AVLTreeRotateLeft_int2

  subroutine AVLTreeRotateRight_int2(root)
    type(psb_treenode_int2), pointer :: root
    type(psb_treenode_int2), pointer :: temp
    if (.not.associated(root)) then 
      return
    endif
    if (.not.associated(root%left)) then 
      return
    endif
    temp       => root%left
    root%left  => temp%right
    temp%right => root
    root       => temp
    
  end subroutine AVLTreeRotateRight_int2


end module psb_avl_mod
