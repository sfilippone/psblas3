/*****************************************************************/
/*                                                               */
/*   avltree.c: balanced AVL tree search and insertion           */
/*              written by: Salvatore Filippone                  */
/*                                                               */
/*   Last updated: Mar  09 2004                                  */
/*                                                               */
/*   Referrences: [1] D. E. Knuth                                */
/*                    The Art of Computer Programming            */
/*                    Vol. 3: Sorting and Searching, sec. 6.2.3  */
/*                    Addison-Wesley                             */
/*                                                               */
/*  General description:                                         */
/*                                                               */
/*  Build and maintain a balanced binary search tree with        */
/*  arbitrary keys. The user is responsible for providing        */
/*  compare functions operating on the keys themselves.          */
/*  Key pointers are stored into nodes that are managed          */
/*  by the subroutine calls; the user should never examine       */
/*  nodes directly.                                              */
/*  The nodes for user items are allocated in batches,           */
/*  and the batches are kept as a doubly linked list.            */
/*                                                               */
/*  Data types:                                                  */
/*   AVLTree: structure containing pointers to the list          */
/*            of node batches and to the root of the binary tree */
/*            structure                                          */
/*                                                               */
/*   AVLNode: binary tree node, containing link pointers         */
/*            a reserved field, and a pointer to user data       */
/*                                                               */
/*                                                               */
/*  User callable functions:                                     */
/*                                                               */
/*   AVLTreePtr GetAVLTree()                                     */
/*       Purpose: allocate a new tree;                           */
/*       Function value: a fresh AVL tree pointer;               */
/*                       returns NULL in case of a memory failure*/
/*                                                               */
/*                                                               */
/*   int AVLTreeReInit(AVLTreePtr Tree)                          */
/*       Purpose: reinitialize an existing AVL Tree, reusing     */
/*                node batches already allocated.                */
/*       Input:  1. Tree                                         */
/*                   A pointer to an existing tree structure     */
/*       Function value:  0 Normal termination                   */
/*                       -1 Invalid input pointer                */
/*                       -3 Memory allocation failure            */
/*                                                               */
/*   AVLNodePtr AVLTreeSearch(AVLTreePtr Tree, void *key,        */
/*                           int (*comp)(void*, void*))          */
/*       Purpose: search an existing AVL Tree for a key          */
/*       Input:  1. Tree                                         */
/*                  A valid pointer to a Tree                    */
/*               2. key                                          */
/*                  The item being searched for                  */
/*               3. comp                                         */
/*                  A comparison function:                       */
/*                    a<b  =>   comp(a,b)<0                      */
/*                    a==b =>   comp(a,b)=0                      */
/*                    a>b  =>   comp(a,b)>0                      */
/*                 The function is always invoked as:            */
/*                    comp(user_key,tree_key);                   */
/*                                                               */
/*                                                               */
/*       Function value:  NULL: input error or item not found    */
/*                      valid pointer: pointer to a node         */
/*                              containing the key               */
/*                                                               */
/*   int AVLTreeInsert(AVLTreePtr Tree, void *key,               */
/*                     int (*comp)(void*,void*),                 */
/*                     void (*update)(void*,void*))              */
/*                                                               */
/*       Purpose: Insert an item into an existing (possibly      */
/*                empty) tree.                                   */
/*                                                               */
/*       Input:  1. Tree                                         */
/*                  The (existing) tree                          */
/*               2. key                                          */
/*                  The (new) item to be inserted                */
/*               3. comp                                         */
/*                  comparison function (as in AVLTreeSearch)    */
/*               4. update                                       */
/*                  A user provided function to be called when   */
/*                  the key is already present in the tree       */
/*                  with the calling sequence:                   */
/*                   update(new_key,existing_key)                */
/*                  enables the user to specify an arbitrary     */
/*                  update procedure.                            */
/*                                                               */
/*                                                               */
/*                                                               */
/*   AVLNodePtr AVLTreeUserInsert(AVLTreePtr Tree, void *key,    */
/*                     int (*comp)(void*,void*))                 */
/*                                                               */
/*       Purpose: Insert an item into an existing (possibly      */
/*                empty) tree; returns a pointer to a node       */
/*                containing the item, even when that node       */
/*                was already existing; does no update           */
/*                                                               */
/*       Input:  1. Tree                                         */
/*                  The (existing) tree                          */
/*               2. key                                          */
/*                  The (new) item to be inserted                */
/*               3. comp                                         */
/*                  comparison function (as in AVLTreeSearch)    */
/*                                                               */
/*       Function value: Valid pointer: pointer to a node        */
/*                            containing the item (possibly      */
/*                            was already there)                 */
/*                       NULL  input error or memory failure     */
/*                                                               */
/*                                                               */
/*   int HowManyKeys(AVLTreePtr Tree)                            */
/*       Purpose: how many keys does Tree contain?               */
/*       Function value:  >=0                                    */
/*                                                               */
/*                                                               */
/*  void AVLTreeInorderTraverse(AVLTreePtr Tree,                 */
/*                void (*func)(     void*, void*), void *data)   */
/*                                                               */
/*       Purpose: visit the nodes of the binary tree in their    */
/*                natural order, performing an arbitrary         */
/*                task upon visit.                               */
/*       Input: 1. Tree                                          */
/*                 A tree pointer                                */
/*              2. func                                          */
/*                 A function performing a user specified        */
/*                 task on each node; the fuction is invoked as  */
/*                  func(    key,data)                           */
/*                 where data is parm. 3                         */
/*              3. data                                          */
/*                 Auxiliary data to be passed to func upon      */
/*                 each visit                                    */
/*                                                               */
/*  int  AVLTreeInorderTraverseWithDelims(AVLTreePtr Tree,       */
/*             void *first, void *last, int (*comp)(void*,void*) */
/*              void (*func)(     void*, void*), void *data)     */
/*                                                               */
/*       Purpose: visit the nodes of the binary tree in their    */
/*                natural order, performing an arbitrary         */
/*                task upon visit, but only on nodes             */
/*                with their key within a specified range.       */
/*                                                               */
/*       Input: 1. Tree                                          */
/*                 A tree pointer                                */
/*              2. first                                         */
/*                 Visit nodes with   first <= node->key         */
/*              3. last                                          */
/*                 Visit nodes with   node->key <= last          */
/*              4. comp                                          */
/*                 comparison function (as in AVLTreeSearch)     */
/*              5. func                                          */
/*                 A function performing a user specified        */
/*                 task on each node; the fuction is invoked as  */
/*                  func(    key,data)                           */
/*                 where data is parm. 3                         */
/*              6. data                                          */
/*                 Auxiliary data to be passed to func upon      */
/*                 each visit                                    */
/*    Function value: total number of nodes visited (>=0)        */
/*                                                               */
/*                                                               */
/*                                                               */
/*  void AVLTreeFree(AVLTreePtr Tree, void (*ffree)(void*))      */
/*       Purpose: free up tree data storage                      */
/*                Does NOT free the Tree pointer itself,         */
/*                rather all the structures that it points to    */
/*       Input: 1. Tree                                          */
/*                 A tree pointer                                */
/*              2. ffree                                         */
/*                 A user specified function invoked on each     */
/*                 key pointer contained in the tree to free     */
/*                 its memory (if necessary). Can be NULL.       */
/*                                                               */
/*                                                               */
/*****************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "avltree.h"

#define  POOLSIZE  4096
#define  MAXSTACK  64
#define  MAX(a,b)  ((a)>=(b) ? (a) : (b))

typedef struct avltvect {
  AVLNode pool[POOLSIZE];
  int avail;
  AVLTVectPtr previous, next;
} AVLTVect;


int HowManyItems(AVLTreePtr Tree)
{
  if (Tree==NULL) {
    return(0);
  } else {
    return(Tree->nnodes);
  }
}


AVLTreePtr GetAVLTree()
{
  AVLTreePtr tree;
  if ((tree=(AVLTreePtr) malloc(sizeof(AVLTree)))!=NULL){
    memset(tree,'\0',sizeof(AVLTree));     
    AVLTreeInit(tree);
  }
  return(tree);
}

int AVLTreeInit(AVLTreePtr Tree)
{
  AVLTVectPtr current;
  if (Tree==NULL) {
    fprintf(stderr,"Cannot initialize a NULL Tree pointer\n");
    return(-1);
  }
   
  if (Tree->first!=NULL) {
    fprintf(stderr,"Cannot initialize a nonempty Tree: call AVLTreeFree first\n");
    return(-2);
  }
  
  if ((current=(AVLTVectPtr)malloc(sizeof(AVLTVect)))==NULL) {
    fprintf(stderr,"Memory allocation failure\n");
    return(-3);
  }
  memset(current,'\0',sizeof(AVLTVect));
  Tree->first=Tree->current=current;
  Tree->nnodes=0;
  Tree->root=NULL;
  return(0);
}

int AVLTreeReInit(AVLTreePtr Tree)
{
  AVLTVectPtr current /* , next */ ;
  if (Tree==NULL) {
    fprintf(stderr,"Cannot ReInitialize a NULL Tree pointer\n");
    return(-1);
  }
   
  if (Tree->first!=NULL) {
    current=Tree->first;
    while (current!=NULL) {
      current->avail=0;     
      memset(current->pool,'\0',POOLSIZE*sizeof(AVLNode));  
      current=current->next;
    }
  } else {  
    if ((current=(AVLTVectPtr)malloc(sizeof(AVLTVect)))==NULL) {
      fprintf(stderr,"Memory allocation failure\n");
      return(-3);
    }
    current->avail=0;
    current->previous=current->next=NULL;
    Tree->first=current;
  }
  Tree->current=Tree->first;
  Tree->nnodes=0;
  Tree->root=NULL;
  return(0);
}




AVLNodePtr AVLTreeSearch(AVLTreePtr Tree, void *key,
                     int (*comp)(void *, void *)) 
{
  AVLNodePtr current; 
  int icmp; 
  if (Tree==NULL) return(NULL);
  current = Tree->root;
  
  while (current != NULL) {
    icmp = (*comp)(key,current->key);
    if (icmp<0) {
      current = current->llink;
    } else if (icmp==0){      
      return(current);
    } else if (icmp>0) {
      current = current->rlink;
    }
  }
  return(current);
}



void  AVLTreeInorderTraverse(AVLTreePtr Tree, void (*func)(void *, void *),
                             void *data)   
{
  int lev;
  AVLNodePtr root;

  AVLNodePtr stack[MAXSTACK+2];
  int choice[MAXSTACK+2];
  root=Tree->root;
  if (root == NULL) return;  
  
  lev=0;
  stack[lev]  = root;
  choice[lev] = -1;
  
  while (lev>=0) {
    if (stack[lev]==NULL) {
      lev--;
    } else {
      if (choice[lev]==-1) {
        stack[lev+1]  = stack[lev]->llink;
        choice[lev+1] = -1;
        choice[lev]  += 1;
        lev++;
      } else if (choice[lev]==0) {
        (*func)(stack[lev]->key,data);    
        stack[lev+1]  = stack[lev]->rlink;
        choice[lev+1] = -1;
        choice[lev]  += 1;
        lev++;
      } else {
        lev--;
      }
    }
  }
}


int AVLTreeInorderTraverseWithDelims(AVLTreePtr Tree, void *first, void *last,
                                     int (*comp)(void*, void*),
                                     void (*func)(void *, void *),
                                     void *data)   
{
  AVLNodePtr root, current;
  int lev, nvisit, icmp;
  AVLNodePtr stack[MAXSTACK+2];
  int choice[MAXSTACK+2];

  root=Tree->root;
  if (root == NULL) return(0);  
  
  nvisit=0;
  lev=0;
  current = root;
  while (current != NULL) {
    stack[lev]  = current;
    icmp = (*comp)(first,current->key);
    if (icmp<=0) {
      choice[lev]=0;
      current = current->llink;
    } else if (icmp>0) {
      current = current->rlink;
      choice[lev]=1;
    }
    lev++;
  } 
  lev--;
  while (lev>=0) {
    if (stack[lev]==NULL) {
      lev--;
    } else {
      if (choice[lev]==-1) {
        stack[lev+1]  = stack[lev]->llink;
        choice[lev+1] = -1;
        choice[lev]  += 1;
        lev++;
      } else if (choice[lev]==0) {
        if (((*comp)(last,stack[lev]->key))<0) {        
          lev--;
        } else {
          (*func)(stack[lev]->key,data);    
          nvisit++;
          stack[lev+1]  = stack[lev]->rlink;
          choice[lev+1] = -1;
          choice[lev]  += 1;
          lev++;
        }
      } else {
        lev--;
      }
    }
  }
  return(nvisit);
}



void  AVLTreePreorderTraverse(AVLTreePtr Tree, void (*func)(void *, void *),
                              void *data)   
{
  AVLNodePtr root;
  int lev;
  AVLNodePtr stack[MAXSTACK+2];
  int choice[MAXSTACK+2];

  root=Tree->root;
  if (root == NULL) return;    
  lev=0;
  stack[lev]  = root;
  choice[lev] = -1;
  
  while (lev>=0) {
    if (stack[lev]==NULL) {
      lev--;
    } else {
      if (choice[lev]==-1) {
        (*func)(stack[lev]->key,data);    
        stack[lev+1]  = stack[lev]->llink;
        choice[lev+1] = -1;
        choice[lev]  += 1;
        lev++;
      } else if (choice[lev]==0) {
        stack[lev+1]  = stack[lev]->rlink;
        choice[lev+1] = -1;
        choice[lev]  += 1;
        lev++;
      } else {
        lev--;
      }
    }
  }
}



void  AVLTreeFree(AVLTreePtr Tree, void (*ffree)(void *))   
{
  AVLTVectPtr current, next;
  int i;
  if (Tree == NULL) return;  
  
  current=Tree->first;
  
  while (current != NULL) {
    next=current->next;
    if (*ffree != NULL) {
      for (i=0; i<current->avail; i++) 
	(*ffree)((current->pool[i]).key);
    }
    free(current);
    current=next;
  }
  Tree->nnodes=0;
  Tree->first=Tree->current=NULL;
  return;
}


AVLNodePtr GetAVLNode(AVLTreePtr Tree)
{
  AVLTVectPtr current, new;
  AVLNodePtr newnode;

  if (Tree==NULL) {
    return(NULL);
  }   
  if ((current=Tree->current)==NULL) {
    return(NULL);
  }

  while  ((current->avail>=POOLSIZE)&&(current->next!=NULL)) 
    current=current->next;
  
  if (current->avail<POOLSIZE) {
    newnode=&(current->pool[current->avail]);
    current->avail += 1;
  } else {    
    if ((new=(AVLTVectPtr)malloc(sizeof(AVLTVect)))==NULL) {
      fprintf(stderr,"Memory allocation failure\n");
      return(NULL);
    }
    memset(new,'\0',sizeof(AVLTVect));   
    newnode=&(new->pool[0]);
    new->avail = 1;
    current->next=new;
    new->previous=current;
    new->next=NULL;
    Tree->current=new;
  }
  return(newnode);
}

int AVLTreeInsert(AVLTreePtr Tree, void *key,int (*comp)(void *, void *),
                  void (*update)(void *, void *)) 
{
  AVLNodePtr root, t, s, p, q, r; 
  int search, bal, icmp;

  if (Tree==NULL) {
    fprintf(stderr,"Fatal error: null tree pointer\n");
    return(-1);
  }

  if ((root = Tree->root) == NULL) {
    if ((t=GetAVLNode(Tree))==NULL) {
      return(-2);
    }
    t->key = key;
    t->rlink=t->llink=NULL;
    t->bal=0;
    Tree->root = t;
    Tree->nnodes=1;
    return(0);
  } 
  t = NULL;
  s = root;
  p = root; 
  search=1;
  while (search) {
    icmp = (*comp)(key,p->key);
    if (icmp<0) {
      if ((q=p->llink)==NULL) {
        if ((q=GetAVLNode(Tree))==NULL) {
          return(-2);
        }
        p->llink=q;
        search=0; 
      } else {
        if (q->bal != 0) {
          t=p;
          s=q;
        }
      }
    } else if (icmp == 0) {
      (*update)(key,p->key);
      return(1);
    } else { 
      if ((q=p->rlink)==NULL) {
        if ((q=GetAVLNode(Tree))==NULL) {
          return(-2);
        }
        p->rlink=q;
        search=0; 
      } else {
        if (q->bal != 0) {
          t=p;
          s=q;
        }
      }
    }
    p=q;
  }  
  q->key=key;
  q->llink=q->rlink=NULL;
  q->bal=0;
  Tree->nnodes += 1;
  
  if ((*comp)(key,s->key)<0) {
    r=p=s->llink;
  } else {
    r=p=s->rlink;
  }

  while (p!=q) {
    if ((*comp)(key,p->key)<0) {
      p->bal=-1;
      p = p->llink;
    } else {
      p->bal=1;
      p=p->rlink;
    }
  }
  
  if ((*comp)(key,s->key)<0) {
    bal=-1;
  } else {
    bal=1;
  }
  
  if (s->bal == 0) {
    s->bal=bal;
    return (0);
  } else if (s->bal == -bal) {
    s->bal=0;
    return (0);
  } else if (s->bal == bal) {
    
    if (r->bal == bal) { 
      /* single rotation */
      p=r;
      if (bal>0) {
        s->rlink=r->llink;
        r->llink=s;
      } else {
        s->llink=r->rlink;
        r->rlink=s;
      }
      s->bal=r->bal=0;
    } else if (r->bal == -bal) { 
      /* double rotation */
      if (bal>0) {
        p=r->llink;
        r->llink=p->rlink;
        p->rlink=r;
        s->rlink=p->llink;
        p->llink=s;     
      } else {
        p=r->rlink;
        r->rlink=p->llink;
        p->llink=r;
        s->llink=p->rlink;
        p->rlink=s;     
      }
      if (p->bal == bal) {
        s->bal=-bal;
        r->bal=0;
      } else if (p->bal==0) {
        s->bal=r->bal=0;
      } else {
        r->bal=bal;
        s->bal=0;
      } 
      p->bal=0;
    }
    if (t==NULL) {
      root=p;
    } else {
      if (t->rlink==s) {
        t->rlink=p;
      } else {
        t->llink=p;
      }
    }
    Tree->root=root;
    return(0);
  }
  return(0);
}

AVLNodePtr AVLTreeUserInsert(AVLTreePtr Tree, void *key,
                             int (*comp)(void *, void *))
{
  AVLNodePtr root, t, s, p, q, r; 
  int search, bal, icmp;

  if (Tree==NULL) {
    fprintf(stderr,"Fatal error: null tree pointer\n");
    return(NULL);
  }

  if ((root = Tree->root) == NULL) {
    if ((t=GetAVLNode(Tree))==NULL) {
      return(NULL);
    }
    t->key = key;
    t->rlink=t->llink=NULL;
    t->bal=0;
    Tree->root = t;
    Tree->nnodes=1;
    return(t);
  } 
  t = NULL;
  s = root;
  p = root; 
  search=1;
  while (search) {
    icmp = (*comp)(key,p->key);
    if (icmp<0) {
      if ((q=p->llink)==(AVLNodePtr) NULL) {
        if ((q=GetAVLNode(Tree))==NULL) {
          return(NULL);
        }
        p->llink=q;
        search=0; 
      } else {
        if (q->bal != 0) {
          t=p;
          s=q;
        }
      }
    } else if (icmp == 0) {
      return(p);
    } else { 
      if ((q=p->rlink)==NULL) {
        if ((q=GetAVLNode(Tree))==NULL) {
          return(NULL);
        }
        p->rlink=q;
        search=0; 
      } else {
        if (q->bal != 0) {
          t=p;
          s=q;
        }
      }
    }
    p=q;
  }  
  q->key=key;
  q->llink=q->rlink=NULL;
  q->bal=0;
  Tree->nnodes += 1;
  
  if ((*comp)(key,s->key)<0) {
    r=p=s->llink;
  } else {
    r=p=s->rlink;
  }

  while (p!=q) {
    if ((*comp)(key,p->key)<0) {
      p->bal=-1;
      p = p->llink;
    } else {
      p->bal=1;
      p=p->rlink;
    }
  }
  
  if ((*comp)(key,s->key)<0) {
    bal=-1;
  } else {
    bal=1;
  }
  
  if (s->bal == 0) {
    s->bal=bal;
    return (q);
  } else if (s->bal == -bal) {
    s->bal=0;
    return (q);
  } else if (s->bal == bal) {
    
    if (r->bal == bal) { 
      /* single rotation */
      p=r;
      if (bal>0) {
        s->rlink=r->llink;
        r->llink=s;
      } else {
        s->llink=r->rlink;
        r->rlink=s;
      }
      s->bal=r->bal=0;
    } else if (r->bal == -bal) { 
      /* double rotation */
      if (bal>0) {
        p=r->llink;
        r->llink=p->rlink;
        p->rlink=r;
        s->rlink=p->llink;
        p->llink=s;     
      } else {
        p=r->rlink;
        r->rlink=p->llink;
        p->llink=r;
        s->llink=p->rlink;
        p->rlink=s;     
      }
      if (p->bal == bal) {
        s->bal=-bal;
        r->bal=0;
      } else if (p->bal==0) {
        s->bal=r->bal=0;
      } else {
        r->bal=bal;
        s->bal=0;
      } 
      p->bal=0;
    }
    if (t==NULL) {
      root=p;
    } else {
      if (t->rlink==s) {
        t->rlink=p;
      } else {
        t->llink=p;
      }
    }
    Tree->root=root;
    return(q);
  }
  return(q);
}


