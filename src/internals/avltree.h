/*   Type definitions for balanced AVL tree search and insertion */
/*   See avltree.c for a full definition of the subroutines      */
/*                                                               */

typedef struct avlnode *AVLNodePtr;
typedef struct avlnode {
  AVLNodePtr llink,rlink;
  int bal;
  void *key;
} AVLNode; 

typedef struct avltvect *AVLTVectPtr;

typedef struct avltree *AVLTreePtr;
typedef struct avltree {
  AVLTVectPtr first, current;
  AVLNodePtr root;
  int nnodes;
} AVLTree;


AVLNodePtr AVLTreeSearch(AVLTreePtr, void *, int (*)(void *, void *));
AVLNodePtr GetAVLNode(AVLTreePtr);
int  AVLTreeInit(AVLTreePtr);
int  AVLTreeReInit(AVLTreePtr);
AVLTreePtr GetAVLTree();
int  AVLTreeInsert(AVLTreePtr, void *, int (*)(void *, void *),
                   void (*)(void *, void *));
AVLNodePtr AVLTreeUserInsert(AVLTreePtr, void *, int (*)(void *, void *));
void AVLTreeInorderTraverse(AVLTreePtr, void (*)(void *, void *), void *);
void AVLTreePreorderTraverse(AVLTreePtr, void (*)(void *, void *), void *);
void AVLTreeFree(AVLTreePtr, void (*)(void *));
int HowManyItems(AVLTreePtr);
int AVLTreeInorderTraverseWithDelims(AVLTreePtr,void*, void*, int (*)(void*,void*),
                                     void (*)(void *, void *), void *);
     


