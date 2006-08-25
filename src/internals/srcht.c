/*
 *             Parallel Sparse BLAS  v2.0
 *   (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
 *                      Alfredo Buttari        University of Rome Tor Vergata
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the PSBLAS group or the names of its contributors may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */
/*****************************************************************/
/*                                                               */
/*   srcht.c : specialized insert/search for (key,val) integer   */
/*             pairs. Written by: Salvatore Filippone            */
/*                                                               */
/*   Last updated: Mar  09 2004                                  */
/*                                                               */
/*   Uses: avltree                                               */
/*                                                               */
/*  Data types:                                                  */
/*                                                               */
/*   KeyType: struct with two integer fields, key and val.       */
/*                                                               */
/*                                                               */
/*  User callable functions:                                     */
/*                                                               */
/*   void InitPairSearchTree(int *iret)                          */
/*       Purpose: initialize a search structure;                 */
/*       Function value: 0: OK                                   */
/*                       -1: failure                             */
/*                                                               */
/*                                                               */
/*   void SearchInsKeyVal(int *key, int *val, int *res,          */
/*                        int *iret)                             */
/*       Purpose: Search for a key, insert it if not present.    */
/*                                                               */
/*       Input:  1. key                                          */
/*                  Key to be searched for.                      */
/*               2. val                                          */
/*                  Value to form a (key,val) pair to be         */
/*                  inserted if key not already present.         */
/*       Output: 3. res                                          */
/*                  The val part of the pair with key; if the    */
/*                  key was freshly inserted then res=val        */
/*       Function value:  0 Normal termination                   */
/*                        1 Key was already present.             */
/*                       -1 Invalid input pointer                */
/*                       -3 Memory allocation failure            */
/*                                                               */
/*                                                               */
/*  void FreePairSearchTree()                                    */
/*       Purpose: free up tree data storage                      */
/*                                                               */
/*                                                               */
/*****************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "avltree.h"

#define  POOLSIZE  4096
#define  CACHESIZE 16
#ifdef Add_
#define InitPairSearchTree initpairsearchtree_
#define FreePairSearchTree freepairsearchtree_
#define SearchInsKeyVal    searchinskeyval_
#define SearchKeyVal       searchkeyval_
#define NPairs             npairs_
#endif
#ifdef AddDouble_
#define InitPairSearchTree initpairsearchtree_
#define FreePairSearchTree freepairsearchtree_
#define SearchInsKeyVal    searchinskeyval_
#define SearchKeyVal       searchkeyval_
#define NPairs             npairs_
#endif
#ifdef NoChange
#define InitPairSearchTree initpairsearchtree
#define FreePairSearchTree freepairsearchtree
#define SearchInsKeyVal    searchinskeyval
#define SearchKeyVal       searchkeyval
#define NPairs             npairs
#endif



typedef struct keypair *KeyPairPtr;
typedef struct keypair {
  int key,val;
} KeyPair;


typedef struct pairvect *PairVectPtr;
typedef struct pairvect {
  KeyPair pool[POOLSIZE];
  int avail;
  PairVectPtr previous, next;
} PairVect;



typedef struct pairtree *PairTreePtr;
typedef struct pairtree {
  int retval;
  int cache[2][CACHESIZE], cpnt;
  PairVectPtr PairPoolRoot,PairPoolCrt;
  AVLTreePtr  tree;
} PairTree;

#ifdef LargeFptr
typedef long long fptr;  /* 32-bit by default */
#else
typedef int fptr;  /* 32-bit by default */
#endif

/* 		 fptr *f_factors,  */
/*     *f_factors = (fptr) LUfactors; */

/* static int retval; */
/* static PairVectPtr PairPoolRoot=NULL,PairPoolCrt=NULL; */
/* static AVLTreePtr  tree=NULL; */

int CompareKeys(void *key1, void *key2)
{
  if (((KeyPairPtr) key1)->key < ((KeyPairPtr) key2)->key){
    return(-1);
  } else if (((KeyPairPtr)key1)->key == ((KeyPairPtr)key2)->key){
    return(0);
  } else {
    return(1);
  }
}

void InitPairSearchTree(fptr *ftree, int *iret)
{
  int i;
  *iret = 0;
  PairTreePtr PTree;
  
  if ((PTree  = malloc(sizeof(PairTree)))==NULL) {
    *iret=-1; return;
  }
  PTree->retval=0;
  for (i=0; i<CACHESIZE; i++) {
    PTree->cache[0][i]=PTree->cache[1][i] = -1;
  }
  PTree->cpnt=0;
  if ((PTree->tree  = GetAVLTree())==NULL) {
    *iret=-1; return;
  }
  if ((PTree->PairPoolRoot=(PairVectPtr)malloc(sizeof(PairVect)))==NULL) {
    *iret=-3; return;
  } else {
    PTree->PairPoolRoot->avail=0;
    PTree->PairPoolRoot->previous=PTree->PairPoolRoot->next=NULL;
    PTree->PairPoolCrt=PTree->PairPoolRoot;
  }
  *ftree = (fptr) PTree;
  return;
}


int  NPairs(fptr *ftree)   
{
  PairTreePtr PTree;
  
  PTree = (PairTreePtr) *ftree;

  return(HowManyItems(PTree->tree));
}

void KeyUpdate( void *key1, void *key2, void *data)
{
  *((int *) data)=((KeyPairPtr) key2)->val;
}


void  FreePairSearchTree(fptr *ftree)   
{
  PairTreePtr PTree;
  PairVectPtr current,next;
  
  PTree = (PairTreePtr) *ftree;

  AVLTreeFree(PTree->tree,NULL);
  
  current=PTree->PairPoolRoot;
  
  while (current != NULL) {
    next=current->next;
    free(current);
    current=next;
  }
  free(PTree->tree);
  free(PTree);
  *ftree = (fptr) NULL;
  return;
}

int  AdvanceKeyPair(PairVectPtr current) 
{
  if (current!=NULL) {
    current->avail +=1;
    return(current->avail);
  }
  return(-1);
}


KeyPairPtr GetKeyPair(PairVectPtr *current)
{
  PairVectPtr new, crt;
  KeyPairPtr newnode;

  crt=*current;
  if (crt==NULL) {
    return(NULL);
  }

  if (crt->avail<POOLSIZE) {
    newnode=&(crt->pool[crt->avail]);
  } else {    
    if ((new=(PairVectPtr)malloc(sizeof(PairVect)))==NULL) {
      fprintf(stderr,"Memory allocation failure\n");
      return(NULL);
    }
    memset(new,'\0',sizeof(PairVect));   
    newnode=&(new->pool[0]);
    crt->next=new;
    new->previous=crt;
    new->next=NULL;
    *current=new;
  }
  return(newnode);
}


/*                                                               */
/*   void SearchInsKeyVal(int *key, int *val, int *res,          */
/*                        int *iret)                             */
/*       Purpose: Search for a key, insert it if not present.    */
/*                                                               */
/*       Input:  1. key                                          */
/*                  Key to be searched for.                      */
/*               2. val                                          */
/*                  Value to form a (key,val) pair to be         */
/*                  inserted if key not already present.         */
/*       Output: 3. res                                          */
/*                  The val part of the pair with key; if the    */
/*                  key was freshly inserted then res=val        */
/*       Function value:  0 Normal termination                   */
/*                       -1 Invalid input pointer                */
/*                       -3 Memory allocation failure            */
/*                                                               */

void SearchInsKeyVal(fptr *ftree, int *key, int *val, int *res, int *iret)
{
  PairTreePtr PTree;
  KeyPairPtr node; 
  int info,i;

  PTree = (PairTreePtr) *ftree;
  node      = GetKeyPair(&(PTree->PairPoolCrt));
  node->key = *key;
  node->val = *val;
  if ((i=PTree->cpnt)<CACHESIZE) {
    PTree->cache[0][i] = *key;
    PTree->cache[1][i] = *val;
    PTree->cpnt=i+1;
  }
  
  info = AVLTreeInsert(PTree->tree,node,CompareKeys,KeyUpdate,&(PTree->retval));
  *iret = info;
  if (info==0) {
    *res = node->val;
    AdvanceKeyPair(PTree->PairPoolCrt);
  } else if (info == 1) {
    *res  = PTree->retval;
  }
  return;
}


void SearchKeyVal(fptr *ftree, int *key, int *res, int *iret)
{
  PairTreePtr PTree;
  KeyPair node; 
  AVLNodePtr noderes;
  KeyPairPtr result;
  int i,sv[2];
  int info;
  
  *iret = 0;
  
  PTree = (PairTreePtr) *ftree;
#if 0
  for (i=0; i<CACHESIZE; i++) {
    if (PTree->cache[0][i] == *key) { 
      *res=PTree->cache[1][i];
      sv[0]=PTree->cache[0][i];
      sv[1]=PTree->cache[1][i];
      PTree->cache[0][i]=PTree->cache[0][0];
      PTree->cache[1][i]=PTree->cache[1][0];
      PTree->cache[0][0]=sv[0];
      PTree->cache[1][0]=sv[1];
      return;
    }
  }
     
#endif
  
  node.key=*key;
  if ((noderes  = AVLTreeSearch(PTree->tree,&node,CompareKeys))==NULL) {
    *res = -1;
    *iret = -1;
  } else {
    result = (KeyPairPtr) noderes->key;
    *res = result->val;
#if 0
    for (i=CACHESIZE-1; i>0; i--) {
      PTree->cache[0][i]=PTree->cache[0][i-1];
      PTree->cache[0][i]=PTree->cache[1][i-1];
    }
    PTree->cache[0][0]=*key;
    PTree->cache[1][0]=result->val;
#endif
  }
#ifdef PROFILE
  *iret = 1;
  *iret = PTree->tree->nsteps;
#endif
  return;
}
