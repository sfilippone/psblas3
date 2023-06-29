/* SPDX-License-Identifier: GPL-2.0-or-later */
/*
  Red Black Trees
  (C) 1999  Andrea Arcangeli <andrea@suse.de>
  

  Userspace GNUC porting,del augmented deps and minimalizing for few OPs only:	Andrea Di Iorio
  linux/include/linux/rbtree.h

  To use rbtrees you'll have to implement your own insert and search cores.
  This will avoid us to use callbacks and to drop drammatically performances.
  I know it's not the cleaner way,  but in C (not in C++) to get
  performances and genericity...

  See Documentation/core-api/rbtree.rst for documentation and samples.
*/
/*
 *              RedBlackTree_linux_userspace
 *    (C) Copyright 2021-2022
 *        Andrea Di Iorio      
 * 
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *    1. Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *    2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions, and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *    3. The name of the RedBlackTree_linux_userspace or the names of its contributors may
 *       not be used to endorse or promote products derived from this
 *       software without specific written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE RedBlackTree_linux_userspace GROUP OR ITS CONTRIBUTORS
 *  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */ 

 /*
 * https://bitbucket.org/andysnake96/redblacktree_linux_userspace
 */ 


#ifndef	_LINUX_RBTREE_H
#define	_LINUX_RBTREE_H

//#include <linux/kernel.h>     //TODO LESS_DEPENDENCIES
//#include <linux/stddef.h>     //TODO LESS_DEPENDENCIES
//#include <linux/rcupdate.h>   //TODO LESS_DEPENDENCIES
#include <stdlib.h>
#include <stdio.h>

#include <macros.h>             ///my usefull macros :)
#include <macrosLinuxMock.h>             ///my usefull macros :)

struct rb_node {
	unsigned long  __rb_parent_color;
	struct rb_node *rb_right;
	struct rb_node *rb_left;
} __attribute__((aligned(sizeof(long))));
    /* The alignment might seem pointless, but allegedly CRIS needs it */

struct rb_root {
	struct rb_node *rb_node;
};
/////MINIMALIZE fulfilling rbtree_augmented deps TODO
///rb_node colors labels
#define	RB_RED		0
#define	RB_BLACK	1
#define __rb_parent(pc)    ((struct rb_node *)(pc & ~3))

#define __rb_color(pc)     ((pc) & 1)
#define __rb_is_black(pc)  __rb_color(pc)
#define __rb_is_red(pc)    (!__rb_color(pc))
#define rb_color(rb)       __rb_color((rb)->__rb_parent_color)
#define rb_is_red(rb)      __rb_is_red((rb)->__rb_parent_color)
#define rb_is_black(rb)    __rb_is_black((rb)->__rb_parent_color)


static inline void rb_set_parent_color(struct rb_node *rb,
				       struct rb_node *p, int color)
{
	rb->__rb_parent_color = (unsigned long)p | color;
}
///rb_set_parent	//TODO GENERIC ADDONLY

static inline void
__rb_change_child(struct rb_node *old, struct rb_node *new,
		  struct rb_node *parent, struct rb_root *root)
{
	if (parent) {
		if (parent->rb_left == old)
			WRITE_ONCE(parent->rb_left, new);
		else
			WRITE_ONCE(parent->rb_right, new);
	} else
		WRITE_ONCE(root->rb_node, new);
}

/////////////////////////////

#define rb_parent(r)   ((struct rb_node *)((r)->__rb_parent_color & ~3))

#define RB_ROOT	(struct rb_root) { NULL, }
#define	rb_entry(ptr, type, member) container_of(ptr, type, member)

#define RB_EMPTY_ROOT(root)  (READ_ONCE((root)->rb_node) == NULL)

/* 'empty' nodes are nodes that are known not to be inserted in an rbtree */
#define RB_EMPTY_NODE(node)  \
	((node)->__rb_parent_color == (unsigned long)(node))
#define RB_CLEAR_NODE(node)  \
	((node)->__rb_parent_color = (unsigned long)(node))


extern void rb_insert_color(struct rb_node *, struct rb_root *);
extern void rb_erase(struct rb_node *, struct rb_root *);


/* Find logical next and previous nodes in a tree */
extern struct rb_node *rb_next(const struct rb_node *);
extern struct rb_node *rb_prev(const struct rb_node *);
extern struct rb_node *rb_first(const struct rb_root *);
extern struct rb_node *rb_last(const struct rb_root *);

/* Postorder iteration - always visit the parent after its children */
extern struct rb_node *rb_first_postorder(const struct rb_root *);
extern struct rb_node *rb_next_postorder(const struct rb_node *);

/* Fast replacement of a single node without remove/rebalance/add/rebalance */
extern void rb_replace_node(struct rb_node *victim, struct rb_node *new,
			    struct rb_root *root);
extern void rb_replace_node_rcu(struct rb_node *victim, struct rb_node *new,
				struct rb_root *root);


static inline void rb_link_node(struct rb_node *node, struct rb_node *parent,
				struct rb_node **rb_link)
{
	node->__rb_parent_color = (unsigned long)parent;
	node->rb_left = node->rb_right = NULL;

	*rb_link = node;
}

/*static inline void rb_link_node_rcu(struct rb_node *node, struct rb_node *parent,
				    struct rb_node **rb_link)
{
	node->__rb_parent_color = (unsigned long)parent;
	node->rb_left = node->rb_right = NULL;

	rcu_assign_pointer(*rb_link, node);
}**/ //TODO LESS_DEPENDENCIES

#define rb_entry_safe(ptr, type, member) \
	({ typeof(ptr) ____ptr = (ptr); \
	   ____ptr ? rb_entry(____ptr, type, member) : NULL; \
	})

/**
 * rbtree_postorder_for_each_entry_safe - iterate in post-order over rb_root of
 * given type allowing the backing memory of @pos to be invalidated
 *
 * @pos:	the 'type *' to use as a loop cursor.
 * @n:		another 'type *' to use as temporary storage
 * @root:	'rb_root *' of the rbtree.
 * @field:	the name of the rb_node field within 'type'.
 *
 * rbtree_postorder_for_each_entry_safe() provides a similar guarantee as
 * list_for_each_entry_safe() and allows the iteration to continue independent
 * of changes to @pos by the body of the loop.
 *
 * Note, however, that it cannot handle other modifications that re-order the
 * rbtree it is iterating over. This includes calling rb_erase() on @pos, as
 * rb_erase() may rebalance the tree, causing us to miss some nodes.
 */
#define rbtree_postorder_for_each_entry_safe(pos, n, root, field) \
	for (pos = rb_entry_safe(rb_first_postorder(root), typeof(*pos), field); \
	     pos && ({ n = rb_entry_safe(rb_next_postorder(&pos->field), \
			typeof(*pos), field); 1; }); \
	     pos = n)

/*
 * Leftmost-cached rbtrees.
 *
 * We do not cache the rightmost node based on footprint
 * size vs number of potential users that could benefit
 * from O(1) rb_last(). Just not worth it, users that want
 * this feature can always implement the logic explicitly.
 * Furthermore, users that want to cache both pointers may
 * find it a bit asymmetric, but that's ok.
 */
struct rb_root_cached {
	struct rb_root rb_root;
	struct rb_node *rb_leftmost;
};

#define RB_ROOT_CACHED (struct rb_root_cached) { {NULL, }, NULL }

/* Same as rb_first(), but O(1) */
#define rb_first_cached(root) (root)->rb_leftmost

static inline void rb_insert_color_cached(struct rb_node *node,
					  struct rb_root_cached *root,
					  bool leftmost)
{
	if (leftmost)
		root->rb_leftmost = node;
	rb_insert_color(node, &root->rb_root);
}

static inline void rb_erase_cached(struct rb_node *node,
				   struct rb_root_cached *root)
{
	if (root->rb_leftmost == node)
		root->rb_leftmost = rb_next(node);
	rb_erase(node, &root->rb_root);
}

static inline void rb_replace_node_cached(struct rb_node *victim,
					  struct rb_node *new,
					  struct rb_root_cached *root)
{
	if (root->rb_leftmost == victim)
		root->rb_leftmost = new;
	rb_replace_node(victim, new, &root->rb_root);
}
///////////////////////////////////////////////////////////////////////////////
/***END OF ORIGINAL LINUX KERNEL HEADER MINIMALIZED
 ***FEW AUX FUNCS TO SUPPORT SpMM-Symbolic
 ***/
#include <string.h>
#include "config.h"
typedef struct{
	idx_t key;
	struct rb_node rb;

	/* following fields used for testing augmented rbtree functionality
	u32 val;
	u32 augmented;	///only for AUGMENTED_TEST
	*/
} rbNode;
//static struct rb_root_cached root = RB_ROOT_CACHED;

typedef struct rb_root_cached rbRoot;


/*
 * return 1 if @node with the given key @key 
 * has been inserted in rbtree rooted at @root; 0 otherwise
 */
static inline int rbInsertStdNewKey(rbRoot *root,rbNode *node, idx_t key)
{
	struct rb_node **new = &root->rb_root.rb_node, *parent = NULL;
	idx_t parentK;

	while (*new) {
		parent = *new;
		parentK = rb_entry(parent, rbNode, rb)->key;
		if	(key < parentK)		new = &parent->rb_left;
		else if	(key > parentK)		new = &parent->rb_right;
		else				return 0; //already in 
		DEBUGCHECKS			assert( *new != parent );
	}
	//insert the node in the correct position, assigning the key
	/*DEBUGCHECKS{	//check for double insertion
		rbNode testNode;
		memset(&testNode,0,sizeof(testNode));
		assert( !memcmp(node,&testNode,sizeof(*node)) );
  		for (struct rb_node* n = rb_first(&root->rb_root); n; n = rb_next(n))
			assert( n != *new );
	}*/
			
	node->key = key;
	rb_link_node(&node->rb, parent, new);
	rb_insert_color(&node->rb, &root->rb_root);
	return 1;
}

static inline int rbInsertCachedNewKey(rbRoot *root,rbNode *node, idx_t key)
{
	struct rb_node **new = &root->rb_root.rb_node, *parent = NULL;
	idx_t parentK;
	bool leftmost = true;

	while (*new) {
		parent = *new;
		parentK = rb_entry(parent, rbNode, rb)->key;
		if (key < parentK)			new = &parent->rb_left;
		else if (key > parentK) {
			new = &parent->rb_right;
			leftmost = false;
		}
		else					return 0;
	}

	//insert the node in the correct position, assigning the key
	/*DEBUGCHECKS{	//check for double insertion
		rbNode testNode;
		memset(&testNode,0,sizeof(testNode));
		assert( !memcmp(node,&testNode,sizeof(*node)) );
  		for (struct rb_node* n = rb_first(&root->rb_root); n; n = rb_next(n))
			assert( n != *new );
	}*/
			
	node->key = key;
	rb_link_node(&node->rb, parent, new);
	rb_insert_color_cached(&node->rb, root, leftmost);
	return 1;
}

static inline int rbInsertNewKey(rbRoot *root,rbNode *node, idx_t key){
	#if	RB_CACHED_INSERT == TRUE
	return rbInsertCachedNewKey(root,node,key);
	#else
	return rbInsertStdNewKey(root,node,key);
	#endif
}

#define rbNodeOrderedVisit(n,root) \
  for (n = rb_first(&root->rb_root); n; n = rb_next(n));

inline void cleanRbNodes(rbRoot* root,rbNode* nodes,idx_t nodesNum){
	memset(nodes,0,nodesNum * sizeof(*nodes));
	memset(root,0,sizeof(*root));
}
#endif	/* _LINUX_RBTREE_H */
