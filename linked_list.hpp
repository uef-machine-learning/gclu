#ifndef _LINKED_LIST_HPP
#define _LINKED_LIST_HPP


struct gItem;

typedef struct linkedListNode {
  void *content;
  linkedListNode *next;
  linkedListNode *prev;
} linkedListNode;

typedef struct linkedList {
  int size;
  linkedListNode *root;
  linkedListNode *end;
} linkedList;


linkedList *initLinkedList();
void ll_free_list(linkedList *ll);
linkedListNode *ll_add_node(linkedList *ll, void *content);
void ll_add_after_node(linkedList *ll, linkedListNode *node, void *content);
void ll_add_node_if_not_exist(linkedList *ll, void *content);
gItem* ll_get_node_if_exist(linkedList *ll, int id);
linkedListNode *ll_pop_first_node(linkedList *ll);
void ll_remove_node(linkedList *ll, int idx);
void ll_remove_node2(linkedList *ll, int idx);
void *ll_get_item(linkedList *ll, int idx);
void ll_write_ints_to_file(linkedList *ll, const char *fn, int N);


#endif
