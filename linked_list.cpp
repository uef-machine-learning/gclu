
#include "linked_list.hpp"




linkedList *initLinkedList() {
  linkedList *ll = (linkedList *)malloc(sizeof(linkedList));
  ll->size = 0;
  ll->root = NULL;
  ll->end = NULL;
  // ll->root = (linkedListNode *)malloc(sizeof(linkedListNode));
  // ll->root->next = NULL;
  // ll->root->prev = NULL;
  return ll;
}

void ll_free_list(linkedList *ll) {
  int i;
  linkedListNode *node = ll->root;
  linkedListNode *last = node;
  for (i = 0; i < ll->size; i++) {
    node = node->next;
    free(last->content);
    free(last);
    last = node;
  }
}

linkedListNode *ll_add_node(linkedList *ll, void *content) {
  int i;
  linkedListNode *node;

  linkedListNode *newnode = (linkedListNode *)malloc(sizeof(linkedListNode));
  newnode->content = content;
  newnode->next = NULL;
  newnode->prev = NULL;

  if (ll->size == 0) {
    ll->root = newnode;
    ll->end = newnode;
  } else {
    node = ll->end;
    node->next = newnode;
    newnode->prev = node;
  }
  ll->end = newnode;
  ll->size++;
  return newnode;
}

void ll_add_after_node(linkedList *ll, linkedListNode *node, void *content) {
  int i;
  // linkedListNode* node = ll->root;
  // for(i=0;i<ll->size;i++) {node = node->next;}
  linkedListNode *newnode = (linkedListNode *)malloc(sizeof(linkedListNode));
  newnode->content = content;
  newnode->next = node->next;
  newnode->prev = node;
  if (newnode->next != NULL) {
    newnode->next->prev = newnode;
  }

  if (node == ll->end) {
    ll->end = newnode;
  }
  node->next = newnode;
  ll->size++;
}

void ll_add_node_if_not_exist(linkedList *ll, void *content) {
  int i;
  linkedListNode *node = ll->root;
  for (i = 0; i < ll->size; i++) {
    if (node->content == content) {
      return;
    }
    node = node->next;
  }
  node->next = (linkedListNode *)malloc(sizeof(linkedListNode));
  node->content = content;
  ll->size++;
}

gItem* ll_get_node_if_exist(linkedList *ll, int id) {
  int i;
  linkedListNode *node = ll->root;
  for (i = 0; i < ll->size; i++) {
    gItem *gi = (gItem *)node->content;
    if (gi->id == id) {
      return gi;
    }
    node = node->next;
  }
  return NULL;
}

linkedListNode *ll_pop_first_node(linkedList *ll) {
  int i;
  linkedListNode *node = ll->root;
  if (ll->size > 1) {
    ll->root = node->next;
  } else {
    ll->root = NULL;
  }
  ll->size--;
  return node;
}

void ll_remove_node(linkedList *ll, int idx) {
  int i;
  linkedListNode *node = ll->root;
  linkedListNode *last = node;
  for (i = 0; i < ll->size; i++) {
    if (*((int *)node->content) == idx) {
      if (i == 0) {
        ll->root = node->next;
      } else {
        last->next = node->next;
      }
      free(node);
      ll->size--;
      break;
    }
    last = node;
    node = node->next;
  }
}

void ll_remove_node2(linkedList *ll, int idx) {
  int i;
  linkedListNode *node = ll->root;
  linkedListNode *last = node;
  for (i = 0; i < ll->size; i++) {
    /*if (*((int*) node->content) == idx) {*/
    if (((gItem *)node->content)->id == idx) {

      if (i == 0) {
        ll->root = node->next;
      } else {
        last->next = node->next;
      }
      free(node);
      ll->size--;
      break;
    }
    last = node;
    node = node->next;
  }
}

void *ll_get_item(linkedList *ll, int idx) {
  int i;
  linkedListNode *node = ll->root;
  for (i = 0; i < ll->size && i < idx; i++) {
    node = node->next;
  }
  return node->content;
}

void ll_write_ints_to_file(linkedList *ll, const char *fn, int N) {
  int i;
  linkedListNode *node = ll->root;
  for (i = 0; i < ll->size; i++) {
    node = node->next;
  }

  FILE *fp;
  fp = fopen(fn, "w");
  for (int i = 0; i < ll->size; i++) {
    for (int j = 0; j < N; j++) {
      fprintf(fp, "%d ", ((int *)ll_get_item(ll, i))[j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

