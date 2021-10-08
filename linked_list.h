#ifndef _LINKED_LIST_H
#define _LINKED_LIST_H


typedef int BOOL;
typedef struct gItem {
    int id;
    double dist;
    BOOL new_item;
    int visited;
} gItem;


typedef struct linkedListNode {
    void * content;
    linkedListNode* next;
} linkedListNode;

typedef struct linkedList {
    int size;
    linkedListNode* root;
} linkedList;

linkedList* initLinkedList() {
    linkedList* ll = (linkedList*) malloc(sizeof(linkedList));
    ll->size=0;
    /*ll->root = NULL;*/
    ll->root = (linkedListNode*) malloc(sizeof(linkedListNode));

    return ll;
}


void ll_free_list(linkedList* ll) {
    int i;
    linkedListNode* node = ll->root;
    linkedListNode* last=node;
    for(i=0;i<ll->size;i++) {node = node->next;free(last);last=node;}
}


void ll_add_node(linkedList* ll,void* content) {
    int i;
    linkedListNode* node = ll->root;
    for(i=0;i<ll->size;i++) {node = node->next;}
    node->next = (linkedListNode*) malloc(sizeof(linkedListNode));
    node->content = content;
    ll->size++;
}

void ll_add_node_if_not_exist(linkedList* ll,void* content) {
    int i;
    linkedListNode* node = ll->root;
    for(i=0;i<ll->size;i++) {
        if(node->content == content) {return;}
        node = node->next;
    }
    node->next = (linkedListNode*) malloc(sizeof(linkedListNode));
    node->content = content;
    ll->size++;
}

gItem* ll_get_node_if_exist(linkedList* ll,int id) {
    int i;
    linkedListNode* node = ll->root;
    for(i=0;i<ll->size;i++) {
        gItem* gi = (gItem*) node->content;
        if( gi->id == id) {return gi;}
        node = node->next;
    }
    return NULL;
}



void ll_remove_node(linkedList* ll,int idx) {
    int i;
    linkedListNode* node = ll->root;
    linkedListNode* last=node;
    for(i=0;i<ll->size;i++) {
        if (*((int*) node->content) == idx) {
            if(i==0) { ll->root=node->next;}
            else {
                last->next = node->next;
            }
            free(node);
            ll->size--;
            break;
        }
        last=node;
        node = node->next;
    }
}

void ll_remove_node2(linkedList* ll,int idx) {
    int i;
    linkedListNode* node = ll->root;
    linkedListNode* last=node;
    for(i=0;i<ll->size;i++) {
        /*if (*((int*) node->content) == idx) {*/
        if (((gItem*) node->content)->id==idx) {

            if(i==0) { ll->root=node->next;}
            else {
                last->next = node->next;
            }
            free(node);
            ll->size--;
            break;
        }
        last=node;
        node = node->next;
    }
}



void* ll_get_item(linkedList* ll, int idx) {
    int i;
    linkedListNode* node = ll->root;
    for(i=0;i<ll->size && i < idx ;i++) {node = node->next;}
    return node->content;
}

void ll_write_ints_to_file(linkedList* ll, const char* fn, int N) {
    int i;
    linkedListNode* node = ll->root;
    for(i=0;i<ll->size;i++) {node = node->next;}

    FILE *fp;
    fp = fopen(fn,"w");
    for(int i=0;i<ll->size;i++) {
        for(int j=0;j<N;j++) {
            fprintf(fp,"%d ", ((int*) ll_get_item(ll,i))[j]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);


}



#endif
