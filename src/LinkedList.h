#ifndef _DJS_LINKEDLIST_H_
#define _DJS_LINKEDLIST_H_

struct Link;

struct LinkedList
{
    struct Link* head; //!< head of the linked list, a dummy sentinel
    struct Link* tail; //!< tail of the linked list, a dummy sentinel
#ifdef SMARTLENGTH
    int length; //!< the number of data items in the linked list, if we are maintaining it.
#endif
};

struct Link
{
    int data; //!< arbitrary data stored in the link
    struct Link* next; //!< the previous link in the chain
    struct Link* prev; //!< the next link in the chain
#ifdef SMARTLENGTH
    struct LinkedList* linkedList; //!< the linked list that this link belongs to.
#endif
};

typedef struct Link Link;

typedef struct LinkedList LinkedList;

int isHead(Link* list);

int isTail(Link* list);

int deleteLink(Link* list);

Link* addAfter(Link* list, int data);

Link* addBefore(Link* list, int data);

void addLinkBefore(Link* list, Link* newLink);

Link* removeLink(Link* list);

LinkedList* createLinkedList(void);

void destroyLinkedList(LinkedList* linkedList);

void copyLinkedList(LinkedList* destination, 
                    LinkedList* source);

int contains(LinkedList* linkedList, int data, int (*comparator)(int,int));

int equal( LinkedList* list1, 
           LinkedList* list2, 
           int (*comparator)(int,int));

void restoreLinksWithReferences(LinkedList* list);

Link* addFirst(LinkedList* linkedList, int data);

Link* addLast(LinkedList* linkedList, int data);

Link* removeFirst(LinkedList* linkedList);

Link* removeLast(LinkedList* linkedList);

void deleteLast(LinkedList* linkedList);

int getFirst(LinkedList* linkedList);

void printListAbbv(LinkedList* linkedList, void (*printFunc)(int));

void printList(LinkedList* linkedList, void (*printFunc)(int));

int length(LinkedList* linkedList);

int isEmpty(LinkedList* linkedList);

#endif
