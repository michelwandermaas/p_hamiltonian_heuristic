#include "BinaryHeap.h"
#include <stdio.h>

void BinaryHeap::Insert(BHNode *node)
{
//printf("1\n");
	int i;
	for(i = ++size; (elements[ i>>1 ]->key - node->key > 0.00001)&&(i>>1 != 0); i /= 2)
	{
		elements[i] = elements[i>>1];
		elements[i]->position = i;
	}
	elements[i] = node;
	elements[i]->position = i;
//printf("2\n");
}

BHNode* BinaryHeap::DeleteMin()
{
	if(size == 0) return NULL;

	BHNode *min = elements[1];
	BHNode *last = elements[size--];

	int i, child;
	for(i = 1; i << 1 <= size; i = child)
	{
		child = i<<1;
		if(child != size && elements[child]->key - elements[child + 1]->key > 0.00001)
			child++;

		if(last->key - elements[child]->key > 0.00001)
		{
			elements[i] = elements[child];
			elements[child]->position = i;
		}
		else
			break;
	}
	elements[i]  = last;
	elements[i]->position = i;

	return min;
}

void BinaryHeap::DecreaseKey(BHNode *node)
{
	int i;
	for(i = node->position; (elements[ i>>1 ]->key - node->key > 0.00001)&&(i>>1 != 0); i /= 2)
	{
		elements[i] = elements[i>>1];
		elements[i]->position = i;
	}
	elements[i] = node;
	elements[i]->position = i;	
}

BinaryHeap::BinaryHeap(int n)
{
	sentinel = new BHNode;
	sentinel->key = MINUS_INFINITY;
	
	elements = new BHNode*[n];
	for(int i = 0; i < n; i++) elements[i] = NULL;

	elements[0] = sentinel;

	size = 0;
}

BinaryHeap::~BinaryHeap()
{
	delete sentinel;
	delete [] elements;
}

void BinaryHeap::Clear()
{
	size = 0;
}
