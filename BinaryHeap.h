#ifndef BINARY_HEAP_H
#define BINARY_HEAP_H

#include <math.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

#define MINUS_INFINITY -10000000000.0

class BHNode
{
public:
	double key;
	int p;
	int position;
};

class BinaryHeap
{
public:
	BinaryHeap(int n);
	~BinaryHeap();

	void Insert(BHNode *node);
	BHNode* DeleteMin();
	void DecreaseKey(BHNode *node);
	void Clear();

	BHNode *sentinel;
	BHNode **elements;
	
	int size;
};


#endif


