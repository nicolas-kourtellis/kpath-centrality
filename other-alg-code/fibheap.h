#ifndef _FIBHEAP_H
#define _FIBHEAP_H

//***************************************************************************
// The Fibonacci heap implementation contained in FIBHEAP.H and FIBHEAP.CPP
// is Copyright (c) 1996 by John Boyer
//
// Once this Fibonacci heap implementation (the software) has been published
// by Dr. Dobb's Journal, permission to use and distribute the software is
// granted provided that this copyright notice remains in the source and
// and the author (John Boyer) is acknowledged in works that use this program.
//
// Every effort has been made to ensure that this implementation is free of
// errors.  Nonetheless, the author (John Boyer) assumes no liability regarding
// your use of this software.
//
// The author would also be very glad to hear from anyone who uses the
// software or has any feedback about it.
// Email: jboyer@gulf.csc.uvic.ca
//***************************************************************************


/************* Change Log ************************
DATE                NAME              DESCRIPTION
11/19/2009          Vincent           Made the key and the vertex poistion into float's because I was unable
to overload the operators in my class that uses FibHeapNode as a base class.


**************************************************/

#include <limits.h>

#define OK      0
#define NOTOK   -1

//======================================================
// Fibonacci Heap Node Class
//======================================================

class FibHeap;

class FibHeapNode
{
  friend class FibHeap;

public:

  FibHeapNode *Left, *Right, *Parent, *Child;
  unsigned long Degree, Mark, NegInfinityFlag;
  //additional variables for my projrct
  unsigned long nhopps, vertexPosition;
  double key;

protected:

  inline int FHN_Cmp(FibHeapNode& RHS);
  inline void FHN_Assign(FibHeapNode& RHS);

public:

  FibHeapNode();
  virtual ~FibHeapNode();

  virtual void operator =(FibHeapNode& RHS);
  virtual int  operator ==(FibHeapNode& RHS);
  virtual int  operator <(FibHeapNode& RHS);

  //virtual void Print();

  //additional functions for my project
  unsigned long Get_nhopps ( ) const { return nhopps; };
  unsigned long Get_vertexPosition ( ) const { return vertexPosition; };
  double Get_key ( ) const { return key; };

  void Set_nhopps ( unsigned long NH ) { nhopps = NH; };
  void Set_vertexPosition ( unsigned long VP ) { vertexPosition = VP; };
  void Set_key ( double K ) { key = K; };
};

//========================================================================
// Fibonacci Heap Class
//========================================================================

class FibHeap
{
  FibHeapNode *MinRoot;
  long NumNodes, NumTrees, NumMarkedNodes;
  int HeapOwnershipFlag;

public:

  FibHeap();
  virtual ~FibHeap();

  // The Standard Heap Operations

  void Insert(FibHeapNode *NewNode);
  void Union(FibHeap *OtherHeap);

  inline FibHeapNode *Minimum();
  FibHeapNode *ExtractMin();

  int DecreaseKey(FibHeapNode *theNode, FibHeapNode& NewKey);
  int Delete(FibHeapNode *theNode);

  // Extra utility functions

  int  GetHeapOwnership() { return HeapOwnershipFlag; };
  void SetHeapOwnership() { HeapOwnershipFlag = 1; };
  void ClearHeapOwnership() { HeapOwnershipFlag = 0; };

  long GetNumNodes() { return NumNodes; };
  long GetNumTrees() { return NumTrees; };
  long GetNumMarkedNodes() { return NumMarkedNodes; };

  //void Print(FibHeapNode *Tree = NULL, FibHeapNode *theParent=NULL);

private:

  // Internal functions that help to implement the Standard Operations

  inline void _Exchange(FibHeapNode*&, FibHeapNode*&);
  void _Consolidate();
  void _Link(FibHeapNode *, FibHeapNode *);
  void _AddToRootList(FibHeapNode *);
  void _Cut(FibHeapNode *, FibHeapNode *);
  void _CascadingCut(FibHeapNode *);
};

#endif
