#ifndef COLLECTIONMASTER_H
#define COLLECTIONMASTER_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "main.h"
#include "NamdTypes.h"

class CollectVectorMsg;

class CollectionMaster : public chare_object
{
public:

  CollectionMaster(void);
  ~CollectionMaster(void);

  void receivePositions(CollectVectorMsg *msg);
  void receiveVelocities(CollectVectorMsg *msg);
  void receiveForces(CollectVectorMsg *msg);

private:

  class CollectVectorInstance
  {
  public:

    CollectVectorInstance(void) : seq(-1) { ; }

    CollectVectorInstance(int s) :
      seq(s), remaining(CNumPes()) { ; }

    // true -> send it and delete it!
    int append(AtomIDList &a, ResizeArray<Vector> &d)
    {
      int size = a.size();
      for( int i = 0; i < size; ++i )
      {
	aid.add(a[i]);
	data.add(d[i]);
      }
      return ( ! --remaining );
    }

    int seq;
    AtomIDList aid;
    ResizeArray<Vector> data;

    operator<(CollectVectorInstance &o) { return (seq < o.seq); }
    operator==(CollectVectorInstance &o) { return (seq == o.seq); }
    void * operator new(size_t size) { return ::operator new(size); }
    void operator delete(void* ptr) { ::operator delete(ptr); }

  private:
    int remaining;

  };

  class CollectVectorSequence
  {
  public:

    CollectVectorInstance* submitData(
	int seq, AtomIDList &i, ResizeArray<Vector> &d)
    {
      CollectVectorInstance *c = data.find(CollectVectorInstance(seq));
      if ( ! c )
      {
	data.add(CollectVectorInstance(seq));
	c = data.find(CollectVectorInstance(seq));
      }
      if ( c->append(i,d) )
      {
	c = new CollectVectorInstance(*c);
	data.del(CollectVectorInstance(seq));
        return c;
      }
      else
      {
        return 0;
      }
    }

    ResizeArray<CollectVectorInstance> data;

    void * operator new(size_t size) { return ::operator new(size); }
    void operator delete(void* ptr) { ::operator delete(ptr); }
  };

  CollectVectorSequence positions;
  CollectVectorSequence velocities;
  CollectVectorSequence forces;

};


class CollectVectorMsg : public comm_object
{
public:

  int seq;
  AtomIDList aid;
  ResizeArray<Vector> data;

  void * pack(int *length);
  void unpack(void *in);

  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }

};



#endif
