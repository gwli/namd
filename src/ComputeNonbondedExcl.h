/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTENONBONDEDEXCL_H
#define COMPUTENONBONDEDEXCL_H

#include "ComputeHomeTuples.h"
#include "ComputeNonbondedUtil.h"
#include "ReductionMgr.h"

class Molecule;
class ExclusionValue;

class NonbondedExclElem : public ComputeNonbondedUtil {
public:
  // ComputeHomeTuples interface
  enum { size = 2 };
  AtomID atomID[size];
  int    localIndex[size];
  TuplePatchElem *p[size];
  void computeForce(BigReal*);

  enum { reductionChecksumLabel = REDUCTION_EXCLUSION_CHECKSUM };
  static void loadTuplesForAtom(void*, AtomID, Molecule*) {};
  static void getMoleculePointers(Molecule*, int*, int***, Exclusion**);
  static void getParameterPointers(Parameters*, const ExclusionValue** v) { *v=0; }

  int hash() const {
    return 0x7FFFFFFF & (atomID[1] << 16 + atomID[0]);
  }

  // Internal data
  Index modified;

  NonbondedExclElem() {
    atomID[0] = -1;
    atomID[1] = -1;
    p[0] = NULL;
    p[1] = NULL;
  }
  NonbondedExclElem(const Exclusion *a, const ExclusionValue *) {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    modified = a->modified;
  }

  NonbondedExclElem(AtomID atom0, AtomID atom1) {
    if (atom0 > atom1) {  // Swap end atoms so lowest is first!
      AtomID tmp = atom1; atom1 = atom0; atom0 = tmp; 
    }
    atomID[0] = atom0;
    atomID[1] = atom1;
  }
  ~NonbondedExclElem() {};

  int operator==(const NonbondedExclElem &a) const {
    return (a.atomID[0] == atomID[0] && a.atomID[1] == atomID[1]);
  }

  int operator<(const NonbondedExclElem &a) const {
    return (atomID[0] < a.atomID[0] ||
            (atomID[0] == a.atomID[0] &&
            (atomID[1] < a.atomID[1]) ));
  }

  private:
};

class ComputeNonbondedExcls : public ComputeHomeTuples<NonbondedExclElem,Exclusion,ExclusionValue>
{
public:

  ComputeNonbondedExcls(ComputeID c) : ComputeHomeTuples<NonbondedExclElem,Exclusion,ExclusionValue>(c) 
  { }

  // void loadTuples(); //overload of the template version
};

#endif

