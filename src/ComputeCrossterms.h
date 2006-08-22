/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTECROSSTERMS_H
#define COMPUTECROSSTERMS_H

#include "ComputeHomeTuples.h"
#include "ComputeSelfTuples.h"
#include "ReductionMgr.h"

class Molecule;
class CrosstermValue;

class CrosstermElem {
public:
    // ComputeHomeTuples interface
    enum { size = 8 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    Real scale;
    void computeForce(BigReal*, BigReal *);
    // The following is evil, but the compiler chokes otherwise. (JCP)
    static void loadTuplesForAtom(void*, AtomID, Molecule*);
    static void getMoleculePointers(Molecule*, int*, int32***, Crossterm**);
    static void getParameterPointers(Parameters*, const CrosstermValue**);

    // pressure profile parameters
    static int pressureProfileSlabs;
    static int pressureProfileAtomTypes;
    static BigReal pressureProfileThickness;
    static BigReal pressureProfileMin;

    // Internal data
    const CrosstermValue *value;

  int hash() const {
    return 0x7FFFFFFF &((atomID[0]<<24) + (atomID[1]<<16) + (atomID[2]<<8) + atomID[3]);
  }
  enum { crosstermEnergyIndex, TENSOR(virialIndex), reductionDataSize };
  enum { reductionChecksumLabel = REDUCTION_CROSSTERM_CHECKSUM };
  static void submitReductionData(BigReal*,SubmitReduction*);

  CrosstermElem() {
	atomID[0] = -1;
	atomID[1] = -1;
	atomID[2] = -1;
	atomID[3] = -1;
	atomID[4] = -1;
	atomID[5] = -1;
	atomID[6] = -1;
	atomID[7] = -1;
	p[0] = NULL;
	p[1] = NULL;
	p[2] = NULL;
	p[3] = NULL;
	p[4] = NULL;
	p[5] = NULL;
	p[6] = NULL;
	p[7] = NULL;
  }
  CrosstermElem(const Crossterm *a, const CrosstermValue *v) {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    atomID[2] = a->atom3;
    atomID[3] = a->atom4;
    atomID[4] = a->atom5;
    atomID[5] = a->atom6;
    atomID[6] = a->atom7;
    atomID[7] = a->atom8;
    value = &v[a->crossterm_type];
  }

  CrosstermElem(AtomID atom0, AtomID atom1, AtomID atom2, AtomID atom3,
                AtomID atom4, AtomID atom5, AtomID atom6, AtomID atom7) {
    atomID[0] = atom0;
    atomID[1] = atom1;
    atomID[2] = atom2;
    atomID[3] = atom3;
    atomID[4] = atom4;
    atomID[5] = atom5;
    atomID[6] = atom6;
    atomID[7] = atom7;
  }
  ~CrosstermElem() {};

  int operator==(const CrosstermElem &a) const {
    return (a.atomID[0] == atomID[0] && a.atomID[1] == atomID[1] &&
        a.atomID[2] == atomID[2] && a.atomID[3] == atomID[3] &&
        a.atomID[4] == atomID[4] && a.atomID[5] == atomID[5] &&
        a.atomID[6] == atomID[6] && a.atomID[7] == atomID[7]);
  }

  int operator<(const CrosstermElem &a) const {
    return  (atomID[0] < a.atomID[0] ||
            (atomID[0] == a.atomID[0] &&
            (atomID[1] < a.atomID[1] ||
            (atomID[1] == a.atomID[1] &&
            (atomID[2] < a.atomID[2] ||
            (atomID[2] == a.atomID[2] &&
            (atomID[3] < a.atomID[3] ||
            (atomID[3] == a.atomID[3] &&
            (atomID[4] < a.atomID[4] ||
            (atomID[4] == a.atomID[4] &&
            (atomID[5] < a.atomID[5] ||
            (atomID[5] == a.atomID[5] &&
            (atomID[6] < a.atomID[6] ||
            (atomID[6] == a.atomID[6] &&
             atomID[7] < a.atomID[7] 
	     ))))))))))))));
  }
};

class ComputeCrossterms : public ComputeHomeTuples<CrosstermElem,Crossterm,CrosstermValue>
{
public:

  ComputeCrossterms(ComputeID c, PatchIDList p) : ComputeHomeTuples<CrosstermElem,Crossterm,CrosstermValue>(c,p) { ; }

};

class ComputeSelfCrossterms : public ComputeSelfTuples<CrosstermElem,Crossterm,CrosstermValue>
{
public:

  ComputeSelfCrossterms(ComputeID c, PatchID p) : ComputeSelfTuples<CrosstermElem,Crossterm,CrosstermValue>(c,p) { ; }

};

#endif
