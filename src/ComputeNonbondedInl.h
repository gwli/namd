/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Common operations for ComputeNonbonded classes
*/

#ifndef COMPUTENONBONDEDINL_H
#define COMPUTENONBONDEDINL_H

#ifndef NAMD_RESTRICT
#define restrict
#endif

#include "ComputeNonbondedUtil.h"
#include "SimParameters.h"
#include "Node.h"
#include "Molecule.h"
#include "LJTable.h"
#include "ReductionMgr.h"
#include "ReserveArray.h"

#include "PressureProfile.h"

inline int pairlist_from_pairlist(BigReal cutoff2,
	BigReal p_i_x, BigReal p_i_y, BigReal p_i_z,
	const CompAtom *p_j,
	const plint *list, int list_size, plint *newlist) {

  plint *nli = newlist;
  if ( list_size > 0 ) {
    int j2 = list[0];
    BigReal p_j_x = p_j[j2].position.x;
    BigReal p_j_y = p_j[j2].position.y;
    BigReal p_j_z = p_j[j2].position.z;
    int g = 0;
    while ( g < list_size ) {
      int j = j2;
      j2 = list[++g];
      BigReal r2 = p_i_x - p_j_x;
      r2 *= r2;
      p_j_x = p_j[j2].position.x;
      BigReal t2 = p_i_y - p_j_y;
      r2 += t2 * t2;
      p_j_y = p_j[j2].position.y;
      t2 = p_i_z - p_j_z;
      r2 += t2 * t2;
      p_j_z = p_j[j2].position.z;
      if ( r2 <= cutoff2 ) { *nli= j; ++nli; }
    }
  }
  return nli - newlist;
}


// clear all
// define interaction type (pair or self)
#define NBPAIR	1
#define NBSELF	2

#endif // COMPUTENONBONDEDINL_H

