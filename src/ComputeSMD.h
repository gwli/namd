/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTESMD_H
#define COMPUTESMD_H

#include "ComputeHomePatch.h"
#include "ReductionMgr.h"

class ComputeSMD : public ComputeHomePatch
{
private:
	int consExp;		//  Exponent for energy function from SimParameters
        BigReal k;              //  Restraint force constant
	int moveAtom;           //  Index of the atom to move
        BigReal moveVel;         // velocity of the restraint movement 
                                // (A/timestep).
        int outputFreq;         // output frequency
        Bool chDirOn;           // is changing direction on?
        Bool chForceOn;         // is changing force on?
	Bool projectForce;	// If true, force applied only in pulling
				// direction.

public:
	ComputeSMD(ComputeID c, PatchID pid); 	//  Constructor
	virtual ~ComputeSMD();			//  Destructor

	virtual void doForce(Position* p, Results* r, AtomProperties* a, Transform* t);

	SubmitReduction *reduction;

};

#endif







