/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef COMPUTEATOM_H
#define COMPUTEATOM_H

#include "main.h"
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "NamdTypes.h"

#include "Templates/Box.h"
#include "Templates/OwnerBox.h"

class Patch;

// Compute object for Atom-wise type
class ComputeHomeAtoms : public Compute
{
  // Array of atoms for which we must report forces.
  AtomArray atom;

  // Array of patches for which we must report all appropo
  // bonded-atom type forces.
  PatchArray patch;
  void (* mapMethod)();
  void mapByAtoms();
  void mapByPatches();
  void mapByHomePatches();

protected:
  // Take list of atoms and registers to appropriate patches
  // which are not necessarily the same as atoms
  // for which forces are computed (generally additional atoms
  // are needed for computation of forces)
  registerDependency(AtomArray &);

public:
  // void argument constructor implies this object will
  // do Atom type computation appropriate atoms in Home Patches.
  ComputeAtoms();

  // Initializes this object to compute appropriate forces
  // for the atoms in AtomArray.
  ComputeAtoms(PatchArray &);
  ComputeAtoms(AtomArray &);
  ~ComputeAtoms(); 

  // Signal to take known atoms and register 
  // dependencies to appropriate patches.
  void mapAtoms();

  // work method
  void virtual doWork(LocalWorkMsg *msg);
};


#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeHomeAtoms.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:35:43 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeHomeAtoms.h,v $
 * Revision 1.777  1997/01/17 19:35:43  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1996/11/01 21:20:45  ari
 * Initial revision
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 * Revision 1.4  1996/07/16 01:54:12  ari
 * *** empty log message ***
 *
 * Revision 1.3  96/07/16  01:10:26  01:10:26  ari (Aritomo Shinozaki)
 * Fixed comments, added methods
 * 
 * Revision 1.2  1996/06/25 21:10:48  gursoy
 * *** empty log message ***
 *
 * Revision 1.1  1996/06/24 14:12:26  gursoy
 * Initial revision
 *
 ***************************************************************************/

