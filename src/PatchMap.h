/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: PatchMap.h
 *
 ***************************************************************************/

#ifndef PATCHMAP_H
#define PATCHMAP_H

#include "NamdTypes.h"
#include "PatchMgr.h"

class Patch;
class PatchMgr;

class PatchMap
{
public:
  static PatchMap *Instance();
  inline static PatchMap *Object() { return _instance; }


  ~PatchMap(void);

  enum { MaxTwoAway = 5*5*5 - 3*3*3 };
  enum { MaxOneAway = 3*3*3 - 1 };
  enum { MaxOneOrTwoAway = MaxOneAway + MaxTwoAway };
  enum ErrCode { OK = 0, ERROR = -1 };

  static void registerPatchMgr(PatchMgr *pmgr) {
    patchMgr = pmgr;
  }

  HomePatchList *homePatchList();
  int numHomePatches(void);

  // numPatches() returns the number of patches being managed 
  // by the map.
  int numPatches(void);

  // xDimension() returns the how many patches span the simulation space
  // along the x axis.
  int xDimension(void);

  // yDim() returns the how many patches span the simulation space
  // along the y axis.
  int yDimension(void);

  // zDim() returns the how many patches span the simulation space
  // along the z axis.
  int zDimension(void);

  // pid(xindex, yindex, zindex) returns the patch id for the given
  // patch coordinates.
  int pid(int xindex, int yindex, int zindex);

  // xIndex(pid) returns the x index for the given patch id.
  int xIndex(int pid);

  // yIndex(pid) returns the y index for the given patch id.
  int yIndex(int pid);

  // zIndex(pid) returns the z index for the given patch id.
  int zIndex(int pid);

  // node(pid) returns the node where the patch currently exists.
  int node(int pid);

  // minX(pid) returns the minimum x coordinate of the region of
  // space the patch is responsible for.
  Coordinate minX(int pid);

  // maxX(pid) returns the maximum x coordinate of the region of
  // space the patch is responsible for.
  Coordinate maxX(int pid);

  // minY(pid) returns the minimum y coordinate of the region of
  // space the patch is responsible for.
  Coordinate minY(int pid);

  // maxY(pid) returns the maximum y coordinate of the region of
  // space the patch is responsible for.
  Coordinate maxY(int pid);
  
  // minZ(pid) returns the minimum z coordinate of the region of
  // space the patch is responsible for.
  Coordinate minZ(int pid);

  // maxZ(pid) returns the maximum z coordinate of the region of
  // space the patch is responsible for.
  Coordinate maxZ(int pid);

  // numCids(pid) returns the number of compute ids which are registered
  // with the patch.  
  int numCids(int pid);
  
  // cid(pid,i) returns the i-th compute id registered
  // with the patch.  
  int cid(int pid, int i);

  // allocatePids(x_dim, y_dim, z_dim) tells the PatchMap to set aside
  // room for x_dim * y_dim * z_dim patches.
  ErrCode allocatePids(int x_dim, int y_dim, int z_dim);

  // PatchID requestPid(int *xi, int *yi, int *zi) tells the
  // PatchMap to return the next unused pid, and the indices of
  // that patch.
  PatchID requestPid(int *xi, int *yi, int *zi);

  // storePatch(pid, node, max_computes, x0, y0, z0, x1, y1, z1)
  // stores the info about the patch into the previously requested
  // pid.
  void storePatch(PatchID pid, int node, int max_computes,
		  Coordinate x0, Coordinate y0, Coordinate z0,
		  Coordinate x1, Coordinate y1, Coordinate z1);

  // newCid(pid,cid) stores a compute id associated with
  // patch id pid.  Error returned when there is no room to store
  // the pid.
  ErrCode newCid(int pid, int cid);

  // oneAwayNeighbors(pid, &n, neighbor_ids) returns the number 
  // and ids of adjacent patches.  The caller is expected to provide
  // sufficient storage for the neighbors.

  int oneAwayNeighbors(int pid, PatchID *neighbor_ids);

  // twoAwayNeighbors(pid, &n, neighbor_ids) returns the number 
  // and ids of all patches exactely two steps distant.  
  // The caller is expected to provide sufficient storage for the neighbors.

  int twoAwayNeighbors(int pid, PatchID *neighbor_ids);

  int oneOrTwoAwayNeighbors(int pid, PatchID *neighbor_ids);

  void printPatchMap(void);

  inline Patch *patch(PatchID pid);

  void registerPatch(PatchID pid, Patch *pptr);
  void unregisterPatch(PatchID pid, Patch *pptr);


protected:
  friend MapDistribMsg;
  void * pack (int *length);
  void unpack (void *in);

  PatchMap(void);

private:
  static PatchMap *_instance;
  static PatchMgr *patchMgr;

  struct PatchData
  {
    int node;
    int xi, yi, zi;
    Coordinate x0, x1, y0, y1, z0, z1;
    int numCids;
    int numCidsAllocated;
    ComputeID *cids;
    Patch *myPatch;
  };
  int curPatch;
  int nPatches;
  PatchData *patchData;
  int xDim, yDim, zDim;

};


//----------------------------------------------------------------------

inline Patch *PatchMap::patch(PatchID pid)
{
  return patchData[pid].myPatch;
}

#endif /* PATCHMAP_H */


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMap.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:36:47 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchMap.h,v $
 * Revision 1.777  1997/01/17 19:36:47  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.10  1997/01/13 21:04:19  jim
 * added numHomePatches()
 *
 * Revision 1.9  1996/12/18 21:07:54  jim
 * added oneOrTwoAwayNeighbors()
 *
 * Revision 1.8  1996/12/12 08:57:17  jim
 * added MapDistribMsg packing / unpacking routines
 *
 * Revision 1.7  1996/11/21 20:39:29  jim
 * small bug fixes
 *
 * Revision 1.6  1996/11/01 21:20:45  ari
 * *** empty log message ***
 *
 * Revision 1.5  1996/10/29 23:35:27  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/10/10 17:23:24  brunner
 * Added patch * in patchmap
 *
 * Revision 1.3  1996/08/23 22:03:52  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/19 21:37:02  brunner
 * Changed Position to Coordinate
 *
 * Revision 1.1  1996/08/16 20:43:53  brunner
 * Initial revision
 *
 * Revision 1.7  1996/08/03 20:08:09  brunner
 * *** empty log message ***
 *
 * Revision 1.6  1996/07/16 20:06:24  brunner
 * *** empty log message ***
 *
 * Revision 1.5  1996/07/16 19:59:00  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/06/11 22:36:35  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/06/10 18:52:26  brunner
 * Initial revision
 *
 *
 ***************************************************************************/
