/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Common operations for ComputeNonbonded classes
*/

#include "ComputeNonbondedInl.h"

// 3 inline functions to handle explicit calculation of separation-shifted
// vdW for FEP and TI, and a shifted electrostatics potential for decoupling

/* ******************************************** */
/* vdW energy, force and lambda2 energy for FEP */
/* ******************************************** */
inline void fep_vdw_forceandenergies (BigReal A, BigReal B, BigReal r2, 
  BigReal myVdwShift, BigReal myVdwShift2, BigReal myVdwShift3,BigReal switchdist2, 
  BigReal cutoff2, BigReal myVdwLambda, BigReal myVdwLambda2, BigReal myVdwLambda3, 
  Bool Fep_WCA_repuOn, Bool Fep_WCA_dispOn, bool Fep_Wham, BigReal WCA_rcut1, 
  BigReal WCA_rcut2, BigReal WCA_rcut3, BigReal switchfactor, 
  Bool vdwForceSwitching, Bool LJcorrection, BigReal* alch_vdw_energy,
  BigReal* alch_vdw_force, BigReal* alch_vdw_energy_2, BigReal* alch_vdw_energy_3,
  BigReal* alch_vdw_energy_2_Left,BigReal* alch_vdw_energy_3_Left,BigReal myVdwShift_r,BigReal myVdwLambda_r) {
  // switching function (this is correct whether switching is active or not)
  const BigReal switchmul = (r2 > switchdist2 ? switchfactor*(cutoff2 - r2) \
			     *(cutoff2 - r2)				\
			     *(cutoff2 - 3.*switchdist2 + 2.*r2) : 1.);
  const BigReal switchmul2 = (r2 > switchdist2 ?			\
			      12.*switchfactor*(cutoff2 - r2)		\
			      *(r2 - switchdist2) : 0.);
  
  *alch_vdw_energy_2_Left = 0.;
  *alch_vdw_energy_3_Left=0;

  if(Fep_WCA_repuOn){
    if(B != 0.0) {// Excluded when some atoms like H, drude, lone pair are involved
      const BigReal Emin = B*B/(4.0*A);
      const BigReal Rmin_SQ = powf(2.0*A/B, 1.f/3);
      const BigReal r2_1 = r2 + (1.-WCA_rcut1)*(1.-WCA_rcut1)*Rmin_SQ;
      const BigReal r2_2 = r2 + (1.-WCA_rcut2)*(1.-WCA_rcut2)*Rmin_SQ;
      const BigReal r2_3 = r2 + (1.-WCA_rcut3)*(1.-WCA_rcut3)*Rmin_SQ;
      const BigReal r6_1 = r2_1*r2_1*r2_1;
      const BigReal r6_2 = r2_2*r2_2*r2_2;
      const BigReal r6_3 = r2_3*r2_3*r2_3;
      const BigReal WCA_rcut1_energy = (r2 <= Rmin_SQ*(1.0 - (1.0-WCA_rcut1) \
						       *(1.0-WCA_rcut1)) ? \
					A/(r6_1*r6_1) - B/r6_1 + Emin : 0.);
      const BigReal WCA_rcut2_energy = (r2 <= Rmin_SQ*(1.0 - (1.0-WCA_rcut2) \
						       *(1.0-WCA_rcut2)) ? \
					A/(r6_2*r6_2) - B/r6_2 + Emin : 0.);
      const BigReal WCA_rcut3_energy = (r2 <= Rmin_SQ*(1.0 - (1.0-WCA_rcut3) \
						       *(1.0-WCA_rcut3)) ? \
					A/(r6_3*r6_3) - B/r6_3 + Emin : 0.);
      const BigReal WCA_rcut1_force = (r2 <= Rmin_SQ*(1.0 - (1.0-WCA_rcut1) \
						      *(1.0-WCA_rcut1)) ? \
				       (12.*(WCA_rcut1_energy)		\
					+ 6.*B/r6_1 - 12.0 * Emin )/r2_1: 0.);
      const BigReal WCA_rcut2_force = (r2 <= Rmin_SQ*(1.0 - (1.0-WCA_rcut2) \
						      *(1.0-WCA_rcut2))? \
				       (12.*(WCA_rcut2_energy)		\
					+ 6.*B/r6_2 - 12.0 * Emin )/r2_2: 0.);
      const BigReal WCA_rcut3_force = (r2 <= Rmin_SQ*(1.0 - (1.0-WCA_rcut3) \
						      *(1.0-WCA_rcut3)) ? \
				       (12.*(WCA_rcut3_energy)		\
					+ 6.*B/r6_3 - 12.0 * Emin )/r2_3: 0.);
      // separation-shifted repulsion force and energy
     //DoubeWide: Update the alch_vdw_energy of the "reverse" calculation
      *alch_vdw_energy = WCA_rcut2_energy; 
      *alch_vdw_force = WCA_rcut2_force;
      if(WCA_rcut1 < WCA_rcut2) {
	*alch_vdw_energy_2_Left = *alch_vdw_energy + WCA_rcut2_energy - WCA_rcut1_energy;
	*alch_vdw_energy_3_Left = *alch_vdw_energy + WCA_rcut2_energy-WCA_rcut1_energy; 
      }
      if(WCA_rcut2 < WCA_rcut3) {
	*alch_vdw_energy_2 = *alch_vdw_energy + WCA_rcut3_energy - WCA_rcut2_energy; 
	*alch_vdw_energy_3=*alch_vdw_energy  + WCA_rcut3_energy-WCA_rcut2_energy;
      }
    }
    else {
      *alch_vdw_energy = 0.0;
      *alch_vdw_force = 0.0;
      *alch_vdw_energy_2_Left = 0.0;
      *alch_vdw_energy_2 = 0.0;
    }
  }
  else if(Fep_WCA_dispOn) {
    // separation-shifted dispersion force and energy
    if(B == 0.0) {	// some atoms like H, drude, lone pair are involved
      *alch_vdw_energy = 0.0;
      *alch_vdw_force = 0.0;
      *alch_vdw_energy_2 = 0.0;
      *alch_vdw_energy_3=0.0;
    }
    else {
      const BigReal Emin = B*B/(4.0*A);
      const BigReal Rmin_SQ = powf(2.0*A/B, 1.f/3);
      const BigReal r2_1 = r2;
      const BigReal r2_2 = r2;
      const BigReal r6_1 = r2_1*r2_1*r2_1;
      const BigReal r6_2 = r2_2*r2_2*r2_2;
      *alch_vdw_energy = r2 > Rmin_SQ? \
                         myVdwLambda*(A/(r6_1*r6_1) - B/r6_1): \
                         A/(r6_1*r6_1) - B/r6_1 + (1.-myVdwLambda)* Emin;
      *alch_vdw_force =  r2 > Rmin_SQ? \
                         (12.*(*alch_vdw_energy) + 6.*myVdwLambda*B/r6_1)/r2_1 * switchmul \
                         + (*alch_vdw_energy) * switchmul2:(12.*(*alch_vdw_energy) \
                         + 6.*B/r6_1 - 12.0*(1.-myVdwLambda)* Emin )/r2_1;
      *alch_vdw_energy *= switchmul;
 
			if(!Fep_Wham){ 
        *alch_vdw_energy_2 = r2 > Rmin_SQ? \
                             myVdwLambda2*switchmul*(A/(r6_1*r6_1) - B/r6_1): \
                             A/(r6_1*r6_1) - B/r6_1 + (1.-myVdwLambda2)* Emin;

	//DoubleWide: Calculate the backward alch_vdw_energy here
        *alch_vdw_energy_3 = r2 > Rmin_SQ? \
                             myVdwLambda3*switchmul*(A/(r6_1*r6_1) - B/r6_1): \
                             A/(r6_1*r6_1) - B/r6_1 + (1.-myVdwLambda3)* Emin;

			}
			else{
				*alch_vdw_energy_2 = r2 > Rmin_SQ? \
					                   switchmul*(A/(r6_1*r6_1) - B/r6_1): - Emin;
				*alch_vdw_energy_3 = r2 > Rmin_SQ? \
					                   switchmul*(A/(r6_1*r6_1) - B/r6_1): - Emin;
        *alch_vdw_energy_2 += *alch_vdw_energy;
	*alch_vdw_energy_3 += *alch_vdw_energy;
			}
    }
  } 
  else {
    //myVdwShift already multplied by relevant (1-vdwLambda)
    const BigReal r2_1 = 1./(r2 + myVdwShift);
    const BigReal r2_2 = 1./(r2 + myVdwShift2);
    const BigReal r2_3 = 1./(r2 + myVdwShift3);
    const BigReal r6_1 = r2_1*r2_1*r2_1;
    const BigReal r6_2 = r2_2*r2_2*r2_2;
    const BigReal r6_3 = r2_3*r2_3*r2_3;
    // separation-shifted vdW force and energy
    const BigReal U1 = A*r6_1*r6_1 - B*r6_1; // NB: unscaled, shorthand only!
    const BigReal U2 = A*r6_2*r6_2 - B*r6_2;
    const BigReal U3 = A*r6_3*r6_3 - B*r6_3;
    *alch_vdw_energy = myVdwLambda*switchmul*U1;
    *alch_vdw_energy_2 = myVdwLambda2*switchmul*U2;
    *alch_vdw_energy_3 = myVdwLambda3*switchmul*U3;
    *alch_vdw_force = (myVdwLambda*(switchmul*(12.*U1 + 6.*B*r6_1)*r2_1 \
                                    + switchmul2*U1));
    // BKR - separation-shifted vdW force switching and potential shifting
    if(vdwForceSwitching){ // add potential shifts term
      const BigReal cutoff6 = cutoff2*cutoff2*cutoff2;
      const BigReal cutoff3 = cutoff2*sqrt(cutoff2);

      const BigReal shifted_switchdist2 = switchdist2 + myVdwShift;
      const BigReal shifted_switchdist = sqrt(shifted_switchdist2);
      const BigReal shifted_switchdist6 = \
          shifted_switchdist2*shifted_switchdist2*shifted_switchdist2;
      const BigReal shifted_switchdist3 = shifted_switchdist2*shifted_switchdist;
     
      const BigReal shifted_switchdist2_2 = switchdist2 + myVdwShift2;
      const BigReal shifted_switchdist_2 = sqrt(shifted_switchdist2_2);
      const BigReal shifted_switchdist6_2 = \
          shifted_switchdist2_2*shifted_switchdist2_2*shifted_switchdist2_2;
      const BigReal shifted_switchdist3_2 = shifted_switchdist2_2*shifted_switchdist_2;

      const BigReal shifted_switchdist2_3 = switchdist2 + myVdwShift3;
      const BigReal shifted_switchdist_3 = sqrt(shifted_switchdist2_3);
      const BigReal shifted_switchdist6_3 = \
          shifted_switchdist2_3*shifted_switchdist2_3*shifted_switchdist2_3;
      const BigReal shifted_switchdist3_3 = shifted_switchdist2_3*shifted_switchdist_3;

      const BigReal v_vdwa = -A / (cutoff6*shifted_switchdist6);
      const BigReal v_vdwb = -B / (cutoff3*shifted_switchdist3);
      const BigReal dU = v_vdwa - v_vdwb; //deltaV2 from Steinbach & Brooks
 
      const BigReal v_vdwa2 = -A / (cutoff6*shifted_switchdist6_2);
      const BigReal v_vdwb2 = -B / (cutoff3*shifted_switchdist3_2);
      const BigReal dU2 = v_vdwa2 - v_vdwb2; //deltaV2 from Steinbach & Brooks

      const BigReal v_vdwa3 = -A / (cutoff6*shifted_switchdist6_3);
      const BigReal v_vdwb3 = -B / (cutoff3*shifted_switchdist3_3);
      const BigReal dU3 = v_vdwa3 - v_vdwb3; //deltaV3 from Steinbach & Brooks

      if(r2 > switchdist2) {
	const BigReal k_vdwa = A*cutoff6 / (cutoff6 - shifted_switchdist6);
	const BigReal k_vdwb = B*cutoff3 / (cutoff3 - shifted_switchdist3);
        const BigReal r_1 = sqrt(r2_1);
        const BigReal tmpa = r6_1 - (1./cutoff6);
        const BigReal tmpb = r2_1*r_1 - (1./cutoff3);

        const BigReal k_vdwa2 = A*cutoff6 / (cutoff6 - shifted_switchdist6_2);
        const BigReal k_vdwb2 = B*cutoff3 / (cutoff3 - shifted_switchdist3_2);
        const BigReal r_2 = sqrt(r2_2);
        const BigReal tmpa2 = r6_2 - (1./cutoff6);
        const BigReal tmpb2 = r2_2*r_2 - (1./cutoff3);

        const BigReal k_vdwa3 = A*cutoff6 / (cutoff6 - shifted_switchdist6_3);
        const BigReal k_vdwb3 = B*cutoff3 / (cutoff3 - shifted_switchdist3_3);
        const BigReal r_3 = sqrt(r2_3);
        const BigReal tmpa3 = r6_2 - (1./cutoff6);
        const BigReal tmpb3 = r2_3*r_3 - (1./cutoff3);

        *alch_vdw_energy = (myVdwLambda*(k_vdwa*tmpa*tmpa - k_vdwb*tmpb*tmpb
                                         - (LJcorrection ? dU : 0.)));
        *alch_vdw_energy_2 = (myVdwLambda2*(k_vdwa2*tmpa2*tmpa2 \
                                            - k_vdwb2*tmpb2*tmpb2
                                            - (LJcorrection ? dU2 : 0.)));
        *alch_vdw_energy_2 = (myVdwLambda3*(k_vdwa3*tmpa3*tmpa3 \
                                            - k_vdwb3*tmpb3*tmpb3
                                            - (LJcorrection ? dU2 : 0.)));
        *alch_vdw_force = (myVdwLambda*r2_1*(12.*k_vdwa*tmpa*r6_1
                                             - 6.*k_vdwb*tmpb*r2_1*r_1));
      }else{
        if(!LJcorrection) {
	  *alch_vdw_energy += myVdwLambda*dU;
	  *alch_vdw_energy_2 += myVdwLambda2*dU2;
          *alch_vdw_energy_3 += myVdwLambda3*dU3;
        }
      }
    }
  }
}

#define FEPFLAG
#define CALCENERGY

#define NBTYPE NBPAIR
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#define NBTYPE NBSELF
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#undef CALCENERGY
#undef FEPFLAG


