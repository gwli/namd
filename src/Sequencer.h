/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef SEQUENCER_H
#define SEQUENCER_H

#include "converse.h"
#include "Priorities.h"
#include "PatchTypes.h"

class HomePatch;
class SimParameters;
class SubmitReduction;
class CollectionMgr;
class ControllerBroadcasts;
class LdbCoordinator;
class Random;

#define SEQUENCER_SOA

class Sequencer
{
    friend class HomePatch;
public:
    Sequencer(HomePatch *p);
    virtual ~Sequencer(void);
    void run(void);             // spawn thread, etc.
    void awaken(void) {
      CthAwakenPrio(thread, CK_QUEUEING_IFIFO, PRIORITY_SIZE, &priority);
    }
    void suspend(void);

protected:
    virtual void algorithm(void);	// subclasses redefine this method

#ifdef SEQUENCER_SOA
    void integrate_SOA(int); // Verlet integrator using SOA data structures
    void rattle1_SOA(BigReal,int);
    void addForceToMomentum_SOA(
        const double scaling,
        double       dt_normal,               // timestep Results::normal = 0
        double       dt_nbond,                // timestep Results::nbond  = 1
        double       dt_slow,                 // timestep Results::slow   = 2
        const float  * __restrict recipMass,
        const double * __restrict f_normal_x, // force    Results::normal = 0
        const double * __restrict f_normal_y,
        const double * __restrict f_normal_z,
        const double * __restrict f_nbond_x,  // force    Results::nbond  = 1
        const double * __restrict f_nbond_y,
        const double * __restrict f_nbond_z,
        const double * __restrict f_slow_x,   // force    Results::slow   = 2
        const double * __restrict f_slow_y,
        const double * __restrict f_slow_z,
        double       * __restrict vel_x,
        double       * __restrict vel_y,
        double       * __restrict vel_z,
        int numAtoms,
        int maxForceNumber
    );
    void addVelocityToPosition_SOA(
        const double dt,  ///< scaled timestep
        const double * __restrict vel_x,
        const double * __restrict vel_y,
        const double * __restrict vel_z,
        double *       __restrict pos_x,
        double *       __restrict pos_y,
        double *       __restrict pos_z,
        int numAtoms      ///< number of atoms
    );
    void submitHalfstep_SOA(
        const int    * __restrict hydrogenGroupSize,
        const float  * __restrict mass,
        const double * __restrict vel_x,
        const double * __restrict vel_y,
        const double * __restrict vel_z,
        int numAtoms
    );
    void submitReductions_SOA(
        const int    * __restrict hydrogenGroupSize,
        const float  * __restrict mass,
        const double * __restrict pos_x,
        const double * __restrict pos_y,
        const double * __restrict pos_z,
        const double * __restrict vel_x,
        const double * __restrict vel_y,
        const double * __restrict vel_z,
        const double * __restrict f_normal_x,
        const double * __restrict f_normal_y,
        const double * __restrict f_normal_z,
        const double * __restrict f_nbond_x,
        const double * __restrict f_nbond_y,
        const double * __restrict f_nbond_z,
        const double * __restrict f_slow_x,
        const double * __restrict f_slow_y,
        const double * __restrict f_slow_z,
        int numAtoms
    );
    void submitCollections_SOA(int step, int zeroVel = 0);
    void maximumMove_SOA(
        const double dt,  ///< scaled timestep
        const double maxvel2, ///< square of bound on velocity
        const double * __restrict vel_x,
        const double * __restrict vel_y,
        const double * __restrict vel_z,
        int numAtoms      ///< number of atoms
    );
    void langevinVelocitiesBBK1_SOA(
        BigReal timestep,
        const float * __restrict langevinParam,
        double      * __restrict vel_x,
        double      * __restrict vel_y,
        double      * __restrict vel_z,
        int numAtoms
    );
    void langevinVelocitiesBBK2_SOA(
        BigReal timestep,
        const float * __restrict langevinParam,
        const float * __restrict langScalVelBBK2,
        const float * __restrict langScalRandBBK2,
        float       * __restrict gaussrand_x,
        float       * __restrict gaussrand_y,
        float       * __restrict gaussrand_z,
        double      * __restrict vel_x,
        double      * __restrict vel_y,
        double      * __restrict vel_z,
        int numAtoms
    );
    void langevinPiston_SOA(
        const int    * __restrict hydrogenGroupSize,
        const float  * __restrict mass,
        double       * __restrict pos_x,
        double       * __restrict pos_y,
        double       * __restrict pos_z,
        double       * __restrict vel_x,
        double       * __restrict vel_y,
        double       * __restrict vel_z,
        int numAtoms,
        int step
    );
    void runComputeObjects_SOA(int migration, int pairlists);
#endif

    void integrate(int); // Verlet integrator
    void minimize(); // CG minimizer
      SubmitReduction *min_reduction;

    void runComputeObjects(int migration = 1, int pairlists = 0, int pressureStep = 0);
    int pairlistsAreValid;
    int pairlistsAge;

    void calcFixVirial(Tensor& fixVirialNormal, Tensor& fixVirialNbond, Tensor& fixVirialSlow,
      Vector& fixForceNormal, Vector& fixForceNbond, Vector& fixForceSlow);

    void submitReductions(int);
    void submitHalfstep(int);
    void submitMinimizeReductions(int, BigReal fmax2);
    void submitCollections(int step, int zeroVel = 0);

    void submitMomentum(int step);
    void correctMomentum(int step, BigReal drifttime);

    void saveForce(const int ftag = Results::normal);
    void addForceToMomentum(BigReal, const int ftag = Results::normal, const int useSaved = 0);
    void addForceToMomentum3(const BigReal timestep1, const int ftag1, const int useSaved1,
        const BigReal timestep2, const int ftag2, const int useSaved2,
        const BigReal timestep3, const int ftag3, const int useSaved3);
    void addVelocityToPosition(BigReal);
    
    void addRotDragToPosition(BigReal);
    void addMovDragToPosition(BigReal);

    void minimizeMoveDownhill(BigReal fmax2);
    void newMinimizeDirection(BigReal);
    void newMinimizePosition(BigReal);
    void quenchVelocities();

    void hardWallDrude(BigReal,int);

    void rattle1(BigReal,int);
    // void rattle2(BigReal,int);

    void maximumMove(BigReal);
    void minimizationQuenchVelocity(void);

    void reloadCharges();

    BigReal adaptTempT;         // adaptive tempering temperature
    void adaptTempUpdate(int); // adaptive tempering temperature update

    void rescaleVelocities(int);
    void rescaleaccelMD(int, int, int); // for accelMD
    int rescaleVelocities_numTemps;
    void reassignVelocities(BigReal,int);
    void reinitVelocities(void);
    void rescaleVelocitiesByFactor(BigReal);
    void tcoupleVelocities(BigReal,int);
    void berendsenPressure(int);
      int berendsenPressure_count;
      int checkpoint_berendsenPressure_count;
    void langevinPiston(int);
      int slowFreq;
    void newtonianVelocities(BigReal, const BigReal, const BigReal, 
                             const BigReal, const int, const int, const int);
    void langevinVelocities(BigReal);
    void langevinVelocitiesBBK1(BigReal);
    void langevinVelocitiesBBK2(BigReal);
    // Multigrator
    void scalePositionsVelocities(const Tensor& posScale, const Tensor& velScale);
    void multigratorPressure(int step, int callNumber);
    void scaleVelocities(const BigReal velScale);
    BigReal calcKineticEnergy();
    void multigratorTemperature(int step, int callNumber);
    SubmitReduction *multigratorReduction;
    int doKineticEnergy;
    int doMomenta;
    // End of Multigrator
    
    void cycleBarrier(int,int);
	void traceBarrier(int);
#ifdef MEASURE_NAMD_WITH_PAPI
	void papiMeasureBarrier(int);
#endif
    void terminate(void);

    Random *random;
    SimParameters *const simParams;	// for convenience
    HomePatch *const patch;		// access methods in patch
    SubmitReduction *reduction;
    SubmitReduction *pressureProfileReduction;

    CollectionMgr *const collection;
    ControllerBroadcasts * broadcast;

    int ldbSteps;
    void rebalanceLoad(int timestep);


private:
    CthThread thread;
    unsigned int priority;
    static void threadRun(Sequencer*);

    LdbCoordinator *ldbCoordinator;
};

#endif
