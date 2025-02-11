//
/// \file exoticphysics/monopole/include/G4MonopoleTransportation.hh
/// \brief Definition of the G4MonopoleTransportation class
//
// ------------------------------------------------------------
//        GEANT 4  include file implementation
// ------------------------------------------------------------
//
// Class description:
//
// G4MonopoleTransportation is a process responsible for the transportation of 
// magnetic monopoles, i.e. the geometrical propagation encountering the 
// geometrical sub-volumes of the detectors.
// It is also tasked with part of updating the "safety".

// =======================================================================
// Created:  3 May 2010, J. Apostolakis, B. Bozsogi
// =======================================================================

#ifndef G4MonopoleTransportation_hh
#define G4MonopoleTransportation_hh 1

#include "G4VProcess.hh"
#include "G4FieldManager.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "MYPropagatorInField.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleChangeForTransport.hh"
#include "MAGCOSMagneticField.hh"

class G4SafetyHelper; 
class G4Monopole;

class G4MonopoleTransportation : public G4VProcess 
{
  // Concrete class that does the geometrical transport 

  public:  // with description

     G4MonopoleTransportation(const G4Monopole* p, G4int verbosityLevel= 1);
     ~G4MonopoleTransportation(); 

     virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                                   G4double  previousStepSize,
                                   G4double  currentMinimumStep, 
                                   G4double& currentSafety,
                                   G4GPILSelection* selection
                            );

     virtual G4VParticleChange* AlongStepDoIt(
                             const G4Track& track,
                             const G4Step& stepData
                            );

     virtual G4VParticleChange* PostStepDoIt(
                             const G4Track& track,
                             const G4Step&  stepData
                            );
       // Responsible for the relocation.

     virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& ,
                             G4double   previousStepSize,
                             G4ForceCondition* pForceCond
                            );
       // Forces the PostStepDoIt action to be called, 
       // but does not limit the step.

     //G4PropagatorInField* GetPropagatorInField();
     //void SetPropagatorInField( G4PropagatorInField* pFieldPropagator);
       // Access/set the assistant class that Propagate in a Field.

     MYPropagatorInField* GetPropagatorInField();
     void SetPropagatorInField( MYPropagatorInField* pFieldPropagator);


     inline G4double GetThresholdWarningEnergy() const; 
     inline G4double GetThresholdImportantEnergy() const; 
     inline G4int GetThresholdTrials() const; 

     inline void SetThresholdWarningEnergy( G4double newEnWarn ); 
     inline void SetThresholdImportantEnergy( G4double newEnImp ); 
     inline void SetThresholdTrials(G4int newMaxTrials ); 

     // Get/Set parameters for killing loopers: 
     //   Above 'important' energy a 'looping' particle in field will 
     //   *NOT* be abandoned, except after fThresholdTrials attempts.
     // Below Warning energy, no verbosity for looping particles is issued

     inline G4double GetMaxEnergyKilled() const; 
     inline G4double GetSumEnergyKilled() const;
     inline void ResetKilledStatistics( G4int report = 1);      
     // Statistics for tracks killed (currently due to looping in field)

     inline void EnableShortStepOptimisation(G4bool optimise=true); 
     // Whether short steps < safety will avoid to call Navigator (if field=0)

  public:  // without description

     virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& ,
                             G4ForceCondition* 
                            ) { return -1.0; };
     // No operation in  AtRestDoIt.

     virtual G4VParticleChange* AtRestDoIt(
                             const G4Track& ,
                             const G4Step&
                            ) {return 0;};
     // No operation in  AtRestDoIt.

  //G4double GetZmagFieldValue() const { return fMagSetup-> GetZmagFieldValue(); }

  virtual void StartTracking(G4Track* aTrack);
       // Reset state for new (potentially resumed) track 

  protected:

     G4bool               DoesGlobalFieldExist();
       // Checks whether a field exists for the "global" field manager.

  private:

     const G4Monopole* fParticleDef;

     MAGCOSMagneticField*  fMagSetup;
    
     G4Navigator*         fLinearNavigator;
     //G4PropagatorInField* fFieldPropagator;
     MYPropagatorInField* fFieldPropagator;

  // The Propagators used to transport the particle

     G4ThreeVector        fTransportEndPosition;
     G4ThreeVector        fTransportEndMomentumDir;
     G4double             fTransportEndKineticEnergy;
     G4ThreeVector        fTransportEndSpin;
     G4bool               fMomentumChanged;
     //  G4bool               fEnergyChanged;
     G4bool               fEndGlobalTimeComputed; 
     G4double             fCandidateEndGlobalTime;
  // The particle's state after this Step, Store for DoIt

     G4bool               fParticleIsLooping;

     G4TouchableHandle    fCurrentTouchableHandle;

     G4bool fGeometryLimitedStep;
       // Flag to determine whether a boundary was reached.

     G4ThreeVector  fPreviousSftOrigin;
     G4double       fPreviousSafety; 
       // Remember last safety origin & value.

     G4ParticleChangeForTransport fParticleChange;
       // New ParticleChange

     G4double endpointDistance;

  // Thresholds for looping particles: 
  // 
     G4double fThreshold_Warning_Energy;     //  Warn above this energy
     G4double fThreshold_Important_Energy;   //  Hesitate above this
     G4int    fThresholdTrials;              //    for this no of trials
       // Above 'important' energy a 'looping' particle in field will 
       //   *NOT* be abandoned, except after fThresholdTrials attempts.
     // G4double fUnimportant_Energy;
       //  Below this energy, no verbosity for looping particles is issued

  // Counter for steps in which particle reports 'looping',
  //   if it is above 'Important' Energy 
     G4int    fNoLooperTrials; 
  // Statistics for tracks abandoned
     G4double fSumEnergyKilled;
     G4double fMaxEnergyKilled;

  // Whether to avoid calling G4Navigator for short step ( < safety)
  //   If using it, the safety estimate for endpoint will likely be smaller.
     G4bool   fShortStepOptimisation; 

     G4SafetyHelper* fpSafetyHelper;  // To pass it the safety value obtained
     G4int noCalls;

};

#include "G4MonopoleTransportation.icc"

#endif  
