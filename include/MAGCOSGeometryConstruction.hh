#ifndef MAGCOSGeometryConstruction_h
#define MAGCOSGeometryConstruction_h 1
// DESCRIPTION
// -----------
//
// This class defines the geometry  for MAGNETOCOSMICS. 
// The geometry id defined by a vacuum sphere representing the Earth and
// placed at the centre of a vacuum world box.
//  
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
//
// void construct()
//     Method where the Geometry is defined. This method is called by the 
//     Geant4 manager.
// 
//void SetMaxStepLength(G4double)
//     Set the maximum tracking step length in space outside the Earth.
//
//void SetMaxTimeOfFlight(G4double)
//     Set the maximum time of flight of a particle when computing its
//     trajectory.
//
//void SetMaxTrackLength()
//     Set the maximum path length for a particle trajectory.
//


#include "globals.hh"
#include "G4Cache.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4Colour.hh"
#include "G4ProductionCuts.hh"


class G4Box;
class G4Sphere;
class G4Orb;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4Element;
class MAGCOSGeometryMessenger;
class G4UserLimits;
class MAGCOSMagneticField;
class MAGCOSSteppingAction;
class G4Region;

//class G4MonopoleFieldSetup;


class MAGCOSGeometryConstruction : public G4VUserDetectorConstruction
{
  public:
     MAGCOSGeometryConstruction();
    ~MAGCOSGeometryConstruction();

  public:
  
     G4VPhysicalVolume* Construct();
     
  public:
  
    void SetMaxStepLength(G4double maxStep);
//    inline G4double GetMaxStepLength(){return MaxStepLength;}
    inline    G4LogicalVolume* GetWorldLogicalVolume() 
                            {return logicWorld;}
    void SetMaxTimeOfFlight(G4double);
    void SetMaxTrackLength(G4double);
    void SetStopAltitude(G4double stop_altitude);
    void SetRemoveEarth(G4bool aBool);
    inline void SetSteppingAction (MAGCOSSteppingAction* aSteppingAction)
     				{theSteppingAction = aSteppingAction;}

    void SetCutsForVolumes(G4ProductionCuts* productionCuts, const G4String& particleName, const G4String& volumeName);
    //virtual void ConstructSDandField();

    MAGCOSMagneticField* GetMonopoleFieldSetup() const { return fMonFieldSetup.Get(); }
 

  private:     
  
  //    physiworld
  //---------------------     
  G4Box*             solidWorld;    //pointer to the solid enveloppe 
  G4LogicalVolume*   logicWorld;    //pointer to the logical enveloppe
  G4VPhysicalVolume* physiWorld;  
  G4UserLimits* theWorldUserLimits;
  G4double maxStepSize;
     
  // Stepping action pointer needed when changing the top of the atmosphere
  MAGCOSSteppingAction* theSteppingAction;
     
  // Earth
  //---------
     
  //G4Sphere*          solidEarth;
 
  G4Material* airMat[150]; 
  G4Orb* atmos[150];
  G4VPhysicalVolume* physAtmos[150];
  G4LogicalVolume* logAtmos[150];
  G4UserLimits* atmosUserLimits[150];  

  G4Material* upperAirMat[100];
  G4Orb* upperAtmos[100];
  G4VPhysicalVolume* upperPhysAtmos[100];
  G4LogicalVolume* upperLogAtmos[100];
  G4UserLimits* upperAtmosUserLimits[100];

  G4Orb*          solidEarth;
  G4LogicalVolume*   logicEarth;  
  G4VPhysicalVolume* physiEarth;
  G4VPhysicalVolume* fPhysCubeSats;
  G4UserLimits* theEarthUserLimits; 
  G4double TopOfAtmosphere;
  
  G4bool RemoveEarth;
     
     
   
  MAGCOSMagneticField* theMagnetosphereMagneticField;   // pointer to the magnetic field 
  
  G4Cache<MAGCOSMagneticField*> fMonFieldSetup;  
  MAGCOSGeometryMessenger*  detectorMessenger;
};

#endif

