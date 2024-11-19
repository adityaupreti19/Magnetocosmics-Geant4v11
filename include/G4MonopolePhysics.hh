//
/// \file exoticphysics/monopole/include/G4MonopolePhysics.hh
/// \brief Definition of the G4MonopolePhysics class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4MonopolePhysics_h
#define G4MonopolePhysics_h 1

//#include "G4VPhysicsConstructor.hh"

#include "G4VUserPhysicsList.hh"
#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "G4ProcessManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MonopolePhysicsMessenger;
class G4Monopole;
class G4VProcess;

//class G4MonopolePhysics : public G4VPhysicsConstructor
class G4MonopolePhysics : public G4VModularPhysicsList
{
public:

  G4MonopolePhysics(const G4String& nam = "Monopole Physics");

  ~G4MonopolePhysics();

  // This method is dummy for physics
  virtual void ConstructParticle();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();
  //virtual void ConstructSM();
  //virtual void ConstructMonopole();
  virtual void SetCuts();

  void SetMagneticCharge(G4double);
  void SetElectricCharge(G4double);
  void SetMonopoleMass(G4double);
  G4VProcess*  GetProcess(const G4String&) const;

private:

  // hide assignment operator
  G4MonopolePhysics & operator=(const G4MonopolePhysics &right);
  G4MonopolePhysics(const G4MonopolePhysics&);

  void AddMyTransportation();
  //void MonopoleTransport();

  G4double    fMagCharge;
  G4double    fElCharge;
  G4double    fMonopoleMass;

  G4MonopolePhysicsMessenger*  fMessenger;
  G4Monopole* fMpl;

protected:
    //virtual void ConstructBosons();
    //virtual void ConstructLeptons();
    //virtual void ConstructBaryons();
    //virtual void ConstructMonopole();
  //virtual void SetCuts();

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

