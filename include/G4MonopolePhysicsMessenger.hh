//
/// \file exoticphysics/monopole/include/G4MonopolePhysicsMessenger.hh
/// \brief Definition of the G4MonopolePhysicsMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4MonopolePhysicsMessenger_h
#define G4MonopolePhysicsMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4MonopolePhysics;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MonopolePhysicsMessenger: public G4UImessenger
{
public:

  G4MonopolePhysicsMessenger(G4MonopolePhysics*);
  ~G4MonopolePhysicsMessenger();
    
  virtual void SetNewValue(G4UIcommand*, G4String);
    
private:

  G4MonopolePhysics*         fPhys;
    
  G4UIdirectory*             fPhysicsDir;    
  G4UIcommand*               fPhysicsCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
