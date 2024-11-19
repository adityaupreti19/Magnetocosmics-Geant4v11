// this file can be modified to provide the magnetic charge, monopole mass 
// through a macro file. the messenger communicates with the macro. 


/// \file exoticphysics/monopole/src/G4MonopolePhysicsMessenger.cc
/// \brief Implementation of the G4MonopolePhysicsMessenger class
//
//
//  12.07.10  S.Burdin (changed the magnetic and electric charge variables 
//            from integer to double)

#include "G4MonopolePhysicsMessenger.hh"

#include "G4MonopolePhysics.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MonopolePhysicsMessenger::G4MonopolePhysicsMessenger(G4MonopolePhysics* p)
  : G4UImessenger(), fPhys(p) 
{
  fPhysicsDir = new G4UIdirectory("/monopole/");
  fPhysicsDir->SetGuidance("histograms control");
   
  fPhysicsCmd = new G4UIcommand("/monopole/setup",this);
  fPhysicsCmd->SetGuidance("Setup monopole");
  //
  G4UIparameter* qmag = new G4UIparameter("qmag",'d',false);
  qmag->SetGuidance("Magnetic charge");
  qmag->SetDefaultValue("1");
  fPhysicsCmd->SetParameter(qmag);

  G4UIparameter* q = new G4UIparameter("qelec",'d',false);
  q->SetGuidance("Electric charge charge");
  q->SetDefaultValue("0");
  fPhysicsCmd->SetParameter(q);
  //    
  G4UIparameter* mass = new G4UIparameter("mass",'d',false);
  mass->SetGuidance("mass");
  mass->SetParameterRange("mass>0.");
  qmag->SetDefaultValue("100");
  fPhysicsCmd->SetParameter(mass);
  //    
  G4UIparameter* unit = new G4UIparameter("unit",'s',false);
  fPhysicsCmd->SetParameter(unit);
  qmag->SetDefaultValue("GeV");
  fPhysicsCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MonopolePhysicsMessenger::~G4MonopolePhysicsMessenger()
{
  delete fPhysicsCmd;
  delete fPhysicsDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MonopolePhysicsMessenger::SetNewValue(G4UIcommand* command, 
                                             G4String newValue)
{ 
  if (command == fPhysicsCmd)
   { G4double q, m; G4double mass; 
     G4String unts;
     std::istringstream is(newValue);
     is >> m >> q >> mass >> unts;
     G4String unit = unts;
     G4double vUnit = G4UIcommand::ValueOf(unit);  
     fPhys->SetMagneticCharge(m);
     fPhys->SetElectricCharge(q);
     fPhys->SetMonopoleMass(mass*vUnit);
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
