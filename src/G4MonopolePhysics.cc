//
/// \file exoticphysics/monopole/src/G4MonopolePhysics.cc
/// \brief Implementation of the G4MonopolePhysics class
//
//  this file defines the monopole magnetic charge, mass, and the physics lists used. 
//--------------------------------------------------------------------------------------
//
// ClassName:   G4MonopolePhysics
//
// Author:      V.Ivanchenko 13.03.2005
//
// Modified:
//
//  12.07.10  S.Burdin (changed the magnetic and electric charge variables 
//            from integer to double)
//----------------------------------------------------------------------------

#include "G4MonopolePhysics.hh"
#include "G4MonopolePhysicsMessenger.hh"

#include "G4Monopole.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "G4StepLimiter.hh"
#include "G4Transportation.hh"
#include "G4MonopoleTransportation.hh"
#include "G4hMultipleScattering.hh"
#include "G4mplIonisation.hh"
#include "G4mplIonisationWithDeltaModel.hh"
#include "G4hhIonisation.hh"
#include "G4hIonisation.hh"

#include "G4PhysicsListHelper.hh"

#include "G4BuilderType.hh"
#include "G4SystemOfUnits.hh"

#include "MYTransportation.hh"
#include "G4ProcessTable.hh"
#include "G4UserSpecialCuts.hh"
#include "G4ParticleWithCuts.hh"
#include "G4PropagatorInField.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4IonConstructor.hh"
#include "G4ios.hh"


#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4IonConstructor.hh"
#include "G4ios.hh"
#include "G4UserSpecialCuts.hh"
//#include "MAGCOSPhysicsMessenger.hh"
#include "G4Transportation.hh"
#include "G4PropagatorInField.hh"
#include "G4MonopoleTransportation.hh"

#include <iomanip>
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MonopolePhysics::G4MonopolePhysics(const G4String& nam)
    :G4VModularPhysicsList(), fMpl(nullptr)

{
  fMagCharge = 1.0;
  fElCharge  = 0.0;
  fMonopoleMass = 1000*GeV;
  //SetPhysicsType(bUnknown);
  fMessenger = new G4MonopolePhysicsMessenger(this);
  defaultCutValue = 3.*mm;

  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV,1*GeV);

  //ConstructBosons();
  //ConstructLeptons();
  //ConstructBaryons();
  //ConstructMonopole();
  G4Alpha::AlphaDefinition();
  G4Deuteron::DeuteronDefinition();
  G4Triton::TritonDefinition();
  G4He3::He3Definition();

  //  generic ion
  G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MonopolePhysics::~G4MonopolePhysics()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MonopolePhysics::ConstructParticle()
{
  if(!fMpl) {
    fMpl = G4Monopole::MonopoleDefinition(fMonopoleMass, fMagCharge, fElCharge);
  } else {
    G4Monopole::Monopole();
  }

  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4Gamma::GammaDefinition();

  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MonopolePhysics::ConstructProcess()
{
  if(verboseLevel > 0) {
    G4cout << "G4MonopolePhysics::ConstructProcess" << G4endl;
  }
  //MYTransportation* theTransportationProcess= new MYTransportation(fMpl);
  auto theParticleIterator=GetParticleIterator(); 
  G4StepLimiter* theStepLimiterProcess = new G4StepLimiter();
  //G4MonopoleTransportation* mplTransport =  new G4MonopoleTransportation(fMpl);  
  while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
	//pmanager->AddProcess(theTransportationProcess);
	//pmanager->SetProcessOrderingToFirst(theTransportationProcess, idxAlongStep);
	//pmanager->SetProcessOrderingToFirst(theTransportationProcess, idxPostStep);
        //pmanager->AddDiscreteProcess(theStepLimiterProcess);
	//AddMyTransportation();
	}

  //G4MonopoleTransportation* mplTransport =  new G4MonopoleTransportation(fMpl);
  G4double magn = fMpl->MagneticCharge();
  G4cout << "magnetic charge is " << magn << G4endl;
  G4double emin = std::min(fMonopoleMass/20000., CLHEP::keV);
  G4double emax = std::max(10.*TeV, fMonopoleMass*100);
  G4int nbin    = G4lrint(10*std::log10(emax/emin));

  if (magn != 0.0){
  //G4MonopoleTransportation* mplTransport =  new G4MonopoleTransportation(fMpl);
  G4cout << "Magnetic monopole processes" << G4endl;
  G4ProcessManager* pmanager = fMpl->GetProcessManager();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  //pmanager->AddProcess(mplTransport);
  pmanager->AddDiscreteProcess(theStepLimiterProcess);
  G4MonopoleTransportation* mplTransport =  new G4MonopoleTransportation(fMpl);
  pmanager->AddProcess(new G4MonopoleTransportation(fMpl), -1, 0, 0);
  //pmanager->AddProcess(mplTransport);
  //pmanager->SetProcessOrdering(mplTransport, idxAlongStep);
  //pmanager->SetProcessOrdering(mplTransport, idxPostStep);  

  G4hIonisation* hhioni = new G4hIonisation();
  hhioni->SetDEDXBinning(nbin);
  hhioni->SetMinKinEnergy(emin);
  hhioni->SetMaxKinEnergy(emax);
  //ph->RegisterProcess(hhioni, fMpl);
  
  G4mplIonisation* mplioni = new G4mplIonisation(magn);
  mplioni->SetDEDXBinning(nbin);
  mplioni->SetMinKinEnergy(emin);
  mplioni->SetMaxKinEnergy(emax);
  ph->RegisterProcess(mplioni, fMpl);
}

  //AddTransportation();
//  AddMyTransportation();


}

//-----------------------------------------------------------------------//


/*void G4MonopolePhysics::AddMyTransportation()
{ MYTransportation* theTransportationProcess= new MYTransportation(fMpl);
  //G4MonopoleTransportation* mplTransport =  new G4MonopoleTransportation(fMpl);
  //theTransportationProcess->GetPropagatorInField()->SetMaxLoopCount(1000);

  // loop over all particles in G4ParticleTable 
  auto theParticleIterator=GetParticleIterator();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        //G4ProcessManager* pmanager_mpl = fMpl->GetProcessManager();
	if (!particle->IsShortLived()) {
                // Add transportation process for all particles other than  "shortlived"
                if ( pmanager == 0) {
                        // Error !! no process manager
                        //G4Exception("G4VUserPhysicsList::AddTransportation : no process manager!");
                }
                else {
                        // add transportation with ordering = ( -1, "first", "first" )
                        
			pmanager ->AddProcess(theTransportationProcess);
                        pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxAlongStep);
                        pmanager ->SetProcessOrderingToFirst(theTransportationProcess, idxPostStep);
                }
        }
	  
	  //if (fMpl){
  	//	G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  	//	G4ProcessManager* pmanager_mpl = fMpl->GetProcessManager();
	//	pmanager_mpl ->AddProcess(theTransportationProcess);
        //        pmanager_mpl ->SetProcessOrderingToFirst(theTransportationProcess, idxAlongStep);
       //         pmanager_mpl ->SetProcessOrderingToFirst(theTransportationProcess, idxPostStep);
//	  }

        else {
                // shortlived particle case
        }
  }
}*/



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void G4MonopolePhysics::SetCuts()
{
  if (verboseLevel >1){
    G4cout << "MAGCOSPhysicsList::SetCuts:";
  }
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  //SetCutValue(1.*m,"gamma");
  SetCutValue(500*m,"e-");
  //SetCutValue(1.*m,"e+");
  //SetCutValue(1*m,"monopole");
  //SetCutValue(1.*m,"proton");
  //SetCutValue(1.*m,"anti_proton");

  //G4Region* region;
  //G4String regName;
  //G4ProductionCuts* cuts;

  //regName = "World";
  //region =  G4RegionStore::GetInstance()->GetRegion(regName);
  //cuts = new G4ProductionCuts;
  //cuts->SetProductionCut(10*m);
  //region->SetProductionCuts(cuts);


  /*G4ProductionCutsTable* theCutsTable = G4ProductionCutsTable::GetProductionCutsTable();
  G4double electronCutValue = 10*m;
  G4double lowEnergy = 1*keV;
  G4double highEnergy = 10*TeV;
  G4ProductionCuts* electronCuts= new G4ProductionCuts();
  electronCuts->SetProductionCut(electronCutValue);
  theCutsTable->SetEnergyRange(lowEnergy, highEnergy);*/

  /*G4double cutValue_low = 10*m;
  G4ProductionCuts* cuts_low = new G4ProductionCuts();
  G4double lowEnergy1 = 1*keV;
  G4double highEnergy1 = 10*TeV;
  //cuts_low->SetEnergyRange(lowEnergy1, highEnergy1);
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowEnergy1, highEnergy1);
  cuts_low->SetProductionCut(cutValue_low, G4ProductionCuts::GetIndex("e-"));
  //cuts_low->SetProductionCut(cutValue_low, G4ProductionCuts::GetIndex("e-"), highEnergy1);*/

  //G4double cutValue_high = 1*mm;
  //G4ProductionCuts* cuts_high = new G4ProductionCuts();
  //G4double lowEnergy2 = 101*MeV;
  //G4double highEnergy2 = 1e10*GeV;
  //cuts_high->SetEnergyRange(lowEnergy2, highEnergy2);
  
  //cuts_high->SetProductionCut(cutValue_high, lowEnergy2, G4Electron::Electron());
  //cuts_high->SetProductionCut(cutValue_high, highEnergy2, G4Electron::Electron());

  //SetCutValueForOthers(defaultCutValue );
}

//---------------------------------------------------------------------------------//


void G4MonopolePhysics::SetMagneticCharge(G4double val)
{
  if ( fMpl ) {
    G4Exception("G4MonopolePhysics", "01", JustWarning, 
                "Cannot set value when monopole is already constructed.");
  } else {
    fMagCharge = val;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MonopolePhysics::SetElectricCharge(G4double val)
{
  if ( fMpl ) {
    G4Exception("G4MonopolePhysics", "01", JustWarning, 
                "Cannot set value when monopole is already constructed.");
  } else {
    fElCharge = val;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MonopolePhysics::SetMonopoleMass(G4double mass)
{
  if ( fMpl ) {
    G4Exception("G4MonopolePhysics", "01", JustWarning, 
                "Cannot set value when monopole is already constructed.");
  } else {
    fMonopoleMass = mass;
  }
}

G4VProcess* G4MonopolePhysics::GetProcess(const G4String& processName) const
{	
  G4ParticleDefinition* particle = G4Monopole::Monopole();
  G4ProcessVector* procList = particle->GetProcessManager()->GetProcessList();
  G4int nbProc = particle->GetProcessManager()->GetProcessListLength();
  G4cout << "Number of Processes for Monopole = " << nbProc << G4endl;
  for (G4int k=0; k<nbProc; k++) {
    G4VProcess* process = (*procList)[k];
    G4cout << "Process names are " << process->GetProcessName() << G4endl;
    if (process->GetProcessName() == processName) return process;
  }
  return nullptr;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

