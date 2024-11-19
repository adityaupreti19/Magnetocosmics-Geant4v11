#include "MAGCOSGeometryConstruction.hh"
//#include "MAGCOSPhysicsList.hh"
#include "G4MonopolePhysics.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "MAGCOSApplicationScenario.hh"
#include "MAGCOSPrimaryGeneratorAction.hh"
#include "MAGCOSEventAction.hh"
#include "MAGCOSTrackingAction.hh"
#include "MAGCOSSteppingAction.hh"
#include "MAGCOSMagneticField.hh"
#include "DurationManager.hh"
#include "G4UImanager.hh"
#include "StackingAction.hh"
#include "G4UIterminal.hh"
//#include "G4UIGAG.hh"
#include "G4UItcsh.hh"
#include "G4UIXm.hh"
//#include "G4UIXaw.hh"
#include "G4RunManager.hh"
//#include "G4UIGainPHPServer.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include <time.h> 

#ifdef G4VIS_USE
#include "MAGCOSVisManager.hh"
#endif

#include "G4ios.hh"
#include "MAGCOSUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "SpaceCoordinateConvertor.hh"
#include "G4ProcessTable.hh"
#include "HistoManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4StepLimiterPhysics.hh"

int main(int argc,char** argv) {

 // G4cout << "geant4 is starting" << G4endl;

  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }
  
  // using a random seed for each simulation. 

  time_t systime = time(NULL);
  long seed = (long) systime;
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(seed);
  G4cout<<"########## The random seed used is ########### => "<<seed<<G4endl;

// Units
// Definition of new units in the unit table should be defined 
// at the beginning before the
// instantiation of the runManager and should be followed by 
// G4UnitDefinition::BuildUnitsTable()  

  new G4UnitDefinition("earth radii","re","Length",re);
  new G4UnitDefinition("earth radii 1","Re","Length",re);
  new G4UnitDefinition("earth radii 2","RE","Length",re);
  new G4UnitDefinition("hour","hour","Time",3600.*s);
  new G4UnitDefinition("minute","minute","Time",60.*s);
  new G4UnitDefinition("day","day","Time",24.*3600.*s);
  new G4UnitDefinition("nanotesla","nT","Magnetic flux density",nT);
  new G4UnitDefinition("gigavolt","GV","Electric potential",1000.*megavolt);
  
  new G4UnitDefinition("/cm3","#/cm3","Density",1./cm3);
  new G4UnitDefinition("/m3","#/m3","Density",1./m3);
  new G4UnitDefinition("/km3","#/km3","Density",1./km3);
  new G4UnitDefinition("/dm3","#/dm3","Density",1.e-3/cm3);
  
  new G4UnitDefinition("km/s","km/s","Speed",km/s);
  new G4UnitDefinition("m/s","m/s","Speed",m/s);
  
  //G4cout << "units defined" << G4endl; 

  G4UnitDefinition::BuildUnitsTable();

  G4double surfTol = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4cout << "surface tolerance is : " << surfTol << G4endl;
  G4double angTol = G4GeometryTolerance::GetInstance()->GetAngularTolerance();
  G4cout << "angular tolerance is : " << angTol << G4endl;
  G4double radTol = G4GeometryTolerance::GetInstance()->GetRadialTolerance();
  G4cout << "radial tolerance is : " << radTol << G4endl;

  G4double newSurfExtent= 1E+12*mm;
  //G4double newAngTol = 1*mm;
  //G4double newRadTol = 1*mm;

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(newSurfExtent);

  G4double radTolNEW = G4GeometryTolerance::GetInstance()->GetRadialTolerance();
  G4cout << "new radial tolerance is : " << radTolNEW << G4endl;


  // Run manager
  G4RunManager * runManager = new G4RunManager;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  //G4cout << "run manager defined" << G4endl;
  // UserInitialization classes
  MAGCOSGeometryConstruction* MAGCOSgeometry = new MAGCOSGeometryConstruction;
  runManager->SetUserInitialization(MAGCOSgeometry);
  //runManager->SetUserInitialization(new MAGCOSPhysicsList);

#ifdef G4VIS_USE

  // Visualization, if you choose to have it!
  G4VisManager* visManager = new MAGCOSVisManager;
  visManager->Initialize();

  
#endif


  G4MonopolePhysics* phys = new G4MonopolePhysics();
  //phys->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(phys);
  
  //Duration manager
  
  DurationManager* theDurationManager;
  theDurationManager = DurationManager::getInstance();
   
  // UserAction classes

  HistoManager* histo = new HistoManager();
  MAGCOSEventAction* event = new MAGCOSEventAction();
  MAGCOSPrimaryGeneratorAction* 
             thePrimaryAction = new MAGCOSPrimaryGeneratorAction();	     
  MAGCOSApplicationScenario* applicationScenario = new MAGCOSApplicationScenario(histo);
  runManager->SetUserAction(applicationScenario); 	     
  runManager->SetUserAction(thePrimaryAction);
  runManager->SetUserAction(new StackingAction);
  runManager->SetUserAction(new MAGCOSTrackingAction);
  runManager->SetUserAction(event);
  MAGCOSSteppingAction* theSteppingAction =new MAGCOSSteppingAction(histo, event);
  MAGCOSgeometry->SetSteppingAction(theSteppingAction);
  runManager->SetUserAction(theSteppingAction);

  //Initialize G4 kernel
  runManager->Initialize(); 
   
  //G4ProcessTable::GetProcessTable()->SetProcessActivation("MYTransportation",false);
  //G4ProcessTable::GetProcessTable()->SetProcessActivation("Transportation",true);
  G4ProcessTable::GetProcessTable()->SetProcessActivation("MonopoleTransportation", true); 
  // User interactions
   
  G4UImanager * UI = G4UImanager::GetUIpointer(); 

  	if (!ui)
        {
                G4String macfileName = argv[1];
                //G4String outfileName = argv[2];

                //G4String command = "/control/execute ";
                //G4String fileName = argv[1];
                //UImanager->ApplyCommand(command+fileName);
                //UImanager->ApplyCommand("/analysis/setFileName "+outfileName);
                UImanager->ApplyCommand("/control/execute "+macfileName);
        }
        else
        {
                UImanager->ApplyCommand("/control/execute init_vis.mac");
                ui->SessionStart();
                delete ui;

        }


#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager; 
  
  return 0;
}



