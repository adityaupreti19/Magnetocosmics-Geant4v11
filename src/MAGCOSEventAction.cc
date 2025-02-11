#include "MAGCOSEventAction.hh"
#include "MAGCOSEventMessenger.hh"
#include "MAGCOSMagneticField.hh"
#include  "MAGCOSUnits.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4strstreambuf.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "SpaceCoordinateConvertor.hh"
#include "MAGCOSSpenvisManager.hh"
#include "DurationManager.hh"
#define pi 3.14159265
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RegionStore.hh"

////////////////////////////////////////////////////////////////////////////////
//
MAGCOSEventAction::MAGCOSEventAction() : G4UserEventAction(), fx(0.), fy(0.), fz(0.), fKE(0.)
{ SetDrawColour(G4Colour(1.,0.,0.));
  DrawingCoordinateSystem="GEO";
  DrawTrajectory=true;
  DrawPoints=false;
  PointSize=1;
  theMessenger = new MAGCOSEventMessenger(this);
  G4cout << "magcos event action " << G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSEventAction::~MAGCOSEventAction()
{}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSEventAction::BeginOfEventAction(const G4Event* evt)
{ DurationManager* theDurationManager = DurationManager::getInstance();
  if (!theDurationManager->CheckDurationAtBeginOfEvent()) {
  	G4RunManager::GetRunManager()->AbortRun(true);
	G4cout<< "The run will be aborted"<<std::endl;
	G4EventManager::GetEventManager()->AbortCurrentEvent();
	
  	
  }

   auto primary = evt->GetPrimaryVertex(0)->GetPrimary(0);
   G4double ini_kenergy = primary->GetKineticEnergy()/GeV;
   G4cout << "initial KE " << ini_kenergy << " GeV" << G4endl;
   //G4double px = primary->GetPx();
   //G4double py = primary->GetPy();
   //G4double pz = primary->GetPz();

   stepLen = 0;

   /*const G4RegionStore* regionStore = G4RegionStore::GetInstance();
   for (G4int i=0; i < regionStore->size(); i++){
	   //G4cout << "changing the production cut at the beginning of event" << G4endl;
	   //G4cout << "# of regions" << regionStore->size() << G4endl;
	   G4Region* region = (*regionStore)[i];
	   G4String regionName = region->GetName();
	   //G4cout << "name of region is " << regionName << G4endl;
	   G4ProductionCuts* cuts = new G4ProductionCuts;
	   cuts->SetProductionCut(100*m, "e-");
	   region->SetProductionCuts(cuts);
   }*/
    
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSEventAction::EndOfEventAction(const G4Event* evt)
{ G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  
  G4cout << "final_KE " << fKE/GeV << " GeV" << G4endl; //G4BestUnit(fKE, "Energy") << G4endl;
  G4cout << "initial_KE " << inKE/GeV << " GeV" << G4endl; 
  G4cout << "Energy lost " << inKE/GeV - fKE/GeV << " GeV" << G4endl;
  G4cout << "final_position_at_the_end_of_event = " << fx << " " << fy << " " << fz << G4endl;
  
  //std::cout << std::setprecision(8);
  
  G4cout << "final radial distance from center of Earth (m) = " << (sqrt(pow(fx, 2) + pow(fy, 2) + pow(fz, 2))/m) << G4endl; 
  G4cout << "distance_travelled_in_ironOre (m)" << " = " << 6.371e06 - (sqrt(pow(fx, 2) + pow(fy, 2) + pow(fz, 2))/m) << G4endl; 

  const G4RegionStore* regionStore = G4RegionStore::GetInstance();
  for (G4int i=0; i < regionStore->size(); i++){
	   //G4cout << "changing the production cut at the end of event" << G4endl;
           //G4cout << "# of regions" << regionStore->size() << G4endl;
           G4Region* region = (*regionStore)[i];
           G4String regionName = region->GetName();
           //G4cout << "name of region is " << regionName << G4endl;
           G4ProductionCuts* cuts = new G4ProductionCuts;
           cuts->SetProductionCut(100*m, "e-");
           region->SetProductionCuts(cuts);
  }
  G4int n_trajectories = 0;
  if(trajectoryContainer){ 
  	n_trajectories = trajectoryContainer->entries(); 
   	//visualisation and registering if needed
   	//---------------------
  	if (DrawTrajectory || DrawPoints){
    		SpaceCoordinateConvertor* theCoordinateConvertor
                        	= SpaceCoordinateConvertor::getInstance();
       		theCoordinateConvertor->SetSystemInAndOut("GEO",DrawingCoordinateSystem);
	  
       		G4int n_point = 
             		(*(evt->GetTrajectoryContainer()))[0]->GetPointEntries();
       
       		G4Polyline pPolyline;
       		G4Polymarker stepPoints;
       		TrajectoryVisAttributes.push_back(new G4VisAttributes(DrawColour));
       		stepPoints.SetMarkerType(G4Polymarker::circles);
       		stepPoints.SetScreenSize(PointSize);
       		stepPoints.SetFillStyle(G4VMarker::filled);
       		stepPoints.SetVisAttributes(TrajectoryVisAttributes.back());
     
     
       		for(G4int i=0; i<n_point; i++){
         		G4ThreeVector pos_geo= 
          			((G4TrajectoryPoint*)
	       				((*(evt->GetTrajectoryContainer()))[0]->GetPoint(i)))
                                                      			->GetPosition();
	  		MAGCOSSpenvisManager::getInstance()->RegisterTrajectoryPoint(pos_geo);
			G4ThreeVector pos= theCoordinateConvertor->Transform(pos_geo,"GEO",DrawingCoordinateSystem);					      
          		//G4cout<<"GEO "<<pos_geo<<std::endl;
			//G4cout<<pos<<std::endl;
			if (DrawTrajectory)  pPolyline.push_back( pos);
	  		if (DrawPoints)  stepPoints.push_back(pos);
	    
         
         	}
	
       
       		pPolyline.SetVisAttributes(TrajectoryVisAttributes.back());
       
       		TrajectoryPolyline.push_back(pPolyline);
       		TrajectoryPoints.push_back(stepPoints);
     	}
 }
}
////////////////////////////////////////////////////////////////////////////////
//
//

void MAGCOSEventAction::finalPos(G4double finalx, G4double finaly, G4double finalz)
{
        fx = finalx;
        fy = finaly;
        fz = finalz;
}

void MAGCOSEventAction::finalKE(G4double finalKE)
{
        fKE = finalKE;
}

void MAGCOSEventAction::initKE(G4double initKE)
{
	inKE = initKE;
} 

/*void MAGCOSEventAction::stepICE(G4double stepLengthIce)
{
	stepLen += stepLengthIce;
}*/


/////////////////////////////////////////////////////////////////////////////

void MAGCOSEventAction::TraceMagnetopauseLine(G4double theta)
{if (DrawTrajectory || DrawPoints) {
	SpaceCoordinateConvertor* theCoordinateConvertor
                        = SpaceCoordinateConvertor::getInstance();
	theCoordinateConvertor->SetSystemInAndOut("GEO",DrawingCoordinateSystem);
	const MAGCOSMagneticField* theField =
                    reinterpret_cast<const MAGCOSMagneticField*>
	                 (G4TransportationManager::GetTransportationManager()->
		                             GetFieldManager()->GetDetectorField());
 	G4ThreeVector soff_pos = theField->FindStandOffPosition(0.001*Re);
	G4double xmin = soff_pos.x();
	G4double xmax = -49.*Re;
	
	
	G4double step =(xmax -xmin) /250.;
	G4Polyline pPolyline1,pPolyline2;
        G4Polymarker stepPoints;
        TrajectoryVisAttributes.push_back(new G4VisAttributes(DrawColour));
        stepPoints.SetMarkerType(G4Polymarker::circles);
        stepPoints.SetScreenSize(PointSize);
        stepPoints.SetFillStyle(G4VMarker::filled);
        stepPoints.SetVisAttributes(TrajectoryVisAttributes.back());
	G4cout<<soff_pos/re<<std::endl;
	G4ThreeVector pos= theCoordinateConvertor->Transform(soff_pos);
	G4cout<<pos/re<<std::endl;
	if (DrawTrajectory){
		pPolyline1.push_back( pos);
		pPolyline2.push_back( pos);
		}
	
	if (DrawPoints){
		stepPoints.push_back(pos);
		}
     
     
        for(G4int i=1; i<251; i++){
		G4double x= xmin + i* step;
		G4ThreeVector pos_gsm=theField
		               ->FindPositionOnMagnetopause(theta,x,.001*Re);
		G4cout<<pos_gsm/re<<std::endl;
		if (pos_gsm.mag() > 0.){
		 	theCoordinateConvertor->SetSystemInAndOut("GSM",DrawingCoordinateSystem);     
			pos= theCoordinateConvertor->Transform(pos_gsm);
			G4cout<<pos/re<<std::endl;
			if (DrawTrajectory){
				pPolyline1.push_back( pos);
			}
	
			if (DrawPoints){
				stepPoints.push_back(pos);
			}
		}
		
			       
		pos_gsm=theField->FindPositionOnMagnetopause(theta +180.*degree,
			                                     x,
			                                     .001*Re);
		if (pos_gsm.mag() > 0.){
		        theCoordinateConvertor->SetSystemInAndOut("GSM",DrawingCoordinateSystem);  
			pos= theCoordinateConvertor->Transform(pos_gsm);
			if (DrawTrajectory){
				pPolyline2.push_back( pos);
			}
	
			if (DrawPoints){
				stepPoints.push_back(pos);
			}
		} 
	}
	pPolyline1.SetVisAttributes(TrajectoryVisAttributes.back());
	pPolyline2.SetVisAttributes(TrajectoryVisAttributes.back());
        TrajectoryPolyline.push_back(pPolyline1);
	TrajectoryPolyline.push_back(pPolyline2);
        TrajectoryPoints.push_back(stepPoints);
	}
 else {G4cout<<"You should first switch on DrawPoint and/or DrawTrajectory"<<std::endl; 
 }		  		 				       
}
////////////////////////////////////////////////////////////////////////////////
// 
void MAGCOSEventAction::DrawTrajectoriesAndFieldLines
                             (G4double , G4double , G4double )
{ unsigned int nline = TrajectoryPolyline.size();
  unsigned int npoints =TrajectoryPoints.size();
     
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (!pVVisManager){ 
   	G4cout<<"For visualising trajectories and magnetic";
       	G4cout<<" field lines you should select a visualisation driver"<<G4endl; 
       	return;
  }
  if (nline ==0){
   	G4cout<<"There is nothing to visualise"<<G4endl;
       	return;
  }
  ((G4VisManager*)pVVisManager)->GetCurrentSceneHandler()-> ClearStore ();
     
  G4Scene* aScene = ((G4VisManager*)pVVisManager)
                            ->GetCurrentSceneHandler()->GetScene();
  if (!aScene){
  	aScene = new G4Scene();
	((G4VisManager*)pVVisManager)
                            ->GetCurrentSceneHandler()->SetScene(aScene);
  } 			    
  
 // aScene->Clear();
        
  //((G4VisManager*)pVVisManager)->GetCurrentScene()-> Clear(); 
  /* ((G4VisManager*)pVVisManager)->GetCurrentViewer()-> SetView ();
  ((G4VisManager*)pVVisManager)->GetCurrentViewer()-> ClearView ();*/
   
  G4UImanager::GetUIpointer () -> ApplyCommand ("/vis/drawVolume");
  G4cout<<"Nlines "<<nline<<std::endl;   
  for (unsigned int i=0;i<nline;i++){
  	pVVisManager->Draw(TrajectoryPolyline[i]); 			     
  }				     
     
  for (unsigned i=0;i<npoints;i++)
        	/*((G4VisManager*)pVVisManager)->GetCurrentSceneHandler()
	                             	     ->AddPrimitive(TrajectoryPoints[i]); */
               pVVisManager->Draw(TrajectoryPoints[i]); 
      
    	/* std::strstream astream;
     	astream<< "/vis/viewer/set/viewpointThetaPhi"<<'\t'<<theta<<'\t'<<phi<<G4endl;
     	G4String cmd1,cmd2,cmd3;
     	astream>>cmd1>>cmd2>>cmd3;			
     	G4String command=cmd1+" "+cmd2+" "+cmd3+" deg";  
     	G4UImanager::GetUIpointer () -> ApplyCommand (command);
     
     
     	G4UImanager::GetUIpointer () -> ApplyCommand 
                               ("/vis/viewer/set/style surface2");
     
     	G4UImanager::GetUIpointer () -> 
                             ApplyCommand ("/vis/viewer/flush ");*/
    	//((G4VisManager*)pVVisManager)->GetCurrentViewer()-> DrawView ();   
   	//  ((G4VisManager*)pVVisManager)->GetCurrentViewer()->ShowView();
      
}

////////////////////////////////////////////////////////////////////////////////
// 
void MAGCOSEventAction::ResetVectorObjectToBeDrawn()
{ TrajectoryVisAttributes.clear();
  TrajectoryPolyline.clear();
  TrajectoryPoints.clear();
}
