//This source file decides what happens in a step
//In the User Stepping Action one can print/save the stepping information 
//or 'kill' a track if it has reached a certain distance, etc. 

#include "MAGCOSSteppingAction.hh"
#include "MAGCOSEventAction.hh"
#include "G4SteppingManager.hh"
#include "MAGCOSMagneticField.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4UserLimits.hh"
#include "SpaceCoordinateConvertor.hh"
#include "DurationManager.hh"
#include "G4EventManager.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4mplIonisationModel.hh"
#include "HistoManager.hh"
#include "G4EmCalculator.hh"


#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RegionStore.hh"


MAGCOSSteppingAction::MAGCOSSteppingAction(HistoManager* histo,  MAGCOSEventAction* event):fHistoManager(histo), fEventAction(event)
{ stop_altitude = 1*m;
  Latitude = 0;
  Longitude = 0;
  Altitude = 0;}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSSteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
  //G4Track* aTrack=aStep->GetTrack();
  //G4int trackID = aTrack->GetTrackID();
  //if (trackID > 1) return;	
  DurationManager* theDurationManager = DurationManager::getInstance();
  if (!theDurationManager->CheckDuration()) {
  	G4cout<< "The event will be aborted"<<std::endl;
	G4EventManager::GetEventManager()->AbortCurrentEvent();
  }
  

  /*if (aStep->GetTrack()->GetCurrentStepNumber() == 1){
	 fStartTime = aStep->GetTrack()->GetGlobalTime();
	 G4cout << "start time is " << fStartTime << G4endl;
  }

  currentTime = aStep->GetTrack()->GetGlobalTime();
  G4cout << "current time is " << G4BestUnit(currentTime, "Time") << G4endl;
  duration = aStep->GetTrack()->GetGlobalTime() - fStartTime;
  G4cout << "elapsed time is " << G4BestUnit(duration, "Time") << G4endl;

  if (duration > 1.0){
	  G4Track* track = aStep->GetTrack();
	  track->SetTrackStatus(fStopAndKill);
	  G4cout << "This track is killed as it took too long to process" << G4endl;
  }*/


  G4double kineticEnergy = aStep->GetTrack()->GetKineticEnergy();
  /*G4double maxStepSize = std::sqrt(kineticEnergy)*1*mm;*/
  G4double stepSize = aStep->GetStepLength();

  //G4cout << "max step size set is " << maxStepSize << G4endl;
  //G4cout << "step size is " << stepSize << G4endl;

  /*if (stepSize > maxStepSize){
	//track->SetTrackStatus(fStopAndKill);
	aStep->GetTrack()->SetStepLength(maxStepSize);
  }*/


  //G4Track* track         = step->GetTrack();
  const G4DynamicParticle* aParticle = aStep->GetTrack()->GetDynamicParticle();
  G4String ParticleName = aParticle->GetParticleDefinition()->GetParticleName();
  
  G4StepPoint* stepPoint = aStep->GetPreStepPoint();
  G4double ParticleEnergy = stepPoint->GetKineticEnergy();
  //G4cout << "Kinetic_energy " << ParticleEnergy << G4endl;
  G4double ParticlePDG = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
  G4double ParticleMass = aStep->GetTrack()->GetParticleDefinition()->GetPDGMass();
  G4double ParticleCharge = aStep->GetTrack()->GetParticleDefinition()->GetPDGCharge();
  
  //G4ThreeVector MomDirection = aStep->GetTrack()->GetMomentumDirection();

  //G4double MomDirX = MomDirection.x();
  //G4double MomDirY = MomDirection.y();
  //G4double MomDirZ = MomDirection.z();

  G4double xMomentumPre = aStep->GetPreStepPoint()->GetMomentum().x();
  G4double yMomentumPre = aStep->GetPreStepPoint()->GetMomentum().y();
  G4double zMomentumPre = aStep->GetPreStepPoint()->GetMomentum().z();

  G4double xMomentumPost = aStep->GetPostStepPoint()->GetMomentum().x();
  G4double yMomentumPost = aStep->GetPostStepPoint()->GetMomentum().y();
  G4double zMomentumPost = aStep->GetPostStepPoint()->GetMomentum().z();

  G4double zTurningPoint = aStep->GetPostStepPoint()->GetPosition().z();
  //G4double kineticEnergy = aStep->GetPostStepPoint()->GetKineticEnergy();
	  
  //if (xMomentumPre * xMomentumPost < 0 || yMomentumPre * yMomentumPost < 0 || zMomentumPre * zMomentumPost < 0) {
  /*if (zMomentumPre * zMomentumPost < 0)	{
	  G4cout << "This is the turning point" << G4endl;
	  G4cout << "z coordinate of turning point is " << zTurningPoint << G4endl;
	  G4cout << "kinetic energy of the particle at turning point " << kineticEnergy << G4endl;
	  G4cout << "Mass of particle : " << ParticleMass << G4endl;
	  
	  fHistoManager->FillNtuple(zTurningPoint, kineticEnergy, ParticleMass);

  }*/

  /*G4cout << "==============================================" << G4endl;
  G4cout << "This particle is a " << ParticleName << G4endl;
  G4cout << "with kinetic energy = " << G4BestUnit(ParticleEnergy, "Energy") << G4endl;
  G4cout << "PDG ID number = " << ParticlePDG << G4endl;
  G4cout << "Mass = " << G4BestUnit(ParticleMass, "Energy")<< G4endl;
  G4cout << "and Electric charge = " << ParticleCharge << G4endl;
  G4cout << "==============================================" << G4endl;*/

  const G4VPhysicalVolume* currentVolume = aStep->GetPostStepPoint()
                ->GetPhysicalVolume();
  G4String volName = currentVolume->GetName();
  //G4cout << "Volume_name_is "<<volName << G4endl;
  G4double volDensity = currentVolume->GetLogicalVolume()->GetMaterial()->GetDensity();
  G4Track* aTrack=aStep->GetTrack();
  G4double StepLength = aTrack->GetStepLength(); 	
  G4int nStep =aTrack->GetCurrentStepNumber();

  if (nStep == 1){
  	G4double iniKE = aStep->GetPreStepPoint()->GetKineticEnergy();
	fEventAction->initKE(iniKE);	
	}


  G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();
  G4double radial_pos = sqrt(pow(position.x(), 2) + pow(position.y(), 2) + pow(position.z(), 2));
  //G4double Altitude, Latitude, Longitude;

  SpaceCoordinateConvertor::getInstance()->ComputeGEOIDCoordinatesFromGEOPosition(position,Altitude,Longitude,Latitude);
  
  //G4cout << "XYZ_position " << position.x()/km << " " << position.y()/km << " " << position.z()/km << G4endl;
  //G4cout << "Radial_position " << radial_pos/km << G4endl;
  //G4cout << "Density_Atmosphere " << volDensity/(g/cm3) << G4endl;

  G4LogicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4StepPoint* postStep = aStep->GetPostStepPoint();

    /*if (volName == "troposphere"){

        //G4LogicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
        G4Region* region = G4RegionStore::GetInstance()->GetRegion(volume->GetRegion()->GetName());
        G4double cutValue = 100*m;

        G4ProductionCuts* cuts = new G4ProductionCuts();
        cuts->SetProductionCut(cutValue, "e-");
        region->SetProductionCuts(cuts);
        //G4cout << "Production cut set to default at change of volumes" << G4endl;
  }*/

  /*if (radial_pos <= 6376.43*km){// && volName == "IceLayer"){
	//if (aStep->GetTrack()->GetKineticEnergy() > 1*MeV){
	G4double KinEnergy = aStep->GetTrack()->GetKineticEnergy();
	//G4cout << "KE of monopole at SLIM height " << G4BestUnit(KinEnergy, "Energy") << G4endl;
	
	//G4cout << "This particle is in the volume" << volName << G4endl;
	G4double stepLengthIce = aTrack->GetStepLength();
	fEventAction->stepICE(stepLengthIce);
	//G4cout << "step length in ICE is : " << G4BestUnit(stepLengthIce, "Length") << G4endl;		
	G4double preKE = preStep->GetKineticEnergy();
	G4double postKE = postStep->GetKineticEnergy();
	G4double EnLoss = preKE - postKE;
	//G4cout << "energy loss in ICE step is : " << G4BestUnit(EnLoss, "Energy") << G4endl;
  }*/


  //if (ParticleName == "monopole" && volName == "IceLayer" && radial_pos < 6369.61*km){

  if (ParticleName == "monopole" && radial_pos < 6371.2*km){
	//G4LogicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();  
	/*G4Region* region = G4RegionStore::GetInstance()->GetRegion(volume->GetRegion()->GetName()); 
 	G4double cutValue = 0.001*m;
	
  	G4ProductionCuts* cuts = new G4ProductionCuts();
  	cuts->SetProductionCut(cutValue, "e-");
  	region->SetProductionCuts(cuts);*/
	//G4double betaMM = preStep->GetBeta();
	//G4cout << "MonopoleBeta_value " << betaMM << G4endl;	

	//G4cout << "Altitude of the end point " << Altitude/km << G4endl;
	//G4cout << "Latitude of the end point " << Latitude/degree << G4endl;
	//G4cout << "Longitude of the end point " << Longitude/degree << G4endl;

	G4cout << "this track reached destination and has been killed" << G4endl;
	aTrack->SetTrackStatus(fStopAndKill);	
  }


  /*if (nStep > 200000){

	G4Region* region = G4RegionStore::GetInstance()->GetRegion(volume->GetRegion()->GetName());
        G4double cutValue = 10*km;

        G4ProductionCuts* cuts = new G4ProductionCuts();
        cuts->SetProductionCut(cutValue, "e-");
        region->SetProductionCuts(cuts);
	aTrack->SetTrackStatus(fStopAndKill); 
        G4cout << "Too many steps, particle killed" << G4endl;

  }*/

  /*if (radial_pos > 6650*km && volName == "Atmosphere"){ aTrack->SetTrackStatus(fStopAndKill);
  			G4cout << "track has crossed atmosphere, its KILLED! " << G4endl;}
  	
  if (radial_pos < 6378*km && volName == "Atmosphere"){ aTrack->SetTrackStatus(fStopAndKill);
                        G4cout << "track is in Earth regime, its KILLED " << G4endl;}*/
  
// -------------------------------------------------------------------------//
// DEDX calc --------------------------------------------------------------//

  G4EmCalculator emCalc;
  //G4mplIonisationModel* calc;
  const G4Material* mat = stepPoint->GetMaterial();
  G4String matName = mat->GetName();
  //const G4DynamicParticle* aParticle = track->GetDynamicParticle();
  G4String particleName = aParticle->GetParticleDefinition()->GetParticleName();
  const G4ParticleDefinition* part = aParticle->GetParticleDefinition();

  //G4double DEDX_perVol;
  //G4cout << "Material is " << matName << " and Particle is " << particleName << " and KE is " << ParticleEnergy << G4endl; 
  //G4double DEDX_perVol = emCalc.GetDEDX(ParticleEnergy, part, mat);
  G4double electronic_DEDX_perVol = emCalc.ComputeElectronicDEDX(kineticEnergy, part, mat);
  //G4double DEDX_perVol = calc->ComputeDEDXPerVolume(mat, part, ParticleEnergy);
  //G4cout << "The_DEDX_perVol_is"  << DEDX_perVol << G4endl;
  //G4cout << "ParticleKE_is " << kineticEnergy << G4endl;
  //G4cout << "TheElectronicDEDXperVol_is " << electronic_DEDX_perVol << G4endl;
  //G4cout << "\n" << G4endl;
  //G4double DEDX_ahlen = calc->ComputeDEDXAhlen(mat, bg2, cutEnergy);

//-l--------------------------------------------------------------------------//


  G4double KE = aStep->GetPostStepPoint()->GetKineticEnergy();
  fEventAction->finalPos(position.x(), position.y(), position.z());
  fEventAction->finalKE(KE);

  //if (nStep > 10 && StepLength == 0.) aTrack->SetTrackStatus(fStopAndKill); 
  /*if (ParticleName == "e-" && StepLength > 100000.*m){
	  G4cout << "This electron has gone too far, time to kill!" <<G4endl;
	  aTrack->SetTrackStatus(fStopAndKill); 
  }*/
  if (currentVolume){		
   	const G4String name = currentVolume->GetName();
    	if (name =="Earth") { 
     		G4ThreeVector position1 = aStep->GetPostStepPoint()->GetPosition();
      		G4double altitude,longitude,latitude;
		G4cout << "This_track_reached_Earth, and killed" << G4endl;
		
    		//  G4cout<<"ok"<<endl;
      		SpaceCoordinateConvertor::getInstance()
              		->ComputeGEOIDCoordinatesFromGEOPosition
	              		(position1,altitude,longitude,latitude);
    		// G4cout<<altitude/km<<std::endl;		      
    		//  G4cout<<"ok1"<<endl;
		
		aTrack->SetTrackStatus(fStopAndKill);
		//G4cout << "The track has been killed" << G4endl;
      		//if (altitude <= stop_altitude )
		                       //aTrack->SetTrackStatus(fStopAndKill);
	   	//		       G4cout << "this track reached Earth, thus killed" << G4endl;	
     	}
    	else {  
     		G4ThreeVector position2 = aStep->GetPostStepPoint()->GetPosition();
      
      		//check if the particle is outside the magnetosphere
      		MAGCOSMagneticField* theField = (MAGCOSMagneticField*)
      		G4TransportationManager::GetTransportationManager()->
                              		GetFieldManager()->GetDetectorField();
      		if (theField->OutsideMagnetosphere(position2)){
           		aTrack->SetTrackStatus(fStopAndKill);
			G4cout << "this track went outside magnetosphere" << G4endl; 
	    		
			return;
		}
	}	
		
	G4UserLimits* theUserLimits = 
                  	currentVolume->GetLogicalVolume()->GetUserLimits();			 
      
      	//stop the particle if the track length is greater than the maximum
      	// allowed track length
      
      	G4double max_track_length =
                      		theUserLimits->GetUserMaxTrackLength(*aTrack);
      
      	if (aTrack->GetTrackLength() >= max_track_length){
           		aTrack->SetTrackStatus(fStopAndKill);
			G4cout << "This track is killed as it went over max track length" << G4endl;
			
	  		return;
	} 
			 
			 
     	//stop the particle if the proper time of the particle is greater than the maximum allowed
     	// proper time of a trajectory
      
      	G4double max_local_time =
                     		theUserLimits->GetUserMaxTime(*aTrack);
      	if (aTrack->GetLocalTime() >= max_local_time){
           	aTrack->SetTrackStatus(fStopAndKill);
	   	G4cout<<"Stop  at local time ="<<aTrack->GetLocalTime()/s<< "second"<<std::endl;
	    	return;
	}   			 					    			 
     	 
  }
}
