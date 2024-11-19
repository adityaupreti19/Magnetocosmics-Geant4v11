#include "MYPropagatorInField.hh"
#include "G4ios.hh"
#include "iomanip"
#include "G4SystemOfUnits.hh"
#include "G4EquationOfMotion.hh" 
#include "G4MonopoleEquation.hh"
#include "MAGCOSMagneticField.hh"
#include"G4TransportationManager.hh"
#include "MAGCOSEquationOfMotion.hh"


const G4double  MYPropagatorInField::fEpsilonMinDefault = 5.0e-7;  
const G4double  MYPropagatorInField::fEpsilonMaxDefault = 0.05;

MYPropagatorInField::MYPropagatorInField(G4Navigator    *theNavigator, 
		                         G4FieldManager *detectorFieldMgr)
  : fDetectorFieldMgr(detectorFieldMgr), 
    fCurrentFieldMgr(detectorFieldMgr), 
    fNavigator(theNavigator),
    End_PointAndTangent(G4ThreeVector(0.,0.,0.),
			G4ThreeVector(0.,0.,0.),0.0,0.0,0.0,0.0,0.0),
    fVerboseLevel(0),
    fEpsilonMin(fEpsilonMinDefault),
    fEpsilonMax(fEpsilonMaxDefault),  
    fmax_loop_count(10000),
    fNoZeroStep(0)
{
     // this->fChordFinder = new G4ChordFinder( (G4MagneticField*)0, 1e-6 );

     fActionThreshold_NoZeroSteps= 2; 
     fSevereActionThreshold_NoZeroSteps= 10; 
     fAbandonThreshold_NoZeroSteps= 50; 
     // fMidPoint_CurveLen_of_LastAttempt= -1;
     fFull_CurveLen_of_LastAttempt= -1; 
     fLast_ProposedStepLength= -1; 

     fLargestAcceptableStep= 1000.0 * meter;
     
     
     
     
     //Motion in G4
     fBSEquationInG4=new BSEquationInG4();
     fBSEquationInG4->SetNavigator(theNavigator);
     
    
     
     
}







///////////////////////////////////////////////////////////////////////////
//
// Compute the next geometric Step


G4double MYPropagatorInField::
         ComputeStep( G4FieldTrack&      pFieldTrack,
	              G4double           CurrentProposedStepLength,
	              G4double&          currentSafety,                // IN/OUT
	              G4VPhysicalVolume* )
{ // Set The Equation of motion for MotionInG4
  
  
  // should be change too complicated
  //------------------
  const G4MagIntegratorStepper* aStepper =fCurrentFieldMgr->GetChordFinder()
	        ->GetIntegrationDriver()->GetStepper();
  		
  G4MagIntegratorStepper* copyStepper = 
  	const_cast<G4MagIntegratorStepper*> (aStepper);   
  G4EquationOfMotion* theEquation  = copyStepper->GetEquationOfMotion(); 
  
  
  //MAGCOSMagneticField* theField = (MAGCOSMagneticField*)
	//  G4TransportationManager::GetTransportationManager()->
         //                     GetFieldManager()->GetDetectorField();


  //G4MonopoleEquation* theEquation = theField->GetEquationOfMotion();
  //MAGCOSEquationOfMotion* theEquation = theField->GetEquationOfMotion();
  //G4cout << "The equation of motion is" << theEquation << G4endl;
  fBSEquationInG4->SetEquationOfMotion(theEquation);
  fBSEquationInG4->SetCrossingDelta(
                     fCurrentFieldMgr->GetDeltaIntersection());	      


  //Compute the integration
  //---------------
  G4double StepLength= 
                 fBSEquationInG4->ComputeTrajectory(pFieldTrack,
		                                CurrentProposedStepLength,
						currentSafety);
     
    
   
 
 
  return StepLength; 




  
}

//////////////////////////////////////////////////////////////
void MYPropagatorInField::SetChargeMomentumMass( 
//			G4double Charge,            // in e+ units
			G4ChargeState particleCharge,
			G4double Momentum,          // in GeV/c 
		        G4double Mass)              // in ? units
   { if  (!fBSEquationInG4->GetEquationOfMotion())
      {const G4MagIntegratorStepper* aStepper =fCurrentFieldMgr->GetChordFinder()
       	      ->GetIntegrationDriver()->GetStepper();
       
       G4MagIntegratorStepper* copyStepper = 
       		const_cast<G4MagIntegratorStepper*> (aStepper);	  
       G4EquationOfMotion* theEquation  = copyStepper->GetEquationOfMotion();

       //MAGCOSMagneticField* theField = (MAGCOSMagneticField*)
	 //      G4TransportationManager::GetTransportationManager()->
           //                   GetFieldManager()->GetDetectorField();

       //G4MonopoleEquation* theEquation  = theField->GetEquationOfMotion(); 
       //MAGCOSEquationOfMotion* theEquation = theField->GetEquationOfMotion();
       fBSEquationInG4->SetEquationOfMotion(theEquation);
      }
     //G4cout << "Now the magnetic charge is set for the Monopole equation" << G4endl; 
     fBSEquationInG4->GetEquationOfMotion()->SetChargeMomentumMass(particleCharge, Momentum, Mass); 
     //G4cout << "Particle charge state is " << particleCharge << G4endl;
     //G4cout << "Mass is " << Mass << G4endl;
   }
