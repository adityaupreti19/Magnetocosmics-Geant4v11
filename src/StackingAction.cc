//here the stacking action is used to track only the primary particles - monopoles
// and kill the secondaries which might be produced, thus, increasing simulation speed. 

#include "StackingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction()
  : G4UserStackingAction(), fKillSecondary(true)
{
  //fStackMessenger = new StackingMessenger(this);
}

StackingAction::~StackingAction()
{
  //delete fStackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  //stack or delete secondaries
  G4ClassificationOfNewTrack status = fUrgent;

  //keep primary particle
  if(aTrack->GetParentID() == 0) { return status; }
  
  //if(fKillSecondary) { status = fKill; }
  //return status;
  else return fKill;
}

