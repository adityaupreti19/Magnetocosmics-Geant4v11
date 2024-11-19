#ifndef StackingAction_h
#define StackingAction_h 1

#include "G4UserStackingAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class StackingAction : public G4UserStackingAction
{
public:
  StackingAction();
  virtual ~StackingAction();

  inline void SetKillStatus(G4bool value) { fKillSecondary = value; };

  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);

private:

  G4bool              fKillSecondary;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

