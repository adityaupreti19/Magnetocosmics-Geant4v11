#ifndef MAGCOSSteppingAction_h
#define MAGCOSSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class MAGCOSEventAction;
class HistoManager;

class MAGCOSSteppingAction : public G4UserSteppingAction
{
  public:
    MAGCOSSteppingAction(HistoManager*, MAGCOSEventAction*);
    virtual ~MAGCOSSteppingAction(){};
    virtual void UserSteppingAction(const G4Step*);
    inline void SetStopAltitude(G4double alt){stop_altitude = alt;}
    
  private:
    HistoManager* fHistoManager;
    MAGCOSEventAction* fEventAction;
    G4double stop_altitude;
    G4double fStartTime;
    G4double currentTime;
    G4double duration;
    G4double Latitude;
    G4double Longitude;
    G4double Altitude;
};

#endif
