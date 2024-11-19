#ifndef MAGCOSEVENTACTION_HH
#define MAGCOSEVENTACTION_HH 

#include "G4UserEventAction.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Timer.hh"
#include "G4EventManager.hh"

class G4Event;
class MAGCOSEventMessenger;
class G4Polyline;
class G4Polymarker;
class G4EventManager;

class MAGCOSEventAction : public G4UserEventAction
{
  public:

   MAGCOSEventAction();
   ~MAGCOSEventAction();

  public:
   void BeginOfEventAction(const G4Event*);
   void EndOfEventAction(const G4Event*);
   void DrawTrajectoriesAndFieldLines(G4double zoom, G4double theta, G4double phi); 
   void ResetVectorObjectToBeDrawn();
   void TraceMagnetopauseLine(G4double theta);
  private:
   G4double fx;
   G4double fy;
   G4double fz;
   G4double fKE;
   G4double inKE;
   G4double stepLen;
   MAGCOSEventMessenger* theMessenger;
   G4Colour DrawColour;
   G4bool DrawTrajectory;
   G4bool DrawPoints;
   G4double PointSize;
   G4String DrawingCoordinateSystem;
   std::vector<G4VisAttributes*> TrajectoryVisAttributes;
   std::vector<G4Polyline>  TrajectoryPolyline;
   std::vector<G4Polymarker>  TrajectoryPoints;

   
   
    
    
 //inline methods
  public:
   inline void SetDrawColour(G4Colour aColour){ DrawColour = aColour;}
   // not yet valid could be in the future ????
  /* inline void SetDrawLineWidth(G4double aVal){
                 TrajectoryVisAttributes.SetLineWidth(aVal);}
   // not yet valid could be in the future ????
   inline void SetDrawLineStyle(G4VisAttributes::LineStyle aStyle){
                 TrajectoryVisAttributes.SetLineStyle(aStyle);}	*/	  		 
   
   void finalPos(G4double finalx, G4double finaly, G4double finalz);
   void finalKE(G4double finalKE);
   void initKE(G4double initKE);
   void stepICE(G4double stepLengthIce);

   inline void SetDrawTrajectory(G4bool aBool){DrawTrajectory=aBool;}
   inline void SetDrawPoints(G4bool aBool){DrawPoints=aBool;}
   inline void SetPointSize(G4double aVal){PointSize=aVal;}
   
   inline G4bool GetDrawTrajectory(){return DrawTrajectory;}
   
   
   
   inline void SetDrawingCoordinateSystem(G4String aCoordinateSystem)
                        {DrawingCoordinateSystem = aCoordinateSystem;}
  /* inline void ResetVectorObjectToBeDrawn()
                        {TrajectoryVisAttributes.clear();
			 TrajectoryPolyline.clear();
                         TrajectoryPoints.clear();}*/
			 
   			 
                  			
   			 
  			
  
};

#endif

    
