/// \brief Implementation of the G4Monopole class
//
//
//---------------------------------------------------------------------------
//
// ClassName:   G4Monopole
//  
// Description: 
//
// Authors:   21.03.05  V.Ivanchenko 
//
// Modified:
//
//  12.07.10  S.Burdin (changed the magnetic and electric charge variables from integer to double)
//----------------------------------------------------------------------------
//

// this file sets the properties of a monopoles, and creates the particle in Geant4
//


#include "G4Monopole.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// ######################################################################
// ###                        Monopole                                ###
// ######################################################################

G4Monopole* G4Monopole::theMonopole = 0;
G4double    G4Monopole::magCharge = 0.0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Monopole::G4Monopole(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable)
 : G4ParticleDefinition( aName, mass, width, charge, iSpin, iParity,
           iConjugation, iIsospin, iIsospin3, gParity, pType,
           lepton, baryon, encoding, stable, lifetime, decaytable )
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Monopole::~G4Monopole()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//     
//    Arguments for constructor are as follows
//               name             mass          width         charge
//             2*spin           parity  C-conjugation
//          2*Isospin       2*Isospin3       G-parity
//               type    lepton number  baryon number   PDG encoding
//             stable         lifetime    decay table 
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Monopole* G4Monopole::MonopoleDefinition(G4double mass, G4double mCharge, 
                                           G4double eCharge)
{    
  if(!theMonopole) {
    magCharge = eplus * mCharge / fine_structure_const * 0.5;
    theMonopole = new G4Monopole(
       "monopole",         mass,       0.0*MeV,       eplus*eCharge, 
                0,               0,             0,          
                0,               0,             0,             
          "boson",               0,             0,           0,
             true,            -1.0,             0);
    
    
    G4cout << "Monopole is created: m(GeV)= " << theMonopole->GetPDGMass()/GeV 
           << " Qel= " << theMonopole->GetPDGCharge()/eplus
           << " Qmag= " << magCharge/eplus
           << G4endl;
  }
  return theMonopole;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Monopole* G4Monopole::Monopole()
{    
  if(!theMonopole) { theMonopole = MonopoleDefinition(); }
  return theMonopole;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4Monopole::MagneticCharge() const
{
  return magCharge;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

