//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4mplIonisationWithDeltaModel
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 06.09.2005
//
// Modifications:
// 12.08.2007 Changing low energy approximation and extrapolation. 
//            Small bug fixing and refactoring (M. Vladymyrov)
// 13.11.2007 Use low-energy asymptotic from [3] (V.Ivanchenko) 
//
//
// -------------------------------------------------------------------
// References
// [1] Steven P. Ahlen: Energy loss of relativistic heavy ionizing particles, 
//     S.P. Ahlen, Rev. Mod. Phys 52(1980), p121
// [2] K.A. Milton arXiv:hep-ex/0602040
// [3] S.P. Ahlen and K. Kinoshita, Phys. Rev. D26 (1982) 2347


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4mplIonisationWithDeltaModel.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4Electron.hh"
#include "G4DynamicParticle.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleTable.hh"
#include "G4Proton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

std::vector<G4double>* G4mplIonisationWithDeltaModel::dedx0 = nullptr;

G4mplIonisationWithDeltaModel::G4mplIonisationWithDeltaModel(G4double mCharge,
                                                             const G4String& nam)
  : G4VEmModel(nam),G4VEmFluctuationModel(nam),
  magCharge(mCharge),
  twoln10(std::log(100.0)),
  betalow(0.01),
  betalim(0.1),
  beta2lim(betalim*betalim),
  bg2lim(beta2lim*(1.0 + beta2lim))
{
  nmpl = G4lrint(std::abs(magCharge) * 2 * fine_structure_const);
  G4cout << "magnetic charge is " << nmpl << G4endl;
  G4cout << "fine structure const is " << fine_structure_const << G4endl;

  //if(nmpl > 6)      { nmpl = 6; }
  //else if(nmpl < 1) { nmpl = 1; }
  if (nmpl < 1) { nmpl = 1; }

  pi_hbarc2_over_mc2 = pi * hbarc * hbarc / (electron_mass_c2);
  G4cout << "hbar c is " << hbarc << G4endl;
  G4cout << "electron mass is " << electron_mass_c2 << G4endl;
  G4cout << "pi hbarc/ mc2 " << pi_hbarc2_over_mc2 << G4endl;

  chargeSquare = magCharge * magCharge;
  dedxlim = 45.*nmpl*nmpl;//*GeV*cm2/g;
  G4cout << "dedx lim is " << dedxlim << G4endl;

  fParticleChange = nullptr;
  theElectron = G4Electron::Electron();
  G4cout << "### Monopole ionisation model with d-electron production, Gmag= " 
         << magCharge/eplus << " and to check if this is working" <<  G4endl;
  monopole = nullptr;
  mass = 0.0;

  G4cout << "Magnetic charge is : "<< nmpl << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4mplIonisationWithDeltaModel::~G4mplIonisationWithDeltaModel()
{
  if(IsMaster()) { delete dedx0; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4mplIonisationWithDeltaModel::SetParticle(const G4ParticleDefinition* p)
{
  monopole = p;
  mass     = monopole->GetPDGMass();
  G4double emin = 
    std::min(LowEnergyLimit(), 0.1*mass*(1./sqrt(1. - betalow*betalow) - 1.)); 
  G4cout << "Low energy limit is " << LowEnergyLimit() << G4endl;
  G4cout << "High energy limit is " << HighEnergyLimit() << G4endl;
  G4cout << "min energy is " << emin <<  G4endl;

  G4double emax = 
    std::max(HighEnergyLimit(),10*mass*(1./sqrt(1. - beta2lim) - 1.)); 
  SetLowEnergyLimit(emin);
  SetHighEnergyLimit(emax);
  G4cout << "changed new Low energy limit is " << LowEnergyLimit() << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4mplIonisationWithDeltaModel::Initialise(const G4ParticleDefinition* p,
                                          const G4DataVector&)
{
  if(!monopole) { SetParticle(p); }
  if(!fParticleChange) { fParticleChange = GetParticleChangeForLoss(); }
  if(IsMaster()) {
    if(!dedx0) { dedx0 = new std::vector<G4double>; }
    G4ProductionCutsTable* theCoupleTable=
      G4ProductionCutsTable::GetProductionCutsTable();
    G4int numOfCouples = theCoupleTable->GetTableSize();
    G4int n = dedx0->size();
    if(n < numOfCouples) { dedx0->resize(numOfCouples); }
    G4Pow* g4calc = G4Pow::GetInstance();

    // initialise vector assuming low conductivity
    for(G4int i=0; i<numOfCouples; ++i) {

      const G4Material* material = 
        theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
      G4double eDensity = material->GetElectronDensity();
      G4double vF2 = 2*electron_Compton_length*g4calc->A13(3.*pi*pi*eDensity);
      (*dedx0)[i] = pi_hbarc2_over_mc2*eDensity*nmpl*nmpl*
        (G4Log(vF2/fine_structure_const) - 0.5)/vF2;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4mplIonisationWithDeltaModel::MinEnergyCut(const G4ParticleDefinition*,
                                            const G4MaterialCutsCouple* couple)
{
  return couple->GetMaterial()->GetIonisation()->GetMeanExcitationEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4mplIonisationWithDeltaModel::ComputeDEDXPerVolume(const G4Material* material,
                                                    const G4ParticleDefinition* p,
                                                    G4double kineticEnergy,
                                                    G4double maxEnergy)
{
  if(!monopole) { SetParticle(p); }
  G4double tmax = MaxSecondaryEnergy(p,kineticEnergy);
  G4double cutEnergy = std::min(tmax, maxEnergy);
  cutEnergy = std::max(LowEnergyLimit(), cutEnergy);
  G4double tau   = kineticEnergy / mass;
  G4double gam   = tau + 1.0;
  G4double bg2   = tau * (tau + 2.0);
  G4double beta2 = bg2 / (gam * gam);
  G4double beta  = sqrt(beta2);
  G4double eDensity = material->GetElectronDensity();
  G4double alpha_M = 1/(4*fine_structure_const);
  mass = monopole->GetPDGMass();
  const G4ElementVector* element = material->GetElementVector();
  const G4Element* firstElement = (*element)[1];
  G4int Z = firstElement->GetZ();
  G4int A = firstElement->GetN();
  //G4cout << "A is " << A << G4endl
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* proton = particleTable->FindParticle("proton");
  G4double protonMass = proton->GetPDGMass();
  //G4cout << "proton Mass is " << protonMass << G4endl;

  G4String matName = material->GetName();
  //G4cout << "Material name is " << matName << G4endl;
  //G4cout << "Atomic number of element is " << Z << G4endl;
  // lowienergy asymptotic formula
  G4double dedx = (*dedx0)[CurrentCouple()->GetIndex()]*beta;
		//+ (6.24*pow(10, -4)*(13.93 + log(beta) - log(nmpl)))*GeV*cm2/g;

  //G4cout << "Material density is : " << CurrentCouple()->GetIndex() << G4endl;

  // above asymptotic
  if(beta > betalow) {

    // high energy
    if(beta >= betalim) {
      dedx = ComputeDEDXAhlen(material, bg2, cutEnergy);
      //G4cout << "Ionization dedx is " << dedx << G4endl;		
	// stochastic energy loss (brem/pair prod implementation)
      if (gam > 1){
		G4double dedx_brem = (16/3)*(Z*eDensity*(pow(nmpl*hbarc/2, 2))*(alpha_M)/(mass))*gam*log(gam);
		//G4cout << "Z is " << Z << G4endl;
		//G4cout << "eDensity is " << eDensity << G4endl;
		//G4cout << "pow(nmpl*hbarc/2, 2) is  " << pow(nmpl*hbarc/2, 2) << G4endl;
		//G4cout << "alpha M is " << alpha_M << G4endl;
		//G4cout << "mass is " << mass << G4endl;
		//G4cout << "gamma is " << gam << G4endl;
		//G4cout << "DEDX_brem " << dedx_brem << G4endl;

		//dedx += dedx_brem;

		G4double gamma_cm = sqrt((gam*mass)/(2*A*protonMass));
		//G4cout << "gamma_cm is " << gamma_cm << G4endl;
		G4double dedx_pair = (16/3e3)*((pow(nmpl*hbarc/2, 2))*Z*fine_structure_const*eDensity/electron_mass_c2)*gamma_cm*log(gam);	
		//G4cout << "DEDX_pair " << dedx_pair << G4endl;

		dedx += dedx_brem + dedx_pair;

	}	
      //G4cout << "Overall dedx including brem/pair is " << dedx << G4endl;	

    } else {
      G4double dedx1 = (*dedx0)[CurrentCouple()->GetIndex()]*betalow;
      G4double dedx2 = ComputeDEDXAhlen(material, bg2lim, cutEnergy);

      // extrapolation between two formula 
      G4double kapa2 = beta - betalow;
      G4double kapa1 = betalim - beta;
      dedx = (kapa1*dedx1 + kapa2*dedx2)/(kapa1 + kapa2);
    }
  }
  
  //G4cout << "DEDX_in_material_is " << G4BestUnit(dedx, "Energy/Length") << G4endl;
  //G4cout << "  " << G4endl;

  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4mplIonisationWithDeltaModel::ComputeDEDXAhlen(const G4Material* material, 
                                                G4double bg2, 
                                                G4double cutEnergy)
{
  G4double eDensity = material->GetElectronDensity();
  //G4cout << "eDensity is " << eDensity*cm3  << G4endl;
  G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
  //G4cout << "MeanExcitationEnergy is " << G4BestUnit(eexc, "Energy") << G4endl;

  /*G4cout << "cDensity is " << material->GetIonisation()->GetCdensity()*cm3 << G4endl;;
  G4cout << "mDensity is " << material->GetIonisation()->GetMdensity()*cm3 << G4endl;;
  G4cout << "aDensity is " << material->GetIonisation()->GetAdensity()*cm3 << G4endl;;
  G4cout << "x0Density is " << material->GetIonisation()->GetX0density()*cm3 << G4endl;;
  G4cout << "x1Density is " << material->GetIonisation()->GetX1density()*cm3 << G4endl;
  G4cout << "D0Density is " << material->GetIonisation()->GetD0density()*cm3 << G4endl;*/
  // Ahlen's formula for nonconductors, [1]p157, f(5.7)
  G4double dedx = 
    0.5*(G4Log(2.0*electron_mass_c2*bg2*cutEnergy/(eexc*eexc)) -1.0);
  //G4cout << "electron mass c2 = " << electron_mass_c2 << G4endl;
  //G4cout << "beta gamma sq = " << bg2 << G4endl;
  //G4cout << "cut Energy is = " << cutEnergy << G4endl;
  //G4cout << "Ionization energy is = " << eexc << G4endl;
  //G4cout << "Ahlen output " << dedx << G4endl;


  // Kazama et al. cross-section correction
  //G4double  k = 0.406;
  //if(nmpl > 1) { k = 0.346; }

  G4double  k;
  if(nmpl == 1) { k = 0.406; }
  else if(nmpl == 2 ) { k = 0.346; }
  else
  	k = 0.300;

  // Bloch correction
  //const G4double B[7] = { 0.0, 0.248, 0.672, 1.022, 1.243, 1.464, 1.685}; 
  //const G4double B[13] = { 0.0, 0.248, 0.672, 1.022, 1.243, 1.464, 1.685, 1.837, 1.969, 2.086, 2.190, 2.285, 2.371};
  const G4double B[21] = { 0.0, 0.248, 0.672, 1.022, 1.243, 1.464, 1.685, 1.837, 1.969, 2.086, 2.190, 2.285, 2.371, 2.451, 2.525, 2.594, 2.658, 2.719, 2.776, 2.830, 2.881};


  if(nmpl < 1) { dedx += 0.5 * k - B[1]; }
  else
    dedx += 0.5 * k - B[nmpl];
  //dedx += 0.5 * k - B[nmpl];
  
  //G4cout << "Magnetic charge is : "<< nmpl << G4endl;
  //G4cout << "Kazama correction factor is : " << k << G4endl;
  //G4cout << "Bloch correction factor is : " << B[nmpl] << G4endl;

  // density effect correction
  G4double x = G4Log(bg2)/twoln10;
  dedx -= material->GetIonisation()->DensityCorrection(x);

  //G4cout << "after density effect correction " << dedx << G4endl;

  // now compute the total ionization loss
  dedx *=  pi_hbarc2_over_mc2 * eDensity * nmpl * nmpl;
  //G4cout << "nmpl = " << nmpl << G4endl;
  //G4cout << " pi_hbarc2_over_mc2 = " << pi_hbarc2_over_mc2 << G4endl;
  //G4double matDensity = material->GetDensity();
  //G4cout << "Material density in g/cm3 " << matDensity/(g/cm3) << G4endl;
  //G4cout << "total ionization loss in MeV cm2/g" << dedx << G4endl;

  //G4cout << "DEDX_ion " << dedx << G4endl;
  dedx = std::max(dedx, 0.0);
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4mplIonisationWithDeltaModel::ComputeCrossSectionPerElectron(
                                           const G4ParticleDefinition* p,
                                           G4double kineticEnergy,
                                           G4double cut,
                                           G4double maxKinEnergy)
{
  if(!monopole) { SetParticle(p); }
  G4double tmax = MaxSecondaryEnergy(p, kineticEnergy);
  G4double maxEnergy = std::min(tmax, maxKinEnergy);
  G4double cutEnergy = std::max(LowEnergyLimit(), cut);
  G4double cross = (cutEnergy < maxEnergy) 
    ? (0.5/cutEnergy - 0.5/maxEnergy)*pi_hbarc2_over_mc2 * nmpl * nmpl : 0.0;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4mplIonisationWithDeltaModel::ComputeCrossSectionPerAtom(
                                          const G4ParticleDefinition* p,
                                          G4double kineticEnergy,
                                          G4double Z, G4double,
                                          G4double cutEnergy,
                                          G4double maxEnergy)
{
  G4double cross = 
    Z*ComputeCrossSectionPerElectron(p,kineticEnergy,cutEnergy,maxEnergy);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4mplIonisationWithDeltaModel::SampleSecondaries(vector<G4DynamicParticle*>* vdp,
                                                 const G4MaterialCutsCouple*,
                                                 const G4DynamicParticle* dp,
                                                 G4double minKinEnergy,
                                                 G4double maxEnergy)
{
  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double tmax = MaxSecondaryEnergy(dp->GetDefinition(),kineticEnergy);

  G4double maxKinEnergy = std::min(maxEnergy,tmax);
  if(minKinEnergy >= maxKinEnergy) { return; }

  //G4cout << "G4mplIonisationWithDeltaModel::SampleSecondaries: E(GeV)= "
  //         << kineticEnergy/GeV << " M(GeV)= " << mass/GeV
  //         << " tmin(MeV)= " << minKinEnergy/MeV << G4endl;

  G4double totEnergy     = kineticEnergy + mass;
  G4double etot2         = totEnergy*totEnergy;
  G4double beta2         = kineticEnergy*(kineticEnergy + 2.0*mass)/etot2;
  
  // sampling without nuclear size effect
  G4double q = G4UniformRand();
  G4double deltaKinEnergy = minKinEnergy*maxKinEnergy
    /(minKinEnergy*(1.0 - q) + maxKinEnergy*q);

  // delta-electron is produced
  G4double totMomentum = totEnergy*sqrt(beta2);
  G4double deltaMomentum =
           sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
  G4double cost = deltaKinEnergy * (totEnergy + electron_mass_c2) /
                                   (deltaMomentum * totMomentum);
  cost = std::min(cost, 1.0);

  G4double sint = sqrt((1.0 - cost)*(1.0 + cost));

  G4double phi = twopi * G4UniformRand() ;

  G4ThreeVector deltaDirection(sint*cos(phi),sint*sin(phi), cost);
  G4ThreeVector direction = dp->GetMomentumDirection();
  deltaDirection.rotateUz(direction);

  // create G4DynamicParticle object for delta ray
  G4DynamicParticle* delta = 
    new G4DynamicParticle(theElectron,deltaDirection,deltaKinEnergy);

  vdp->push_back(delta);

  // Change kinematics of primary particle
  kineticEnergy       -= deltaKinEnergy;
  G4ThreeVector finalP = direction*totMomentum - deltaDirection*deltaMomentum;
  finalP               = finalP.unit();

  fParticleChange->SetProposedKineticEnergy(kineticEnergy);
  fParticleChange->SetProposedMomentumDirection(finalP);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/*G4double G4mplIonisationWithDeltaModel::SampleFluctuations(
                                       const G4MaterialCutsCouple* couple,
                                       const G4DynamicParticle* dp,
                                       G4double tmax,
                                       G4double length,
                                       G4double meanLoss)
{
  G4double siga = Dispersion(couple->GetMaterial(),dp,tmax,length);
  G4double loss = meanLoss;
  siga = sqrt(siga);
  G4double twomeanLoss = meanLoss + meanLoss;

  if(twomeanLoss < siga) {
    G4double x;
    do {
      loss = twomeanLoss*G4UniformRand();
      x = (loss - meanLoss)/siga;
      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while (1.0 - 0.5*x*x < G4UniformRand());
  } else {
    do {
      loss = G4RandGauss::shoot(meanLoss,siga);
      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while (0.0 > loss || loss > twomeanLoss);
  }
  return loss;
}*/

G4double G4mplIonisationWithDeltaModel::SampleFluctuations(
                                       const G4MaterialCutsCouple* couple,
                                       const G4DynamicParticle* dp,
                                       const G4double tcut,
                                       const G4double tmax,
                                       const G4double length,
                                       const G4double meanLoss)
{
  G4double siga = Dispersion(couple->GetMaterial(),dp,tcut,tmax,length);
  G4double loss = meanLoss;
  siga = std::sqrt(siga);
  G4double twomeanLoss = meanLoss + meanLoss;

  if(twomeanLoss < siga) {
    G4double x;
    do {
      loss = twomeanLoss*G4UniformRand();
      x = (loss - meanLoss)/siga;
      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while (1.0 - 0.5*x*x < G4UniformRand());
  } else {
    do {
      loss = G4RandGauss::shoot(meanLoss,siga);
      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while (0.0 > loss || loss > twomeanLoss);
  }
  return loss;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/*G4double 
G4mplIonisationWithDeltaModel::Dispersion(const G4Material* material,
                                          const G4DynamicParticle* dp,
                                          G4double tmax,
                                          G4double length)
{
  G4double siga = 0.0;
  G4double tau   = dp->GetKineticEnergy()/mass;
  if(tau > 0.0) { 
    G4double electronDensity = material->GetElectronDensity();
    G4double gam   = tau + 1.0;
    G4double invbeta2 = (gam*gam)/(tau * (tau+2.0));
    siga  = (invbeta2 - 0.5) * twopi_mc2_rcl2 * tmax * length
      * electronDensity * chargeSquare;
  }
  return siga;
}*/

G4double
G4mplIonisationWithDeltaModel::Dispersion(const G4Material* material,
                                          const G4DynamicParticle* dp,
                                          const G4double tcut,
                                          const G4double tmax,
                                          const G4double length)
{
  G4double siga = 0.0;
  G4double tau   = dp->GetKineticEnergy()/mass;
  if(tau > 0.0) {
    const G4double beta = dp->GetBeta();
    siga  = (tmax/(beta*beta) - 0.5*tcut) * twopi_mc2_rcl2 * length
      * material->GetElectronDensity() * chargeSquare;
  }
  return siga;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4mplIonisationWithDeltaModel::MaxSecondaryEnergy(const G4ParticleDefinition*,
                                                  G4double kinEnergy)
{
  G4double tau = kinEnergy/mass;
  return 2.0*electron_mass_c2*tau*(tau + 2.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
