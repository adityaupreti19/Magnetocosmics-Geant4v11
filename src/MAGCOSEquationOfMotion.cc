//this source file sets the equation of motion
//here I've modified the source file to include the equation of motion of a monopole. 

#include "MAGCOSEquationOfMotion.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>


////////////////////////////////////////////////////////////////////////////////
//
MAGCOSEquationOfMotion::MAGCOSEquationOfMotion( G4MagneticField* MagField )
                              : G4Mag_EqRhs( MagField ) 
{ G4cout << "G4MonopoleEquation::G4MonopoleEquation" << G4endl;
  reverse_time=false;
  direction=1.;
  //Equation=&MAGCOSEquationOfMotion::MonopoleMotion;
  Equation=&MAGCOSEquationOfMotion::LorentzMotion;  
  		//initially set as Lorentz motion, then changed later
}
////////////////////////////////////////////////////////////////////////////////
//


void MAGCOSEquationOfMotion::SetChargeMomentumMass( G4ChargeState particleChargeState,
                                           G4double      ,           // momentum,
                                           G4double particleMass)
{
   G4double particleMagneticCharge= particleChargeState.MagneticCharge();
   G4double particleElectricCharge= particleChargeState.GetCharge();

  //   fElCharge = particleElectricCharge;
  fElCharge =eplus* particleElectricCharge*c_light;

  fMagCharge =  eplus*particleMagneticCharge*c_light ;

  //G4cout << " G4MonopoleEquation: ElectricCharge=" << particleElectricCharge
  //             << "; MagneticCharge=" << particleMagneticCharge
  //           << G4endl;

  fMassCof = particleMass*particleMass ;
}

void MAGCOSEquationOfMotion::EvaluateRhsGivenB( const G4double y[],
			                        const G4double B[],
				                      G4double dydx[] ) const
{ 
  (this->*Equation)( y,  B, dydx);
  return ;
}
/////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSEquationOfMotion::MotionAlongBfieldLine( const G4double*,
			                            const G4double B[3],
				                          G4double dydx[] ) const
{
  //G4cout << "Motion along the B field line is being used" << G4endl;	
  G4double Bmag=direction*std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
  dydx[0] = B[0]/Bmag;       
  dydx[1] = B[1]/Bmag;       
  dydx[2] = B[2]/Bmag;
   
  dydx[3]=0.;
  dydx[4]=0.;
  dydx[5]=0.; 
  return ;
}

/////////////////////////////////////////////////////////////////////////////


//equation of motion of a monopole. 

void MAGCOSEquationOfMotion::MonopoleMotion(const G4double y[],
						const G4double B[], 
						G4double dydx[] ) const
{
  //G4cout << "The monopole motion equation is being used" << G4endl;	
  G4double pSquared = y[3]*y[3] + y[4]*y[4] + y[5]*y[5] ;
  G4double Energy   = std::sqrt( pSquared + fMassCof );

  G4double pModuleInverse  = 1.0/std::sqrt(pSquared);

  G4double inverse_velocity = Energy * pModuleInverse / c_light;

  G4double cofEl     = fElCharge * pModuleInverse ;
  G4double cofMag = fMagCharge * Energy * pModuleInverse;

  //G4cout << "cofMag is : " << cofMag << G4endl;

  dydx[0] = y[3]*pModuleInverse ;
  dydx[1] = y[4]*pModuleInverse ;
  dydx[2] = y[5]*pModuleInverse ;

  // G4double magCharge = twopi * hbar_Planck / (eplus * mu0);    
  // magnetic charge in SI units A*m convention
  //  see http://en.wikipedia.org/wiki/Magnetic_monopole   
  //   G4cout  << "Magnetic charge:  " << magCharge << G4endl;   
  // dp/ds = dp/dt * dt/ds = dp/dt / v = Force / velocity
  // dydx[3] = fMagCharge * Field[0]  * inverse_velocity  * c_light;    
  // multiplied by c_light to convert to MeV/mm
  //     dydx[4] = fMagCharge * Field[1]  * inverse_velocity  * c_light; 
  //     dydx[5] = fMagCharge * Field[2]  * inverse_velocity  * c_light; 

  dydx[3] = cofMag * B[0] + cofEl * (y[4]*B[2] - y[5]*B[1]);
  dydx[4] = cofMag * B[1] + cofEl * (y[5]*B[0] - y[3]*B[2]);
  dydx[5] = cofMag * B[2] + cofEl * (y[3]*B[1] - y[4]*B[0]);

  //        G4cout << std::setprecision(5)<< "E=" << Energy
  //               << "; p="<< 1/pModuleInverse
  //               << "; mC="<< magCharge
  //               <<"; x=" << y[0]
  //               <<"; y=" << y[1]
  //               <<"; z=" << y[2]
  //               <<"; dydx[3]=" << dydx[3]
  //               <<"; dydx[4]=" << dydx[4]
  //               <<"; dydx[5]=" << dydx[5]
  //               << G4endl;

  dydx[6] = 0.;//not used

  // Lab Time of flight
  dydx[7] = inverse_velocity;
  return;
}

/////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSEquationOfMotion::LorentzMotion( const G4double y[],
			                    const G4double B[3],
				            G4double dydx[] ) const
{
  //G4cout << "The lorentz motion is being used" << G4endl;	
  G4double momentum_mag_square = y[3]*y[3] +  y[4]*y[4] +  y[5]*y[5];
  G4double inv_momentum_magnitude = 
            direction * 1.0 / sqrt( momentum_mag_square );

  G4double cof =   FCof()*inv_momentum_magnitude;

  dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
  dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
  dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

  dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
  dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
  dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)
   

  return ;
}

////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSEquationOfMotion::SetReverseTimeMode(G4bool abool)
{ reverse_time=abool;
  direction=1.;
  if (reverse_time) direction= -1.;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSEquationOfMotion::SetEquationType(G4String aString)
{ if (aString=="LORENTZ_MOTION")  
          Equation=&MAGCOSEquationOfMotion::LorentzMotion;
  if (aString=="BFIELD_LINE")
           Equation=&MAGCOSEquationOfMotion::MotionAlongBfieldLine;
  if (aString=="MONOPOLE")
           Equation=&MAGCOSEquationOfMotion::MonopoleMotion;
}




