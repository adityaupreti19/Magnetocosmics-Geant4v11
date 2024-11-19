//constructs the geometry of the Earth in the simulation
//this file should be edited to increase the extent of the 'World' volume
//also sets the altitude at which the track must be 'killed' when it reaches Earth. 

#include "MAGCOSGeometryConstruction.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ParticleTable.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include "G4UserLimits.hh"
#include "MAGCOSGeometryMessenger.hh"
#include "MAGCOSMagneticField.hh"
#include "MAGCOSUnits.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "MAGCOSSteppingAction.hh"

/////monopole/////
//#include "G4MonopoleFieldSetup.hh"
#include "G4AutoDelete.hh"
#include "G4NistManager.hh"

#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"

MAGCOSGeometryConstruction::MAGCOSGeometryConstruction()
{ 
  detectorMessenger = new MAGCOSGeometryMessenger(this);
  theMagnetosphereMagneticField=new MAGCOSMagneticField();
  TopOfAtmosphere =1*m;
  RemoveEarth =false;
}
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSGeometryConstruction::~MAGCOSGeometryConstruction()
{ delete detectorMessenger;
  delete theMagnetosphereMagneticField;
}
////////////////////////////////////////////////////////////////////////////////
//
G4VPhysicalVolume* MAGCOSGeometryConstruction::Construct()
{ 
  // Clean old geometry, if any
  //
  
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();


  G4double z, a,density, fractionmass;
  G4String name, symbol;
  G4int nel, ncomponents, natoms;
 
  
  //Vacuum
     
    a = 1.*g/mole; 
    density = 1.e-10*g/cm3;
    G4Element* elN = new G4Element( "Nitrogen", "N", 7. , 14.00674*g/mole );
    G4Element* elO = new G4Element( "Oxygen", "O", 8. , 15.9994*g/mole );
    G4Element* elC = new G4Element(name = "Carbon", symbol = "C", z = 6, a = 12.01*g/mole);
    G4Element* elH = new G4Element(name = "Hydrogen", symbol = "H", z = 1, a = 1.01*g/mole);
    G4Element* elHe = new G4Element(name = "Helium", symbol="He", z=2, a =4.0026*g/mole);
    G4Element* elAr = new G4Element(name = "Argon", symbol = "Ar", z = 18, a = 39.948*g/mole);
  
    G4Material* Vacuum =new G4Material(name="Vacuum",density,
                     nel=2,kStateSolid, 293.0*kelvin, 1.e-6*atmosphere );
    Vacuum->AddElement(elN, .7);
    Vacuum->AddElement(elO, .3);


    G4double density_sandstone = 2.65 * g/cm3; // typical density for sandstone
    G4Material* sandstone = new G4Material("Sandstone", density_sandstone, 4);
    sandstone->AddElement(G4NistManager::Instance()->FindOrBuildElement("Si"), 1); // Silicon
    sandstone->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), 2);  // Oxygen
    sandstone->AddElement(G4NistManager::Instance()->FindOrBuildElement("Al"), 1); // Aluminum
    sandstone->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 4);  // Hydrogen
 
    G4double density_Hematite = 5.3 * g/cm3;
    G4Material* hematite = new G4Material("Hematite", density_Hematite, 2);
    hematite->AddElement(G4NistManager::Instance()->FindOrBuildElement("Fe"), 2);
    hematite->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), 3);

    auto nistManager = G4NistManager::Instance();
    nistManager->FindOrBuildMaterial("G4_AIR");
    G4Material* Air = G4Material::GetMaterial("G4_AIR");
    nistManager->FindOrBuildMaterial("G4_WATER");
    G4Material* Water = G4Material::GetMaterial("G4_WATER");
    nistManager->FindOrBuildMaterial("G4_Al");
    G4Material* Aluminum = G4Material::GetMaterial("G4_Al"); 
 
    G4Element* elCa = new G4Element(name = "Calcium", symbol = "Ca", z = 20, a = 40.078*g/mole);
    G4Element* elMg = new G4Element(name = "Magnesium", symbol = "Mg", z = 12, a = 24.305*g/mole);

    G4Material* CaCO3 = new G4Material(name="CaCO3", density = 2.71*g/cm3, ncomponents = 3);
    CaCO3 -> AddElement(elO, natoms = 3);
    CaCO3 -> AddElement(elC, natoms = 1);
    CaCO3 -> AddElement(elCa, natoms = 1);

    G4Material* MgCO3 = new G4Material(name="MgCO3", density = 2.96*g/cm3, ncomponents = 3);
    MgCO3 -> AddElement(elO, natoms = 3);
    MgCO3 -> AddElement(elC, natoms = 1);
    MgCO3 -> AddElement(elMg, natoms = 1);

    G4Material* limestone = new G4Material("Limestone", density = 2.711*g/cm3, ncomponents=2);
    limestone->AddMaterial(CaCO3, fractionmass=90*perCent);
    limestone->AddMaterial(MgCO3, fractionmass=10*perCent);


    ///////////// ICE /////////////////////////////////////////////////

    density = 0.917*g/cm3;
    G4double temp = 273.15 *kelvin;
    G4double pressure = 1.0 * atmosphere;
    G4double meanExcitationEnergy = 75.0 * eV;

    G4Material* ice = new G4Material("ice", density, 2);
    ice->AddElement(elH, 2);
    ice->AddElement(elO, 1);

    G4MaterialPropertiesTable* iceProp = new G4MaterialPropertiesTable();
    //iceProp->AddConstProperty("TEMPERATURE", temp);
    //iceProp->AddConstProperty("PRESSURE", pressure);
    ice->SetMaterialPropertiesTable(iceProp);
    ice->GetIonisation()->SetMeanExcitationEnergy(meanExcitationEnergy);   


 //Visualisation attributes
   G4VisAttributes * VisAttEarth = new G4VisAttributes(G4Colour(0.,0.,1.));
	VisAttEarth->SetVisibility(true);		  
   G4VisAttributes * VisAttWorld = new G4VisAttributes();
        VisAttWorld->SetVisibility(false);
    
  
  /*G4double cutValue = 0.000001*m;
  G4ProductionCuts* productionCuts = new G4ProductionCuts();
  productionCuts->SetProductionCut(cutValue);

  SetCutsForVolumes(productionCuts, "e-", "World"); */
  

  /*G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(10*eV, 100*TeV);
  const G4double electronCut = G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCut("e-");
  SetCutValue(10*m, "e-");*/
   
  
  
  //------------------------------ 
  // World
  //------------------------------ 
   G4double MagSizeX=200*re;
   G4double MagSizeY=200*re;
   G4double MagSizeZ=200*re;
   solidWorld= 
          new G4Box("World",MagSizeX/2.,MagSizeY/2.,MagSizeZ/2.);
   logicWorld= new G4LogicalVolume( solidWorld, Vacuum,
                                                  "World", 0, 0, 0);
   maxStepSize=1.*km;
   theWorldUserLimits=new G4UserLimits(maxStepSize);
   theWorldUserLimits->SetUserMaxTrackLength(200.*re);
   logicWorld->SetUserLimits(theWorldUserLimits);
   logicWorld->SetVisAttributes(VisAttWorld);
   /*physiWorld= new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
				 "WorldPV",       // its name
                                 logicWorld,      // its logical volume
                                 0,               // its mother  volume
                                 false,           // no boolean operations
                                 0);              // no field specific to volume*/

    physiWorld= new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
                                 //"WorldPV",       // its name
                                 logicWorld,      // its logical volume
                                 "World",
                                 0,              // its mother  volume
                                 false,           // no boolean operations
                                 0);              // no field specific to volume

 


  G4VisAttributes * VisAttSpace = new G4VisAttributes();
  VisAttSpace->SetVisibility(true);    

  G4double radSpace = 20*re;
  G4Orb* Space = new G4Orb("Space", radSpace);
  G4LogicalVolume* logSpace = new G4LogicalVolume(Space, Vacuum, "Space");
  G4UserLimits* theSpaceUserLimits=new G4UserLimits(1*km);//100*km);
  logSpace->SetUserLimits(theSpaceUserLimits);
  theSpaceUserLimits->SetUserMaxTrackLength(200*re);
  //logSpace->SetUserLimits(theSpaceUserLimits);
  logSpace->SetVisAttributes(VisAttSpace);

  G4VPhysicalVolume* fPhysSpace = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                  logSpace, "Space", logicWorld, false, 0);

///----------------upper atmosphere beyond 100km above sea level---//


 G4double upperDensityAtmos[] = {1.647e-12, 1.982e-13, 3.982e-14, 1.011e-14, 2.94e-15, 9.352e-16, 3.247e-16, 1.266e-16, 5.78e-17, 3.148e-17, 1.991e-17, 1.391e-17, 1.028e-17, 7.834e-18, 6.073e-18, 4.758e-18, 3.757e-18, 2.986e-18};
 G4double upper_N2_percent[] = {91.05991296, 92.13538036, 86.70235757, 60.76053178, 22.0004288, 4.78738325, 0.90476058, 0.1692613, 0.03217586, 0.00624428, 0.00123691, 0.00024982, 5.142e-05, 1.078e-05, 2.3e-06, 5e-07, 1.1e-07, 2e-08};
 G4double upper_O2_percent[] = {8.57102884, 6.30886841, 4.45350211, 2.32947051, 0.6304389, 0.10291625, 0.01465077, 0.00207443, 0.00029952, 4.432e-05, 6.72e-06, 1.04e-06, 1.7e-07, 3e-08, 0.0, 0.0, 0.0, 0.0};
 G4double upper_He_percent[] = {0.1880897, 1.46806073, 8.66647753, 36.11418469, 75.3221013, 92.00962385, 95.11788506, 94.95061081, 94.03155571, 92.8124917, 91.36117826, 89.67224622, 87.72247264, 85.50299492, 82.98834891, 80.18856625, 77.10136415, 73.72184482};
 G4double upper_Ar_percent[] = {0.17626625, 0.06480003, 0.0243646, 0.0070032, 0.00105604, 9.705e-05, 7.84e-06, 6.4e-07, 5e-08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
 G4double upper_H_percent[] = {0.00470224, 0.02289046, 0.15329819, 0.78880982, 2.04597495, 3.0999796, 3.96269575, 4.87805282, 5.93596886, 7.18121969, 8.63757812, 10.32750291, 12.27747577, 14.49699427, 17.01164878, 19.81143325, 22.89863574, 26.27815516};
 G4double upperTempAtmos[] = {613.5, 712.4, 729.1, 731.9, 732.4, 732.5, 732.6, 732.6, 732.6, 732.6, 732.6, 732.6, 732.6, 732.6, 732.6, 732.6, 732.6, 732.6};
 G4String upperVolName[] = {"air150", "air200", "air250", "air300", "air350", "air400", "air450", "air500", "air550", "air600", "air650", "air700", "air750", "air800", "air850", "air900", "air950", "air1000"};



for (G4int i=0; i<18; i++)
{
        upperAirMat[i] = new G4Material(name=upperVolName[i], upperDensityAtmos[i]*g/cm3, ncomponents=5, kStateGas, upperTempAtmos[i]*kelvin);
        upperAirMat[i]->AddElement(elN, upper_N2_percent[i]*perCent);
        upperAirMat[i]->AddElement(elO, upper_O2_percent[i]*perCent);
        upperAirMat[i]->AddElement(elHe, upper_He_percent[i]*perCent);
        upperAirMat[i]->AddElement(elAr, upper_Ar_percent[i]*perCent);
        upperAirMat[i]->AddElement(elH, upper_H_percent[i]*perCent);
}

        G4int l = 17;
        upperAtmos[l] = new G4Orb("upperAtmos", (6371.2+1000)*km);
        upperLogAtmos[l] = new G4LogicalVolume(upperAtmos[l], upperAirMat[l], "upperAtmos");
        upperAtmosUserLimits[l] = new G4UserLimits(100*m);
        upperLogAtmos[l]->SetUserLimits(upperAtmosUserLimits[l]);
        upperAtmosUserLimits[l]->SetUserMaxTrackLength(200*re);
        upperPhysAtmos[l] = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                        upperLogAtmos[l], "upperAtmos", logSpace, false, 0);

  for (G4int i = 16; i>=0; i--)
    {
        upperAtmos[i] = new G4Orb("upperAtmos", (6371.2+i*50+151)*km);
        upperLogAtmos[i] = new G4LogicalVolume(upperAtmos[i], upperAirMat[i], "upperAtmos");
        upperAtmosUserLimits[i] = new G4UserLimits(100*m);
        upperLogAtmos[i]->SetUserLimits(upperAtmosUserLimits[i]);
        upperAtmosUserLimits[i]->SetUserMaxTrackLength(200*re);
        upperPhysAtmos[i] = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                        upperLogAtmos[i], "upperAtmos", upperLogAtmos[i+1], false, 0);

    }




  //--------------------density in atmosphere-----------------------//
  
  G4double densityAtmos[] = {0.001295, 0.001155, 0.001036, 0.0009324, 0.0008389, 0.0007528, 0.0006724, 0.000597, 0.0005262, 0.0004601, 0.0003989, 0.000343, 0.0002932, 0.0002499, 0.0002129, 0.0001815, 0.0001551, 0.0001327, 0.0001136, 9.709e-05, 8.284e-05, 7.05e-05, 5.987e-05, 5.074e-05, 4.294e-05, 3.63e-05, 3.066e-05, 2.588e-05, 2.184e-05, 1.844e-05, 1.557e-05, 1.315e-05, 1.112e-05, 9.411e-06, 7.985e-06, 6.792e-06, 5.794e-06, 4.955e-06, 4.25e-06, 3.656e-06, 3.154e-06, 2.728e-06, 2.367e-06, 2.059e-06, 1.796e-06, 1.571e-06, 1.378e-06, 1.211e-06, 1.066e-06, 9.401e-07, 8.295e-07, 7.324e-07, 6.467e-07, 5.709e-07, 5.039e-07, 4.443e-07, 3.915e-07, 3.446e-07, 3.029e-07, 2.66e-07, 2.332e-07, 2.042e-07, 1.786e-07, 1.56e-07, 1.36e-07, 1.184e-07, 1.029e-07, 8.931e-08, 7.741e-08, 6.7e-08, 5.791e-08, 4.998e-08, 4.309e-08, 3.722e-08, 3.206e-08, 2.762e-08, 2.379e-08, 2.05e-08, 1.765e-08, 1.519e-08, 1.307e-08, 1.124e-08, 9.654e-09, 8.286e-09, 7.105e-09, 6.085e-09, 5.206e-09, 4.449e-09, 3.796e-09, 3.235e-09, 2.753e-09, 2.34e-09, 1.985e-09, 1.68e-09, 1.419e-09, 1.195e-09, 1.005e-09, 8.422e-10, 7.042e-10, 5.874e-10, 4.888e-10, 4.06e-10, 3.368e-10, 2.791e-10, 2.314e-10, 1.919e-10, 1.593e-10, 1.324e-10, 1.103e-10, 9.207e-11, 7.706e-11, 6.466e-11, 5.443e-11, 4.597e-11, 3.898e-11, 3.318e-11, 2.837e-11, 2.437e-11, 2.105e-11, 1.828e-11, 1.598e-11, 1.407e-11, 1.249e-11, 1.119e-11, 1.015e-11, 9.22e-12, 8.406e-12, 7.688e-12, 7.051e-12, 6.484e-12, 5.977e-12, 5.522e-12, 5.113e-12, 4.744e-12, 4.409e-12, 4.106e-12, 3.83e-12, 3.577e-12, 3.347e-12, 3.135e-12, 2.941e-12, 2.763e-12, 2.598e-12, 2.446e-12, 2.305e-12, 2.174e-12, 2.053e-12, 1.941e-12, 1.836e-12, 1.738e-12};


  G4double N2_percent[] = {78.114, 78.109, 78.105, 78.115, 78.104, 78.119, 78.112, 78.112, 78.107, 78.114, 78.115, 78.106, 78.105, 78.106, 78.113, 78.113, 78.113, 78.109, 78.111, 78.108, 78.105, 78.112, 78.114, 78.114, 78.113, 78.108, 78.108, 78.115, 78.11, 78.113, 78.109, 78.107, 78.109, 78.117, 78.111, 78.109, 78.113, 78.108, 78.112, 78.111, 78.112, 78.104, 78.111, 78.112, 78.112, 78.11, 78.111, 78.114, 78.107, 78.111, 78.107, 78.107, 78.108, 78.107, 78.109, 78.114, 78.11, 78.114, 78.107, 78.107, 78.111, 78.112, 78.11, 78.117, 78.134, 78.146, 78.165, 78.178, 78.195, 78.196, 78.22, 78.233, 78.247, 78.259, 78.278, 78.296, 78.308, 78.333, 78.355, 78.38, 78.403, 78.436, 78.469, 78.5, 78.549, 78.589, 78.635, 78.695, 78.756, 78.822, 78.897, 78.98, 79.073, 79.176, 79.289, 79.417, 79.553, 79.711, 79.876, 80.073, 80.279, 80.50373336, 80.75540355, 81.01205289, 81.29266487, 81.5889176, 81.89185701, 82.21272807, 82.53829063, 82.87455928, 83.20949577, 83.54552075, 83.88417594, 84.21836316, 84.54529461, 84.87621866, 85.19344416, 85.5001667, 85.80117195, 86.09105431, 86.37107414, 86.64122859, 86.90284569, 87.1595817, 87.39964859, 87.63335074, 87.86112326, 88.07994621, 88.28276221, 88.48127751, 88.67707076, 88.86410195, 89.03512635, 89.20197338, 89.36306687, 89.51371941, 89.65776925, 89.79496965, 89.92659217, 90.0495704, 90.16603743, 90.27844931, 90.38299748, 90.48528862, 90.58024645, 90.67021217, 90.75519542, 90.83592928, 90.91362346, 90.98980071};

  G4double O2_percent[] = {20.952, 20.957, 20.96, 20.951, 20.961, 20.947, 20.954, 20.953, 20.958, 20.951, 20.95, 20.959, 20.96, 20.959, 20.953, 20.953, 20.952, 20.956, 20.954, 20.957, 20.96, 20.953, 20.951, 20.952, 20.952, 20.957, 20.957, 20.95, 20.955, 20.953, 20.957, 20.958, 20.956, 20.949, 20.954, 20.956, 20.953, 20.957, 20.953, 20.954, 20.953, 20.962, 20.954, 20.953, 20.954, 20.956, 20.954, 20.951, 20.958, 20.954, 20.958, 20.958, 20.958, 20.958, 20.956, 20.951, 20.955, 20.951, 20.959, 20.959, 20.955, 20.954, 20.955, 20.948, 20.931, 20.919, 20.9, 20.887, 20.871, 20.869, 20.845, 20.832, 20.819, 20.807, 20.787, 20.77, 20.758, 20.734, 20.712, 20.687, 20.665, 20.633, 20.6, 20.57, 20.523, 20.484, 20.439, 20.381, 20.322, 20.259, 20.187, 20.108, 20.02, 19.922, 19.815, 19.694, 19.567, 19.419, 19.265, 19.082, 18.89, 18.68250614, 18.44913027, 18.21256446, 17.95347964, 17.68021384, 17.40083255, 17.10437912, 16.8030972, 16.49153705, 16.18110821, 15.86884748, 15.55322556, 15.24103099, 14.93485, 14.62359812, 14.32466107, 14.03468426, 13.7491708, 13.47357588, 13.20608613, 12.94755315, 12.69610033, 12.44854377, 12.21621186, 11.98969788, 11.76830875, 11.55544248, 11.35754514, 11.16360795, 10.97184368, 10.78817926, 10.62004366, 10.455621, 10.29641044, 10.14718232, 10.00421231, 9.86746944, 9.73594279, 9.61270057, 9.49537965, 9.38192119, 9.27577849, 9.1715527, 9.07409922, 8.98137555, 8.89309704, 8.80862662, 8.72688262, 8.64622183};

  G4double He_percent[] = {0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524, 0.000525, 0.000525, 0.000527, 0.000528, 0.000529, 0.00053, 0.000531, 0.000532, 0.000533, 0.000534, 0.000535, 0.000537, 0.000539, 0.000541, 0.000543, 0.000546, 0.00055, 0.000554, 0.000558, 0.000564, 0.00057, 0.000577, 0.000586, 0.000596, 0.000608, 0.000622, 0.000639, 0.000659, 0.000682, 0.00071, 0.000743, 0.000783, 0.000831, 0.00089, 0.00096, 0.001047, 0.001151, 0.001277, 0.00143188, 0.00161982, 0.00184977, 0.00212764, 0.00246493, 0.00287334, 0.00336387, 0.00395196, 0.00465231, 0.00548274, 0.00645877, 0.00759594, 0.00891242, 0.01042064, 0.01213156, 0.01405626, 0.01619691, 0.01855486, 0.02111986, 0.02388407, 0.02683271, 0.02994446, 0.0331681, 0.03648963, 0.03997779, 0.04362967, 0.04742653, 0.0514231, 0.05555392, 0.05984118, 0.06431133, 0.06893335, 0.07373719, 0.07871069, 0.08386445, 0.08920963, 0.09475032, 0.10048622, 0.10642222, 0.11260481, 0.11897842, 0.12560509, 0.13243893, 0.13958092, 0.14691087, 0.15457267, 0.16252126, 0.17072093, 0.17920463};

  G4double Ar_percent[] = {0.934394, 0.934225, 0.934575, 0.93418, 0.934607, 0.93436, 0.934054, 0.934544, 0.934216, 0.934356, 0.934271, 0.934331, 0.934277, 0.934312, 0.934191, 0.934174, 0.934198, 0.93427, 0.93429, 0.934519, 0.934824, 0.934485, 0.934379, 0.934389, 0.934353, 0.934255, 0.934223, 0.934561, 0.934455, 0.934121, 0.934217, 0.934581, 0.934542, 0.934336, 0.933967, 0.934614, 0.934001, 0.934306, 0.93434, 0.934388, 0.93427, 0.934285, 0.934206, 0.934221, 0.93413, 0.934317, 0.934196, 0.934278, 0.934313, 0.934469, 0.934618, 0.93466, 0.934172, 0.93426, 0.934375, 0.934407, 0.934306, 0.934436, 0.934328, 0.934281, 0.934197, 0.93428, 0.934469, 0.934321, 0.934006, 0.934104, 0.933876, 0.934052, 0.933994, 0.934052, 0.934388, 0.934032, 0.934008, 0.933863, 0.933854, 0.93343, 0.933314, 0.932866, 0.932428, 0.932086, 0.931792, 0.930971, 0.929989, 0.929361, 0.928117, 0.926709, 0.925501, 0.923255, 0.921016, 0.918275, 0.91516, 0.911241, 0.906605, 0.901256, 0.89515, 0.887415, 0.879195, 0.868841, 0.857751, 0.844363, 0.828944, 0.81208734, 0.79357742, 0.77323207, 0.75139037, 0.72802419, 0.7040104, 0.67904882, 0.65412057, 0.62864605, 0.60323545, 0.57841669, 0.55416133, 0.53076122, 0.50840606, 0.48692146, 0.46660256, 0.44760689, 0.42964559, 0.41268011, 0.39727287, 0.38259086, 0.36920489, 0.35669676, 0.34554115, 0.33476546, 0.32463147, 0.31478113, 0.30576767, 0.29696254, 0.28855085, 0.28061798, 0.27301233, 0.26568971, 0.25873849, 0.2520665, 0.2455467, 0.23945325, 0.23352683, 0.22776004, 0.22233285, 0.21691019, 0.21177931, 0.20678048, 0.20203065, 0.19735616, 0.19288139, 0.18856219, 0.18430296, 0.18019048};

 G4double H_percent[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2e-06, 8e-06, 2.1e-05, 3.8e-05, 5.6e-05, 7.1e-05, 8.2e-05, 9.1e-05, 9.7e-05, 0.000103, 0.000107, 0.000112, 0.000116, 0.000122, 0.000128, 0.000135, 0.000143, 0.000154, 0.000166, 0.00018, 0.000197, 0.000218, 0.00024128, 0.00026894, 0.00030081, 0.00033748, 0.00037944, 0.00042671, 0.00048012, 0.00053964, 0.0006053, 0.00067783, 0.00075631, 0.00084123, 0.00093222, 0.0010287, 0.00113019, 0.00123595, 0.00134524, 0.0014568, 0.00156984, 0.00168279, 0.00179468, 0.00190464, 0.00200966, 0.00210878, 0.00220812, 0.00230686, 0.00240365, 0.00250188, 0.00259809, 0.00269353, 0.00278947, 0.00288432, 0.00297872, 0.00307352, 0.00316732, 0.00326211, 0.00335733, 0.00345198, 0.00354677, 0.00364527, 0.00374089, 0.00383963, 0.00393927, 0.00404276, 0.00414525, 0.00425348, 0.00436064, 0.00447003, 0.00458235};


G4double tempAtmos[] = {272.1, 268.9, 263.6, 256.9, 249.6, 242.1, 234.8, 228.3, 222.6, 218.1, 214.9, 213.2, 212.4, 212.3, 212.4, 212.2, 211.5, 210.4, 209.1, 207.8, 206.7, 206.0, 205.8, 205.9, 206.3, 207.2, 208.4, 209.9, 211.7, 213.9, 216.5, 219.4, 222.7, 226.2, 229.9, 233.6, 237.2, 240.8, 244.2, 247.5, 250.6, 253.3, 255.8, 257.8, 259.4, 260.5, 261.0, 261.0, 260.6, 259.7, 258.5, 257.0, 255.2, 253.2, 251.2, 249.0, 246.8, 244.6, 242.4, 240.2, 238.1, 236.0, 234.0, 232.1, 230.2, 228.5, 226.9, 225.4, 224.0, 222.8, 221.7, 220.8, 220.0, 219.4, 218.6, 217.7, 216.7, 215.5, 214.2, 212.8, 211.3, 209.7, 208.0, 206.3, 204.5, 202.7, 200.8, 198.8, 196.9, 195.0, 193.0, 191.0, 189.2, 187.4, 185.9, 184.6, 183.7, 183.1, 183.0, 183.3, 184.3, 185.8, 188.1, 191.0, 194.5, 198.8, 203.7, 209.5, 216.0, 223.4, 231.7, 241.0, 251.4, 263.0, 275.6, 289.3, 304.2, 320.0, 336.7, 354.0, 371.5, 388.6, 404.7, 419.0, 429.5, 440.3, 450.6, 460.6, 470.2, 479.5, 488.5, 497.1, 505.4, 513.5, 521.2, 528.7, 535.9, 542.8, 549.5, 556.0, 562.2, 568.2, 574.0, 579.6, 585.0, 590.2, 595.2, 600.0, 604.7, 609.2};

G4String volName[] = {"air0", "air1", "air2", "air3", "air4", "air5", "air6", "air7", "air8", "air9", "air10", "air11", "air12", "air13", "air14", "air15", "air16", "air17", "air18", "air19", "air20", "air21", "air22", "air23", "air24", "air25", "air26", "air27", "air28", "air29", "air30", "air31", "air32", "air33", "air34", "air35", "air36", "air37", "air38", "air39", "air40", "air41", "air42", "air43", "air44", "air45", "air46", "air47", "air48", "air49", "air50", "air51", "air52", "air53", "air54", "air55", "air56", "air57", "air58", "air59", "air60", "air61", "air62", "air63", "air64", "air65", "air66", "air67", "air68", "air69", "air70", "air71", "air72", "air73", "air74", "air75", "air76", "air77", "air78", "air79", "air80", "air81", "air82", "air83", "air84", "air85", "air86", "air87", "air88", "air89", "air90", "air91", "air92", "air93", "air94", "air95", "air96", "air97", "air98", "air99", "air100", "air101", "air102", "air103", "air104", "air105", "air106", "air107", "air108", "air109", "air110", "air111", "air112", "air113", "air114", "air115", "air116", "air117", "air118", "air119", "air120", "air121", "air122", "air123", "air124", "air125", "air126", "air127", "air128", "air129", "air130", "air131", "air132", "air133", "air134", "air135", "air136", "air137", "air138", "air139", "air140", "air141", "air142", "air143", "air144", "air145", "air146", "air147", "air148", "air149"};

//G4cout << "did it reach here before the material? " << G4endl; 

  for (G4int i=0; i<150; i++)
    {     
        airMat[i] = new G4Material(name=volName[i], densityAtmos[i]*g/cm3, ncomponents=5, kStateGas, tempAtmos[i]*kelvin);
	//G4cout << "inside the loop " << i << G4endl;  
        airMat[i]->AddElement(elN, N2_percent[i]*perCent);
        airMat[i]->AddElement(elO, O2_percent[i]*perCent);
	airMat[i]->AddElement(elHe, He_percent[i]*perCent);
	airMat[i]->AddElement(elAr, Ar_percent[i]*perCent);
	airMat[i]->AddElement(elH, H_percent[i]*perCent);
    }

//G4cout << "were the materials defined?" << G4endl; 
 
  	G4int k = 149;
	atmos[k] = new G4Orb("atmos", (6371.2+k+1)*km);
        logAtmos[k] = new G4LogicalVolume(atmos[k], airMat[k], "atmos");
        atmosUserLimits[k] = new G4UserLimits(100*m);
        logAtmos[k]->SetUserLimits(atmosUserLimits[k]);
        atmosUserLimits[k]->SetUserMaxTrackLength(200*re);
        physAtmos[k] = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                        logAtmos[k], "atmos", upperLogAtmos[0], false, 0);

	
  for (G4int i = 148; i>=0; i--)
    {
	atmos[i] = new G4Orb("atmos", (6371.2+i+1)*km);
	logAtmos[i] = new G4LogicalVolume(atmos[i], airMat[i], "atmos");
	atmosUserLimits[i] = new G4UserLimits(100*m);
	logAtmos[i]->SetUserLimits(atmosUserLimits[i]);
	atmosUserLimits[i]->SetUserMaxTrackLength(200*re);
	physAtmos[i] = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
			logAtmos[i], "atmos", logAtmos[i+1], false, 0); 
	
    }	
  //-----------------layers of atmosphere-----------------------//

  /*G4int natoms;  
  density = 1.225*mg/cm3;
  G4Material* Air1 = new G4Material(name="Air1", density, ncomponents=3);
  Air1->AddElement(elN, natoms = 2, 78.0*perCent);
  Air1->AddElement(elO, natoms = 2, 21.0*perCent);
  Air1->AddElement(elAr, 1.0*perCent);

  density = 0.2*mg/cm3;
  G4Material* Air2 = new G4Material(name="Air2", density, ncomponents=3);
  Air2->AddElement(elN, 78.0*perCent);
  Air2->AddElement(elO, 21.0*perCent);
  Air2->AddElement(elAr, 1.0*perCent);


  density = 0.001*mg/cm3;
  G4Material* Air3 = new G4Material(name="Air3", density, ncomponents=3);
  Air3->AddElement(elN, 78.0*perCent);
  Air3->AddElement(elO, 21.0*perCent);
  Air3->AddElement(elAr, 1.0*perCent);


  density = 1E-06*mg/cm3;
  G4Material* Air4 = new G4Material(name="Air4", density, ncomponents=3);
  Air4->AddElement(elN, 78.0*perCent);
  Air4->AddElement(elO, 21.0*perCent);
  Air4->AddElement(elAr, 1.0*perCent);*/


  //G4double threshold = 1.0*GeV;

  /*G4Orb* thermosphere = new G4Orb("thermosphere", 6600*km);
  G4LogicalVolume* logthermo = new G4LogicalVolume(thermosphere, Air4, "thermosphere");
  G4UserLimits* thermoUserLimits = new G4UserLimits(100*m);
  logthermo->SetUserLimits(thermoUserLimits);
  //thermoUserLimits->SetUserMinEkine(threshold);
  thermoUserLimits->SetUserMaxTrackLength(200*re);
  G4VPhysicalVolume* fPhysThermo = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                  logthermo,  "thermosphere", logSpace, false, 0);


  G4Orb* mesosphere = new G4Orb("mesosphere", 6463*km);
  G4LogicalVolume* logmeso = new G4LogicalVolume(mesosphere, Air3, "mesosphere");
  G4UserLimits* mesoUserLimits = new G4UserLimits(100*m);
  logmeso->SetUserLimits(mesoUserLimits);
  //mesoUserLimits->SetUserMinEkine(threshold);
  mesoUserLimits->SetUserMaxTrackLength(200*re);
  G4VPhysicalVolume* fPhysMeso = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                  logmeso,  "mesosphere", logthermo, false, 0);


  G4Orb* stratosphere = new G4Orb("stratosphere", 6429*km);
  G4LogicalVolume* logstrato = new G4LogicalVolume(stratosphere, Air2, "stratosphere");
  G4UserLimits* stratoUserLimits = new G4UserLimits(100*m);
  logstrato->SetUserLimits(stratoUserLimits);
  //stratoUserLimits->SetUserMinEkine(threshold);
  stratoUserLimits->SetUserMaxTrackLength(200*re);
  G4VPhysicalVolume* fPhysStrato = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                  logstrato,  "stratosphere", logmeso, false, 0);

  G4Orb* troposphere = new G4Orb("troposphere", 6390*km);
  G4LogicalVolume* logtropo = new G4LogicalVolume(troposphere, Air1, "troposphere");
  G4UserLimits* tropoUserLimits = new G4UserLimits(100*m);
  logtropo->SetUserLimits(tropoUserLimits);
  //tropoUserLimits->SetUserMinEkine(threshold);
  tropoUserLimits->SetUserMaxTrackLength(200*re);
  G4VPhysicalVolume* fPhysTropo = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                  logtropo,  "troposphere", logstrato, false, 0);*/

  /*G4Orb* sandstone_orb = new G4Orb("Sandstone", 6371.00*km);
  G4LogicalVolume* sandLV = new G4LogicalVolume(sandstone_orb, sandstone, "Sandstone");
  G4UserLimits* sandUserLimits = new G4UserLimits(0.01*m);
  sandLV->SetUserLimits(sandUserLimits);
  G4VPhysicalVolume* sandPV = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), 
		sandLV, "Sandstone", logAtmos[0], false, 0);*/

  /*G4Orb* hematite_orb = new G4Orb("Hematite", 6371.000*km);
  G4LogicalVolume* hemaLV = new G4LogicalVolume(hematite_orb, hematite, "Hematite");
  G4UserLimits* hemaUserLimits = new G4UserLimits(10*m);
  hemaLV->SetUserLimits(hemaUserLimits);
  G4VPhysicalVolume* hemaPV = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                hemaLV, "Hematite", logAtmos[0], false, 0);*/


  /*G4Orb* Al_orb = new G4Orb("aluminum", 6371.2*km);
  G4LogicalVolume* aluLV = new G4LogicalVolume(Al_orb, Aluminum, "aluminum");
  G4UserLimits* aluUserLimits = new G4UserLimits(10*m);
  aluLV->SetUserLimits(aluUserLimits);
  G4VPhysicalVolume* aluPV = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                aluLV, "aluminum", logAtmos[0], false, 0);
*/

  /*G4Orb* iceLayer = new G4Orb("IceLayer", 6371.2*km); //6370.900*km);
  //G4LogicalVolume* iceLV = new G4LogicalVolume(iceLayer, Water, "IceLayer");
  G4LogicalVolume* iceLV = new G4LogicalVolume(iceLayer, ice, "IceLayer");  
  G4UserLimits* iceUserLimits = new G4UserLimits(100*m);
  iceLV->SetUserLimits(iceUserLimits);
  G4VPhysicalVolume* icePV = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), 
		  iceLV, "IceLayer", logAtmos[0], false, 0);

*/


  /*G4Orb* Atmos = new G4Orb("Atmosphere", 6600*km);
  //G4Orb* Atmos = new G4Orb("Atmosphere", 6480*km);
  G4LogicalVolume* logAtmos = new G4LogicalVolume(Atmos, Air, "Atmosphere");
  G4UserLimits* theAtmosUserLimits=new G4UserLimits(100*m);//100*km);
  logAtmos->SetUserLimits(theAtmosUserLimits);
  theAtmosUserLimits->SetUserMaxTrackLength(2*re);

  G4VPhysicalVolume* fPhysAtmos = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                  logAtmos, "Atmosphere", logSpace, false, 0);

  G4VisAttributes * VisAtmos = new G4VisAttributes(G4Colour(1,0,0.));
  VisAtmos->SetVisibility(true);
  logAtmos->SetVisAttributes(VisAtmos);*/

 
  //----------------------
  //Earth
  //-----------------------


  /*G4VisAttributes * VisAttRock = new G4VisAttributes();
  VisAttRock->SetVisibility(true);

  G4double radRock = 1.4*km;
  G4Orb* rock = new G4Orb("Rock", 6371.2*km);
  G4LogicalVolume* logRock = new G4LogicalVolume(rock, limestone, "Rock");
  G4UserLimits* theRockUserLimits=new G4UserLimits(100*m);//100*km);
  logRock->SetUserLimits(theRockUserLimits);
  //theRockUserLimits->SetUserMaxTrackLength(200*re);
  logRock->SetUserLimits(theRockUserLimits);
  logRock->SetVisAttributes(VisAttRock);
  
  G4VPhysicalVolume* fPhysRock = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                        logRock, "Rock", logAtmos[0], false, 0);
*/
   
  /* solidEarth=new G4Sphere("Earth",0.
                      ,6378.16*km+top_altitude,0.,360.*deg,0.,180.*deg);*/
   if (!RemoveEarth){		      
   	solidEarth=new G4Orb("Earth",6371.2*km/*+TopOfAtmosphere*/);
	logicEarth=new G4LogicalVolume(solidEarth,Vacuum,"Earth",0,0,0);
   	theEarthUserLimits=new G4UserLimits(1.*km);
   	logicEarth->SetUserLimits(theEarthUserLimits);
	theEarthUserLimits->SetUserMaxTrackLength(200.*re);
   	logicEarth->SetVisAttributes(VisAttEarth);
   	/*physiEarth=new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(0.,0.,0.), // at (0,0,0)
				 "Earth",       // its name
                                 logicEarth,      // its logical volume
                                 fPhysAtmos,               // its mother  volume
                                 false,           // no boolean operations
                                 0);  */

	physiEarth=new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(0.,0.,0.), // at (0,0,0)
                                // "Earth",       // its name
                                 logicEarth,      // its logical volume
                                 "Earth",
                                 logAtmos[0], //iceLV   // its mother  volume
                                 false,           // no boolean operations
                                 0);

				 
  /*G4VisAttributes * VisAttCubeSats = new G4VisAttributes(G4Colour(1.,0.,0.));
  VisAttCubeSats->SetVisibility(true);  

  G4Material* Glass = new G4Material(name="Glass", density, ncomponents=2);
  Glass->AddElement(elC, 91.533 * perCent);
  Glass->AddElement(elH, 8.467 * perCent);

  G4double fCubeSatsSize = 100*km;   //cubesats thickness
  G4Box* SolidCubeSats = new G4Box("CubeSats", 50*km, 50*km, fCubeSatsSize/2);
  G4LogicalVolume* fLogCubeSats = new G4LogicalVolume(SolidCubeSats, Glass, "Glass");
  G4ThreeVector cubesatspos = G4ThreeVector(0, 0, 6378.137*km+TopOfAtmosphere+100*km);
  fPhysCubeSats = new G4PVPlacement(0, cubesatspos, fLogCubeSats, "Glass", logicWorld, false, 0);
  fLogCubeSats->SetVisAttributes(VisAttCubeSats);*/
	  
  
  //G4ProductionCutsTable* theCutsTable = G4ProductionCutsTable::GetProductionCutsTable();
  /*G4Region* region1 = new G4Region("Atmospheric1");
  G4Region* region2 = new G4Region("Atmospheric2");
  G4Region* region3 = new G4Region("Atmospheric3");
  G4Region* region4 = new G4Region("Atmospheric4");


  region1->AddRootLogicalVolume(logthermo);
  region2->AddRootLogicalVolume(logmeso);
  region3->AddRootLogicalVolume(logstrato);
  region4->AddRootLogicalVolume(logtropo);*/
  
  /*G4double cutValue1 = 100*m;
  G4double cutValue2 = 100*m;
  G4double cutValue3 = 100*m;
  G4double cutValue4 = 100*m;

  //G4double lowEnergy = 1*eV;
  //G4double highEnergy = 1*GeV;
  //theCutsTable->SetEnergyRange(lowEnergy, highEnergy);

  G4ProductionCuts* cuts1 = new G4ProductionCuts();  
  G4ProductionCuts* cuts2 = new G4ProductionCuts();
  G4ProductionCuts* cuts3 = new G4ProductionCuts();
  G4ProductionCuts* cuts4 = new G4ProductionCuts();


  cuts1->SetProductionCut(cutValue1, "e-");
  cuts2->SetProductionCut(cutValue2, "e-");
  cuts3->SetProductionCut(cutValue3, "e-");
  cuts4->SetProductionCut(cutValue4, "e-");

  //cuts->SetProductionCut(lowEnergy, "e-");
  //cuts->SetProductionCut(highEnergy, "e-");
  //theCutsTable->SetEnergyRange(lowEnergy, highEnergy);

  region1->SetProductionCuts(cuts1);
  region2->SetProductionCuts(cuts2);
  region3->SetProductionCuts(cuts3);
  region4->SetProductionCuts(cuts4);*/

  /*G4double cutValue = 0.001*mm;
  G4double lowEnergy = 1*eV;
  G4double highEnergy = 1*GeV;
  theCutsTable->SetEnergyRange(lowEnergy, highEnergy);
  G4ProductionCuts* cuts = new G4ProductionCuts();
  cuts->SetProductionCut(cutValue, "e-");
  region1->SetProductionCuts(cuts);
  region2->SetProductionCuts(cuts);
  region3->SetProductionCuts(cuts);
  region4->SetProductionCuts(cuts);*/




  }

 return physiWorld;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSGeometryConstruction::SetStopAltitude(G4double stop_altitude)
{
 if (stop_altitude <0.) {
 	G4cout<<"The top altitude below which particle are stopped ";
	G4cout<<"should be >0."<<std::endl; 
 	return;
 }
 else {
 	TopOfAtmosphere = stop_altitude;
	G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
	theSteppingAction->SetStopAltitude(stop_altitude);
	
 }	 
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSGeometryConstruction::SetRemoveEarth(G4bool aBool)
{RemoveEarth = aBool;
 G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSGeometryConstruction::SetMaxStepLength(G4double maxStep)
{ theWorldUserLimits->SetMaxAllowedStep(maxStep);
       maxStepSize=maxStep;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSGeometryConstruction::SetMaxTrackLength
                                                   (G4double maxTrackLength)
{ theWorldUserLimits->SetUserMaxTrackLength(maxTrackLength);
  theEarthUserLimits->SetUserMaxTrackLength(maxTrackLength);
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSGeometryConstruction::SetMaxTimeOfFlight
                                                   (G4double maxTime)
{ theWorldUserLimits->SetUserMaxTime(maxTime);
  theEarthUserLimits->SetUserMaxTime(maxTime);
}

//////////////////////////////////////////////////////////////////


void MAGCOSGeometryConstruction::SetCutsForVolumes(G4ProductionCuts* productionCuts, const G4String& particleName, const G4String& volumeName)
{
	G4LogicalVolume* volume = G4LogicalVolumeStore::GetInstance()->GetVolume(volumeName);
	G4cout << "does it even call this cut function" << G4endl;
	if (volume){
		G4Region* region = new G4Region(volumeName);
		region->AddRootLogicalVolume(volume);
		G4ProductionCuts* cuts = new G4ProductionCuts();
		cuts->SetProductionCut(productionCuts->GetProductionCut(particleName));
		region->SetProductionCuts(cuts);
	}
}


//###########################################################################//
///*void MAGCOSGeometryConstruction::ConstructSDandField()

//  // Define magnetic field
//    if ( ! fMonFieldSetup.Get() ) {
//        G4MonopoleFieldSetup* fieldSetup
//              = new G4MonopoleFieldSetup();
//                  G4AutoDelete::Register(fieldSetup); // Kernel will delete the F01FieldSetup
//                      fMonFieldSetup.Put(fieldSetup);
//                        }
//                          fMonFieldSetup.Get()->ConstructMagField(); // add field value
//                          }*/
//                          //##########################################################################//
//

