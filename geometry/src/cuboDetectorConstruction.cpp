#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4Transform3D.hh"
#include "G4Box.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4GenericMessenger.hh"

#include "../include/cuboDetectorConstruction.h"
#include "geometry/pluginmanagers/GGSGeoPluginMacros.h"

GeometryPlugin(cuboDetectorConstruction);
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
cuboDetectorConstruction::cuboDetectorConstruction()
{
 _messenger = new G4GenericMessenger(this, "/GGS/geometry/cuboDetectorConstruction/");
 _messenger->DeclareProperty("scintillationYield", _scintillationYield, "Number of photons per MeV");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

cuboDetectorConstruction::~cuboDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* cuboDetectorConstruction::Construct()
{

//      ------------- Dimension -------------
  G4double abs_width = 0;
  G4double _World_Size_X = 10*cm;
  G4double _World_Size_Y = _World_Size_X;
  G4double _World_Size_Z = 10*cm;
  G4double _Cube_Size = 3.6*cm;
  G4double _Absorber_Size_X = 9*cm;
  G4double _Absorber_Size_Y = _Absorber_Size_X;
  G4double _Absorber_Size_Z = abs_width*1.86*cm;
  G4double abs_cub_gap = 0.3*cm;
  G4double _Diode_Size_X = 9.2*mm;
  G4double _Diode_Size_Y = 0.125*mm;
  G4double _Diode_Size_Z = 9.2*mm;
  G4double _Resin_Size_X = 15*mm;
  G4double _Resin_Size_Y = 0.1*mm;
  G4double _Resin_Size_Z = 15*mm;
  G4double Case_Size_X = 15.*mm;
  G4double Case_Size_Y = 2.*mm;
  G4double Case_Size_Z = 15.*mm;

//	------------- Materials -------------

  G4double a, z, density;
  G4int nelements;

// Air
// 
  G4Element * N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element * O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  G4Element * Ar = new G4Element("Argon"  , "Ar", z=18 , a=39.948*g/mole);


  G4Material * Air = new G4Material("Air", density=1.29*mg/cm3, nelements=3);
  Air->AddElement(N, 78.*perCent);
  Air->AddElement(O, 21.*perCent);
  Air->AddElement(Ar, 1.*perCent);

// CsI
// 
  G4Element * Cs = new G4Element("Cesium", "Cs", z=55 , a=132.9054519*g/mole);
  G4Element * I = new G4Element("Iodine", "I", z=53 , a=126.90447*g/mole);
  G4Material * CsImaterial = new G4Material("CsImaterial", density= 4.51*g/cm3, 
			       nelements=2);
  CsImaterial->AddElement(Cs, 1);
  CsImaterial->AddElement(I, 1);

// Fake CsI
// 
  G4Material * FakeCsImaterial = new G4Material("FakeCsImaterial", density= 4.51*g/cm3, 
				   nelements=2);
  FakeCsImaterial->AddElement(Cs, 1);
  FakeCsImaterial->AddElement(I, 1);

// Silicon
// 
  G4Material * SiliconMaterial = new G4Material("SiliconMaterial", z=14, 
					       a=28.09*g/mole, 
					       density=2.329*g/cm3);

// Optical Resin
//
  G4Material * ResinMaterial = new G4Material("ResinMaterial", z=6, a=12*g/mole, 
					     density=1*g/cm3);

// Ceramic
//
  G4Material * CeramicMaterial = new G4Material("ResinMaterial", z=6, a=12*g/mole, 
					     density=2*g/cm3);

//
// ------------ Generate & Add Material Properties Table ------------
//
  const G4int nEntries = 9;
  const G4int nEntries_2 = 12;
  G4double PhotonEnergy[nEntries] =
    { 1.65*eV, 1.77*eV, 1.91*eV, 2.07*eV, 2.25*eV, 2.48*eV, 2.76*eV, 3.10*eV,
      3.87*eV };
  //{  750*nm,  700*nm,  650*nm,  600*nm,  550*nm,  500*nm,  450*nm,  400*nm,
  //   320*nm }

  G4double PhotonEnergy_2[nEntries_2] =
    { 1.65*eV, 1.77*eV, 1.91*eV, 2.07*eV, 2.25*eV, 2.48*eV, 2.76*eV, 3.10*eV,
      3.54*eV, 3.87*eV, 4.13*eV, 4.28*eV};
  //{  750*nm,  700*nm,  650*nm,  600*nm,  550*nm,  500*nm,  450*nm,  400*nm,
  //   350*nm,  320*nm,  300*nm,  290*nm}

//
// CsI
//	      
  G4double RefractiveIndex_CsI[nEntries_2] =
    {   1.769,   1.773,   1.778,   1.785,   1.794,   1.806,   1.824,   1.850, 
        1.894,   1.937,   1.979,   2.006 };

  G4double Absorption_CsI[nEntries] =
    { 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 
      0.*cm };

  //G4double Absorption_CsI[nEntries] =
  //  { 463.9*cm, 431.2*cm, 422.6*cm, 452.2*cm, 466.4*cm, 462.1*cm, 516.9*cm, 
  //    751.7*cm, 0.2*cm  };

  G4double ScintilFast[nEntries] =
    { 0.00, 0.18, 0.39, 0.70, 1.0, 0.75, 0.28, 0.12, 0.00 };

  G4MaterialPropertiesTable * csiMPT = new G4MaterialPropertiesTable();

  csiMPT->AddProperty("RINDEX", PhotonEnergy_2, RefractiveIndex_CsI, nEntries_2)
        ->SetSpline(true);
  csiMPT->AddProperty("ABSLENGTH", PhotonEnergy, Absorption_CsI, nEntries)
        ->SetSpline(true);
  csiMPT->AddProperty("FASTCOMPONENT",PhotonEnergy, ScintilFast, nEntries)
        ->SetSpline(true);

  //csiMPT->AddConstProperty("SCINTILLATIONYIELD", 54000./MeV);
  //csiMPT->AddConstProperty("SCINTILLATIONYIELD", 0./MeV);
  csiMPT->AddConstProperty("SCINTILLATIONYIELD", 54./MeV);
  csiMPT->AddConstProperty("RESOLUTIONSCALE", 1);
  csiMPT->AddConstProperty("FASTTIMECONSTANT", 1000.*ns);
  csiMPT->AddConstProperty("YIELDRATIO", 1.);

  CsImaterial->SetMaterialPropertiesTable(csiMPT);

  // Set the Birks Constant for the scintillator
  //CsImaterial->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

//
// Air
//
  G4double RefractiveIndex_Air[nEntries] =
    { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 };

  G4MaterialPropertiesTable* airMPT = new G4MaterialPropertiesTable();
  airMPT->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_Air, nEntries);
  
  Air->SetMaterialPropertiesTable(airMPT);

//
// Optical Resin
//
  G4double RefractiveIndex_Resin[nEntries] =
    { 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55 };

  G4MaterialPropertiesTable * resinMPT = new G4MaterialPropertiesTable();
  resinMPT->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_Resin, nEntries);
  
  ResinMaterial->SetMaterialPropertiesTable(resinMPT);

//
// Plexiglass
//
//   G4double RefractiveIndex_Plexiglass[nEntries] =
//     { 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49 };

//   plexiMPT = new G4MaterialPropertiesTable();
//   plexiMPT->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_Plexiglass, 
// 			nEntries);
//   PlexiglassMaterial->SetMaterialPropertiesTable(plexiMPT);

//
// Silicon
//
  G4double RefractiveIndex_Si[nEntries] =
    { 3.733, 3.783, 3.851, 3.947, 4.084, 4.297, 4.674, 5.57, 5.023 };

  G4double Absorption_Si[nEntries] =
    { 
      7.0e-3*mm,
      5.0e-3*mm,
      4.0e-3*mm,
      2.4e-3*mm,
      1.5e-3*mm,
      0.8e-3*mm,
      0.3e-3*mm,
      0.1e-3*mm,
      0.08e-3*mm 
    };

  G4MaterialPropertiesTable * SiMPT = new G4MaterialPropertiesTable();
  SiMPT->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_Si, nEntries)
       ->SetSpline(true);
  SiMPT->AddProperty("ABSLENGTH", PhotonEnergy, Absorption_Si, nEntries);

  SiliconMaterial->SetMaterialPropertiesTable(SiMPT);


//
//	------------- Volumes --------------

// The experimental Hall
//
  G4Box * _Solid_World = new G4Box("World", _World_Size_X / 2., 
			   _World_Size_Y / 2., _World_Size_Z / 2.);

  G4LogicalVolume * _Logic_World = new G4LogicalVolume(_Solid_World, Air, "World");

  G4PVPlacement * _Physical_World
    = new G4PVPlacement(0, G4ThreeVector(), _Logic_World ,"World", 0, false, 0);

// experimental Hall - PD Side
  G4Box * _Solid_PDSideWorld = new G4Box("PDSideWorld", _World_Size_X / 2., 
				 (_World_Size_Y - _Cube_Size) / 4.,
				 _World_Size_Z / 2.);

  G4LogicalVolume * _Logic_PDSideWorld
    = new G4LogicalVolume(_Solid_PDSideWorld, Air, "PDSideWorld");

  G4PVPlacement * _Physical_PDSideWorld
    = new G4PVPlacement(0, G4ThreeVector(0., (_World_Size_Y+_Cube_Size)/4., 0.), 
			_Logic_PDSideWorld ,"PDSideWorld", 
			_Logic_World, false, 0);

  
// The CsI cube
//	
  G4Box * _Solid_Cube = new G4Box("CsI_Cube", _Cube_Size / 2., 
				 _Cube_Size / 2., _Cube_Size / 2.);

  G4LogicalVolume * _Logic_Cube  = new G4LogicalVolume(_Solid_Cube, CsImaterial, 
						      "CsI_Cube");

  G4double cube_Z = -(_Absorber_Size_Z + _Cube_Size/2. + abs_cub_gap);

  G4PVPlacement * _Physical_Cube = 
    new G4PVPlacement(0,G4ThreeVector(
				      0., 
				      0., 
				      cube_Z
				      ),
		      _Logic_Cube , "CsI_Cube", _Logic_World, false, 0);
  
// Optical Resin

  G4Box * _Solid_Resin = new G4Box("Optical_Resin", _Resin_Size_X / 2., 
			   _Resin_Size_Y / 2., _Resin_Size_Z / 2.);

  G4LogicalVolume * _Logic_Resin  = new G4LogicalVolume(_Solid_Resin, 
				      ResinMaterial, 
				      "Optical_Resin");

  G4PVPlacement * _Physical_Resin = 
    new G4PVPlacement(0, G4ThreeVector(
				       0., 
				       (_Cube_Size - _World_Size_Y)/4. + _Resin_Size_Y /2., 
				       cube_Z
                                       //-Case_Size_Z/2
				       ), 
		      _Logic_Resin , "Optical_Resin", _Logic_PDSideWorld, 
		      false, 0);

// Photodiode case
//
  G4Box * _Solid_Case = new G4Box("Photodiode_case", Case_Size_X / 2., Case_Size_Y / 2., 
	 		  Case_Size_Z / 2.);

  G4LogicalVolume * _Logic_Case  = new G4LogicalVolume(_Solid_Case, 
				     CeramicMaterial, 
				     "Photodiode_case");
  G4PVPlacement * _Physical_Case = 
    new G4PVPlacement(0, G4ThreeVector(
				       0.,
				       _Cube_Size/4. - _World_Size_Y/4. 
				       + _Resin_Size_Y + Case_Size_Y / 2.,
				       //-Case_Size_Z/2
                                       cube_Z
				       ),
		      _Logic_Case , "Photodiode_case", _Logic_PDSideWorld, 
		      false, 0
		      );

// The Photodiode
//
  G4Box * _Solid_Diode = new G4Box("Photodiode", _Diode_Size_X / 2., 
			   _Diode_Size_Y / 2., _Diode_Size_Z / 2.);

  G4LogicalVolume * _Logic_Diode  = new G4LogicalVolume(_Solid_Diode, 
				      SiliconMaterial, 
				      "Photodiode");
  G4PVPlacement * _Physical_Diode = 
    new G4PVPlacement(0, G4ThreeVector(0.,
				       - (Case_Size_Y - _Diode_Size_Y) / 2.,
				       0.
				       ),
		      _Logic_Diode , "Photodiode", _Logic_Case, false, 0);


// Defining Sensitive Detector
  //G4String SDname, SDname_Edep, SDname_Diode;
  //sd = new cuboSD(SDname = "cuboSD");
  //sd_Edep = new cuboSD_Edep(SDname_Edep = "cuboSD_Edep");
  //sd_Diode = new cuboSD_Diode(SDname_Diode = "cuboSD_Diode");

 //  The SD has to be registered to both SDManager and LogicalVolume.
  //G4SDManager* SDman = G4SDManager::GetSDMpointer();
  //SDman->AddNewDetector(sd);
  //SDman->AddNewDetector(sd_Edep);
  //SDman->AddNewDetector(sd_Diode);
  //_Logic_Resin->SetSensitiveDetector(sd);
  //_Logic_Cube->SetSensitiveDetector(sd_Edep);
  //_Logic_Diode->SetSensitiveDetector(sd_Diode);

//	------------- Surfaces --------------
//

  // CsI - Air Surface
  G4OpticalSurface * Op_CsI_Surface = new G4OpticalSurface("CsI_Surface");
  Op_CsI_Surface->SetType(dielectric_dielectric);
  Op_CsI_Surface->SetFinish(polishedfrontpainted);
  //Op_CsI_Surface->SetFinish(polishedbackpainted);
  Op_CsI_Surface->SetModel(unified);
  //Op_CsI_Surface->SetPolish(1.);
  //Op_CsI_Surface->SetSigmaAlpha(0.3);

  G4LogicalBorderSurface * CsI_Air_Surf = 
    new G4LogicalBorderSurface("CsI_Surface1", _Physical_Cube, _Physical_World,
			       Op_CsI_Surface);
/*
  G4LogicalBorderSurface * Air_CsI_Surf = 
    new G4LogicalBorderSurface("CsI_Surface2", _Physical_World, _Physical_Cube,
			       Op_CsI_Surface);
*/                               

  const G4int num = 2;
  G4double Ephoton[num] = {1.65*eV, 4.28*eV};

  G4double Reflectivity[num] = {0., 0.};
  G4double Efficiency[num]   = {0., 0.};
  G4double cubeReflectivity = 1.;

  G4double specularlobe[num] = {1., 1.};
  G4double specularspike[num] = {0., 0.};
  G4double backscatter[num] = {0., 0.};

  G4double backpaint[num] = {1., 1.}; // refractive index

  G4MaterialPropertiesTable * csiSPT = new G4MaterialPropertiesTable();

  csiSPT->AddProperty("REFLECTIVITY", Ephoton, Reflectivity, num);
  csiSPT->AddProperty("EFFICIENCY",   Ephoton, Efficiency,   num);
  csiSPT->AddProperty("SPECULARLOBECONSTANT", Ephoton, specularlobe, num);
  csiSPT->AddProperty("SPECULARSPIKECONSTANT", Ephoton, specularspike, num);
  csiSPT->AddProperty("BACKSCATTERCONSTANT", Ephoton, backscatter, num);
  csiSPT->AddProperty("RINDEX", Ephoton, backpaint, num);

  Op_CsI_Surface->SetMaterialPropertiesTable(csiSPT);

  // CsI - Air Surface (Diode Side)

  G4OpticalSurface * Op_CsI_Surface_DSide = new G4OpticalSurface("CsI_Surface_DSide");
  Op_CsI_Surface_DSide->SetType(dielectric_dielectric);
  Op_CsI_Surface_DSide->SetFinish(polishedfrontpainted);
  Op_CsI_Surface_DSide->SetModel(unified);
  Op_CsI_Surface_DSide->SetSigmaAlpha(0.);
  Op_CsI_Surface_DSide->SetPolish(1.);

  G4LogicalBorderSurface * CsI_AirD_Surf = 
    new G4LogicalBorderSurface("CsI_Surface1_DSide", _Physical_Cube,
			       _Physical_PDSideWorld, Op_CsI_Surface_DSide);

  G4LogicalBorderSurface * AirD_CsI_Surf = 
    new G4LogicalBorderSurface("CsI_Surface2_DSide", _Physical_PDSideWorld,
			     _Physical_Cube, Op_CsI_Surface_DSide);

  const G4int num_DSide = 2;
  G4double Ephoton_DSide[num_DSide] = {1.65*eV, 4.28*eV};

  G4double Reflectivity_DSide[num_DSide] = {0., 0.};
  G4double Efficiency_DSide[num_DSide]   = {0., 0.};

  G4double specularlobe_DSide[num_DSide] = {1., 1.};
  G4double specularspike_DSide[num_DSide] = {0., 0.};
  G4double backscatter_DSide[num_DSide] = {0., 0.};

  G4double backpaint_DSide[num_DSide] = {1., 1.}; // refractive index

  G4MaterialPropertiesTable * csiSPT_DSide = new G4MaterialPropertiesTable();

  csiSPT_DSide->AddProperty("REFLECTIVITY", 
		      Ephoton_DSide, Reflectivity_DSide, num_DSide);
  csiSPT_DSide->AddProperty("EFFICIENCY",
		      Ephoton_DSide, Efficiency_DSide,   num_DSide);
  csiSPT_DSide->AddProperty("SPECULARLOBECONSTANT",
		      Ephoton_DSide, specularlobe_DSide, num_DSide);
  csiSPT_DSide->AddProperty("SPECULARSPIKECONSTANT",
		      Ephoton_DSide, specularspike_DSide, num_DSide);
  csiSPT_DSide->AddProperty("BACKSCATTERCONSTANT", 
		      Ephoton_DSide, backscatter_DSide, num_DSide);
  csiSPT_DSide->AddProperty("RINDEX", Ephoton_DSide, backpaint_DSide, num_DSide);

  Op_CsI_Surface_DSide->SetMaterialPropertiesTable(csiSPT_DSide);

  // CsI - Resin Surface
  G4OpticalSurface * Op_Resin_Surface = new G4OpticalSurface("Resin_Surface");
  Op_Resin_Surface->SetType(dielectric_dielectric);
  Op_Resin_Surface->SetFinish(polished);
  Op_Resin_Surface->SetModel(unified);
  Op_Resin_Surface->SetSigmaAlpha(0.);
  Op_Resin_Surface->SetPolish(1.);

  G4LogicalBorderSurface * CsI_Res_Surf = 
    new G4LogicalBorderSurface("CsI_to_Resin", _Physical_Cube, 
  			       _Physical_Resin, Op_Resin_Surface);

  G4LogicalBorderSurface * Res_CsI_Surf = 
     new G4LogicalBorderSurface("Resin_to_CsI", _Physical_Resin, 
   			       _Physical_Cube, Op_Resin_Surface);

  G4double Reflectivity_Resin[num] = {1., 1.};
  G4double Efficiency_Resin[num]   = {0., 0.};

  G4double specularlobe_Resin[num] = {1., 1.};
  G4double specularspike_Resin[num] = {0., 0.};
  G4double backscatter_Resin[num] = {0., 0.};

  G4MaterialPropertiesTable *resinSPT = new G4MaterialPropertiesTable();

  resinSPT->AddProperty("REFLECTIVITY", Ephoton, Reflectivity_Resin, num);
  resinSPT->AddProperty("EFFICIENCY",   Ephoton, Efficiency_Resin,   num);
  resinSPT->AddProperty("SPECULARLOBECONSTANT",
  			Ephoton, specularlobe_Resin, num);
  resinSPT->AddProperty("SPECULARSPIKECONSTANT",
  			Ephoton, specularspike_Resin, num);
  resinSPT->AddProperty("BACKSCATTERCONSTANT", 
  			Ephoton, backscatter_Resin, num);

  Op_Resin_Surface->SetMaterialPropertiesTable(resinSPT);


  // Resin - Photodiode Surface
  G4OpticalSurface * Op_Diode_Surface = new G4OpticalSurface("Diode_Surface");
  Op_Diode_Surface->SetType(dielectric_metal);
  Op_Diode_Surface->SetFinish(polished);
  Op_Diode_Surface->SetModel(glisur);
  Op_Diode_Surface->SetPolish(1.);

  G4LogicalBorderSurface * Res_Dio_Surf = 
    new G4LogicalBorderSurface("Resin_to_Diode", _Physical_Resin, 
			       _Physical_Diode, Op_Diode_Surface);
 
  G4LogicalBorderSurface * Dio_Res_Surf =
    new G4LogicalBorderSurface("Diode_to_Resin", _Physical_Diode,
                               _Physical_Resin, Op_Diode_Surface);
                              

  const G4int num_Diode = 11;
  G4double Ephoton_Diode[num_Diode] =
    { 1.65*eV, 1.77*eV, 1.91*eV, 2.07*eV, 2.25*eV, 2.48*eV, 2.76*eV, 
      3.10*eV, 3.87*eV, 4.13*eV, 4.28*eV};
  //{  750*nm,  700*nm,  650*nm,  600*nm,  550*nm,  500*nm,  450*nm,  
  //   400*nm,  320*nm }
  G4double Reflectivity_Diode[num_Diode] = 
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  G4double Efficiency_Diode[num_Diode] = 
    {   0.794,   0.788,   0.763,   0.723,   0.721,   0.670,   0.634, 
	0.558,   0.388, 0.3, 0.3};

  G4MaterialPropertiesTable * diodeSPT = new G4MaterialPropertiesTable();

  diodeSPT->AddProperty("REFLECTIVITY", Ephoton_Diode, 
			Reflectivity_Diode, num_Diode);
  diodeSPT->AddProperty("EFFICIENCY",   Ephoton_Diode, 
			Efficiency_Diode,   num_Diode);

  Op_Diode_Surface->SetMaterialPropertiesTable(diodeSPT);


//always return the physical World
  return _Physical_World;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void cuboDetectorConstruction::SetReflectivity(G4double ref, 
					       G4OpticalSurface *surface, G4double& cubeReflectivity ){
  G4double ephot[2] = {1.65*eV, 4.28*eV};
  G4double reflect[2];
  reflect[0]=ref;
  reflect[1]=ref;
  G4MaterialPropertiesTable *mpt = surface->GetMaterialPropertiesTable();
  mpt->RemoveProperty("REFLECTIVITY");
  mpt->AddProperty("REFLECTIVITY", ephot, reflect, 2);
  surface->SetMaterialPropertiesTable(mpt);

  cubeReflectivity = ref;
}

