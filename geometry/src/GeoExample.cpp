/*
 * GeoExample.cpp
 *
 *  Created on: 29 Aug 2013
 *      Author: Nicola Mori
 */

#include "../include/GeoExample.h"

#include "geometry/GGSMaterials.h"
#include "geometry/pluginmanagers/GGSGeoPluginMacros.h"

#include "G4GenericMessenger.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

GeometryPlugin(GeoExample);

GeoExample::GeoExample() :
    _world(NULL), _size1(5 * cm), _size2(20 * cm), _distance(10 * cm) {
  G4cout << "[GeoExample::GeoExample] Constructor" << G4endl;
  _messenger = new G4GenericMessenger(this, "/GGS/examples/geometry/GeoExample/");
  _messenger->DeclareProperty("size1", _size1, "Set the size of upstream cube.").SetUnit("cm");
  _messenger->DeclareProperty("size2", _size2, "Set the size of downstream cube.").SetUnit("cm");
  _messenger->DeclareProperty("distance", _distance, "Set the distance between cubes.").SetUnit("cm");
}

GeoExample::~GeoExample() {
  G4cout << "[GeoExample::~GeoExample] Destructor" << G4endl;
}

G4VPhysicalVolume* GeoExample::Construct() {
  G4cout << "[GeoExample::Construct] GeoExample construction" << G4endl;

  // Run geometry configuration script before building
  if (_geoDataCard != "") {
    G4UImanager::GetUIpointer()->ApplyCommand(G4String("/control/execute " + _geoDataCard));
  }
  // Delete the messenger so that the commands for configuring the geometry won't be available anymore
  delete _messenger;
  _messenger = NULL;

  // Create materials
  GGSMaterials *materials = new GGSMaterials;

  // Create the world volume
  float transvSize = fmax(_size1, _size2)/2.;
  G4Box *worldSolid = new G4Box("worldSolid", transvSize*1.1, transvSize*1.1, _size1 + _size2 + _distance);
  G4LogicalVolume *worldLogic = new G4LogicalVolume(worldSolid, G4Material::GetMaterial("AirMaterial"), "worldLogic");
  _world = new G4PVPlacement(NULL, G4ThreeVector(0, 0, 0), worldLogic, "worlPhys", NULL, false, 0);

  // Create small cube
  G4Box *cube1Solid = new G4Box("cube1Solid", _size1 / 2., _size1 / 2., _size1 / 2.);
  G4LogicalVolume *cube1Logic = new G4LogicalVolume(cube1Solid, G4Material::GetMaterial("PVTMaterial"), "cube1Logic");
  G4PVPlacement *cube1Phys = new G4PVPlacement(NULL, G4ThreeVector(0, 0, -_size1 / 2.), cube1Logic, "cube1Phys",
      worldLogic, false, 0);

  // Create big cube
  G4Box *cube2Solid = new G4Box("cube2Solid", _size2 / 2., _size2 / 2., _size2 / 2.);
  G4LogicalVolume *cube2Logic = new G4LogicalVolume(cube2Solid, G4Material::GetMaterial("BGOMaterial"), "cube2Logic");
  G4PVPlacement *cube2Phys = new G4PVPlacement(NULL, G4ThreeVector(0, 0, -_size1 - _distance - _size2/2.), cube2Logic, "cube2Phys",
      worldLogic, false, 0);

  return _world;

}

G4VPhysicalVolume* GeoExample::GetVolume() {
  return _world;
}

void GeoExample::HelpUsage() const {
  G4cout << "[GeoExample::HelpUsage] GGS example geometry" << G4endl;
  G4cout << "    This is a sample geometry consisting of two cubes placed along the Z axis." << G4endl;

}

// Acceptance check: particle must hit on the upstream face of the small cube
bool GeoExample::IsInsideAcceptance(const G4ThreeVector &generationPosition, const G4ThreeVector &direction) const {

  float deltaZ = generationPosition[2]; // Upstream face is at Z=0
  float impactX = generationPosition[0] + direction[0] / direction[2] * deltaZ;
  float impactY = generationPosition[1] + direction[1] / direction[2] * deltaZ;
  if (fabs(impactX) > _size1 / 2. || fabs(impactY) > _size1 / 2.)
    return false;

  return true;
}
