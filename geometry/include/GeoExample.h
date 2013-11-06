/*
 * GeoExample.h
 *
 *  Created on: 29 Aug 2013
 *      Author: Nicola Mori
 */

#ifndef GEOEXAMPLE_H_
#define GEOEXAMPLE_H_

#include "geometry/GGSMaterials.h"
#include "geometry/GGSVGeometryConstruction.h"
#include "geometry/pluginmanagers/GGSGeoPluginMacros.h"

class G4GenericMessenger;

class GeoExample: public GGSVGeometryConstruction{
public:

  GeoExample();
  ~GeoExample();

  G4VPhysicalVolume* Construct();
  G4VPhysicalVolume* GetVolume();

  void HelpUsage() const;
  bool IsInsideAcceptance(const G4ThreeVector &generationPosition, const G4ThreeVector &direction) const;

private:

  G4VPhysicalVolume *_world;
  float _size1, _size2, _distance;
  G4GenericMessenger *_messenger;

};


#endif /* GEOEXAMPLE_H_ */
