#ifndef cuboPhysicsList_h
#define cuboPhysicsList_h 1

#include "globals.hh"
#include "G4VUserPhysicsList.hh"
#include "G4Decay.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hIonisation.hh"

class G4Cerenkov;
class G4Scintillation;
class G4OpAbsorption;
class G4OpRayleigh;
class G4OpMieHG;
class G4OpBoundaryProcess;

class cuboPhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class cuboPhysicsList : public G4VUserPhysicsList
{
public:

  cuboPhysicsList();
  ~cuboPhysicsList();

public:

  void ConstructParticle();
  void ConstructProcess();

  void SetCuts();

  //these methods Construct particles
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();
  void ConstructIons();

  //these methods Construct physics processes and register them
  void ConstructGeneral();
  void ConstructEM();
  void ConstructOp();
    
  //for the Messenger 
  void SetVerbose(G4int);
  void SetNbOfPhotonsCerenkov(G4int);
    
private:

  G4Decay* theDecayProcess;

  G4GammaConversion* theGammaConversionProcess;
  G4ComptonScattering* theComptonScatteringProcess;
  G4PhotoElectricEffect* thePhotoElectricEffectProcess;

  G4eMultipleScattering* TheeMultipleScatteringProcess;
  G4eIonisation* TheeIonisationProcess;
  G4eBremsstrahlung* TheeBremsstrahlungProcess;

  G4eMultipleScattering* TheeMultipleScatteringProcess_p;
  G4eIonisation* TheeIonisationProcess_p;
  G4eBremsstrahlung* TheeBremsstrahlungProcess_p;
  G4eplusAnnihilation* TheeplusAnnihilationProcess;

  G4MuMultipleScattering* TheMuMultipleScatteringProcess;
  G4MuIonisation* TheMuIonisationProcess;
  G4MuBremsstrahlung* TheMuBremsstrahlungProcess;
  G4MuPairProduction* TheMuPairProductionProcess;

  G4hMultipleScattering* ThehMultipleScatteringProcess;
  G4hIonisation* ThehIonisationProcess;

  G4Cerenkov*          theCerenkovProcess;
  G4Scintillation*     theScintillationProcess;
  G4OpAbsorption*      theAbsorptionProcess;
  G4OpRayleigh*        theRayleighScatteringProcess;
  G4OpBoundaryProcess* theBoundaryProcess;
  
  cuboPhysicsListMessenger* pMessenger;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* cuboPhysicsList_h */
