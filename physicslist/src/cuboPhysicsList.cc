#include "globals.hh"
#include "../include/cuboPhysicsList.hh"
#include "../include/cuboPhysicsListMessenger.hh"
#include "montecarlo/pluginmanagers/GGSMCPluginMacros.h"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4IonConstructor.hh"

#include "G4ProcessManager.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListPlugin(cuboPhysicsList);

cuboPhysicsList::cuboPhysicsList() :  G4VUserPhysicsList()
{
  theCerenkovProcess           = NULL;
  theScintillationProcess      = NULL;
  theAbsorptionProcess         = NULL;
  theRayleighScatteringProcess = NULL;
  theBoundaryProcess           = NULL;
  
  pMessenger = new cuboPhysicsListMessenger(this);  
  SetVerboseLevel(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

cuboPhysicsList::~cuboPhysicsList()
{
  delete pMessenger;
  delete theDecayProcess;
  delete theGammaConversionProcess;
  delete theComptonScatteringProcess;
  delete thePhotoElectricEffectProcess;
  delete TheeMultipleScatteringProcess;
  delete TheeIonisationProcess;
  delete TheeBremsstrahlungProcess;
  delete TheeMultipleScatteringProcess_p;
  delete TheeIonisationProcess_p;
  delete TheeBremsstrahlungProcess_p;
  delete TheeplusAnnihilationProcess;
  delete TheMuMultipleScatteringProcess;
  delete TheMuIonisationProcess;
  delete TheMuBremsstrahlungProcess;
  delete TheMuPairProductionProcess;
  delete ThehMultipleScatteringProcess;
  delete ThehIonisationProcess;
  delete theCerenkovProcess;
  delete theScintillationProcess;
  delete theAbsorptionProcess;
  delete theRayleighScatteringProcess;
  delete theBoundaryProcess;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void cuboPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  ConstructIons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void cuboPhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void cuboPhysicsList::ConstructLeptons()
{
  // leptons
  //  e+/-
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  // mu+/-
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  // nu_e
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  // nu_mu
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void cuboPhysicsList::ConstructMesons()
{
  //  mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();

  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void cuboPhysicsList::ConstructBaryons()
{
  //  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();

  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void cuboPhysicsList::ConstructIons()
{
  //  ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void cuboPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructGeneral();
  ConstructEM();
  ConstructOp();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void cuboPhysicsList::ConstructGeneral()
{
  // Add Decay Process
  theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) {
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void cuboPhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
    // gamma
      // Construct processes for gamma
      theGammaConversionProcess = new G4GammaConversion();
      theComptonScatteringProcess = new G4ComptonScattering();
      thePhotoElectricEffectProcess = new G4PhotoElectricEffect();
      pmanager->AddDiscreteProcess(theGammaConversionProcess);
      pmanager->AddDiscreteProcess(theComptonScatteringProcess);
      pmanager->AddDiscreteProcess(thePhotoElectricEffectProcess);

    } else if (particleName == "e-") {
    //electron
      // Construct processes for electron
      TheeMultipleScatteringProcess = new G4eMultipleScattering();
      TheeIonisationProcess = new G4eIonisation();
      TheeBremsstrahlungProcess = new G4eBremsstrahlung();
      pmanager->AddProcess(TheeMultipleScatteringProcess,-1, 1, 1);
      pmanager->AddProcess(TheeIonisationProcess,       -1, 2, 2);
      pmanager->AddProcess(TheeBremsstrahlungProcess,   -1, 3, 3);

    } else if (particleName == "e+") {
    //positron
      // Construct processes for positron
      TheeMultipleScatteringProcess_p = new G4eMultipleScattering();
      TheeIonisationProcess_p = new G4eIonisation();
      TheeBremsstrahlungProcess_p = new G4eBremsstrahlung();
      TheeplusAnnihilationProcess = new G4eplusAnnihilation();
      pmanager->AddProcess(TheeMultipleScatteringProcess_p,-1, 1, 1);
      pmanager->AddProcess(TheeIonisationProcess_p,       -1, 2, 2);
      pmanager->AddProcess(TheeBremsstrahlungProcess_p,   -1, 3, 3);
      pmanager->AddProcess(TheeplusAnnihilationProcess,  0,-1, 4);

    } else if( particleName == "mu+" ||
               particleName == "mu-"    ) {
    //muon
     // Construct processes for muon
      TheMuMultipleScatteringProcess = new G4MuMultipleScattering();
      TheMuIonisationProcess = new G4MuIonisation();
      TheMuBremsstrahlungProcess = new G4MuBremsstrahlung();
      TheMuPairProductionProcess = new G4MuPairProduction();
      pmanager->AddProcess(TheMuMultipleScatteringProcess,-1, 1, 1);
      pmanager->AddProcess(TheMuIonisationProcess,      -1, 2, 2);
      pmanager->AddProcess(TheMuBremsstrahlungProcess,  -1, 3, 3);
      pmanager->AddProcess(TheMuPairProductionProcess,  -1, 4, 4);

    } else {
      if ((particle->GetPDGCharge() != 0.0) &&
          (particle->GetParticleName() != "chargedgeantino")) {
	// all others charged particles except geantino
	ThehMultipleScatteringProcess = new G4hMultipleScattering();
	ThehIonisationProcess = new G4hIonisation();
	pmanager->AddProcess(ThehMultipleScatteringProcess,-1,1,1);
	pmanager->AddProcess(ThehIonisationProcess,       -1,2,2);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void cuboPhysicsList::ConstructOp()
{
  theCerenkovProcess           = new G4Cerenkov("Cerenkov");
  theScintillationProcess      = new G4Scintillation("Scintillation");
  theAbsorptionProcess         = new G4OpAbsorption();
  theRayleighScatteringProcess = new G4OpRayleigh();
  theBoundaryProcess           = new G4OpBoundaryProcess();

//  theCerenkovProcess->DumpPhysicsTable();
//  theScintillationProcess->DumpPhysicsTable();
//  theRayleighScatteringProcess->DumpPhysicsTable();

  SetVerbose(0);
  
  theCerenkovProcess->SetMaxNumPhotonsPerStep(20);
  theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
  theCerenkovProcess->SetTrackSecondariesFirst(true);
  
  theScintillationProcess->SetScintillationYieldFactor(1.);
  theScintillationProcess->SetTrackSecondariesFirst(true);

  // Use Birks Correction in the Scintillation process
  G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  theScintillationProcess->AddSaturation(emSaturation);

  G4OpticalSurfaceModel themodel = unified;
  //theBoundaryProcess->SetModel(themodel);

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    if (theCerenkovProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(theCerenkovProcess);
      pmanager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
    }
    if (theScintillationProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(theScintillationProcess);
      pmanager->SetProcessOrderingToLast(theScintillationProcess, idxAtRest);
      pmanager->SetProcessOrderingToLast(theScintillationProcess, idxPostStep);
    }
    if (particleName == "opticalphoton") {
      G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
      pmanager->AddDiscreteProcess(theAbsorptionProcess);
      pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
      pmanager->AddDiscreteProcess(theBoundaryProcess);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void cuboPhysicsList::SetVerbose(G4int verbose)
{
  theCerenkovProcess->SetVerboseLevel(verbose);
  theScintillationProcess->SetVerboseLevel(verbose);
  theAbsorptionProcess->SetVerboseLevel(verbose);
  theRayleighScatteringProcess->SetVerboseLevel(verbose);
  theBoundaryProcess->SetVerboseLevel(verbose);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void cuboPhysicsList::SetNbOfPhotonsCerenkov(G4int MaxNumber)
{  
  theCerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumber);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void cuboPhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types
  // 
  SetCutsWithDefault();
  
  if (verboseLevel>0) DumpCutValuesTable();   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
