#include "../include/GiulioOptical.hh"

#include "G4OpticalPhysics.hh"
#include "G4DataQuestionaire.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "HadronPhysicsFTFP_BERT.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4NeutronTrackingCut.hh"

#include "montecarlo/pluginmanagers/GGSMCPluginMacros.h"

PhysicsListPlugin(GiulioOptical);

GiulioOptical::GiulioOptical():G4VModularPhysicsList(){
  // default cut value  (1.0mm)
  // defaultCutValue = 1.0*CLHEP::mm;
  G4DataQuestionaire it(photon);
  G4cout << "<<< Geant4 Physics List simulation engine: Giulio 12.0"<<G4endl;
  G4cout <<G4endl;
  this->defaultCutValue = 0.7*CLHEP::mm;
  this->SetVerboseLevel(1);  
   // EM Physics
  this->RegisterPhysics( new G4EmStandardPhysics());
  this->ConstructParticle();

  G4cout << "Particle Iterator " << theParticleIterator << G4endl;
  G4cout << "Particle table " << theParticleTable << G4endl;




  // Synchroton Radiation & GN Physics
  this->RegisterPhysics( new G4EmExtraPhysics() );

  // Decays
  this->RegisterPhysics( new G4DecayPhysics() );

   // Hadron Elastic scattering
  this->RegisterPhysics( new G4HadronElasticPhysics() );

   // Hadron Physics
  this->RegisterPhysics(  new HadronPhysicsFTFP_BERT());

  // Stopping Physics
  this->RegisterPhysics( new G4StoppingPhysics() );

  // Ion Physics
  this->RegisterPhysics( new G4IonPhysics());

  // Neutron tracking cut
  this->RegisterPhysics( new G4NeutronTrackingCut());
  
  opticalPhysics = new G4OpticalPhysics();

  opticalPhysics->SetWLSTimeProfile("delta");

  opticalPhysics->SetScintillationYieldFactor(1.);
  opticalPhysics->SetScintillationExcitationRatio(1.0);

  opticalPhysics->SetMaxNumPhotonsPerStep(100);
  opticalPhysics->SetMaxBetaChangePerStep(10.0);

  opticalPhysics->SetTrackSecondariesFirst(kCerenkov,true);
  opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);
  this->RegisterPhysics( opticalPhysics );
  this->SetVerboseLevel(500);

  G4int iphys =0;
  while (this->GetPhysics(iphys)){
    G4cout << "Registered " << this->GetPhysics(iphys)->GetPhysicsName() << G4endl;
    iphys++;
  }

}

void GiulioOptical::SetCuts()
{
  
  this->SetCutsWithDefault();
   G4cout << "Particle Iterator in setcuts " << theParticleIterator << G4endl;
     G4cout << "Particle table in setcuts " << theParticleTable << G4endl;

}

