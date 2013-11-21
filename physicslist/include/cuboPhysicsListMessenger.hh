#ifndef cuboPhysicsListMessenger_h
#define cuboPhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class cuboPhysicsList;
class G4UIdirectory;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class cuboPhysicsListMessenger: public G4UImessenger
{
  public:  
    cuboPhysicsListMessenger(cuboPhysicsList* );
   ~cuboPhysicsListMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:  
    cuboPhysicsList*     pPhysicsList;
    
    G4UIdirectory*        cuboDir;
    G4UIdirectory*        physDir;
    G4UIcmdWithAnInteger* verboseCmd;
    G4UIcmdWithAnInteger* cerenkovCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

