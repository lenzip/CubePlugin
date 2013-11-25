#ifndef GiulioOptical_h
#define GiulioOptical_h 1


#include "G4VModularPhysicsList.hh"
class G4OpticalPhysics;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class GiulioOptical : public G4VModularPhysicsList 
{
public:

  GiulioOptical();
  ~GiulioOptical(){}
  void SetCuts();

private:
  G4OpticalPhysics* opticalPhysics;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /* FTFP_BERT_Optical_h */
