#ifndef DRsimSteppingAction_h
#define DRsimSteppingAction_h 1

#include "DRsimInterface.h"
#include "DRsimEventAction.hh"
#include "DRsimDetectorConstruction.hh"

#include "G4UserSteppingAction.hh"
#include "G4LogicalVolume.hh"
#include "G4Step.hh"
#include "G4SteppingManager.hh"

using namespace std;

class DRsimSteppingAction : public G4UserSteppingAction {
public:
  DRsimSteppingAction(DRsimEventAction* eventAction);
  virtual ~DRsimSteppingAction();
  virtual void UserSteppingAction(const G4Step*);

private:
  const DRsimDetectorConstruction* DetectorConstruction;
  int fNumBarrel;
  DRsimEventAction* fEventAction;
  DRsimInterface::DRsimEdepData fEdep;
  DRsimInterface::DRsimLeakageData fLeak;

  G4VPhysicalVolume* GetMotherTower(G4TouchableHandle touchable) { return touchable->GetVolume(touchable->GetHistoryDepth()-1); }
  G4int GetTowerIPhi(G4VPhysicalVolume* motherTower_) { return motherTower_->GetCopyNo(); }
  G4int GetTowerITheta(G4String towerName);
};

#endif
