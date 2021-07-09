#ifndef DRsimDetectorConstruction_h
#define DRsimDetectorConstruction_h 1

#include "DRsimMagneticField.hh"
#include "DRsimMaterials.hh"
#include "DRsimSiPMHit.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4Trap.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VisAttributes.hh"
#include "G4GenericMessenger.hh"
#include "G4FieldManager.hh"
#include "G4ThreeVector.hh"
#include "G4Region.hh"

#include "DRparamBarrel.hh"
#include "DRparamEndcap.hh"

using namespace std;

class DRsimMagneticField;
class TString;

class DRsimDetectorConstruction : public G4VUserDetectorConstruction {
public:
  DRsimDetectorConstruction();
  virtual ~DRsimDetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();

  G4bool HasBarrel() const { return fDoBarrel; }
  G4bool HasEndcap() const { return fDoEndcap; }

  G4int GetNumBarrel() const { return mNumBarrel; }
  G4int GetNumEndcap() const { return mNumEndcap; }
  G4int GetNumZRot() const { return mNumZRot; }

private:
  void DefineCommands();
  void DefineMaterials();
  G4Material* FindMaterial(G4String matName) { return fMaterials->GetMaterial(matName); }
  G4OpticalSurface* FindSurface(G4String surfName) { return fMaterials->GetOpticalSurface(surfName); }

  void initGeoParam(TString funcFormula, G4double thres, G4double inner, G4double innerHalf, std::vector<G4double>& paramVec, G4int& param);

  void implementTowers(DRparamBase* paramBase, std::vector<G4LogicalVolume*>& PMTcathLogical, std::vector<DRsimInterface::DRsimTowerProperty>& towerProps);
  void implementFibers(DRparamBase* paramBase, G4double dTheta, G4LogicalVolume* towerLogical);

  G4bool checkOverlaps;
  G4GenericMessenger* fMessenger;
  DRsimMaterials* fMaterials;

  static G4ThreadLocal DRsimMagneticField* fMagneticField;
  static G4ThreadLocal G4FieldManager* fFieldMgr;

  G4VisAttributes* fVisAttrOrange;
  G4VisAttributes* fVisAttrBlue;
  G4VisAttributes* fVisAttrGray;
  G4VisAttributes* fVisAttrGreen;
  G4VisAttributes* fVisAttrWhite;

  G4Region* fScintRegion;
  G4Region* fCerenRegion;

  G4double fPMTT;
  G4double fFilterT;
  G4double fSiPMT;

  G4bool fDoBarrel;
  G4bool fDoEndcap;
  G4bool fDoFiber;
  G4String mTowerMaterial;

  G4int mNumBarrel;
  G4int mNumEndcap;
  G4int mNumZRot;
  G4double mTowerH;
  G4double mInnerX;
  G4double mInnerZ;

  std::unique_ptr<DRparamBarrel> pDimB;
  std::unique_ptr<DRparamEndcap> pDimE;
  G4bool fIsEndcap = false;

  char name[20];
  G4Trap* mTowerSolid;
  G4Trap* pmtcath;
  G4Tubs* fiber;
  G4Tubs* fiber_S;
  G4Tubs* fiber_C;
  G4Tubs* fiberS;
  G4Tubs* fiberC;
  G4VSolid* intersect;
  G4VSolid* intersect_;

  std::vector<G4LogicalVolume*> mPMTcathLogicalBR;
  std::vector<G4LogicalVolume*> mPMTcathLogicalBL;
  std::vector<G4LogicalVolume*> mPMTcathLogicalER;
  std::vector<G4LogicalVolume*> mPMTcathLogicalEL;

  DRsimInterface::hitXY fTowerXY;
  std::vector<DRsimInterface::DRsimTowerProperty> fTowerBL;
  std::vector<DRsimInterface::DRsimTowerProperty> fTowerBR;
  std::vector<DRsimInterface::DRsimTowerProperty> fTowerEL;
  std::vector<DRsimInterface::DRsimTowerProperty> fTowerER;

  G4double clad_C_rMin;
  G4double clad_C_rMax;
  G4double clad_C_Sphi;
  G4double clad_C_Dphi;

  G4double core_C_rMin;
  G4double core_C_rMax;
  G4double core_C_Sphi;
  G4double core_C_Dphi;

  G4double clad_S_rMin;
  G4double clad_S_rMax;
  G4double clad_S_Sphi;
  G4double clad_S_Dphi;

  G4double core_S_rMin;
  G4double core_S_rMax;
  G4double core_S_Sphi;
  G4double core_S_Dphi;

  std::vector<G4float> fFiberX;
  std::vector<G4float> fFiberY;

  std::vector<G4double> fDThetaBarrel = { // kept for backward compatibility
    0.02222,0.02220,0.02217,0.02214,0.02209,0.02203,0.02196,0.02188,0.02179,0.02169,
    0.02158,0.02146,0.02133,0.02119,0.02105,0.02089,0.02073,0.02056,0.02039,0.02020,
    0.02002,0.01982,0.01962,0.01941,0.01920,0.01898,0.01876,0.01854,0.01831,0.01808,
    0.01785,0.01761,0.01738,0.01714,0.01689,0.01665,0.01641,0.01616,0.01592,0.01567,
    0.01543,0.01518,0.01494,0.01470,0.01445,0.01421,0.01397,0.01373,0.01350,0.01326,
    0.01303,0.01280
  }; //apply the significance digit
  std::vector<G4double> fDThetaEndcap;
  G4double mTransition;
  G4double mFullTheta;

  G4LogicalVolume* worldLogical;

  G4String setTowerName(bool rbool, G4String BorE, int i) {
    if (rbool) return "T" + BorE + "R" + std::to_string(i);
    else return "T" + BorE + "L" + std::to_string(i);
  }
};

#endif
