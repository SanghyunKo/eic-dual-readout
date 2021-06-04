#ifndef DRparamBase_h
#define DRparamBase_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include <vector>
#include <cmath>

class DRparamBase {
public:
  DRparamBase();
  virtual ~DRparamBase();

  void SetIsRHS(bool isRHS) { fIsRHS = isRHS; }
  void SetInnerX(double innerX) { fInnerX = innerX; }
  void SetTowerH(double towerH) { fTowerH = towerH; }
  void SetNumZRot(int num) { fNumZRot = num; fPhiZRot = 2*M_PI/(double)num; }
  void SetDeltaTheta(double theta) { fDeltaTheta = theta; }
  void SetThetaOfCenter(double theta) { fThetaOfCenter = theta; }
  void SetSipmHeight(double SipmHeight) { fSipmHeight = SipmHeight; }

  bool GetIsRHS() { return fIsRHS; }
  double GetCurrentInnerR() { return fCurrentInnerR; }
  double GetSipmHeight() { return fSipmHeight; }

  G4RotationMatrix GetRotationMat(int numPhi);
  G4ThreeVector GetTowerPos(int numPhi);
  G4ThreeVector GetSipmLayerPos(int numPhi);
  void GetPt(G4ThreeVector* pt);
  void GetPtSipm(G4ThreeVector* pt);

  G4ThreeVector GetV1() { return fV1; }
  G4ThreeVector GetV2() { return fV2; }
  G4ThreeVector GetV3() { return fV3; }
  G4ThreeVector GetV4() { return fV4; }

  int signedTowerNo(int unsignedTowerNo) { return fIsRHS ? unsignedTowerNo : -unsignedTowerNo-1; }
  int unsignedTowerNo(int signedTowerNo) { return signedTowerNo >= 0 ? signedTowerNo : -signedTowerNo-1; }

  virtual void SetDeltaThetaByTowerNo(int , int ) {}
  virtual void SetThetaOfCenterByTowerNo(int , int ) {}
  void SetIsRHSByTowerNo(int signedTowerNo) { fIsRHS = ( signedTowerNo >=0 ? true : false ); }

  int GetTotTowerNum() { return fTotNum; }
  void SetTotTowerNum(int totNum) { fTotNum = totNum; }

  int GetCurrentTowerNum() { return fCurrentTowerNum; }
  void SetCurrentTowerNum(int numEta) { fCurrentTowerNum = numEta; }

  virtual void init() {};
  void filled() { fFilled = true; }
  void finalized() { fFinalized = true; }
  bool IsFinalized() { return fFinalized; }

protected:
  bool fIsRHS;
  double fPhiZRot;
  double fInnerX;
  double fTowerH;
  int fNumZRot;
  double fDeltaTheta;
  double fThetaOfCenter;
  double fCurrentInnerR;
  G4ThreeVector fCurrentCenter;
  G4ThreeVector fV1;
  G4ThreeVector fV2;
  G4ThreeVector fV3;
  G4ThreeVector fV4;
  G4ThreeVector fV2sipm;
  G4ThreeVector fV4sipm;
  double fSipmHeight;

  double fCurrentInnerHalf;
  double fCurrentOuterHalf;
  double fCurrentOuterHalfSipm;

  int fTotNum;
  int fCurrentTowerNum;
  std::vector<double> fDeltaThetaVec;
  std::vector<double> fThetaOfCenterVec;
  bool fFilled;
  bool fFinalized;
};

#endif
