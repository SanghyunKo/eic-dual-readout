#include "DRparamBase.hh"

#include <stdexcept>

DRparamBase::DRparamBase() {
  fIsRHS = 0;
  fPhiZRot = 0.;
  fInnerX = 0.;
  fTowerH = 0.;
  fNumZRot = 0;
  fDeltaTheta = 0.;
  fThetaOfCenter = 0.;
  fCurrentInnerR = 0.;
  fPhiZRot = 0;
  fCurrentCenter = G4ThreeVector();
  fV1 = G4ThreeVector();
  fV2 = G4ThreeVector();
  fV3 = G4ThreeVector();
  fV4 = G4ThreeVector();
  fSipmHeight = 0.;
  fCurrentInnerHalf = 0.;
  fCurrentOuterHalf = 0.;
  fCurrentTowerNum = 0;
  fFilled = false;
  fFinalized = false;
}

DRparamBase::~DRparamBase() {}

G4RotationMatrix DRparamBase::GetRotationMat(int numPhi) {
  double numPhi_ = (double)numPhi;
  double xRot = fIsRHS ? -fThetaOfCenter : fThetaOfCenter;
  double zRot = fIsRHS ? -M_PI/2. : M_PI/2.;
  G4RotationMatrix rot;
  rot.rotateZ(zRot);
  rot.rotateY(M_PI/2.+xRot);
  rot.rotateZ(numPhi_*fPhiZRot);

  return rot;
}

G4ThreeVector DRparamBase::GetTowerPos(int numPhi) {
  double numPhi_ = (double)numPhi;
  double x = std::cos(numPhi_*fPhiZRot)*fCurrentCenter.x();
  double y = std::sin(numPhi_*fPhiZRot)*fCurrentCenter.x();
  double z = fIsRHS ? fCurrentCenter.z() : -fCurrentCenter.z();
  G4ThreeVector pos = G4ThreeVector(x,y,z);

  return pos;
}

G4ThreeVector DRparamBase::GetSipmLayerPos(int numPhi) {
  double numPhi_ = (double)numPhi;
  double x = std::cos(numPhi_*fPhiZRot)*fCurrentCenter.x()*(fCurrentCenter.mag()+fTowerH/2.+fSipmHeight/2.)/fCurrentCenter.mag();
  double y = std::sin(numPhi_*fPhiZRot)*fCurrentCenter.x()*(fCurrentCenter.mag()+fTowerH/2.+fSipmHeight/2.)/fCurrentCenter.mag();
  double z_abs = fCurrentCenter.z()*(fCurrentCenter.mag()+fTowerH/2.+fSipmHeight/2.)/fCurrentCenter.mag();
  double z = fIsRHS ? z_abs : -z_abs;
  G4ThreeVector pos = G4ThreeVector(x,y,z);

  return pos;
}

void DRparamBase::GetPt(G4ThreeVector* pt) {
  pt[0]=G4ThreeVector(-(fV3.x()*std::tan(fPhiZRot/2.)),-fCurrentInnerHalf,-fTowerH/2.);
  pt[1]=G4ThreeVector((fV3.x()*std::tan(fPhiZRot/2.)),-fCurrentInnerHalf,-fTowerH/2.);
  pt[2]=G4ThreeVector(-(fV1.x()*std::tan(fPhiZRot/2.)),fCurrentInnerHalf,-fTowerH/2.);
  pt[3]=G4ThreeVector((fV1.x()*std::tan(fPhiZRot/2.)),fCurrentInnerHalf,-fTowerH/2.);
  pt[4]=G4ThreeVector(-(fV4.x()*std::tan(fPhiZRot/2.)),-fCurrentOuterHalf,fTowerH/2.);
  pt[5]=G4ThreeVector((fV4.x()*std::tan(fPhiZRot/2.)),-fCurrentOuterHalf,fTowerH/2.);
  pt[6]=G4ThreeVector(-(fV2.x()*std::tan(fPhiZRot/2.)),fCurrentOuterHalf,fTowerH/2.);
  pt[7]=G4ThreeVector((fV2.x()*std::tan(fPhiZRot/2.)),fCurrentOuterHalf,fTowerH/2.);
}

void DRparamBase::GetPtSipm(G4ThreeVector* pt) {
  pt[0] = G4ThreeVector(-(fV4.x()*tan(fPhiZRot/2.)),-fCurrentOuterHalf,-fSipmHeight/2.);
  pt[1] = G4ThreeVector((fV4.x()*tan(fPhiZRot/2.)),-fCurrentOuterHalf,-fSipmHeight/2.);
  pt[2] = G4ThreeVector(-(fV2.x()*tan(fPhiZRot/2.)),fCurrentOuterHalf,-fSipmHeight/2.);
  pt[3] = G4ThreeVector((fV2.x()*tan(fPhiZRot/2.)),fCurrentOuterHalf,-fSipmHeight/2.);
  pt[4] = G4ThreeVector(-(fV4.x()*tan(fPhiZRot/2.)),-fCurrentOuterHalf,fSipmHeight/2.);
  pt[5] = G4ThreeVector((fV4.x()*tan(fPhiZRot/2.)),-fCurrentOuterHalf,fSipmHeight/2.);
  pt[6] = G4ThreeVector(-(fV2.x()*tan(fPhiZRot/2.)),fCurrentOuterHalf,fSipmHeight/2.);
  pt[7] = G4ThreeVector((fV2.x()*tan(fPhiZRot/2.)),fCurrentOuterHalf,fSipmHeight/2.);
}
