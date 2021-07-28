#include "DRsimInterface.h"
#include "TF1.h"

#include "DRsimDetectorConstruction.hh"
#include "DRsimCellParameterisation.hh"
#include "DRsimFilterParameterisation.hh"
#include "DRsimSiPMSD.hh"
#include "FastOpTransportModel.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "G4IntersectionSolid.hh"
#include "G4SDManager.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GeometryManager.hh"
#include "G4Transform3D.hh"
#include "G4Colour.hh"

#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

G4ThreadLocal DRsimMagneticField* DRsimDetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* DRsimDetectorConstruction::fFieldMgr = 0;

DRsimDetectorConstruction::DRsimDetectorConstruction()
: G4VUserDetectorConstruction(), fMessenger(0), fMaterials(NULL) {
  DefineCommands();
  DefineMaterials();

  clad_C_rMin = 0.49*CLHEP::millimeter;
  clad_C_rMax = 0.50*CLHEP::millimeter;
  clad_C_Sphi = 0.;
  clad_C_Dphi = 2.*M_PI;

  core_C_rMin = 0.*CLHEP::millimeter;
  core_C_rMax = 0.49*CLHEP::millimeter;
  core_C_Sphi = 0.;
  core_C_Dphi = 2.*M_PI;

  clad_S_rMin = 0.485*CLHEP::millimeter;
  clad_S_rMax = 0.50*CLHEP::millimeter;
  clad_S_Sphi = 0.;
  clad_S_Dphi = 2.*M_PI;

  core_S_rMin = 0.*CLHEP::millimeter;
  core_S_rMax = 0.485*CLHEP::millimeter;
  core_S_Sphi = 0.;
  core_S_Dphi = 2.*M_PI;
  fPMTT = 0.3*CLHEP::millimeter;
  fFilterT = 0.01*CLHEP::millimeter;
  fSiPMT = 0.01*CLHEP::millimeter;

  fVisAttrOrange = new G4VisAttributes(G4Colour(1.0,0.5,0.,1.0));
  fVisAttrOrange->SetVisibility(true);
  fVisAttrBlue = new G4VisAttributes(G4Colour(0.,0.,1.0,1.0));
  fVisAttrBlue->SetVisibility(true);
  fVisAttrGray = new G4VisAttributes(G4Colour(0.3,0.3,0.3,0.3));
  fVisAttrGray->SetVisibility(true);
  fVisAttrGreen = new G4VisAttributes(G4Colour(0.3,0.7,0.3));
  fVisAttrGreen->SetVisibility(true);
  fVisAttrWhite = new G4VisAttributes(G4Colour(1.,1.,1.));
  fVisAttrWhite->SetVisibility(true);

  fDoBarrel = false;
  fDoEndcap = true;
  fDoFiber = true;
  mTowerMaterial = "Copper"; //"Tungsten"

  mNumBarrel = 52;
  mNumEndcap = 40;
  mNumZRot = 283;
  mTowerH = 1250.*CLHEP::millimeter; // tower height
  mInnerX = 1.8*CLHEP::meter;
  mInnerZ = 2.556*CLHEP::meter;
  fDThetaEndcap = std::vector<G4double>(mNumEndcap,fDThetaBarrel.at(52-1)); // default (for Copper)
  mTransition = 0.95717;
  mFullTheta = 0.;

  if (mTowerMaterial!="Copper") {
    fDThetaBarrel.clear();
    fDThetaEndcap.clear();

    // modify according to Moliere_radius of the material
    G4double innerHalf = 28.5365/2.*CLHEP::millimeter; // Tungsten, 2*innerHalf = Moliere_radius * std::sqrt(M_PI)
    mNumZRot = 165; // 2*M_PI*mInnerZ*std::tan(0.285) = 2*innerHalf*mNumZRot

    TString funcFormula = "[1]*TMath::Tan(x/2.)/TMath::Cos([0]+x/2.)";
    G4double thres = 0.95;
    G4double inner = mInnerX;
    initGeoParam(funcFormula,thres,inner,innerHalf,fDThetaBarrel,mNumBarrel);

    mTransition = mFullTheta;
    mInnerZ = mInnerX*std::tan(mTransition);

    funcFormula = "[1]*TMath::Tan(x/2.)/TMath::Sin([0]+x/2.)";
    thres = 1.47;
    inner = mInnerZ;
    initGeoParam(funcFormula,thres,inner,innerHalf,fDThetaEndcap,mNumEndcap);
  }
}

DRsimDetectorConstruction::~DRsimDetectorConstruction() {
  delete fMessenger;
  delete fMaterials;

  delete fVisAttrOrange;
  delete fVisAttrBlue;
  delete fVisAttrGray;
  delete fVisAttrGreen;
}

void DRsimDetectorConstruction::DefineMaterials() {
  fMaterials = DRsimMaterials::GetInstance();
}

G4VPhysicalVolume* DRsimDetectorConstruction::Construct() {
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  checkOverlaps = false;

  G4VSolid* worldSolid = new G4Box("worldBox",10.*CLHEP::meter,10.*CLHEP::meter,10.*CLHEP::meter);
  worldLogical = new G4LogicalVolume(worldSolid,FindMaterial("G4_Galactic"),"worldLogical");
  G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,false,0,checkOverlaps);

  fiber = new G4Tubs("fiber",0,clad_C_rMax,mTowerH/2.,0*CLHEP::degree,360.*CLHEP::degree);// S is the same
  fiberC = new G4Tubs("fiberC",0,core_C_rMax,mTowerH/2.,0*CLHEP::degree,360.*CLHEP::degree);
  fiberS = new G4Tubs("fiberS",0,core_S_rMax,mTowerH/2.,0*CLHEP::degree,360.*CLHEP::degree);

  fCerenRegion = new G4Region("cerenRegion");
  fScintRegion = new G4Region("scintRegion");

  // barrel
  if (fDoBarrel) {
    pDimB = std::make_unique<DRparamBarrel>();
    pDimB->SetInnerX(mInnerX);
    pDimB->SetTowerH(mTowerH);
    pDimB->SetNumZRot(mNumZRot);
    pDimB->SetSipmHeight(fPMTT+fFilterT);

    pDimB->SetIsRHS(true);
    mFullTheta = 0.;
    implementTowers(pDimB.get(),mPMTcathLogicalBR,fTowerBR);

    pDimB->SetIsRHS(false);
    mFullTheta = 0.;
    implementTowers(pDimB.get(),mPMTcathLogicalBL,fTowerBL);
  }

  fIsEndcap = true;

  // endcap
  if (fDoEndcap) {
    pDimE = std::make_unique<DRparamEndcap>();
    pDimE->SetInnerX(mInnerZ);
    pDimE->SetTowerH(mTowerH);
    pDimE->SetNumZRot(mNumZRot);
    pDimE->SetSipmHeight(fPMTT+fFilterT);

    mFullTheta = mTransition;
    pDimE->SetIsRHS(true);

    implementTowers(pDimE.get(),mPMTcathLogicalER,fTowerER);

    // endcap L
    mFullTheta = mTransition;
    pDimE->SetIsRHS(false);

    implementTowers(pDimE.get(),mPMTcathLogicalEL,fTowerEL);
  }

  return worldPhysical;
}

void DRsimDetectorConstruction::ConstructSDandField() {
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String SiPMName = "SiPMSD";

  if (!fDoFiber) return;

  // Not a memory leak - SDs are deleted by G4SDManager. Deleting them manually will cause double delete!
  if (fDoBarrel) {
    for (int i = 0; i < mNumBarrel; i++) {
      DRsimSiPMSD* SiPMSDBR = new DRsimSiPMSD("BR"+std::to_string(i),"BRC"+std::to_string(i),fTowerBR.at(i));
      SDman->AddNewDetector(SiPMSDBR);
      mPMTcathLogicalBR.at(i)->SetSensitiveDetector(SiPMSDBR);
    }

    for (int i = 0; i < mNumBarrel; i++) {
      DRsimSiPMSD* SiPMSDBL = new DRsimSiPMSD("BL"+std::to_string(i),"BLC"+std::to_string(i),fTowerBL.at(i));
      SDman->AddNewDetector(SiPMSDBL);
      mPMTcathLogicalBL.at(i)->SetSensitiveDetector(SiPMSDBL);
    }
  }

  if (fDoEndcap) {
    for (int i = 0; i < mNumEndcap; i++) {
      DRsimSiPMSD* SiPMSDER = new DRsimSiPMSD("ER"+std::to_string(i),"ERC"+std::to_string(i),fTowerER.at(i));
      SDman->AddNewDetector(SiPMSDER);
      mPMTcathLogicalER.at(i)->SetSensitiveDetector(SiPMSDER);
    }

    for (int i = 0; i < mNumEndcap; i++) {
      DRsimSiPMSD* SiPMSDEL = new DRsimSiPMSD("EL"+std::to_string(i),"ELC"+std::to_string(i),fTowerEL.at(i));
      SDman->AddNewDetector(SiPMSDEL);
      mPMTcathLogicalEL.at(i)->SetSensitiveDetector(SiPMSDEL);
    }
  }

  FastOpTransportModel* cerenModel = new FastOpTransportModel("fastOpTransportCeren",fCerenRegion);
  FastOpTransportModel* scintModel = new FastOpTransportModel("fastOpTransportScint",fScintRegion);
  cerenModel->SetFiberLength(mTowerH);
  cerenModel->SetCoreMaterial(FindMaterial("PMMA"));
  scintModel->SetFiberLength(mTowerH);
  scintModel->SetCoreMaterial(FindMaterial("Polystyrene"));
}

void DRsimDetectorConstruction::initGeoParam(TString funcFormula, G4double thres, G4double inner, G4double innerHalf, std::vector<G4double>& paramVec, G4int& param) {
  std::unique_ptr<TF1> calcTowerInnerHalf = std::make_unique<TF1>("calcTowerInnerHalf",funcFormula,0.,0.1); // FIXME need to sole without introducing ROOT, causing warnings!

  while (mFullTheta < thres) {
    calcTowerInnerHalf->SetParameter(1,inner);
    calcTowerInnerHalf->SetParameter(0,mFullTheta);
    G4double dTheta = calcTowerInnerHalf->GetX(innerHalf,0.,0.1);
    mFullTheta += dTheta;
    paramVec.push_back(dTheta);
  }

  param = paramVec.size();
}

void DRsimDetectorConstruction::implementTowers(DRparamBase* paramBase, std::vector<G4LogicalVolume*>& PMTcathLogical, std::vector<DRsimInterface::DRsimTowerProperty>& towerProps) {
  G4String moduleType = (fIsEndcap ? "E" : "B");

  for (int i = 0; i < (fIsEndcap ? mNumEndcap : mNumBarrel); i++) {
    double dTheta = (fIsEndcap ? fDThetaEndcap.at(i) : fDThetaBarrel.at(i)); // default geometry
    double towerTheta = mFullTheta + dTheta/2.;
    G4ThreeVector pt[8];
    paramBase->SetDeltaTheta( dTheta );
    paramBase->SetThetaOfCenter(towerTheta);
    paramBase->init();
    paramBase->GetPt(pt);
    G4String towerName = setTowerName(paramBase->GetIsRHS(), moduleType, i);

    mTowerSolid = new G4Trap("Tower"+moduleType,pt);
    G4LogicalVolume* towerLogical = new G4LogicalVolume(mTowerSolid,FindMaterial(mTowerMaterial),towerName);

    paramBase->GetPtSipm(pt);
    G4Trap* pmtg = new G4Trap("PMTG"+moduleType,pt);
    G4LogicalVolume* PMTGLogical = new G4LogicalVolume(pmtg,FindMaterial("G4_AIR"),towerName);

    for(int j=0;j<mNumZRot;j++){
      auto rot = paramBase->GetRotationMat(j);
      new G4PVPlacement(G4Transform3D(rot,paramBase->GetTowerPos(j)),towerLogical,towerName,worldLogical,false,j,checkOverlaps);
      new G4PVPlacement(G4Transform3D(rot,paramBase->GetSipmLayerPos(j)),PMTGLogical,towerName,worldLogical,false,j,checkOverlaps);
    }

    if (fDoFiber)
      implementFibers(paramBase,dTheta,towerLogical);

    int iTheta = paramBase->signedTowerNo( fIsEndcap ? mNumBarrel+i : i );
    float signedTowerTheta = paramBase->GetIsRHS() ? towerTheta : -towerTheta;
    DRsimInterface::DRsimTowerProperty towerProp;
    towerProp.towerXY = fTowerXY;
    towerProp.towerTheta = std::make_pair(iTheta,signedTowerTheta);
    towerProp.innerR = paramBase->GetCurrentInnerR();
    towerProp.towerH = mTowerH;
    towerProp.dTheta = dTheta;
    towerProps.push_back(towerProp);

    if (fDoFiber) {
      G4VSolid* SiPMlayerSolid = new G4Box("SiPMlayerSolid",fTowerXY.first*1.5/2.*CLHEP::millimeter,fTowerXY.second*1.5/2.*CLHEP::millimeter,fPMTT/2.);
      G4LogicalVolume* SiPMlayerLogical = new G4LogicalVolume(SiPMlayerSolid,FindMaterial("G4_POLYVINYL_CHLORIDE"),"SiPMlayerLogical");
      new G4PVPlacement(0,G4ThreeVector(0.,0.,fFilterT/2.),SiPMlayerLogical,"SiPMlayerPhysical",PMTGLogical,false,0,checkOverlaps);

      G4VSolid* filterlayerSolid = new G4Box("filterlayerSolid",fTowerXY.first*1.5/2.*CLHEP::millimeter,fTowerXY.second*1.5/2.*CLHEP::millimeter,fFilterT/2.);
      G4LogicalVolume* filterlayerLogical = new G4LogicalVolume(filterlayerSolid,FindMaterial("G4_POLYVINYL_CHLORIDE"),"filterlayerLogical");
      new G4PVPlacement(0,G4ThreeVector(0.,0.,-fPMTT/2.),filterlayerLogical,"filterlayerPhysical",PMTGLogical,false,0,checkOverlaps);

      G4VSolid* PMTcellSolid = new G4Box("PMTcellSolid",1.2/2.*CLHEP::millimeter,1.2/2.*CLHEP::millimeter,fPMTT/2.);
      G4LogicalVolume* PMTcellLogical = new G4LogicalVolume(PMTcellSolid,FindMaterial("Glass"),"PMTcellLogical");

      DRsimCellParameterisation* PMTcellParam = new DRsimCellParameterisation(fTowerXY.first,fTowerXY.second);
      G4PVParameterised* PMTcellPhysical = new G4PVParameterised("PMTcellPhysical",PMTcellLogical,SiPMlayerLogical,kXAxis,fTowerXY.first*fTowerXY.second,PMTcellParam);

      G4VSolid* PMTcathSolid = new G4Box("PMTcathSolid",1.2/2.*CLHEP::millimeter,1.2/2.*CLHEP::millimeter,fSiPMT/2.*CLHEP::millimeter);
      G4LogicalVolume* pmtCathLogic = new G4LogicalVolume(PMTcathSolid,FindMaterial("Silicon"),"PMTcathLogical");
      pmtCathLogic->SetVisAttributes(fVisAttrGreen);
      PMTcathLogical.push_back(pmtCathLogic);

      new G4PVPlacement(0,G4ThreeVector(0.,0.,(fPMTT-fSiPMT)/2.*CLHEP::millimeter),pmtCathLogic,"PMTcathPhysical",PMTcellLogical,false,0,checkOverlaps);
      new G4LogicalSkinSurface("Photocath_surf",pmtCathLogic,FindSurface("SiPMSurf"));

      G4VSolid* filterSolid = new G4Box("filterSolid",1.2/2.*CLHEP::millimeter,1.2/2.*CLHEP::millimeter,fFilterT/2.);
      G4LogicalVolume* PMTfilterLogical = new G4LogicalVolume(filterSolid,FindMaterial("Glass"),"PMTfilterLogical");

      DRsimFilterParameterisation* filterParam = new DRsimFilterParameterisation(fTowerXY.first,fTowerXY.second,FindMaterial("Glass"),FindMaterial("Gelatin"));
      filterParam->SetFilterVis(fVisAttrOrange);
      filterParam->SetGlassVis(fVisAttrWhite);

      G4PVParameterised* filterPhysical = new G4PVParameterised("filterPhysical",PMTfilterLogical,filterlayerLogical,kXAxis,fTowerXY.first*fTowerXY.second,filterParam);
      new G4LogicalBorderSurface("filterSurf",filterPhysical,PMTcellPhysical,FindSurface("FilterSurf"));
    }

    mFullTheta += dTheta;
  }
}

void DRsimDetectorConstruction::DefineCommands() {}

void DRsimDetectorConstruction::implementFibers(DRparamBase* paramBase, G4double dTheta, G4LogicalVolume* towerLogical) {
  fFiberX.clear();
  fFiberY.clear();

  G4ThreeVector v4 = paramBase->GetV4();

  G4double outerSide_half = (paramBase->GetCurrentInnerR()+mTowerH)*std::tan(dTheta/2.);
  G4double phiUnit = 2*M_PI/static_cast<G4double>(mNumZRot);

  int numx = (int)(((v4.getX()*std::tan(phiUnit/2.)*2)-1.2*CLHEP::millimeter)/(1.5*CLHEP::millimeter)) + 1;
  int numy = (int)((outerSide_half*2-1.2*CLHEP::millimeter)/(1.5*CLHEP::millimeter)) + 1;
  fTowerXY = std::make_pair(numx,numy);

  for (int j = 0; j < numy; j++) {
    for (int k = 0; k < numx; k++) {
      G4float fX = -1.5*CLHEP::millimeter*(numx/2) + k*1.5*CLHEP::millimeter + ( numx%2==0 ? 0.75*CLHEP::millimeter : 0 );
      G4float fY = -1.5*CLHEP::millimeter*(numy/2) + j*1.5*CLHEP::millimeter + ( numy%2==0 ? 0.75*CLHEP::millimeter : 0 );
      fFiberX.push_back(fX);
      fFiberY.push_back(fY);
    }
  }

  for (unsigned int j = 0; j<fFiberX.size();j++) {
    G4int column = j % numx;
    G4int row = j / numx;

    if ( DRsimInterface::IsCerenkov(column,row) ) { //c fibre
      intersect = new G4IntersectionSolid("fiber_",fiber,mTowerSolid,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
      G4LogicalVolume* cladLogical = new G4LogicalVolume(intersect,FindMaterial("FluorinatedPolymer"),name);
      new G4PVPlacement(0,G4ThreeVector(fFiberX.at(j),fFiberY.at(j),0),cladLogical,name,towerLogical,false,j,checkOverlaps);

      intersect_ = new G4IntersectionSolid("fiber_",fiberC,mTowerSolid,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
      G4LogicalVolume* coreLogical = new G4LogicalVolume(intersect_,FindMaterial("PMMA"),name);
      new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),coreLogical,name,cladLogical,false,j,checkOverlaps);

      fCerenRegion->AddRootLogicalVolume(cladLogical);
      fCerenRegion->AddRootLogicalVolume(coreLogical);
      cladLogical->SetRegion(fCerenRegion);
      coreLogical->SetRegion(fCerenRegion);

      cladLogical->SetVisAttributes(fVisAttrGray);
      coreLogical->SetVisAttributes(fVisAttrBlue);
    } else { // s fibre
      intersect = new G4IntersectionSolid("fiber_",fiber,mTowerSolid,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
      G4LogicalVolume* cladLogical = new G4LogicalVolume(intersect,FindMaterial("PMMA"),name);
      new G4PVPlacement(0,G4ThreeVector(fFiberX.at(j),fFiberY.at(j),0),cladLogical,name,towerLogical,false,j,checkOverlaps);

      intersect_ = new G4IntersectionSolid("fiber_",fiberS,mTowerSolid,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
      G4LogicalVolume* coreLogical = new G4LogicalVolume(intersect_,FindMaterial("Polystyrene"),name);
      new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),coreLogical,name,cladLogical,false,j,checkOverlaps);

      fScintRegion->AddRootLogicalVolume(cladLogical);
      fScintRegion->AddRootLogicalVolume(coreLogical);
      cladLogical->SetRegion(fScintRegion);
      coreLogical->SetRegion(fScintRegion);

      cladLogical->SetVisAttributes(fVisAttrGray);
      coreLogical->SetVisAttributes(fVisAttrOrange);
    }
  }
}
