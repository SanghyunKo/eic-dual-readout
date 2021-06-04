#include "DRsimInterface.h"
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
#include "G4SystemOfUnits.hh"

using namespace std;

G4ThreadLocal DRsimMagneticField* DRsimDetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* DRsimDetectorConstruction::fFieldMgr = 0;

DRsimDetectorConstruction::DRsimDetectorConstruction()
: G4VUserDetectorConstruction(), fMessenger(0), fMaterials(NULL) {
  DefineCommands();
  DefineMaterials();

  clad_C_rMin = 0.49*mm;
  clad_C_rMax = 0.50*mm;
  clad_C_Sphi = 0.;
  clad_C_Dphi = 2.*M_PI;

  core_C_rMin = 0.*mm;
  core_C_rMax = 0.49*mm;
  core_C_Sphi = 0.;
  core_C_Dphi = 2.*M_PI;

  clad_S_rMin = 0.485*mm;
  clad_S_rMax = 0.50*mm;
  clad_S_Sphi = 0.;
  clad_S_Dphi = 2.*M_PI;

  core_S_rMin = 0.*mm;
  core_S_rMax = 0.485*mm;
  core_S_Sphi = 0.;
  core_S_Dphi = 2.*M_PI;
  fPMTT = 0.3*mm;
  fFilterT = 0.01*mm;
  fSiPMT = 0.01*mm;

  theta_unit=0;
  phi_unit=0;
  fDThetaEndcap = fDThetaBarrel[52-1];

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
  mTowerMaterial = "Copper"; // "Tungsten"

  mNumBarrel = 52;
  mNumEndcap = 40;
  mNumZRot = 283;
  mTowerH = 1250.*mm; // tower height
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

  G4VSolid* worldSolid = new G4Box("worldBox",10.*m,10.*m,10.*m);
  worldLogical = new G4LogicalVolume(worldSolid,FindMaterial("G4_Galactic"),"worldLogical");
  G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,false,0,checkOverlaps);

  innerR = 1800.;
  fulltheta = 0.;
  phi_unit = 2*M_PI/(G4double)mNumZRot;

  fiber = new G4Tubs("fiber",0,clad_C_rMax,mTowerH/2.,0*deg,360.*deg);// S is the same
  fiberC = new G4Tubs("fiberC",0,core_C_rMax,mTowerH/2.,0*deg,360.*deg);
  fiberS = new G4Tubs("fiberS",0,core_S_rMax,mTowerH/2.,0*deg,360.*deg);

  fCerenRegion = new G4Region("cerenRegion");
  fScintRegion = new G4Region("scintRegion");

  // barrel
  if (fDoBarrel) {
    pDimB = std::make_unique<DRparamBarrel>();
    pDimB->SetInnerX(innerR);
    pDimB->SetTowerH(mTowerH);
    pDimB->SetNumZRot(mNumZRot);
    pDimB->SetSipmHeight(fPMTT+fFilterT);

    pDimB->SetIsRHS(true);
    fulltheta = 0.;
    Barrel(towerLogicalBR,PMTGLogicalBR,PMTfilterLogicalBR,PMTcellLogicalBR,PMTcathLogicalBR,fiberLogical_BR,fiberLogical_BR_,fTowerBR);

    pDimB->SetIsRHS(false);
    fulltheta = 0.;
    Barrel(towerLogicalBL,PMTGLogicalBL,PMTfilterLogicalBL,PMTcellLogicalBL,PMTcathLogicalBL,fiberLogical_BL,fiberLogical_BL_,fTowerBL);
  }

  // endcap
  if (fDoEndcap) {
    pDimE = std::make_unique<DRparamEndcap>();
    pDimE->SetInnerX(2.556*m);
    pDimE->SetTowerH(mTowerH);
    pDimE->SetNumZRot(mNumZRot);
    pDimE->SetDeltaTheta(fDThetaEndcap);
    pDimE->SetSipmHeight(fPMTT+fFilterT);

    fulltheta = 0.95717;
    pDimE->SetIsRHS(true);

    Endcap(towerLogicalER,PMTGLogicalER,PMTfilterLogicalER,PMTcellLogicalER,PMTcathLogicalER,fiberLogical_ER,fiberLogical_ER_,fTowerER);

    // endcap L
    fulltheta = 0.95717;
    pDimE->SetIsRHS(false);

    Endcap(towerLogicalEL,PMTGLogicalEL,PMTfilterLogicalEL,PMTcellLogicalEL,PMTcathLogicalEL,fiberLogical_EL,fiberLogical_EL_,fTowerEL);
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
      PMTcathLogicalBR[i]->SetSensitiveDetector(SiPMSDBR);
    }

    for (int i = 0; i < mNumBarrel; i++) {
      DRsimSiPMSD* SiPMSDBL = new DRsimSiPMSD("BL"+std::to_string(i),"BLC"+std::to_string(i),fTowerBL.at(i));
      SDman->AddNewDetector(SiPMSDBL);
      PMTcathLogicalBL[i]->SetSensitiveDetector(SiPMSDBL);
    }
  }

  if (fDoEndcap) {
    for (int i = 0; i < mNumEndcap; i++) {
      DRsimSiPMSD* SiPMSDER = new DRsimSiPMSD("ER"+std::to_string(i),"ERC"+std::to_string(i),fTowerER.at(i));
      SDman->AddNewDetector(SiPMSDER);
      PMTcathLogicalER[i]->SetSensitiveDetector(SiPMSDER);
    }

    for (int i = 0; i < mNumEndcap; i++) {
      DRsimSiPMSD* SiPMSDEL = new DRsimSiPMSD("EL"+std::to_string(i),"ELC"+std::to_string(i),fTowerEL.at(i));
      SDman->AddNewDetector(SiPMSDEL);
      PMTcathLogicalEL[i]->SetSensitiveDetector(SiPMSDEL);
    }
  }

  FastOpTransportModel* cerenModel = new FastOpTransportModel("fastOpTransportCeren",fCerenRegion);
  FastOpTransportModel* scintModel = new FastOpTransportModel("fastOpTransportScint",fScintRegion);
  cerenModel->SetFiberLength(mTowerH);
  cerenModel->SetCoreMaterial(FindMaterial("PMMA"));
  scintModel->SetFiberLength(mTowerH);
  scintModel->SetCoreMaterial(FindMaterial("Polystyrene"));
}

void DRsimDetectorConstruction::Barrel(G4LogicalVolume* towerLogical[], G4LogicalVolume* PMTGLogical[], G4LogicalVolume* PMTfilterLogical[], G4LogicalVolume* PMTcellLogical[],
  G4LogicalVolume* PMTcathLogical[], std::vector<G4LogicalVolume*> fiberLogical[], std::vector<G4LogicalVolume*> fiberLogical_[], std::vector<DRsimInterface::DRsimTowerProperty>& towerProps) {

  for (int i = 0; i < mNumBarrel; i++) {
    float towerTheta = fulltheta + fDThetaBarrel[i]/2.;
    G4ThreeVector pt[8];
    pDimB->SetDeltaTheta(fDThetaBarrel[i]);
    pDimB->SetThetaOfCenter(towerTheta);
    pDimB->init();
    pDimB->GetPt(pt);
    towerName = setTowerName(pDimB->GetIsRHS(), "B", i);

    tower = new G4Trap("TowerB",pt);
    towerLogical[i] = new G4LogicalVolume(tower,FindMaterial(mTowerMaterial),towerName);

    pDimB->GetPtSipm(pt);
    pmtg = new G4Trap("PMTGB",pt);
    PMTGLogical[i] = new G4LogicalVolume(pmtg,FindMaterial("G4_AIR"),towerName);

    for(int j=0;j<mNumZRot;j++){
      auto rot = pDimB->GetRotationMat(j);
      new G4PVPlacement(G4Transform3D(rot,pDimB->GetTowerPos(j)),towerLogical[i],towerName,worldLogical,false,j,checkOverlaps);
      new G4PVPlacement(G4Transform3D(rot,pDimB->GetSipmLayerPos(j)),PMTGLogical[i],towerName,worldLogical,false,j,checkOverlaps);
    }

    if (fDoFiber)
      fiberBarrel(i,fDThetaBarrel[i],towerLogical,fiberLogical,fiberLogical_);

    int iTheta = pDimB->GetIsRHS() ? i : -i-1;
    float signedTowerTheta = pDimB->GetIsRHS() ? towerTheta : -towerTheta;
    DRsimInterface::DRsimTowerProperty towerProp;
    towerProp.towerXY = fTowerXY;
    towerProp.towerTheta = std::make_pair(iTheta,signedTowerTheta);
    towerProp.innerR = pDimB->GetCurrentInnerR();
    towerProp.towerH = mTowerH;
    towerProp.dTheta = fDThetaBarrel[i];
    towerProps.push_back(towerProp);

    if (fDoFiber) {
      G4VSolid* SiPMlayerSolid = new G4Box("SiPMlayerSolid",fTowerXY.first*1.5/2.*mm,fTowerXY.second*1.5/2.*mm,fPMTT/2.);
      G4LogicalVolume* SiPMlayerLogical = new G4LogicalVolume(SiPMlayerSolid,FindMaterial("G4_POLYVINYL_CHLORIDE"),"SiPMlayerLogical");
      new G4PVPlacement(0,G4ThreeVector(0.,0.,fFilterT/2.),SiPMlayerLogical,"SiPMlayerPhysical",PMTGLogical[i],false,0,checkOverlaps);

      G4VSolid* filterlayerSolid = new G4Box("filterlayerSolid",fTowerXY.first*1.5/2.*mm,fTowerXY.second*1.5/2.*mm,fFilterT/2.);
      G4LogicalVolume* filterlayerLogical = new G4LogicalVolume(filterlayerSolid,FindMaterial("G4_POLYVINYL_CHLORIDE"),"filterlayerLogical");
      new G4PVPlacement(0,G4ThreeVector(0.,0.,-fPMTT/2.),filterlayerLogical,"filterlayerPhysical",PMTGLogical[i],false,0,checkOverlaps);

      G4VSolid* PMTcellSolid = new G4Box("PMTcellSolid",1.2/2.*mm,1.2/2.*mm,fPMTT/2.);
      PMTcellLogical[i] = new G4LogicalVolume(PMTcellSolid,FindMaterial("Glass"),"PMTcellLogical");

      DRsimCellParameterisation* PMTcellParam = new DRsimCellParameterisation(fTowerXY.first,fTowerXY.second);
      G4PVParameterised* PMTcellPhysical = new G4PVParameterised("PMTcellPhysical",PMTcellLogical[i],SiPMlayerLogical,kXAxis,fTowerXY.first*fTowerXY.second,PMTcellParam);

      G4VSolid* PMTcathSolid = new G4Box("PMTcathSolid",1.2/2.*mm,1.2/2.*mm,fSiPMT/2.*mm);
      PMTcathLogical[i] = new G4LogicalVolume(PMTcathSolid,FindMaterial("Silicon"),"PMTcathLogical");
      PMTcathLogical[i]->SetVisAttributes(fVisAttrGreen);
      new G4PVPlacement(0,G4ThreeVector(0.,0.,(fPMTT-fSiPMT)/2.*mm),PMTcathLogical[i],"PMTcathPhysical",PMTcellLogical[i],false,0,checkOverlaps);
      new G4LogicalSkinSurface("Photocath_surf",PMTcathLogical[i],FindSurface("SiPMSurf"));

      G4VSolid* filterSolid = new G4Box("filterSolid",1.2/2.*mm,1.2/2.*mm,fFilterT/2.);
      PMTfilterLogical[i] = new G4LogicalVolume(filterSolid,FindMaterial("Glass"),"PMTfilterLogical");

      DRsimFilterParameterisation* filterParam = new DRsimFilterParameterisation(fTowerXY.first,fTowerXY.second,FindMaterial("Glass"),FindMaterial("Gelatin"));
      filterParam->SetFilterVis(fVisAttrOrange);
      filterParam->SetGlassVis(fVisAttrWhite);

      G4PVParameterised* filterPhysical = new G4PVParameterised("filterPhysical",PMTfilterLogical[i],filterlayerLogical,kXAxis,fTowerXY.first*fTowerXY.second,filterParam);
      new G4LogicalBorderSurface("filterSurf",filterPhysical,PMTcellPhysical,FindSurface("FilterSurf"));
    }

    fulltheta = fulltheta+fDThetaBarrel[i];
  }
}

void DRsimDetectorConstruction::Endcap(G4LogicalVolume* towerLogical[], G4LogicalVolume* PMTGLogical[], G4LogicalVolume* PMTfilterLogical[], G4LogicalVolume* PMTcellLogical[],
  G4LogicalVolume* PMTcathLogical[], std::vector<G4LogicalVolume*> fiberLogical[], std::vector<G4LogicalVolume*> fiberLogical_[], std::vector<DRsimInterface::DRsimTowerProperty>& towerProps) {

  for(int i=0;i<mNumEndcap;i++) {
    float towerTheta = fulltheta + fDThetaEndcap/2.;
    G4ThreeVector pt[8];
    pDimE->SetThetaOfCenter(towerTheta);
    pDimE->init();
    pDimE->GetPt(pt);
    towerName = setTowerName(pDimE->GetIsRHS(), "E", i);

    tower = new G4Trap("TowerE",pt);
    towerLogical[i] = new G4LogicalVolume(tower,FindMaterial(mTowerMaterial),towerName);

    pDimE->GetPtSipm(pt);
    pmtg = new G4Trap("PMTGE",pt);
    PMTGLogical[i] = new G4LogicalVolume(pmtg,FindMaterial("G4_AIR"),towerName);

    for(int j=0;j<mNumZRot;j++){
      auto rot = pDimE->GetRotationMat(j);
      new G4PVPlacement(G4Transform3D(rot,pDimE->GetTowerPos(j)),towerLogical[i],towerName,worldLogical,false,j,checkOverlaps);
      new G4PVPlacement(G4Transform3D(rot,pDimE->GetSipmLayerPos(j)),PMTGLogical[i],towerName,worldLogical,false,j,checkOverlaps);
    }

    if (fDoFiber)
      fiberEndcap(i,fDThetaEndcap,towerLogical,fiberLogical,fiberLogical_);

    int iTheta = pDimE->GetIsRHS() ? i+mNumBarrel : -i-mNumBarrel-1;
    float signedTowerTheta = pDimE->GetIsRHS() ? towerTheta : -towerTheta;
    DRsimInterface::DRsimTowerProperty towerProp;
    towerProp.towerXY = fTowerXY;
    towerProp.towerTheta = std::make_pair(iTheta,signedTowerTheta);
    towerProp.innerR = pDimE->GetCurrentInnerR();
    towerProp.towerH = mTowerH;
    towerProp.dTheta = fDThetaEndcap;
    towerProps.push_back(towerProp);

    if (fDoFiber) {
      G4VSolid* SiPMlayerSolid = new G4Box("SiPMlayerSolid",fTowerXY.first*1.5/2.*mm,fTowerXY.second*1.5/2.*mm,fPMTT/2.);
      G4LogicalVolume* SiPMlayerLogical = new G4LogicalVolume(SiPMlayerSolid,FindMaterial("G4_POLYVINYL_CHLORIDE"),"SiPMlayerLogical");
      new G4PVPlacement(0,G4ThreeVector(0.,0.,fFilterT/2.),SiPMlayerLogical,"SiPMlayerPhysical",PMTGLogical[i],false,0,checkOverlaps);

      G4VSolid* filterlayerSolid = new G4Box("filterlayerSolid",fTowerXY.first*1.5/2.*mm,fTowerXY.second*1.5/2.*mm,fFilterT/2.);
      G4LogicalVolume* filterlayerLogical = new G4LogicalVolume(filterlayerSolid,FindMaterial("G4_POLYVINYL_CHLORIDE"),"filterlayerLogical");
      new G4PVPlacement(0,G4ThreeVector(0.,0.,-fPMTT/2.),filterlayerLogical,"filterlayerPhysical",PMTGLogical[i],false,0,checkOverlaps);

      G4VSolid* PMTcellSolid = new G4Box("PMTcellSolid",1.2/2.*mm,1.2/2.*mm,fPMTT/2.);
      PMTcellLogical[i] = new G4LogicalVolume(PMTcellSolid,FindMaterial("Glass"),"PMTcellLogical");

      DRsimCellParameterisation* PMTcellParam = new DRsimCellParameterisation(fTowerXY.first,fTowerXY.second);
      G4PVParameterised* PMTcellPhysical = new G4PVParameterised("PMTcellPhysical",PMTcellLogical[i],SiPMlayerLogical,kXAxis,fTowerXY.first*fTowerXY.second,PMTcellParam);

      G4VSolid* PMTcathSolid = new G4Box("PMTcathSolid",1.2/2.*mm,1.2/2.*mm,fSiPMT/2.*mm);
      PMTcathLogical[i] = new G4LogicalVolume(PMTcathSolid,FindMaterial("Silicon"),"PMTcathLogical");
      PMTcathLogical[i]->SetVisAttributes(fVisAttrGreen);
      new G4PVPlacement(0,G4ThreeVector(0.,0.,(fPMTT-fSiPMT)/2.*mm),PMTcathLogical[i],"PMTcathPhysical",PMTcellLogical[i],false,0,checkOverlaps);
      new G4LogicalSkinSurface("Photocath_surf",PMTcathLogical[i],FindSurface("SiPMSurf"));

      G4VSolid* filterSolid = new G4Box("filterSolid",1.2/2.*mm,1.2/2.*mm,fFilterT/2.);
      PMTfilterLogical[i] = new G4LogicalVolume(filterSolid,FindMaterial("Glass"),"PMTfilterLogical");

      DRsimFilterParameterisation* filterParam = new DRsimFilterParameterisation(fTowerXY.first,fTowerXY.second,FindMaterial("Glass"),FindMaterial("Gelatin"));
      filterParam->SetFilterVis(fVisAttrOrange);
      filterParam->SetGlassVis(fVisAttrWhite);

      G4PVParameterised* filterPhysical = new G4PVParameterised("filterPhysical",PMTfilterLogical[i],filterlayerLogical,kXAxis,fTowerXY.first*fTowerXY.second,filterParam);
      new G4LogicalBorderSurface("filterSurf",filterPhysical,PMTcellPhysical,FindSurface("FilterSurf"));
    }

    fulltheta = fulltheta+fDThetaEndcap;
  }
}

void DRsimDetectorConstruction::DefineCommands() {}

void DRsimDetectorConstruction::fiberBarrel(G4int i, G4double deltatheta_,G4LogicalVolume* towerLogical[], std::vector<G4LogicalVolume*> fiberLogical[], std::vector<G4LogicalVolume*> fiberLogical_[]) {

  fFiberX.clear();
  fFiberY.clear();

  v1 = pDimB->GetV1();
  v2 = pDimB->GetV2();
  v3 = pDimB->GetV3();
  v4 = pDimB->GetV4();

  innerSide_half = pDimB->GetCurrentInnerR()*tan(deltatheta_/2.);
  outerSide_half = (pDimB->GetCurrentInnerR()+mTowerH)*tan(deltatheta_/2.);

  int numx = (int)(((v4.getX()*tan(phi_unit/2.)*2)-1.2*mm)/(1.5*mm)) + 1;
  int numy = (int)((outerSide_half*2-1.2*mm)/(1.5*mm)) + 1;
  fTowerXY = std::make_pair(numx,numy);

  for (int j = 0; j < numy; j++) {
    for (int k = 0; k < numx; k++) {
      G4float fX = -1.5*mm*(numx/2) + k*1.5*mm + ( numx%2==0 ? 0.75*mm : 0 );
      G4float fY = -1.5*mm*(numy/2) + j*1.5*mm + ( numy%2==0 ? 0.75*mm : 0 );
      fFiberX.push_back(fX);
      fFiberY.push_back(fY);
    }
  }

  for (unsigned int j = 0; j<fFiberX.size();j++) {
    G4int column = j % numx;
    G4int row = j / numx;

    if ( DRsimInterface::IsCerenkov(column,row) ) { //c fibre
      intersect = new G4IntersectionSolid("fiber_",fiber,tower,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
      G4LogicalVolume* cladLogical = new G4LogicalVolume(intersect,FindMaterial("FluorinatedPolymer"),name);
      fiberLogical[i].push_back(cladLogical);
      new G4PVPlacement(0,G4ThreeVector(fFiberX.at(j),fFiberY.at(j),0),fiberLogical[i].at(j),name,towerLogical[i],false,j,checkOverlaps);

      intersect_ = new G4IntersectionSolid("fiber_",fiberC,tower,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
      G4LogicalVolume* coreLogical = new G4LogicalVolume(intersect_,FindMaterial("PMMA"),name);
      fiberLogical_[i].push_back(coreLogical);
      new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberLogical_[i].at(j),name,fiberLogical[i].at(j),false,j,checkOverlaps);

      fCerenRegion->AddRootLogicalVolume(cladLogical);
      fCerenRegion->AddRootLogicalVolume(coreLogical);
      cladLogical->SetRegion(fCerenRegion);
      coreLogical->SetRegion(fCerenRegion);

      fiberLogical[i].at(j)->SetVisAttributes(fVisAttrGray);
      fiberLogical_[i].at(j)->SetVisAttributes(fVisAttrBlue);
    } else { // s fibre
      intersect = new G4IntersectionSolid("fiber_",fiber,tower,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
      G4LogicalVolume* cladLogical = new G4LogicalVolume(intersect,FindMaterial("PMMA"),name);
      fiberLogical[i].push_back(cladLogical);
      new G4PVPlacement(0,G4ThreeVector(fFiberX.at(j),fFiberY.at(j),0),fiberLogical[i].at(j),name,towerLogical[i],false,j,checkOverlaps);

      intersect_ = new G4IntersectionSolid("fiber_",fiberS,tower,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
      G4LogicalVolume* coreLogical = new G4LogicalVolume(intersect_,FindMaterial("Polystyrene"),name);
      fiberLogical_[i].push_back(coreLogical);
      new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberLogical_[i].at(j),name,fiberLogical[i].at(j),false,j,checkOverlaps);

      fScintRegion->AddRootLogicalVolume(cladLogical);
      fScintRegion->AddRootLogicalVolume(coreLogical);
      cladLogical->SetRegion(fScintRegion);
      coreLogical->SetRegion(fScintRegion);

      fiberLogical[i].at(j)->SetVisAttributes(fVisAttrGray);
      fiberLogical_[i].at(j)->SetVisAttributes(fVisAttrOrange);
    }
  }
}

void DRsimDetectorConstruction::fiberEndcap(G4int i, G4double deltatheta_, G4LogicalVolume* towerLogical[], std::vector<G4LogicalVolume*> fiberLogical[], std::vector<G4LogicalVolume*> fiberLogical_[]) {

  fFiberX.clear();
  fFiberY.clear();

  v1 = pDimE->GetV1();
  v2 = pDimE->GetV2();
  v3 = pDimE->GetV3();
  v4 = pDimE->GetV4();

  innerSide_half = pDimE->GetCurrentInnerR()*tan(deltatheta_/2.);
  outerSide_half = (pDimE->GetCurrentInnerR()+mTowerH)*tan(deltatheta_/2.);

  int numx = (int)(((v4.getX()*tan(phi_unit/2.)*2)-1.2*mm)/(1.5*mm)) + 1;
  int numy = (int)((outerSide_half*2-1.2*mm)/(1.5*mm)) + 1;
  fTowerXY = std::make_pair(numx,numy);

  for (int j = 0; j < numy; j++) {
    for (int k = 0; k < numx; k++) {
      G4float fX = -1.5*mm*(numx/2) + k*1.5*mm + ( numx%2==0 ? 0.75*mm : 0 );
      G4float fY = -1.5*mm*(numy/2) + j*1.5*mm + ( numy%2==0 ? 0.75*mm : 0 );
      fFiberX.push_back(fX);
      fFiberY.push_back(fY);
    }
  }

  for (unsigned int j = 0; j<fFiberX.size();j++) {
    G4int column = j % numx;
    G4int row = j / numx;

    if ( DRsimInterface::IsCerenkov(column,row) ) { //c fibre
      intersect = new G4IntersectionSolid("fiber_",fiber,tower,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
      G4LogicalVolume* cladLogical = new G4LogicalVolume(intersect,FindMaterial("FluorinatedPolymer"),name);
      fiberLogical[i].push_back(cladLogical);
      new G4PVPlacement(0,G4ThreeVector(fFiberX.at(j),fFiberY.at(j),0),fiberLogical[i].at(j),name,towerLogical[i],false,j,checkOverlaps);

      intersect_ = new G4IntersectionSolid("fiber_",fiberC,tower,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
      G4LogicalVolume* coreLogical = new G4LogicalVolume(intersect_,FindMaterial("PMMA"),name);
      fiberLogical_[i].push_back(coreLogical);
      new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberLogical_[i].at(j),name,fiberLogical[i].at(j),false,j,checkOverlaps);

      fCerenRegion->AddRootLogicalVolume(cladLogical);
      fCerenRegion->AddRootLogicalVolume(coreLogical);
      cladLogical->SetRegion(fCerenRegion);
      coreLogical->SetRegion(fCerenRegion);

      fiberLogical[i].at(j)->SetVisAttributes(fVisAttrGray);
      fiberLogical_[i].at(j)->SetVisAttributes(fVisAttrBlue);
    } else { // s fibre
      intersect = new G4IntersectionSolid("fiber_",fiber,tower,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
      G4LogicalVolume* cladLogical = new G4LogicalVolume(intersect,FindMaterial("PMMA"),name);
      fiberLogical[i].push_back(cladLogical);
      new G4PVPlacement(0,G4ThreeVector(fFiberX.at(j),fFiberY.at(j),0),fiberLogical[i].at(j),name,towerLogical[i],false,j,checkOverlaps);

      intersect_ = new G4IntersectionSolid("fiber_",fiberS,tower,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
      G4LogicalVolume* coreLogical = new G4LogicalVolume(intersect_,FindMaterial("Polystyrene"),name);
      fiberLogical_[i].push_back(coreLogical);
      new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberLogical_[i].at(j),name,fiberLogical[i].at(j),false,j,checkOverlaps);

      fScintRegion->AddRootLogicalVolume(cladLogical);
      fScintRegion->AddRootLogicalVolume(coreLogical);
      cladLogical->SetRegion(fScintRegion);
      coreLogical->SetRegion(fScintRegion);

      fiberLogical[i].at(j)->SetVisAttributes(fVisAttrGray);
      fiberLogical_[i].at(j)->SetVisAttributes(fVisAttrOrange);
    }
  }
}
