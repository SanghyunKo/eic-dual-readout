#ifndef DRparamEndcap_h
#define DRparamEndcap_h 1

#include "DRparamBase.hh"

#include "G4ThreeVector.hh"

#include <vector>
#include <cmath>

class DRparamEndcap : public DRparamBase {
public:
  DRparamEndcap();
  virtual ~DRparamEndcap();

  virtual void SetDeltaThetaByTowerNo(int signedTowerNo, int BEtrans) override;
  virtual void SetThetaOfCenterByTowerNo(int signedTowerNo, int BEtrans) override;

  virtual void init() override;
};

#endif
