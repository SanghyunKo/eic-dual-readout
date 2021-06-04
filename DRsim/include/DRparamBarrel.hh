#ifndef DRparamBarrel_h
#define DRparamBarrel_h 1

#include "DRparamBase.hh"

#include "G4ThreeVector.hh"

#include <vector>
#include <cmath>

class DRparamBarrel : public DRparamBase {
public:
  DRparamBarrel();
  virtual ~DRparamBarrel();

  virtual void SetDeltaThetaByTowerNo(int signedTowerNo, int) override;
  virtual void SetThetaOfCenterByTowerNo(int signedTowerNo, int) override;

  virtual void init() override;
};

#endif
