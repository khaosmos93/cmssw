#include "RecoMuon/TrackerSeedGenerator/interface/SeedMvaEstimator.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CommonTools/MVAUtils/interface/GBRForestTools.h"

SeedMvaEstimator::SeedMvaEstimator(const edm::FileInPath& weightsfile,
                                   const std::vector<double>& scale_mean,
                                   const std::vector<double>& scale_std) {
  gbrForest_ = createGBRForest(weightsfile);
  scale_mean_ = scale_mean;
  scale_std_ = scale_std;
}

SeedMvaEstimator::~SeedMvaEstimator() {}

namespace {
  enum inputIndexes {
    kTsosErr0,       // 0
    kTsosErr2,       // 1
    kTsosErr5,       // 2
    kTsosDxdz,       // 3
    kTsosDydz,       // 4
    kTsosQbp,        // 5
    kDRdRL1SeedP,    // 6
    kDPhidRL1SeedP,  // 7
    kLastL1,         // 8

    kDRdRL2SeedP = 8,  // 8
    kDPhidRL2SeedP,    // 9
    kLastL2,           // 10
  };
}  // namespace

void SeedMvaEstimator::getL1MuonVariables(const GlobalVector& global_p,
                                          const l1t::MuonBxCollection& l1Muons,
                                          int minL1Qual,
                                          float& dR2dRL1SeedP,
                                          float& dPhidRL1SeedP) const {
  for (int ibx = l1Muons.getFirstBX(); ibx <= l1Muons.getLastBX(); ++ibx) {
    if (ibx != 0)
      continue;  // -- only take when ibx == 0 -- //

    for (auto it = l1Muons.begin(ibx); it != l1Muons.end(ibx); it++) {
      if (it->hwQual() < minL1Qual)
        continue;

      float dR2tmp = reco::deltaR2(it->etaAtVtx(), it->phiAtVtx(), global_p.eta(), global_p.phi());
      if (dR2tmp < dR2dRL1SeedP) {
        dR2dRL1SeedP = dR2tmp;
        dPhidRL1SeedP = reco::deltaPhi(it->phiAtVtx(), global_p.phi());
      }
    }
  }
}

void SeedMvaEstimator::getL2MuonVariables(const GlobalVector& global_p,
                                          const reco::RecoChargedCandidateCollection& l2Muons,
                                          float& dR2dRL2SeedP,
                                          float& dPhidRL2SeedP) const {
  for (auto it = l2Muons.begin(); it != l2Muons.end(); it++) {
    float dR2tmp = reco::deltaR2(*it, global_p);
    if (dR2tmp < dR2dRL2SeedP) {
      dR2dRL2SeedP = dR2tmp;
      dPhidRL2SeedP = reco::deltaPhi(it->phi(), global_p.phi());
    }
  }
}

double SeedMvaEstimator::computeMva(const TrajectorySeed& seed,
                                    const GlobalVector& global_p,
                                    const l1t::MuonBxCollection& l1Muons,
                                    int minL1Qual,
                                    const reco::RecoChargedCandidateCollection& l2Muons,
                                    bool isFromL1) const {
  static constexpr float initDRdPhi(99999.);
  const int kLast = isFromL1 ? kLastL1 : kLastL2;
  float varL1[kLastL1]{};
  float varL2[kLastL2]{};
  float* var = isFromL1 ? varL1 : varL2;

  var[kTsosErr0] = seed.startingState().error(0);
  var[kTsosErr2] = seed.startingState().error(2);
  var[kTsosErr5] = seed.startingState().error(5);
  var[kTsosDxdz] = seed.startingState().parameters().dxdz();
  var[kTsosDydz] = seed.startingState().parameters().dydz();
  var[kTsosQbp] = seed.startingState().parameters().qbp();

  float dR2dRL1SeedP = initDRdPhi;
  float dPhidRL1SeedP = initDRdPhi;
  getL1MuonVariables(global_p, l1Muons, minL1Qual, dR2dRL1SeedP, dPhidRL1SeedP);

  var[kDRdRL1SeedP] = std::sqrt(dR2dRL1SeedP);
  var[kDPhidRL1SeedP] = dPhidRL1SeedP;

  if (!isFromL1) {
    float dR2dRL2SeedP = initDRdPhi;
    float dPhidRL2SeedP = initDRdPhi;
    getL2MuonVariables(global_p, l2Muons, dR2dRL2SeedP, dPhidRL2SeedP);

    var[kDRdRL2SeedP] = std::sqrt(dR2dRL2SeedP);
    var[kDPhidRL2SeedP] = dPhidRL2SeedP;
  }

  for (int iv = 0; iv < kLast; ++iv) {
    var[iv] = (var[iv] - scale_mean_.at(iv)) / scale_std_.at(iv);
  }

  return gbrForest_->GetResponse(var);
}
