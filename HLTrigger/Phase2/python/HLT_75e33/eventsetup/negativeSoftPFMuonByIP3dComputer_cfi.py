import FWCore.ParameterSet.Config as cms

negativeSoftPFMuonByIP3dComputer = cms.ESProducer("LeptonTaggerByIPESProducer",
    ipSign = cms.string('negative'),
    use3d = cms.bool(True)
)
