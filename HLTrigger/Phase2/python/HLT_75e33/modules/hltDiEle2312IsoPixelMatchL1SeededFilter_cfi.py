import FWCore.ParameterSet.Config as cms

hltDiEle2312IsoPixelMatchL1SeededFilter = cms.EDFilter("HLTElectronPixelMatchFilter",
    candTag = cms.InputTag("hltDiEG2312IsoHcalIsoL1SeededFilter"),
    l1EGCand = cms.InputTag("hltEgammaCandidatesL1Seeded"),
    l1PixelSeedsTag = cms.InputTag("hltEgammaElectronPixelSeedsL1Seeded"),
    ncandcut = cms.int32(2),
    npixelmatchcut = cms.double(1.0),
    pixelVeto = cms.bool(False),
    s2_threshold = cms.double(0.4),
    s_a_phi1B = cms.double(0.0069),
    s_a_phi1F = cms.double(0.0076),
    s_a_phi1I = cms.double(0.0088),
    s_a_phi2B = cms.double(0.00037),
    s_a_phi2F = cms.double(0.00906),
    s_a_phi2I = cms.double(0.0007),
    s_a_rF = cms.double(0.04),
    s_a_rI = cms.double(0.027),
    s_a_zB = cms.double(0.012),
    saveTags = cms.bool(True),
    tanhSO10BarrelThres = cms.double(0.35),
    tanhSO10ForwardThres = cms.double(1.0),
    tanhSO10InterThres = cms.double(1.0),
    useS = cms.bool(False)
)
