import FWCore.ParameterSet.Config as cms

muonSeededTracksInOutSelector = cms.EDProducer("MultiTrackSelector",
    beamspot = cms.InputTag("offlineBeamSpot"),
    src = cms.InputTag("muonSeededTracksInOut"),
    trackSelectors = cms.VPSet(
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(False),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(10.0),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.55, 4.0),
            d0_par2 = cms.vdouble(0.55, 4.0),
            dz_par1 = cms.vdouble(0.65, 4.0),
            dz_par2 = cms.vdouble(0.45, 4.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(4),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(7),
            minNumber3DLayers = cms.uint32(0),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(5),
            nSigmaZ = cms.double(4.0),
            name = cms.string('muonSeededTracksInOutLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.01),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(False),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(1.0),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.3, 4.0),
            d0_par2 = cms.vdouble(0.4, 4.0),
            dz_par1 = cms.vdouble(0.35, 4.0),
            dz_par2 = cms.vdouble(0.4, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(3),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(10),
            minNumber3DLayers = cms.uint32(2),
            minNumberLayers = cms.uint32(5),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(6),
            nSigmaZ = cms.double(4.0),
            name = cms.string('muonSeededTracksInOutTight'),
            preFilterName = cms.string('muonSeededTracksInOutLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.01),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(False),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.4),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.3, 4.0),
            d0_par2 = cms.vdouble(0.4, 4.0),
            dz_par1 = cms.vdouble(0.35, 4.0),
            dz_par2 = cms.vdouble(0.4, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(2),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(2),
            minNumberLayers = cms.uint32(5),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(7),
            nSigmaZ = cms.double(4.0),
            name = cms.string('muonSeededTracksInOutHighPurity'),
            preFilterName = cms.string('muonSeededTracksInOutTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        )
    ),
    useVertices = cms.bool(True),
    useVtxError = cms.bool(False),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)
