import FWCore.ParameterSet.Config as cms

RecoParticleFlowAOD = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*',
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*',
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*',
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*',
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*',
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*',
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*',
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*',
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*',
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*',
        'keep recoCaloClusters_particleFlowEGamma_*_*',
        'keep recoSuperClusters_particleFlowEGamma_*_*',
        'keep recoCaloClusters_particleFlowSuperClusterECAL_*_*',
        'keep recoSuperClusters_particleFlowSuperClusterECAL_*_*',
        'keep recoConversions_particleFlowEGamma_*_*',
        'keep recoPFCandidates_particleFlow_*_*',
        'keep recoPFCandidates_particleFlowTmp_AddedMuonsAndHadrons_*',
        'keep recoPFCandidates_particleFlowTmp_CleanedCosmicsMuons_*',
        'keep recoPFCandidates_particleFlowTmp_CleanedFakeMuons_*',
        'keep recoPFCandidates_particleFlowTmp_CleanedHF_*',
        'keep recoPFCandidates_particleFlowTmp_CleanedPunchThroughMuons_*',
        'keep recoPFCandidates_particleFlowTmp_CleanedPunchThroughNeutralHadrons_*',
        'keep recoPFCandidates_particleFlowTmp_CleanedTrackerAndGlobalMuons_*',
        'keep *_particleFlow_electrons_*',
        'keep *_particleFlow_photons_*',
        'keep *_particleFlow_muons_*',
        'keep recoCaloClusters_pfElectronTranslator_*_*',
        'keep recoPreshowerClusters_pfElectronTranslator_*_*',
        'keep recoSuperClusters_pfElectronTranslator_*_*',
        'keep recoCaloClusters_pfPhotonTranslator_*_*',
        'keep recoPreshowerClusters_pfPhotonTranslator_*_*',
        'keep recoSuperClusters_pfPhotonTranslator_*_*',
        'keep recoPhotons_pfPhotonTranslator_*_*',
        'keep recoPhotonCores_pfPhotonTranslator_*_*',
        'keep recoConversions_pfPhotonTranslator_*_*',
        'keep *_particleFlowPtrs_*_*',
        'keep *_particleFlowTmpPtrs_*_*',
        'keep *_chargedHadronPFTrackIsolation_*_*',
        'keep recoPFRecHits_particleFlowRecHitHGC_Cleaned_*',
        'keep recoPFClusters_particleFlowClusterHGCal__*',
        'keep recoPFClusters_particleFlowClusterHGCalFromMultiCl__*',
        'keep recoSuperClusters_simPFProducer_*_*',
        'keep *_ecalBarrelClusterFastTimer_*_*'
    )
)