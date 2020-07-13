import FWCore.ParameterSet.Config as cms

def getHLTTDRProductsToDrop():
    return [
"drop *_hltGtStage2Digis_*_*",
"drop *_simBmtfDigis_*_*",
"drop *_simCaloStage2Digis_*_*",
"drop *_simCaloStage2Layer1Digis_*_*",
"drop *_simEmtfDigis_*_*",
"drop *_simGmtStage2Digis_*_*",
"drop *_simGtStage2Digis_*_*",
"drop *_simOmtfDigis_*_*",
"drop l1tHGCalTriggerCellBXVector_hgcalVFEProducer_HGCalVFEProcessorSums_*",
"drop recoHGCalMultiClusters_ticlMultiClustersFromTrackstersMerge__*",
"drop recoHGCalMultiClusters_ticlMultiClustersFromTrackstersTrk__*",
"drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTClusteredmNewDetSetVector_TTClustersFromPhase2TrackerDigis_ClusterInclusive_HLT",
"drop recoHGCalMultiClusters_ticlMultiClustersFromTrackstersMIP__*",
"drop recoPFClusters_particleFlowClusterHGCalFromMultiCl__*",
"drop l1tHGCalClusterBXVector_hgcalBackEndLayer1Producer_HGCalBackendLayer1Processor2DClustering_*",
"drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTClusterAssociationMap_TTClusterAssociatorFromPixelDigis_ClusterInclusive_HLT",
"drop l1tHGCalTowerMapBXVector_hgcalTowerMapProducer_HGCalTowerMapProcessor_*",
"drop recoGsfTrackExtras_electronGsfTracks__*",
"drop recoCaloClusters_particleFlowSuperClusterHGCal__*",
"drop l1tHGCalTriggerCellBXVector_hgcalConcentratorProducer_HGCalConcentratorProcessorSelection_*",
"drop recoSuperClusters_particleFlowSuperClusterHGCal__*",
"drop recoGsfTrackExtras_electronGsfTracksFromMultiCl__*",
"drop recoHGCalMultiClusters_hgcalMultiClusters__*",
"drop recoGsfTrackExtras_electronGsfTracks__*",
"drop recoCaloClusters_particleFlowSuperClusterHGCalFromMultiCl__*",
"drop recoSuperClusters_particleFlowSuperClusterHGCalFromMultiCl__*",
"drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTClusteredmNewDetSetVector_TTStubsFromPhase2TrackerDigis_ClusterAccepted_HLT",
"drop recoElectronSeeds_electronMergedSeedsFromMultiCl__*",
"drop recoTrackExtras_electronGsfTracks__*",
"drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTClusterAssociationMap_TTClusterAssociatorFromPixelDigis_ClusterAccepted_HLT",
"drop recoHGCalMultiClusters_ticlMultiClustersFromTrackstersHAD__*",
"drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTStubedmNewDetSetVector_TTStubsFromPhase2TrackerDigis_StubAccepted_HLT",
"drop recoTrackExtras_electronGsfTracksFromMultiCl__*",
"drop recoGsfElectrons_ecalDrivenGsfElectronsFromMultiCl__*",
"drop l1tHGCalMulticlusterBXVector_hgcalBackEndLayer2Producer_HGCalBackendLayer2Processor3DClustering_*",
"drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTStubAssociationMap_TTStubAssociatorFromPixelDigis_StubAccepted_HLT",
"drop recoGsfTracks_electronGsfTracks__*",
"drop CaloTowersSorted_towerMaker__*",
"drop recoGsfTracks_electronGsfTracksFromMultiCl__*",
"drop l1tHGCalTowerBXVector_hgcalTowerProducer_HGCalTowerProcessor_*",
"drop TrackingRecHitsOwned_electronGsfTracks__*",
"drop Phase2TrackerDigiedmDetSetVectorPhase2TrackerDigiPhase2TrackerDigiedmrefhelperFindForDetSetVectoredmRefTTTracks_TTTracksFromExtendedTrackletEmulation_Level1TTTracks_*",
"drop recoHGCalMultiClusters_ticlMultiClustersFromTrackstersEM__*",
]

def dropInputProducts(source):
    if not hasattr(source,"inputCommands"):
        source.inputCommands = cms.untracked.vstring("keep *")
    source.inputCommands.extend(getHLTTDRProductsToDrop())
    source.dropDescendantsOfDroppedBranches = cms.untracked.bool(False)

def dropOutputProducts(output):
    if hasattr(output,"outputCommands"):
        output.outputCommands = cms.untracked.vstring("drop *")
    output.outputCommands.extend(getHLTTDRProductsToDrop())        


