import FWCore.ParameterSet.Config as cms

RECOEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( (
        'drop *', 
        'keep DetIdedmEDCollection_siStripDigis_*_*', 
        'keep DetIdedmEDCollection_siPixelDigis_*_*', 
        'keep PixelFEDChanneledmNewDetSetVector_siPixelDigis_*_*', 
        'keep *_siPixelClusters_*_*', 
        'keep *_siStripClusters_*_*', 
        'keep ClusterSummary_clusterSummaryProducer_*_*', 
        'keep *_siPhase2Clusters_*_*', 
        'keep *_dt1DRecHits_*_*', 
        'keep *_dt1DCosmicRecHits_*_*', 
        'keep *_csc2DRecHits_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep *_gemRecHits_*_*', 
        'keep *_gemSegments_*_*', 
        'keep *_me0RecHits_*_*', 
        'keep *_me0Segments_*_*', 
        'keep *_hbhereco_*_*', 
        'keep *_hbheprereco_*_*', 
        'keep *_hfprereco_*_*', 
        'keep *_hfreco_*_*', 
        'keep *_horeco_*_*', 
        'keep HBHERecHitsSorted_hbherecoMB_*_*', 
        'keep HORecHitsSorted_horecoMB_*_*', 
        'keep HFRecHitsSorted_hfrecoMB_*_*', 
        'keep ZDCDataFramesSorted_hcalDigis_*_*', 
        'keep ZDCDataFramesSorted_castorDigis_*_*', 
        'keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_ZDC_*', 
        'keep ZDCRecHitsSorted_zdcreco_*_*', 
        'keep *_castorreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_HGCalRecHit_*_*', 
        'keep recoCaloClusters_hgcalLayerClusters_*_*', 
        'keep *_hgcalLayerClusters_timeLayerCluster_*', 
        'keep *_hgcalLayerClusters_InitialLayerClustersMask_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep *_ecalCompactTrigPrim_*_*', 
        'keep *_ecalTPSkim_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*', 
        'keep *_mix_EBTimeDigi_*', 
        'keep *_mix_EETimeDigi_*', 
        'keep *_ecalDetailedTimeRecHit_*_*', 
        'keep *_hgcalMultiClusters_*_*', 
        'keep *_iterHGCalMultiClusters_*_*', 
        'keep *_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep *_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep *_particleFlowSuperClusterECAL_*_*', 
        'keep *_particleFlowSuperClusterOOTECAL_*_*', 
        'drop recoClusterShapes_*_*_*', 
        'drop recoBasicClustersToOnerecoClusterShapesAssociation_*_*_*', 
        'drop recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*', 
        'drop recoSuperClusters_multi5x5SuperClusters_multi5x5BarrelSuperClusters_*', 
        'keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep recoCaloClusters_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_hybridSuperClusters_uncleanOnlyHybridSuperClusters_*', 
        'keep recoCaloClusters_multi5x5SuperClusters_multi5x5EndcapBasicClusters_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterOOTECAL_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterOOTECAL_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterHGCal__*', 
        'keep recoCaloClusters_particleFlowSuperClusterHGCal__*', 
        'keep recoSuperClusters_particleFlowSuperClusterHGCalFromMultiCl__*', 
        'keep recoCaloClusters_particleFlowSuperClusterHGCalFromMultiCl__*', 
        'keep *_particleFlowSuperClusterHGCal_*_*', 
        'keep *_particleFlowSuperClusterHGCalFromMultiCl_*_*', 
        'keep *_CkfElectronCandidates_*_*', 
        'keep *_GsfGlobalElectronTest_*_*', 
        'keep *_electronMergedSeeds_*_*', 
        'keep recoGsfTrackExtras_electronGsfTracks_*_*', 
        'keep recoTrackExtras_electronGsfTracks_*_*', 
        'keep TrackingRecHitsOwned_electronGsfTracks_*_*', 
        'keep recoTracks_GsfGlobalElectronTest_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoGsfTracks_electronGsfTracksFromMultiCl_*_*', 
        'keep recoGsfTracks_electronGsfTracksFromMultiCl_*_*', 
        'keep recoGsfTrackExtras_electronGsfTracksFromMultiCl_*_*', 
        'keep recoTrackExtras_electronGsfTracksFromMultiCl_*_*', 
        'keep TrackingRecHitsOwned_electronGsfTracksFromMultiCl_*_*', 
        'keep *_electronMergedSeedsFromMultiCl_*_*', 
        'keep recoTrackExtras_generalTracks_*_*', 
        'keep TrackingRecHitsOwned_generalTracks_*_*', 
        'keep TrackingRecHitsOwned_extraFromSeeds_*_*', 
        'keep uints_extraFromSeeds_*_*', 
        'keep recoTrackExtras_beamhaloTracks_*_*', 
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*', 
        'keep recoTrackExtras_conversionStepTracks_*_*', 
        'keep TrackingRecHitsOwned_conversionStepTracks_*_*', 
        'keep *_ctfPixelLess_*_*', 
        'keep *_dedxTruncated40_*_*', 
        'keep recoTracks_cosmicDCTracks_*_*', 
        'keep recoTrackExtras_cosmicDCTracks_*_*', 
        'keep TrackingRecHitsOwned_cosmicDCTracks_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTracks_ctfPixelLess_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_towerMaker_*_*', 
        'keep *_ak4JetTracksAssociatorAtCaloFace_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep recoCaloJets_ak4CaloJets_*_*', 
        'keep *_ak4CaloJets_rho_*', 
        'keep *_ak4CaloJets_sigma_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep floatedmValueMap_puppi_*_*', 
        'keep *_ak4PFJetsPuppi_*_*', 
        'keep *_ak8PFJetsPuppi_*_*', 
        'keep *_ak8PFJetsPuppiSoftDrop_*_*', 
        'keep recoPFJets_ak4PFJets_*_*', 
        'keep *_ak4PFJets_rho_*', 
        'keep *_ak4PFJets_sigma_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep recoBasicJets_ak5CastorJets_*_*', 
        'keep *_ak5CastorJets_rho_*', 
        'keep *_ak5CastorJets_sigma_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep recoBasicJets_ak7CastorJets_*_*', 
        'keep *_ak7CastorJets_rho_*', 
        'keep *_ak7CastorJets_sigma_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'keep *_ak8PFJetsPuppiSoftDropMass_*_*', 
        'keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep recoEcalHaloData_EcalHaloData_*_*', 
        'keep recoHcalHaloData_HcalHaloData_*_*', 
        'keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep floatedmValueMap_puppiNoLep_*_*', 
        'keep recoPFMETs_pfMetPuppi_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep *_MuonSeed_*_*', 
        'keep *_ancientMuonSeed_*_*', 
        'keep *_displacedMuonSeeds_*_*', 
        'keep TrackingRecHitsOwned_globalMuons_*_*', 
        'keep TrackingRecHitsOwned_tevMuons_*_*', 
        'keep *_CosmicMuonSeed_*_*', 
        'keep recoTrackExtras_cosmicMuons_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons_*_*', 
        'keep recoTrackExtras_cosmicMuons1Leg_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons1Leg_*_*', 
        'keep recoTracks_cosmicsVetoTracks_*_*', 
        'keep recoMuons_muons_*_*', 
        'keep booledmValueMap_muons_*_*', 
        'keep doubleedmValueMap_muons_muPFMean*_*', 
        'keep doubleedmValueMap_muons_muPFSum*_*', 
        'keep *_muons_muonShowerInformation_*', 
        'keep recoMuonTimeExtraedmValueMap_muons_*_*', 
        'keep recoMuonCosmicCompatibilityedmValueMap_muons_*_*', 
        'keep uintedmValueMap_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*', 
        'keep *_softPFMuonsTagInfos_*_*', 
        'keep *_softPFElectronsTagInfos_*_*', 
        'keep *_pfImpactParameterTagInfos_*_*', 
        'keep *_pfSecondaryVertexTagInfos_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderTagInfos_*_*', 
        'keep *_pfGhostTrackVertexTagInfos_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderCvsLTagInfos_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*', 
        'keep *_pixelClusterTagInfos_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseElectronRejection_*_*', 
        'keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauBasicDiscriminators_*_*', 
        'keep *_hpsPFTauBasicDiscriminatorsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByMuonRejection3_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6ElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*', 
        'keep *_offlinePrimaryVertices4D__*', 
        'keep *_offlinePrimaryVertices4DWithBS__*', 
        'keep *_trackTimeValueMapProducer_*_*', 
        'keep *_offlinePrimaryVertices4DnoPID__*', 
        'keep *_offlinePrimaryVertices4DnoPIDWithBS__*', 
        'keep *_tofPID_*_*', 
        'keep *_gedPhotonCore_*_*', 
        'keep *_gedPhotons_*_*', 
        'keep recoPhotons_mustachePhotons_*_*', 
        'keep recoPhotonCores_mustachePhotonCore_*_*', 
        'keep recoTrackExtras_ckfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep recoPhotonCores_gedPhotonCore_*_*', 
        'keep recoPhotons_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoConversions_conversions_*_*', 
        'drop recoConversions_conversions_uncleanedConversions_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'keep *_gsfTracksOpenConversions_*_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep *_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*', 
        'keep recoCaloClusters_lowPtGsfElectronSuperClusters_*_*', 
        'keep recoGsfElectrons_lowPtGsfElectrons_*_*', 
        'keep recoGsfElectronCores_lowPtGsfElectronCores_*_*', 
        'keep recoGsfTracks_lowPtGsfEleGsfTracks_*_*', 
        'keep *_lowPtGsfToTrackLinks_*_*', 
        'keep recoSuperClusters_lowPtGsfElectronSuperClusters_*_*', 
        'keep floatedmValueMap_lowPtGsfElectronSeedValueMaps_*_*', 
        'keep floatedmValueMap_lowPtGsfElectronID_*_*', 
        'keep *_ecalDrivenGsfElectronCores_*_*', 
        'keep *_ecalDrivenGsfElectrons_*_*', 
        'keep *_ecalDrivenGsfElectronCoresFromMultiCl_*_*', 
        'keep *_ecalDrivenGsfElectronsFromMultiCl_*_*', 
        'keep *_photonCoreFromMultiCl_*_*', 
        'keep *_photonsFromMultiCl_*_*', 
        'keep *_pixelTracks_*_*', 
        'keep *_pixelVertices_*_*', 
        'keep recoPFClusters_particleFlowClusterECAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHCAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHO_*_*', 
        'keep recoPFClusters_particleFlowClusterHF_*_*', 
        'keep recoPFClusters_particleFlowClusterPS_*_*', 
        'keep recoPFBlocks_particleFlowBlock_*_*', 
        'keep recoPFCandidates_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlowTmp_electrons_*', 
        'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_pfPhotonTranslator_*_*', 
        'keep *_trackerDrivenElectronSeeds_preid_*', 
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
        'keep *_ecalBarrelClusterFastTimer_*_*', 
        'keep *_particleFlowSuperClusterHGCalFromMultiCl_*_*', 
        'keep recoPFBlocks_simPFProducer_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep L1MuGMTReadoutCollection_gtDigis_*_*', 
        'keep L1GctEmCand*_gctDigis_*_*', 
        'keep L1GctJetCand*_gctDigis_*_*', 
        'keep L1GctEtHad*_gctDigis_*_*', 
        'keep L1GctEtMiss*_gctDigis_*_*', 
        'keep L1GctEtTotal*_gctDigis_*_*', 
        'keep L1GctHtMiss*_gctDigis_*_*', 
        'keep L1GctJetCounts*_gctDigis_*_*', 
        'keep L1GctHFRingEtSums*_gctDigis_*_*', 
        'keep L1GctHFBitCounts*_gctDigis_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'keep *_gtStage2Digis_*_*', 
        'keep *_gmtStage2Digis_*_*', 
        'keep *_caloStage2Digis_*_*', 
        'drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep *_hltScoutingTrackPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep CTPPSRecord_onlineMetaDataDigis_*_*', 
        'keep DCSRecord_onlineMetaDataDigis_*_*', 
        'keep OnlineLuminosityRecord_onlineMetaDataDigis_*_*', 
        'keep recoBeamSpot_onlineMetaDataDigis_*_*', 
        'keep *_tcdsDigis_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoSingleTauDiscriminatorContaineredmValueMap_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep TotemTriggerCounters_totemTriggerRawToDigi_*_*', 
        'keep TotemFEDInfos_totemRPRawToDigi_*_*', 
        'keep TotemRPDigiedmDetSetVector_totemRPRawToDigi_*_*', 
        'keep TotemVFATStatusedmDetSetVector_totemRPRawToDigi_*_*', 
        'keep TotemRPClusteredmDetSetVector_totemRPClusterProducer_*_*', 
        'keep TotemRPRecHitedmDetSetVector_totemRPRecHitProducer_*_*', 
        'keep TotemRPUVPatternedmDetSetVector_totemRPUVPatternFinder_*_*', 
        'keep TotemRPLocalTrackedmDetSetVector_totemRPLocalTrackFitter_*_*', 
        'keep TotemFEDInfos_ctppsDiamondRawToDigi_*_*', 
        'keep CTPPSDiamondDigiedmDetSetVector_ctppsDiamondRawToDigi_*_*', 
        'keep TotemVFATStatusedmDetSetVector_ctppsDiamondRawToDigi_*_*', 
        'keep CTPPSDiamondRecHitedmDetSetVector_ctppsDiamondRecHits_*_*', 
        'keep CTPPSDiamondLocalTrackedmDetSetVector_ctppsDiamondLocalTracks_*_*', 
        'keep TotemTimingDigiedmDetSetVector_totemTimingRawToDigi_*_*', 
        'keep TotemTimingRecHitedmDetSetVector_totemTimingRecHits_*_*', 
        'keep TotemTimingLocalTrackedmDetSetVector_totemTimingLocalTracks_*_*', 
        'keep CTPPSPixelDigiedmDetSetVector_ctppsPixelDigis_*_*', 
        'keep CTPPSPixelDataErroredmDetSetVector_ctppsPixelDigis_*_*', 
        'keep CTPPSPixelClusteredmDetSetVector_ctppsPixelClusters_*_*', 
        'keep CTPPSPixelRecHitedmDetSetVector_ctppsPixelRecHits_*_*', 
        'keep CTPPSPixelLocalTrackedmDetSetVector_ctppsPixelLocalTracks_*_*', 
        'keep CTPPSLocalTrackLites_ctppsLocalTrackLiteProducer_*_*', 
        'keep recoForwardProtons_ctppsProtons_*_*', 
        'keep *_ticlTrackstersTrkEM_*_*', 
        'keep *_ticlTrackstersEM_*_*', 
        'keep *_ticlTrackstersHAD_*_*', 
        'keep *_ticlTrackstersTrk_*_*', 
        'keep *_ticlTrackstersMIP_*_*', 
        'keep *_ticlTrackstersMerge_*_*', 
        'keep *_ticlTrackstersHFNoseEM_*_*', 
        'keep *_ticlTrackstersHFNoseMIP_*_*', 
        'keep *_ticlTrackstersHFNoseMerge_*_*', 
        'keep *_pfTICL_*_*', 
        'keep *_ticlMultiClustersFromTrackstersEM_*_*', 
        'keep *_ticlMultiClustersFromTrackstersHAD_*_*', 
        'keep *_ticlMultiClustersFromTrackstersTrk_*_*', 
        'keep *_ticlMultiClustersFromTrackstersTrkEM_*_*', 
        'keep *_ticlMultiClustersFromTrackstersMIP_*_*', 
        'keep *_ticlMultiClustersFromTrackstersMerge_*_*'
     ) ),
    splitLevel = cms.untracked.int32(0)
)