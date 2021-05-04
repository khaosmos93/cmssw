import FWCore.ParameterSet.Config as cms

ticlTrackstersEML1 = cms.EDProducer("TrackstersProducer",
    algo_verbosity = cms.int32(0),
    detector = cms.string('HGCAL'),
    eid_graph_path = cms.string('RecoHGCal/TICL/data/tf_models/energy_id_v0.pb'),
    eid_input_name = cms.string('input'),
    eid_min_cluster_energy = cms.double(1),
    eid_n_clusters = cms.int32(10),
    eid_n_layers = cms.int32(50),
    eid_output_name_energy = cms.string('output/regressed_energy'),
    eid_output_name_id = cms.string('output/id_probabilities'),
    energy_em_over_total_threshold = cms.double(0.9),
    etaLimitIncreaseWindow = cms.double(2.1),
    filter_on_categories = cms.vint32(0, 1),
    filtered_mask = cms.InputTag("filteredLayerClustersEML1Seeded","EM"),
    itername = cms.string('EM'),
    layer_clusters = cms.InputTag("hgcalLayerClustersL1Seeded"),
    layer_clusters_hfnose_tiles = cms.InputTag("ticlLayerTileHFNose"),
    layer_clusters_tiles = cms.InputTag("ticlLayerTileProducerL1Seeded"),
    max_delta_time = cms.double(3),
    max_longitudinal_sigmaPCA = cms.double(10),
    max_missing_layers_in_trackster = cms.int32(1),
    max_out_in_hops = cms.int32(1),
    mightGet = cms.optional.untracked.vstring,
    min_cos_pointing = cms.double(0.9),
    min_cos_theta = cms.double(0.97),
    min_layers_per_trackster = cms.int32(10),
    oneTracksterPerTrackSeed = cms.bool(False),
    original_mask = cms.InputTag("hgcalLayerClustersL1Seeded","InitialLayerClustersMask"),
    out_in_dfs = cms.bool(True),
    pid_threshold = cms.double(0.5),
    promoteEmptyRegionToTrackster = cms.bool(False),
    root_doublet_max_distance_from_seed_squared = cms.double(9999),
    seeding_regions = cms.InputTag("ticlSeedingL1"),
    shower_start_max_layer = cms.int32(5),
    skip_layers = cms.int32(2),
    time_layerclusters = cms.InputTag("hgcalLayerClustersL1Seeded","timeLayerCluster")
)
