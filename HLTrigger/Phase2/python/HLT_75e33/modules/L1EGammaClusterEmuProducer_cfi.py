import FWCore.ParameterSet.Config as cms

L1EGammaClusterEmuProducer = cms.EDProducer("L1EGCrystalClusterEmulatorProducer",
    calib = cms.PSet(
        etaBins = cms.vdouble(
            0.087, 0.174, 0.261, 0.348, 0.435, 
            0.522, 0.609, 0.696, 0.783, 0.87, 
            0.957, 1.044, 1.131, 1.218, 1.305, 
            1.392, 1.479
        ),
        ptBins = cms.vdouble(
            12, 20, 30, 40, 55, 
            90, 1000000.0
        ),
        scale = cms.vdouble(
            1.298, 1.287, 1.309, 1.298, 1.309, 
            1.309, 1.309, 1.298, 1.309, 1.298, 
            1.309, 1.309, 1.309, 1.32, 1.309, 
            1.32, 1.309, 1.1742, 1.1639, 1.1639, 
            1.1639, 1.1639, 1.1639, 1.1639, 1.1742, 
            1.1742, 1.1639, 1.1639, 1.1742, 1.1639, 
            1.1639, 1.1742, 1.1742, 1.1536, 1.11, 
            1.11, 1.11, 1.11, 1.11, 1.11, 
            1.11, 1.11, 1.11, 1.11, 1.11, 
            1.11, 1.11, 1.11, 1.11, 1.11, 
            1.1, 1.09, 1.09, 1.09, 1.09, 
            1.09, 1.09, 1.09, 1.09, 1.09, 
            1.09, 1.09, 1.09, 1.09, 1.09, 
            1.09, 1.09, 1.09, 1.07, 1.07, 
            1.07, 1.07, 1.07, 1.07, 1.07, 
            1.08, 1.07, 1.07, 1.08, 1.08, 
            1.07, 1.08, 1.08, 1.08, 1.08, 
            1.06, 1.06, 1.06, 1.06, 1.05, 
            1.05, 1.06, 1.06, 1.06, 1.06, 
            1.06, 1.06, 1.06, 1.06, 1.06, 
            1.06, 1.06, 1.04, 1.04, 1.04, 
            1.04, 1.05, 1.04, 1.05, 1.05, 
            1.05, 1.05, 1.05, 1.05, 1.05, 
            1.05, 1.05, 1.05, 1.05
        )
    ),
    ecalTPEB = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","HLT"),
    hcalTP = cms.InputTag("simHcalTriggerPrimitiveDigis","","HLT")
)

