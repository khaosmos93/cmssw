import FWCore.ParameterSet.Config as cms

from ..modules.hltPrePhoton100EBTightIDTightIsoOpenL1Seeded_cfi import *
from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTEndSequence_cfi import *
from ..sequences.HLTPhoton100EBTightIDTightIsoOpenL1SeededSequence_cfi import *

HLT_Photon100EB_TightID_TightIso_Open_L1Seeded_v1 = cms.Path(HLTBeginSequence+hltPrePhoton100EBTightIDTightIsoOpenL1Seeded+HLTPhoton100EBTightIDTightIsoOpenL1SeededSequence+HLTEndSequence)
