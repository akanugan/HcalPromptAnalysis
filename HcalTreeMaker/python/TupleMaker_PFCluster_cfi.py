import FWCore.ParameterSet.Config as cms

tuplePFClusterHGCal = cms.EDProducer("TupleMaker_PFCluster",
  source    = cms.untracked.InputTag('particleFlowClusterHGCal', ''), # HGCalFromMultiCL, HO, HF, PS, HCAL, ECAL

  Prefix    = cms.untracked.string  ("PFClusterHGCal"),
  Suffix    = cms.untracked.string  ("")
)

tuplePFClusterHGCalFromMultiCL = cms.EDProducer("TupleMaker_PFCluster",
  source    = cms.untracked.InputTag('particleFlowClusterHGCalFromMultiCl', ''), # HGCalFromMultiCL, HO, HF, PS, HCAL, ECAL

  Prefix    = cms.untracked.string  ("PFClusterHGCalFromMultiCL"),
  Suffix    = cms.untracked.string  ("")
)

tuplePFClusterHO = cms.EDProducer("TupleMaker_PFCluster",
  source    = cms.untracked.InputTag('particleFlowClusterHO', ''), # HGCalFromMultiCL, HO, HF, PS, HCAL, ECAL

  Prefix    = cms.untracked.string  ("PFClusterHO"),
  Suffix    = cms.untracked.string  ("")
)

tuplePFClusterHCAL = cms.EDProducer("TupleMaker_PFCluster",
  source    = cms.untracked.InputTag('particleFlowClusterHCAL', ''), # HGCalFromMultiCL, HO, HF, PS, HCAL, ECAL

  Prefix    = cms.untracked.string  ("PFClusterHCAL"),
  Suffix    = cms.untracked.string  ("")
)

tuplePFClusterECAL = cms.EDProducer("TupleMaker_PFCluster",
  source    = cms.untracked.InputTag('particleFlowClusterECAL', ''), # HGCalFromMultiCL, HO, HF, PS, HCAL, ECAL

  Prefix    = cms.untracked.string  ("PFClusterECAL"),
  Suffix    = cms.untracked.string  ("")
)

tuplePFClusterHF = cms.EDProducer("TupleMaker_PFCluster",
  source    = cms.untracked.InputTag('particleFlowClusterHF', ''), # HGCalFromMultiCL, HO, HF, PS, HCAL, ECAL

  Prefix    = cms.untracked.string  ("PFClusterHF"),
  Suffix    = cms.untracked.string  ("")
)

#tuplePFClusterHBHE = cms.EDProducer("TupleMaker_PFCluster",
#  source    = cms.untracked.InputTag('particleFlowClusterHBHE', ''), # HGCalFromMultiCL, HO, HF, PS, HCAL, ECAL
#
#  Prefix    = cms.untracked.string  ("PFClusterHBHE"),
#  Suffix    = cms.untracked.string  ("")
#)

tuplePFClusterPS = cms.EDProducer("TupleMaker_PFCluster",
  source    = cms.untracked.InputTag('particleFlowClusterPS', ''), # HGCalFromMultiCL, HO, HF, PS, HCAL, ECAL

  Prefix    = cms.untracked.string  ("PFClusterPS"),
  Suffix    = cms.untracked.string  ("")
)
