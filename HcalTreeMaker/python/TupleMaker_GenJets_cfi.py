import FWCore.ParameterSet.Config as cms

tupleGenJets = cms.EDProducer("TupleMaker_GenJets",
  source    = cms.untracked.InputTag('ak4GenJets', ''),
  #PackedCandidate = cms.untracked.bool(False),
  Prefix    = cms.untracked.string  ("GenJets"),
  Suffix    = cms.untracked.string  ("")
)
