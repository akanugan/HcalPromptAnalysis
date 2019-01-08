import FWCore.ParameterSet.Config as cms

tuplePFJets = cms.EDProducer("TupleMaker_PFJets",
  source    = cms.untracked.InputTag('ak4PFJets', ''),
  #PackedCandidate = cms.untracked.bool(False),
  Prefix    = cms.untracked.string  ("PFJets"),
  Suffix    = cms.untracked.string  ("")
)
