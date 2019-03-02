import FWCore.ParameterSet.Config as cms

tuplePFSimParticles = cms.EDProducer("TupleMaker_PFSimParticles",
  source    = cms.untracked.InputTag('trueParticles', ''),
  PackedCandidate = cms.untracked.bool(False),
  Prefix    = cms.untracked.string  ("PFSim"),
  Suffix    = cms.untracked.string  ("")
)


