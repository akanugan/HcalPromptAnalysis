import FWCore.ParameterSet.Config as cms

tuplePFSimParticles = cms.EDProducer("TupleMaker_PFSimParticles",
  source    = cms.untracked.InputTag('trueParticles', ''),
  Prefix    = cms.untracked.string  (""),
  Suffix    = cms.untracked.string  ("")
)

