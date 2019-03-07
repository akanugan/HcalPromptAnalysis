import FWCore.ParameterSet.Config as cms

tuplePFSimParticles = cms.EDProducer("TupleMaker_PFSimParticles",
  source    = cms.untracked.InputTag('particleFlowSimParticle'),
  Prefix    = cms.untracked.string  ("PFSimPar"),
  Suffix    = cms.untracked.string  ("")
)

