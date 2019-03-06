import FWCore.ParameterSet.Config as cms

tuplePFSimParticles = cms.EDProducer("TupleMaker_PFSimParticles",
  source    = cms.untracked.InputTag('particleFlowSimParticle'),
  Prefix    = cms.untracked.string  ("PFSim"),
  Suffix    = cms.untracked.string  ("")
)

###-- Dump config ------------------------------------------------------------
file = open('allDump_cfg.py','w')
file.write(str(process.dumpPython()))
file.close()
