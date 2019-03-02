import FWCore.ParameterSet.Config as cms

hcalTupleSimTracks = cms.EDProducer("HcalTupleMaker_SimTracks",
  Source = cms.untracked.InputTag("g4SimHits",""),
  Source_SimVtx = cms.untracked.InputTag("g4SimHits",""),
  Prefix = cms.untracked.string  ("SimTracks"),
  Suffix = cms.untracked.string  ("")
)
