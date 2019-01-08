import FWCore.ParameterSet.Config as cms

tupleGenMet = cms.EDProducer("TupleMaker_GenMet",
  source    = cms.untracked.InputTag('genMetTrue', ''),
  #PackedCandidate = cms.untracked.bool(False),
  Prefix    = cms.untracked.string  ("Gen"),
  Suffix    = cms.untracked.string  ("")
)

#tuplePackedPFCandidates = tuplePFCandidates.clone(
#  source    = cms.untracked.InputTag('packedPFCandidates', ''),
#  PackedCandidate = cms.untracked.bool(True),
#)
