import FWCore.ParameterSet.Config as cms

tuplePFMet = cms.EDProducer("TupleMaker_PFMet",
  source    = cms.untracked.InputTag('pfMet', ''),
  #PackedCandidate = cms.untracked.bool(False),
  Prefix    = cms.untracked.string  ("PF"),
  Suffix    = cms.untracked.string  ("")
)

#tuplePackedPFCandidates = tuplePFCandidates.clone(
#  source    = cms.untracked.InputTag('packedPFCandidates', ''),
#  PackedCandidate = cms.untracked.bool(True),
#)
