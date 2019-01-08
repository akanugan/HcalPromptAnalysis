#ifndef TupleMaker_PFMet_h
#define TupleMaker_PFMet_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
//#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

class TupleMaker_PFMet : public edm::EDProducer {
 public:
  explicit TupleMaker_PFMet(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag;
  const std::string     prefix,suffix;

  edm::EDGetTokenT<reco::PFMETCollection>  pfMETsToken_;
  //edm::EDGetTokenT<std::vector<pat::PackedCandidate> > pflowPackedToken_;

  bool debug;
  //bool bool_PackedCandidate;
  
};

#endif
