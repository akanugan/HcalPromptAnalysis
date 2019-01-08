#include "HcalPromptAnalysis/HcalTreeMaker/interface/TupleMaker_GenMet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/ESHandle.h"

TupleMaker_GenMet::TupleMaker_GenMet(const edm::ParameterSet& iConfig):
  inputTag    (iConfig.getUntrackedParameter<edm::InputTag>("source")),
  prefix      (iConfig.getUntrackedParameter<std::string>  ("Prefix")),
  suffix      (iConfig.getUntrackedParameter<std::string>  ("Suffix"))
  //bool_SlimmedMET (iConfig.getUntrackedParameter<bool>("SlimmedMET"))
{

  produces < double > (prefix + "MET" + suffix);
  produces < double > (prefix + "METPhi" + suffix);

  //if(bool_PackedCandidate){
  //pflowPackedToken_ = consumes<std::vector<pat::PackedCandidate> >(inputTag);
  //}else{
  genMETsToken_ = consumes<reco::GenMETCollection> (inputTag);
  //}

  debug=false;}


void TupleMaker_GenMet::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::unique_ptr<double>              MET                  ( new double        ());
  std::unique_ptr<double>              METPhi               ( new double        ());

  //
  //-----
  //
  
  //if(bool_PackedCandidate){
  //  edm::Handle<std::vector<pat::PackedCandidate> > packedParticleFlow;
  //  iEvent.getByToken(pflowPackedToken_, packedParticleFlow);
  //
  //  for (unsigned int i = 0; i < packedParticleFlow->size(); i++) {
  //    const pat::PackedCandidate& c = packedParticleFlow->at(i);
  //

  //
 
  //    MET->push_back(c.MET());
  //  }
  //    
  //
  //-----
  //
  //} else {
  edm::Handle<reco::GenMETCollection> genMETs;
  iEvent.getByToken(genMETsToken_, genMETs);
  if (!genMETs.isValid()) return;

  reco::MET met;
  met = genMETs->front();

  const double MetPt = met.pt();
  const double MetPhi = met.phi();


  *MET    =  MetPt;
  *METPhi =  MetPhi;
 
  //}
    
  //
  //-----
  //



  iEvent.put(move( MET                ) , prefix + "MET"        + suffix );
  iEvent.put(move( METPhi             ) , prefix + "METPhi"        + suffix );

}


