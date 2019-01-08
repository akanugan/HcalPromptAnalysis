#include "HcalPromptAnalysis/HcalTreeMaker/interface/TupleMaker_PFCluster.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TLorentzVector.h"

TupleMaker_PFCluster::TupleMaker_PFCluster(const edm::ParameterSet& iConfig):
  inputTag    (iConfig.getUntrackedParameter<edm::InputTag>("source")),
  prefix      (iConfig.getUntrackedParameter<std::string>  ("Prefix")),
  suffix      (iConfig.getUntrackedParameter<std::string>  ("Suffix"))
{
   
  produces< std::vector< double > >(prefix + "PosX" + suffix );
  produces< std::vector< double > >(prefix + "PosY" + suffix );
  produces< std::vector< double > >(prefix + "PosZ" + suffix );
  //produces< std::vector< double > >(prefix + "Layer" + suffix );
  produces< std::vector< double > >(prefix + "Energy" + suffix );
  produces< std::vector< double > >(prefix + "Eta" + suffix );
  produces< std::vector< double > >(prefix + "Phi" + suffix );
  produces< std::vector< double > >(prefix + "Size" + suffix );
  produces< std::vector< float  > >(prefix + "Time" + suffix );
  produces< std::vector< float  > >(prefix + "TimeE" + suffix );  
  produces< std::vector< double > >(prefix + "Depth" + suffix );
  
  
  PFClusterTok_ = consumes<std::vector<reco::PFCluster>>(inputTag);

}

void TupleMaker_PFCluster::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::unique_ptr<std::vector<double> > posX   ( new std::vector<double>       ());
  std::unique_ptr<std::vector<double> > posY   ( new std::vector<double>       ());
  std::unique_ptr<std::vector<double> > posZ   ( new std::vector<double>       ());
  std::unique_ptr<std::vector<double> > energy ( new std::vector<double>       ());
  std::unique_ptr<std::vector<double> > eta    ( new std::vector<double>       ());
  std::unique_ptr<std::vector<double> > phi    ( new std::vector<double>       ());
  std::unique_ptr<std::vector<double> > size   ( new std::vector<double>       ());
  std::unique_ptr<std::vector<float> > time   ( new std::vector<float>       ());
  std::unique_ptr<std::vector<float> > timeE   ( new std::vector<float>       ());
  std::unique_ptr<std::vector<double> > depth   ( new std::vector<double>       ());
  
  edm::Handle< std::vector<reco::PFCluster> > pfClusters;
  iEvent.getByToken(PFClusterTok_, pfClusters);

  for(unsigned int iPart = 0 ; iPart < pfClusters->size(); ++iPart){
    //const reco::PFCluster& cluster = (*pfClusters)[iPart];
    
    posX->push_back(pfClusters->at(iPart).x());
    posY->push_back(pfClusters->at(iPart).y());
    posZ->push_back(pfClusters->at(iPart).z());
    energy->push_back(pfClusters->at(iPart).energy());
    eta->push_back(pfClusters->at(iPart).eta());
    phi->push_back(pfClusters->at(iPart).phi());
    size->push_back(pfClusters->at(iPart).size());
    time->push_back(pfClusters->at(iPart).time());
    timeE->push_back(pfClusters->at(iPart).timeError());
    depth->push_back(pfClusters->at(iPart).depth());

  }

  iEvent.put(move( posX             ) , prefix + "PosX"     + suffix );
  iEvent.put(move( posY             ) , prefix + "PosY"     + suffix );
  iEvent.put(move( posZ             ) , prefix + "PosZ"     + suffix );
  iEvent.put(move( energy           ) , prefix + "Energy"   + suffix );
  iEvent.put(move( eta              ) , prefix + "Eta"      + suffix );
  iEvent.put(move( phi              ) , prefix + "Phi"      + suffix );
  iEvent.put(move( size             ) , prefix + "Size"     + suffix );
  iEvent.put(move( time             ) , prefix + "Time"     + suffix );
  iEvent.put(move( timeE             ) , prefix + "TimeE"     + suffix );
  iEvent.put(move( depth             ) , prefix + "Depth"     + suffix );
		      
}
		      
    
