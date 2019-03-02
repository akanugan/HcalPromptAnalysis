#ifndef HGCalTupleMaker_SimTracks_h
#define HGCalTupleMaker_SimTracks_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/transform.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/transform.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

//
// Class definition
// 
class HGCalTupleMaker_SimTracks : public edm::EDProducer {
 protected:

  edm::EDGetTokenT<edm::SimTrackContainer> SimTrackContainer_;
  edm::EDGetTokenT<std::vector<SimVertex>> SimVertexContainer_;

  const edm::InputTag   m_inputTag;
  const edm::InputTag   m_inputTag_SimVtx;
  const std::string     m_prefix;
  const std::string     m_suffix;

  bool debug=false;

  void produce( edm::Event & iEvent, const edm::EventSetup & iSetup ) { 

    if (debug) std::cout << "HGCalTupleMaker_SimTracks: produce starts" << std::endl;
    
    //-----------------------------------------------------
    // Prepare to put things into event
    //-----------------------------------------------------
    
    loadAlgo();

    //-----------------------------------------------------
    // edm::Handles
    //-----------------------------------------------------

    edm::Handle<edm::SimTrackContainer> simTracks;
    edm::Handle<std::vector<SimVertex>> simVertices;

    iEvent.getByToken(SimTrackContainer_, simTracks);
    iEvent.getByToken(SimVertexContainer_,simVertices);

    //
    // Loop over SimTracks
    //
    const edm::SimTrackContainer& trks = *(simTracks.product());
    edm::SimTrackContainer::const_iterator trksiter;

    if (debug) std::cout << "SimTrack Size: " << trks.size() << std::endl;

    for (trksiter = trks.begin(); trksiter != trks.end(); trksiter++){ 

      double simPt=trksiter->momentum().pt();
      double simEta=trksiter->momentum().eta();
      double simPhi=trksiter->momentum().phi();
      int    simCharge=trksiter->charge();
      int    simPID=trksiter->type();
      int    vertexId=trksiter->vertIndex();

      double simR=-1.;
      double simZ= 0.;
      if ( vertexId >= 0 ) {
	simR=(*simVertices)[vertexId].position().pt();
	simZ=(*simVertices)[vertexId].position().z();
      }

      if (debug) std::cout << "SimTrack pt, eta, phi: " << simPt << " " << simEta << " " << simPhi << std::endl;
      if (debug) std::cout << simCharge << " " << simPID << std::endl;
      if (debug) std::cout << simR      << " " << simZ   << std::endl;

      v_pt     -> push_back ( simPt );
      v_eta    -> push_back ( simEta );
      v_phi    -> push_back ( simPhi );
      
      v_charge -> push_back ( simCharge );
      v_pid    -> push_back ( simPID );
      v_r      -> push_back ( simR   );
      v_z      -> push_back ( simZ   );

    }

    //-----------------------------------------------------
    // Put things into the event
    //-----------------------------------------------------
    
    dumpAlgo(iEvent);

  }
  
 public:
  
 HGCalTupleMaker_SimTracks(const edm::ParameterSet& iConfig) :
    m_inputTag        (iConfig.getUntrackedParameter<edm::InputTag>("Source")),
    m_inputTag_SimVtx (iConfig.getUntrackedParameter<edm::InputTag>("Source_SimVtx")),
    m_prefix          (iConfig.getUntrackedParameter<std::string>  ("Prefix")),
    m_suffix          (iConfig.getUntrackedParameter<std::string>  ("Suffix")) {

    SimTrackContainer_ = consumes<edm::SimTrackContainer>(m_inputTag);
    SimVertexContainer_ = consumes<std::vector<SimVertex>>(m_inputTag_SimVtx);

    produces<std::vector<float> > ( m_prefix + "Pt"     + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Eta"    + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Phi"    + m_suffix );

    produces<std::vector<int> >   ( m_prefix + "Charge" + m_suffix );
    produces<std::vector<int> >   ( m_prefix + "PID"    + m_suffix );
    produces<std::vector<float> > ( m_prefix + "R"      + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Z"      + m_suffix );

  }

  std::unique_ptr<std::vector<float> > v_pt;
  std::unique_ptr<std::vector<float> > v_eta;
  std::unique_ptr<std::vector<float> > v_phi;

  std::unique_ptr<std::vector<int>   > v_charge;
  std::unique_ptr<std::vector<int>   > v_pid;
  std::unique_ptr<std::vector<float> > v_r;
  std::unique_ptr<std::vector<float> > v_z;

 private:

  void beginRun(const edm::Run&, const edm::EventSetup& iSetup){
  }
  
 protected:

  void loadAlgo(){
    v_pt     = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_eta    = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_phi    = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());

    v_charge = std::unique_ptr<std::vector<int>   > ( new std::vector<int> ());
    v_pid    = std::unique_ptr<std::vector<int>   > ( new std::vector<int> ());
    v_r      = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_z      = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());

  }
  
  void dumpAlgo( edm::Event & iEvent ){

    iEvent.put( move(v_pt     ), m_prefix + "Pt"     + m_suffix );
    iEvent.put( move(v_eta    ), m_prefix + "Eta"    + m_suffix );
    iEvent.put( move(v_phi    ), m_prefix + "Phi"    + m_suffix );

    iEvent.put( move(v_charge ), m_prefix + "Charge" + m_suffix );
    iEvent.put( move(v_pid    ), m_prefix + "PID"    + m_suffix );
    iEvent.put( move(v_r      ), m_prefix + "R"      + m_suffix );
    iEvent.put( move(v_z      ), m_prefix + "Z"      + m_suffix );

  }  
    
};

#endif
