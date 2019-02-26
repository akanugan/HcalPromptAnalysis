#include "HcalPromptAnalysis/HcalTreeMaker/interface/TupleMaker_PFSimParticles.h"
//#include "DataFormats/Common/interface/Handle.h"
//#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h" 
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h" 
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h" 
#include "DataFormats/HcalRecHit/interface/HcalRecHitDefs.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"


#include <TROOT.h>
#include <TVector3.h>


TupleMaker_PFSimParticles::TupleMaker_PFSimParticles(const edm::ParameterSet& iConfig):
  inputTag    (iConfig.getUntrackedParameter<edm::InputTag>("source")),
  prefix      (iConfig.getUntrackedParameter<std::string>  ("Prefix")),
  suffix      (iConfig.getUntrackedParameter<std::string>  ("Suffix"))
  stime       (iConfig.getUntrackedParameter<bool>("PFSimParticleCollection"))
{

  produces< std::vector< double > >(prefix + "Pt"  + suffix );
  produces< std::vector< double > >(prefix + "Eta" + suffix );
  produces< std::vector< double > >(prefix + "Phi" + suffix );
  produces< std::vector< double > >(prefix + "M"   + suffix );

  produces< std::vector< int > >(prefix + "PdgId" + suffix );
  produces< std::vector< int > >(prefix + "Status" + suffix );

  // To be changed.
  //pflowToken_ = consumes<std::vector<reco::PFCandidate> >(inputTag);
  if(stime){
 tokenPFSimParticles_ = consumes<std::vector<reco::PFSimParticleCollection> >(inputTag);
  }

  debug=false;
  
}

void TupleMaker_PFSimParticles::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::unique_ptr<std::vector<double> >            pt                ( new std::vector<double>           ());
  std::unique_ptr<std::vector<double> >            eta               ( new std::vector<double>           ());
  std::unique_ptr<std::vector<double> >            phi               ( new std::vector<double>           ());
  std::unique_ptr<std::vector<double> >            mass              ( new std::vector<double>           ());

  std::unique_ptr<std::vector<int   > >            pdgid             ( new std::vector<int>              ());
  std::unique_ptr<std::vector<int   > >            status            ( new std::vector<int>              ());

  //
  //-----
  //
  if(stime ) {
    edm::Handle<std::vector<reco::<PFSimParticleCollection> > trueParticles;
    iEvent.getByToken(tokenPFSimParticles_, trueParticles);
 
    for (unsigned int i = 0; i < trueParticles->size(); i++) {
      const reco::PFSimParticleCollection& c = trueParticles->at(i);
    //reco::PFTrajectoryPoint::LayerType ecalEntrance = reco::PFTrajectoryPoint::ECALEntrance;
    //const reco::PFTrajectoryPoint& tpatecal = ((*trueParticles)[0]).extrapolatedPoint( ecalEntrance );
    // eta_ = tpatecal.positionREP().Eta();
    // phi_ = tpatecal.positionREP().Phi();
    // true_ = std::sqrt(tpatecal.momentum().Vect().Mag2());

    // pt->push_back(c.pt());
    // eta->push_back(c.eta());
    // phi->push_back(c.phi());
    // mass->push_back(c.mass());

    // pdgid->push_back(c.pdgId());
    // status->push_back(c.status());
      pt->push_back(c.pt());
      eta->push_back(c.eta());
      phi->push_back(c.phi());
      mass->push_back(c.mass());

      pdgid->push_back(c.pdgId());
      status->push_back(c.status());

  }
  
  //
  //-----
  //
  // /*
  // if(bool_PackedCandidate){
  //   edm::Handle<std::vector<pat::PackedCandidate> > packedParticleFlow;
  //   iEvent.getByToken(pflowPackedToken_, packedParticleFlow);

  //   for (unsigned int i = 0; i < packedParticleFlow->size(); i++) {
  //     const pat::PackedCandidate& c = packedParticleFlow->at(i);
      
  //     if (debug){
  // 	//if (c.pdgId()==130){ // K0L neutral hadron
  // 	if (fabs(c.pdgId())==211){ // charged hadron
  // 	  std::cout << "pt,eta: " << c.pt() << " " << c.eta() << std::endl;
  // 	  //std::cout << c.ecalEnergy() << std::endl;
  // 	  std::cout << c.rawCaloFraction() << std::endl;
  // 	  std::cout << c.hcalFraction() << std::endl;
  // 	  /*
  // 	    std::cout << c.hoEnergy() << std::endl;
  // 	    std::cout << c.hcalDepthEnergyFraction(1) << " "
  // 	    << c.hcalDepthEnergyFraction(2) << " "
  // 	    << c.hcalDepthEnergyFraction(3) << std::endl;
  // 	  */
  // 	}
  //     }

  //     pt->push_back(c.pt());
  //     eta->push_back(c.eta());
  //     phi->push_back(c.phi());
  //     mass->push_back(c.mass());

  //     pdgid->push_back(c.pdgId());
  //     status->push_back(c.status());

  //   }
      
  // //
  // //-----
  // //
  // } else {
  //   edm::Handle<std::vector<reco::PFCandidate> > particleFlow;
  //   iEvent.getByToken(pflowToken_, particleFlow);

  //   for (unsigned int i = 0; i < particleFlow->size(); i++) {
  //     const reco::PFCandidate& c = particleFlow->at(i);

  //     double track_Pt=0.;
  //     const reco::Track * tr=c.bestTrack();
  //     if(tr!=nullptr) track_Pt=tr->pt();

  //     if (debug){
  // 	std::cout << "pt,eta: " << c.pt() << " " << c.eta() << std::endl;
  // 	std::cout << c.ecalEnergy() << std::endl;
  // 	std::cout << c.hcalEnergy() << std::endl;
  // 	std::cout << c.hoEnergy() << std::endl;
  // 	std::cout << c.hcalDepthEnergyFraction(1) << " "
  // 		  << c.hcalDepthEnergyFraction(2) << " "
  // 		  << c.hcalDepthEnergyFraction(3) << std::endl;
  // 	std::cout << track_Pt << std::endl;
  //     }
      
  //     pt->push_back(c.pt());
  //     eta->push_back(c.eta());
  //     phi->push_back(c.phi());
  //     mass->push_back(c.mass());

  //     pdgid->push_back(c.particleId());
  //     status->push_back(c.status());

  //     if (c.energy()>0.){
  // 	ecalEnergyFrac->push_back(c.ecalEnergy()/c.energy());
  // 	hcalEnergyFrac->push_back(c.hcalEnergy()/c.energy());
  // 	hoEnergyFrac->push_back(c.hoEnergy()/c.energy());
  //     } else {
  // 	ecalEnergyFrac->push_back(0.);
  // 	hcalEnergyFrac->push_back(0.);
  // 	hoEnergyFrac->push_back(0.);
  //     }
  //     trackPt->push_back(track_Pt);

  //     hcalFrac1->push_back(c.hcalDepthEnergyFraction(1));
  //     hcalFrac2->push_back(c.hcalDepthEnergyFraction(2));
  //     hcalFrac3->push_back(c.hcalDepthEnergyFraction(3));
  //     hcalFrac4->push_back(c.hcalDepthEnergyFraction(4));
  //     hcalFrac5->push_back(c.hcalDepthEnergyFraction(5));
  //     hcalFrac6->push_back(c.hcalDepthEnergyFraction(6));
  //     hcalFrac7->push_back(c.hcalDepthEnergyFraction(7));

      
  //   }
  //   */

  //
  //-----
  //
  }

  iEvent.put(move( pt              ) , prefix + "Pt"            + suffix );
  iEvent.put(move( eta             ) , prefix + "Eta"           + suffix );
  iEvent.put(move( phi             ) , prefix + "Phi"           + suffix );
  iEvent.put(move( mass            ) , prefix + "M"             + suffix );

  iEvent.put(move( pdgid           ) , prefix + "PdgId"         + suffix );
  iEvent.put(move( status          ) , prefix + "Status"        + suffix );

// iEvent.put(move( ecalEnergyFrac  ) , prefix + "EcalEnergyFrac"  + suffix );
//iEvent.put(move( hcalEnergyFrac  ) , prefix + "HcalEnergyFrac"  + suffix );
// iEvent.put(move( hoEnergyFrac    ) , prefix + "HOEnergyFrac"    + suffix );
// iEvent.put(move( trackPt         ) , prefix + "TrackPt"         + suffix );
  		            		            
// iEvent.put(move( hcalFrac1       ) , prefix + "HcalFrac1"  + suffix );
// iEvent.put(move( hcalFrac2       ) , prefix + "HcalFrac2"  + suffix );
// iEvent.put(move( hcalFrac3       ) , prefix + "HcalFrac3"  + suffix );
// iEvent.put(move( hcalFrac4       ) , prefix + "HcalFrac4"  + suffix );
// iEvent.put(move( hcalFrac5       ) , prefix + "HcalFrac5"  + suffix );
// iEvent.put(move( hcalFrac6       ) , prefix + "HcalFrac6"  + suffix );
// iEvent.put(move( hcalFrac7       ) , prefix + "HcalFrac7"  + suffix );
  
}


