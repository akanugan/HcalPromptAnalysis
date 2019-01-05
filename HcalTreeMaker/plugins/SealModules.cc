#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_Tree.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_Event.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_HcalDigis.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_QIE10Digis.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_QIE11Digis.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_HcalRecHits.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_HcalSimHits.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_GenParticles.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_HcalTriggerPrimitives.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/TupleMaker_PFCandidates.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/TupleMaker_GenJets.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/TupleMaker_GenMet.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/TupleMaker_PFJets.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/TupleMaker_PFMet.h"
#include "HcalPromptAnalysis/HcalTreeMaker/interface/TupleMaker_PFCluster.h"


DEFINE_FWK_MODULE(HcalTupleMaker_Tree);
DEFINE_FWK_MODULE(HcalTupleMaker_Event);
DEFINE_FWK_MODULE(HcalTupleMaker_GenParticles);
DEFINE_FWK_MODULE(HcalTupleMaker_HBHERecHits);
DEFINE_FWK_MODULE(HcalTupleMaker_HORecHits);
DEFINE_FWK_MODULE(HcalTupleMaker_HFRecHits);
DEFINE_FWK_MODULE(HcalTupleMaker_HcalSimHits);
DEFINE_FWK_MODULE(HcalTupleMaker_HBHEDigis);
DEFINE_FWK_MODULE(HcalTupleMaker_HODigis);
DEFINE_FWK_MODULE(HcalTupleMaker_HFDigis);
DEFINE_FWK_MODULE(HcalTupleMaker_QIE10Digis);
DEFINE_FWK_MODULE(HcalTupleMaker_QIE11Digis);
DEFINE_FWK_MODULE(HcalTupleMaker_HcalTriggerPrimitives);
DEFINE_FWK_MODULE(TupleMaker_PFCandidates);
DEFINE_FWK_MODULE(TupleMaker_GenJets);
DEFINE_FWK_MODULE(TupleMaker_GenMet);
DEFINE_FWK_MODULE(TupleMaker_PFJets);
DEFINE_FWK_MODULE(TupleMaker_PFMet);
DEFINE_FWK_MODULE(TupleMaker_PFCluster);
