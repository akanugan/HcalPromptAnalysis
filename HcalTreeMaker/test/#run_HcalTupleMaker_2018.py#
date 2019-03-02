#------------------------------------------------------------------------------------
# Imports
#------------------------------------------------------------------------------------
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
import FWCore.ParameterSet.VarParsing as VarParsing

#------------------------------------------------------------------------------------
# Declare the process and input variables
#------------------------------------------------------------------------------------
#process = cms.Process('NOISE',eras.Run2_50ns)#for 50ns 13 TeV data
#process = cms.Process('NOISE',eras.Run2_25ns)#for 25ns 13 TeV data
options = VarParsing.VarParsing ('analysis')
process = cms.Process("Trees",eras.Run2_2018)

options.register ('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "no of skipped events")
#
# Dataset e.g.
# dasgoclient --query 'dataset dataset=/RelValTTbar_13/CMSSW_10_2_0_pre3-*realistic*/GEN-SIM-RECO'                 
# dasgoclient --query 'file dataset=/RelValTTbar_13/CMSSW_10_2_0_pre3-101X_upgrade2018_realistic_v7-v1/GEN-SIM-RECO'
#
# TTbar sample
options.inputFiles = '/store/relval/CMSSW_10_2_0_pre3/RelValTTbar_13/GEN-SIM-RECO/101X_upgrade2018_realistic_v7-v1/20000/D8D612AA-0153-E811-9080-0CC47A4D7614.root','/store/relval/CMSSW_10_2_0_pre3/RelValTTbar_13/GEN-SIM-RECO/101X_upgrade2018_realistic_v7-v1/20000/0C125E68-0153-E811-93FB-0025905A60A8.root','/store/relval/CMSSW_10_2_0_pre3/RelValTTbar_13/GEN-SIM-RECO/101X_upgrade2018_realistic_v7-v1/20000/C2B32CAF-0253-E811-9221-0025905A607E.root'
#options.secondaryInputFiles = '/store/relval/CMSSW_10_2_0_pre3/RelValTTbar_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/D8EFF5FA-F552-E811-AC44-0CC47A7C3638.root','/store/relval/CMSSW_10_2_0_pre3/RelValTTbar_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/B28C2E00-F652-E811-B1E4-0CC47A78A3EC.root','/store/relval/CMSSW_10_2_0_pre3/RelValTTbar_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/BA2C73FF-F552-E811-9923-0CC47A4D767A.root','/store/relval/CMSSW_10_2_0_pre3/RelValTTbar_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/82A8DE13-F652-E811-87CD-0CC47A7C35A8.root','/store/relval/CMSSW_10_2_0_pre3/RelValTTbar_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/EACBA822-F652-E811-AEA8-0CC47A4C8E26.root','/store/relval/CMSSW_10_2_0_pre3/RelValTTbar_13/GEN-SIM-DIGI-RAW/101X_upgrade2018_realistic_v7-v1/20000/0C0BAB5E-F752-E811-958A-0025905A48D0.root'
options.outputFile = 'relval_ttbar_2018.root'
#
options.maxEvents = 10 # -1 means all events
#options.skipEvents = 0 # default is 0.

#------------------------------------------------------------------------------------
# Get and parse the command line arguments
#------------------------------------------------------------------------------------
options.parseArguments()
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(options.secondaryInputFiles),
    skipEvents = cms.untracked.uint32(options.skipEvents) # default is 0.
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string(options.outputFile)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    Rethrow = cms.untracked.vstring("ProductNotFound"), # make this exception fatal
    fileMode  =  cms.untracked.string('NOMERGE') # no ordering needed, but calls endRun/beginRun etc. at file boundaries
)

#------------------------------------------------------------------------------------
# import of standard configurations
#------------------------------------------------------------------------------------
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#KH
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#------------------------------------------------------------------------------------
# Set up our analyzer
#------------------------------------------------------------------------------------
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_Tree_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_Event_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_GenParticles_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HBHERecHits_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HBHEDigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HFDigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HODigis_cfi")
process.load("HcalPromptAnalysis.HcalTreeMaker.HcalTupleMaker_HcalTriggerPrimitives_cfi")

#------------------------------------------------------------------------------------
# Specify Global Tag
#------------------------------------------------------------------------------------
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')

#------------------------------------------------------------------------------------
# HcalTupleMaker sequence definition
#------------------------------------------------------------------------------------
process.tuple_step = cms.Sequence(
    # Make HCAL tuples: Event, run, ls number
    process.hcalTupleEvent*
    # Make HCAL tuples: digi info
    #process.hcalTupleHBHEDigis*
    #process.hcalTupleHODigis*
    #process.hcalTupleHFDigis*
    #process.hcalTupleTriggerPrimitives*
    # Make HCAL tuples: reco info
    process.hcalTupleHBHERecHits*
    # Make HCAL tuples: gen info
    process.hcalTupleGenParticles*
    #
    process.hcalTupleTree
)


#-----------------------------------------------------------------------------------
# Path and EndPath definitions
#-----------------------------------------------------------------------------------
process.preparation = cms.Path(
    process.tuple_step
)
