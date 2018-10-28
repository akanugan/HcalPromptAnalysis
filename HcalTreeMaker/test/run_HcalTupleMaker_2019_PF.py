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
process = cms.Process("Trees",eras.Run3)

##
## Setup command line options
##
options.register ('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "no of skipped events")
options.register ('isMINIAOD', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "MINIAODSIM input file(s)?")

##
## Default
##
options.maxEvents = 10 # -1 means all events
#options.skipEvents = 0 # default is 0.

##
## get and parse the command line arguments
##
options.parseArguments()
print("isMINIAOD: ", options.isMINIAOD)
print("maxEvents: ", options.maxEvents)

#
# Dataset e.g.
# dasgoclient --query 'dataset dataset=/RelValTTbar_13/CMSSW_10_2_0_pre3-*realistic*/GEN-SIM-RECO'                 
# dasgoclient --query 'file dataset=/RelValTTbar_13/CMSSW_10_2_0_pre3-101X_upgrade2018_realistic_v7-v1/GEN-SIM-RECO'
#
# TTbar sample
#
# MINIAODSIM
if options.isMINIAOD: 
    options.inputFiles = '/store/relval/CMSSW_10_3_0_pre4/RelValTTbar_13/MINIAODSIM/PUpmx25ns_103X_upgrade2018_realistic_v4-v1/20000/17E57223-3406-7340-B867-2DDC36E7C371.root'
    options.outputFile = 'trees_relval_ttbar_2019_pmx25ns_miniaodsim.root'
# GEN-SIM-RECO
else:
    options.inputFiles = '/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_1.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_10.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_11.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_12.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_13.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_14.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_15.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_16.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_17.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_18.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_19.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_2.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_20.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_21.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_22.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_23.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_24.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_25.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_26.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_27.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_28.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_29.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_3.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_30.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_31.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_32.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_33.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_34.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_35.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_36.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_37.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_38.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_39.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_4.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_40.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_41.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_42.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_43.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_44.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_45.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_46.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_47.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_48.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_49.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_5.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_50.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_6.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_7.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_8.root','/store/user/hatake/crab_outputs/SingleK0L/CMSSW_10_3_0_pre6_Step3_v1/181023_023516/0000/step3_9.root'
    options.outputFile = 'trees_relval_singleK0L_2019.root'
#

#
#
#
print("maxEvents: ", options.maxEvents)
print("inputFiles: ", options.inputFiles)
print("outputFile: ", options.outputFile)

#------------------------------------------------------------------------------------
# Get and parse the command line arguments
#------------------------------------------------------------------------------------
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
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('Configuration.StandardSequences.PATMC_cff')
#process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
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

process.load("HcalPromptAnalysis.HcalTreeMaker.TupleMaker_PFCandidates_cfi")
#process.tuplePackedPFCandidates.debug = cms.untracked.bool(True)

#------------------------------------------------------------------------------------
# Specify Global Tag
#------------------------------------------------------------------------------------
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2019_realistic', '')

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
    #process.hcalTupleHBHERecHits*
    # Make HCAL tuples: gen info
    process.hcalTupleGenParticles*
    #
    process.tuplePFCandidates*
    #
    process.hcalTupleTree
)

#
# in case we are using MINIAOD files
#
if options.isMINIAOD: 
    process.tuple_step = cms.Sequence(
        # Make HCAL tuples: Event, run, ls number
        process.hcalTupleEvent*
        # Make HCAL tuples: gen info
        #process.hcalTupleGenParticles*
        #
        process.tuplePackedPFCandidates*
        #
        process.hcalTupleTree
    )

#-----------------------------------------------------------------------------------
# Path and EndPath definitions
#-----------------------------------------------------------------------------------
process.preparation = cms.Path(
    process.tuple_step
)

