from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = ''
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles=['step3.root','step3_inMINIAODSIM.root']

config.JobType.maxMemoryMB = 3000

config.section_("Data")
# MC example
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/'
config.Data.inputDataset = '/SingleK0L/hatake-CMSSW_10_3_0_pre6_Step2_v1-846a642314dcc4ef3ea99a82b4cfab59/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
#KH config.Data.totalUnits = 100
# MC example ends
# Data example
#config.Data.inputDataset = '/ZeroBias/Run2018A-v1/RAW'
#config.Data.inputDataset = '/HLTPhysics/Run2018A-v1/RAW'
#config.Data.runRange = '315361-315690'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 10
# Data example ends

config.Data.publication = True
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/' # Parameter Data.publishDbsUrl has been renamed to Data.publishDBS
config.Data.outputDatasetTag = 'CMSSW_10_3_0_pre6_Step3_v1' # <== Check!!!

config.Data.outLFNDirBase = '/store/user/hatake/crab_outputs'  # Data.outLFN has been renamed to Data.outLFNDirBase
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.storageSite = 'T3_US_Baylor'
#KH (this whitelisting below is not really necessary. we can use any T2/T3 for running jobs. we can still send output to Baylor)
config.Site.whitelist = ['T3_US_Baylor','T2_US_*']
