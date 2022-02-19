from WMCore.Configuration import Configuration 

config = Configuration()
config.section_("General")
config.General.requestName = "DoubleMuon_Run2016F"
config.General.transferLogs = False 
config.General.workArea = "crab2016pre"

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "PSet.py"
config.JobType.scriptExe = "./WWG_crab_script.sh" 
config.JobType.inputFiles = ["../../../scripts/haddnano.py","../WWG_selector/WWG_postproc.py","../WWG_selector/WWG_Module.py","../WWG_selector/WWG_keep_and_drop.txt","../WWG_selector/WWG_output_branch.txt","../WWG_selector/DAS_filesearch.py"] #hadd nano will not be needed once nano tools are in cmssw 
config.JobType.scriptArgs = ["isdata=data","year=2016pre","era=F"] 
config.JobType.sendPythonFolder  = True
config.JobType.allowUndistributedCMSSW = True 

config.section_("Data")
config.Data.inputDataset = "/DoubleMuon/Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD" 
config.Data.inputDBS = "global"
# config.Data.splitting = "FileBased"
# config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.ignoreLocality = True
config.Data.allowNonValidInputDataset = True
config.Data.outputDatasetTag = "DoubleMuon_Run2016F" 

config.section_("Site")
config.Site.storageSite = "T3_CH_CERNBOX"
config.Site.whitelist = ["T2_US_MIT","T2_US_Wisconsin","T2_US_Purdue","T2_US_UCSD","T2_US_Caltech","T2_US_Nebraska"] 