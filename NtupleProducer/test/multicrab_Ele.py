name = 'Run3Winter21_112X_v1'
# name = 'RunIIAutumn18_102X_v2'

dataset = {
   'DYToLL_M50_Run3': '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM',
   'DYEE_FlatPU20to70' : '/ZToEE_TuneCUETP8M1_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU20to70_for_DNN_112X_mcRun3_2021_realistic_v16_ext1-v1/MINIAODSIM',
   'DYEE_NoPURAW' : '/ZToEE_TuneCUETP8M1_14TeV-pythia8/Run3Winter21DRMiniAOD-NoPURAW_for_DNN_112X_mcRun3_2021_realistic_v16_ext1-v1/MINIAODSIM',
   'TT_FlatPU20to70':  '/TT_TuneCP5_14TeV-powheg-pythia8/Run3Winter21DRMiniAOD-FlatPU20to70_for_DNN_112X_mcRun3_2021_realistic_v16_ext1-v1/MINIAODSIM',
   'TT_FlatPU30to80':  '/TT_TuneCP5_14TeV-powheg-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM',
   # 'ZprimeEE' : '/ZprimeToEE_M-6000_TuneCP5_14TeV-pythia8/Run3Winter20DRMiniAOD-FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/MINIAODSIM',
   # 'TTSemilep' : '/TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8/Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v2/MINIAODSIM',
   'QCD_Pt15to7000' : '/QCD_Pt-15to7000_TuneCUETP8M1_Flat_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU20to70_for_DNN_112X_mcRun3_2021_realistic_v16_ext1-v1/MINIAODSIM',
   'DYJetsToLL_2018': '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
   }


#nevents = -1
#lumisPerJob = {
#   'Run2017B':        100,
#   'Run2017C':        100,
#   'Run2017D':        100,
#   'Run2017E':        100,
#   'Run2017F':        100,
#   }

listOfSamples = [
   # Run3Winter
    # 'DYEE_FlatPU20to70',
    # 'DYEE_NoPURAW',
    # # 'ZprimeEE',
    # # 'TTSemilep',
    # 'TT_FlatPU20to70',
    # 'TT_FlatPU30to80',
    # 'QCD_Pt15to7000',
    'DYToLL_M50_Run3',

    # RunIIAutumn
    # 'DYJetsToLL_2018'
   ]


if __name__ == '__main__':

   from CRABClient.UserUtilities import config
   config = config()

   from CRABAPI.RawCommand import crabCommand
   from multiprocessing import Process

   def submit(config):
       res = crabCommand('submit', config = config)

   config.General.workArea = 'crab_'+name
   config.General.transferLogs = False
   config.General.transferOutputs = True
   config.JobType.allowUndistributedCMSSW = True
   config.JobType.pluginName = 'Analysis'
   config.JobType.psetName = 'runNtupler.py'
   config.JobType.outputFiles = ['TnP_ntuple.root']

   config.Data.inputDBS = 'global'
   config.Data.splitting = 'Automatic'
#   config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
   config.Data.publication = False
   config.Data.totalUnits = -1
   config.Data.outLFNDirBase = '/store/group/phys_egamma/Run3TriggerStudies/rasharma/Ntuples/' + name
   config.Data.allowNonValidInputDataset = True

   config.Site.storageSite = 'T2_CH_CERN'
 #  config.Site.blacklist = ['T2_BR_SPRACE', 'T2_US_Wisconsin', 'T1_RU_JINR', 'T2_RU_JINR', 'T2_EE_Estonia']

   listOfSamples.reverse()
   for sample in listOfSamples:

      config.General.requestName = sample
      config.Data.inputDataset = dataset[sample]
#      config.Data.unitsPerJob = lumisPerJob[sample]
      config.Data.outputDatasetTag = sample
      p = Process(target=submit, args=(config,))
      p.start()
      p.join()
