import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Ntupler")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# NOTE: the pick the right global tag!
#    for Spring15 50ns MC: global tag is 'auto:run2_mc_50'
#    for Spring15 25ns MC: global tag is 'auto:run2_mc'
#    for Run 2 data: global tag is 'auto:run2_data'
#  as a rule, find the "auto" global tag in $CMSSW_RELEASE_BASE/src/Configuration/AlCa/python/autoCond.py
#  This auto global tag will look up the "proper" global tag
#  that is typically found in the DAS under the Configs for given dataset
#  (although it can be "overridden" by requirements of a given release)
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc','')

#
# Define input data to read
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

#import FWCore.Utilities.FileUtils as FileUtils
#
#inputFilesAOD = cms.untracked.vstring()
#
#
#inputFilesAOD = cms.untracked.vstring( FileUtils.loadListFromFile (os.environ['CMSSW_BASE']+'/src/Efficiency/Analyzer/test/'+'inputFiles.txt') )
#
inputFilesAOD = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/data/Run2017C/SingleElectron/AOD/12Sep2017-v1/70000/80EF56E5-69A6-E711-AB37-48FD8E2824D7.root',
    )

inputFilesMiniAOD = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/data/Run2017B/SingleElectron/MINIAOD/31Mar2018-v1/90000/FC89D712-AF37-E811-AD13-008CFAC93F84.root'
#'root://cms-xrd-global.cern.ch//store/mc/Run3Summer19MiniAOD/DYToEE_M-50_NNPDF31_TuneCP5_14TeV-powheg-pythia8/MINIAODSIM/2021Scenario_106X_mcRun3_2021_realistic_v3-v2/270000/087644E5-95DB-E24C-B7A9-81646685DBA9.root'
#'root://cms-xrd-global.cern.ch//store/mc/Run3Summer19MiniAOD/DYToEE_M-50_NNPDF31_TuneCP5_14TeV-powheg-pythia8/MINIAODSIM/2023Scenario_106X_mcRun3_2023_realistic_v3-v2/260000/007DAFDA-8CA6-504C-A54C-4479B4E75DFF.root'
#'file:DY2023.root'
#'root://cms-xrd-global.cern.ch//store/mc/Run3Winter20DRPremixMiniAOD/DYJetsToEE_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/110X_mcRun3_2021_realistic_v6-v1/60000/02F993F1-A40E-7048-90EE-E0941D52577D.root'
# 'root://cms-xrd-global.cern.ch//store/mc/Run3Winter21DRMiniAOD/ZToEE_TuneCUETP8M1_14TeV-pythia8/MINIAODSIM/FlatPU20to70_for_DNN_112X_mcRun3_2021_realistic_v16_ext1-v1/100000/221bde1b-774c-486f-a666-45aa584a243f.root'
    'root://cms-xrd-global.cern.ch//store/mc/Run3Winter21DRMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/MINIAODSIM/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/70000/aed45160-b48b-4b42-8d56-77160031d2a7.root'
# 'root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/80000/E53E0D52-39FF-6F42-A8B2-EEF28EEB4C43.root'
    )
#
# You can list here either AOD or miniAOD files, but not both types mixed
#
useAOD = False
if useAOD == True :
    inputFiles = inputFilesAOD
    outputFile = "electron_ntuple.root"
    pileupProductName = "addPileupInfo"
    print("AOD input files are used")
else :
    inputFiles = inputFilesMiniAOD
    outputFile = "TnP_ntuple.root"
    pileupProductName = "slimmedAddPileupInfo"
    print("MiniAOD input files are used")
process.source = cms.Source ("PoolSource", fileNames = inputFiles )

#
# Set up electron ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
		 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
		 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',
		 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff',
		 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
]

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.minPtLeptons = cms.EDFilter("PtMinCandViewSelector",
#    src = cms.InputTag("slimmedMuons"),
    src = cms.InputTag("slimmedElectrons"),
    ptMin = cms.double(9)
  )

process.atLeastTwoLeptons = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("minPtLeptons"),
                                       minNumber = cms.uint32(2)
                                       )

#
# Configure the ntupler module
#
process.ntupler = cms.EDAnalyzer('Ntupler',
                                 # The module automatically detects AOD vs miniAOD, so we configure both
                                 #
                                 # Common to all formats objects
                                 #
                                 pileup   = cms.InputTag( pileupProductName ),
                                 rho      = cms.InputTag("fixedGridRhoFastjetAll"),
                                 beamSpot = cms.InputTag('offlineBeamSpot'),
                                 genEventInfoProduct = cms.InputTag('generator'),
                                 #
                                 # Objects specific to AOD format
                                 #
                                 electrons    = cms.InputTag("gedGsfElectrons"),
                                 genParticles = cms.InputTag("genParticles"),
                                 vertices     = cms.InputTag("offlinePrimaryVertices"),
                                 conversions  = cms.InputTag('allConversions'),
                                 triggerResultTag     = cms.InputTag("TriggerResults", "", "HLT"),
                                 triggerSummaryTag = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                                 l1EGTag      = cms.InputTag("caloStage2Digis","EGamma","RECO"),

				 pathsToSave  = cms.vstring( "HLT_Ele32_WPTight_Gsf_v",
								),
                                 filterToMatch= cms.vstring(
                                                          "hltEGL1SingleEGOrFilter",
                                                          "hltEG32L1SingleEGOrEtFilter",
                                                          "hltEle32WPTightClusterShapeFilter",
                                                          "hltEle32WPTightHEFilter",
                                                          "hltEle32WPTightEcalIsoFilter",
                                                          "hltEle32WPTightHcalIsoFilter",
                                                          "hltEle32WPTightPixelMatchFilter",
                                                          "hltEle32WPTightPMS2Filter",
                                                          "hltEle32WPTightGsfOneOEMinusOneOPFilter",
                                                          "hltEle32WPTightGsfMissingHitsFilter",
                                                          "hltEle32WPTightGsfDetaFilter",
                                                          "hltEle32WPTightGsfDphiFilter",
                                                          "hltEle32WPTightGsfTrackIsoFilter"
										),

				HLTprocess = cms.string("HLT"),

				#
                                 # Objects specific to MiniAOD format
                                #
                                 electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
                                 genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                 verticesMiniAOD     = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 conversionsMiniAOD  = cms.InputTag('reducedEgamma:reducedConversions'),
                                 trigger     = cms.InputTag("TriggerResults", "", "HLT"),
			         prescale = cms.InputTag("patTrigger"),
                                 objects = cms.InputTag('slimmedPatTrigger'),
                                 # Effective areas for computing PU correction for isolations
                                 effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
                                 # ID decisions (common to all formats)
                                 #
                                 # all IDs listed below are available given the content of "my_id_modules" defined above.
                                 # only one is exercised for this example.
                                 #

                                  eleIdMapLoose = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose"),
                                  eleIdMapMedium = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium"),
                                  eleIdMapTight = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight"),
				  eleMVA90noIso =  cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp90'),
				  eleMVA80noIso =  cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp80'),
				  eleMVA90Iso    =  cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-wp90'),
				  eleMVA80Iso    =  cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-wp80'),
				  eleMVALoosenoIso =  cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wpLoose'),
				  eleMVALooseIso =  cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-wpLoose'),
                                  egInputTag = cms.InputTag("caloStage2Digis","EGamma","RECO"),
                                 isMC = cms.bool(True),
                                 doEle = cms.bool(True)

                                 )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )


#process.p = cms.Path(process.egmGsfElectronIDSequence * process.ntupler)
process.p = cms.Path(process.egmGsfElectronIDSequence * process.minPtLeptons * process.atLeastTwoLeptons * process.ntupler)

