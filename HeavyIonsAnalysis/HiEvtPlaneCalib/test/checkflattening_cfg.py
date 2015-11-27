import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("FlatCalib")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
process.load("RecoHI.HiEvtPlaneAlgos.hiEvtPlaneFlat_cfi")
process.load("HeavyIonsAnalysis.HiEvtPlaneCalib/checkflattening_cfi")
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load('GeneratorInterface.HiGenCommon.HeavyIon_cff')
process.load("CondCore.DBCommon.CondDBCommon_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_mcRun2_HeavyIon_v5', '')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")

process.MessageLogger.cerr.FwkReport.reportEvery=1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.GlobalTag.toGet.extend([
   cms.PSet(record = cms.string("HeavyIonRcd"),
      tag = cms.string("CentralityTable_HFtowers200_HydjetDrum5_v750x02_mc"),
      connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
      label = cms.untracked.string("HFtowers")
   ),
    cms.PSet(record = cms.string("HeavyIonRPRcd"),
                 tag = cms.string("HeavyIonRPRcd_Hydjet_75x_v03_mc"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
                 )
])

#process.GlobalTag.toGet = cms.VPSet(
#    cms.PSet(record = cms.string("HeavyIonRPRcd"),
#                 tag = cms.string("HeavyIonRPRcd_Hydjet_74x_v03_mc"),
#                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
#                 )
#        )

#process.CondDBCommon.connect = "sqlite_file:HeavyIonRPRcd_Hydjet_74x_v03_mc.db"
#process.PoolDBESSource = cms.ESSource("PoolDBESSource",
#                                       process.CondDBCommon,
#                                       toGet = cms.VPSet(cms.PSet(record = cms.string('HeavyIonRPRcd'),
#                                                                  tag = cms.string('HeavyIonRPRcd_Hydjet_74x_v03_mc')
#                                                                  )
#                                                         )
#                                      )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_1.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_10.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_100.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_101.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_102.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_103.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_104.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_105.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_106.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_107.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_108.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_109.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_11.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_110.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_111.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_112.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_113.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_114.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_115.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_116.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_117.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_118.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_119.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_12.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_120.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_121.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_122.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_123.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_124.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_125.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_126.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_127.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_128.root',
       'root://xrootd.unl.edu//store/user/tuos/HIAOD2015/round3/June01/HydjetMB/Hydjet_Quenched_MinBias_5020GeV/HydjetMB_AOD_750pre5_round3v01/150601_223002/0000/step2_RAW2DIGI_L1Reco_MB_AODSIM_129.root'
 ),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            inputCommands=cms.untracked.vstring(
        'keep *',
        'drop *_hiEvtPlane_*_*'
        ),
                            dropDescendantsOfDroppedBranches=cms.untracked.bool(False)
                            )



process.TFileService = cms.Service("TFileService",
    fileName = cms.string("check.root")
)

process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.hiEvtPlane.loadDB_ = cms.untracked.bool(True)
process.hiEvtPlaneFlat.genFlatPsi_ = cms.untracked.bool(True)
process.p = cms.Path(process.centralityBin*process.hiEvtPlane*process.hiEvtPlaneFlat*process.checkflattening)


                        

