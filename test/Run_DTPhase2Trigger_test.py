### FOR DATA

#import FWCore.ParameterSet.Config as cms

#process = cms.Process("DTNT")
## process = cms.Process("DTNTandRPC")
##process = cms.Process("RECLUSTERIZATION")

 
#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load('Configuration/StandardSequences/Services_cff')

#process.load('Configuration/StandardSequences/GeometryRecoDB_cff')  ##  solve STA problem
#process.load('Configuration/EventContent/EventContent_cff')
#process.load("Geometry.DTGeometryBuilder.dtGeometryDB_cfi")
#process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

#process.load("EventFilter.DTRawToDigi.dtunpackerDDUGlobal_cfi")
#process.dtunpacker.readOutParameters.debug = True
#process.dtunpacker.readOutParameters.rosParameters.debug = True

## for DTTF (Not used from 2016)
#process.load("EventFilter.DTTFRawToDigi.dttfunpacker_cfi")
#process.load("EventFilter.DTTFRawToDigi.dttfpacker_cfi")

## for TWINMUX (Start to use in 2016)
#process.load("EventFilter.L1TXRawToDigi.twinMuxStage2Digis_cfi")

##for RAW
#process.load("Configuration.StandardSequences.Reconstruction_cff")

#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")  #DB v2, at least since GR_E_V42

#process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v2'

## for the emulator
#process.load("L1TriggerConfig.DTTPGConfigProducers.L1DTTPGConfigFromDB_cff")
#process.load("L1Trigger.DTTrigger.dtTriggerPrimitiveDigis_cfi")
#process.dtTriggerPrimitiveDigis.debug = True
#process.L1DTConfigFromDB.debug = True

#process.load('EventFilter.ScalersRawToDigi.ScalersRawToDigi_cfi')
#process.load('RecoLuminosity.LumiProducer.lumiProducer_cfi')

## process.load('EventFilter.L1TRawToDigi.l1tRawtoDigiBMTF_cfi')
#process.load('EventFilter.L1TRawToDigi.bmtfDigis_cfi')

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.source = cms.Source("PoolSource",

  #fileNames = cms.untracked.vstring
  #(
  #'file:/pool/cienciasrw/userstorage/vrbouza/ntuples/DT_ntuples/simulaciones/MuonGuns/CMSSW_10_1_1/src/UserCode/producers/Pt2_100_Eta1_24_B2B/step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM_5.root',
  #),
  #secondaryFileNames = cms.untracked.vstring(
  #)
#)

##this is to select collisions
#process.primaryVertexFilter = cms.EDFilter("VertexSelector",
   #src = cms.InputTag("offlinePrimaryVertices"),
   #cut = cms.string("!isFake && ndof > 4"), # && abs(z) <= 15 && position.Rho <= 2" # tracksSize() > 3 for the older cut
   #filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
#)

#process.noscraping = cms.EDFilter("FilterOutScraping",
                                  #applyfilter = cms.untracked.bool(True),
                                  #debugOn = cms.untracked.bool(False),
                                  #numtrack = cms.untracked.uint32(10),
                                  #thresh = cms.untracked.double(0.25)
#)

#process.DTMuonSelection = cms.EDFilter("DTMuonSelection",
                                 #src = cms.InputTag('muons'),
                                 #Muons = cms.InputTag('muons'),
                                 #SAMuons = cms.InputTag('standAloneMuons'),
                                 #dtSegmentLabel = cms.InputTag('dt4DSegments'),
                                 #etaMin = cms.double(-1.25),
                                 #etaMax = cms.double(1.25),
                                 #ptMin = cms.double(0.),#3.),
                                 #tightness = cms.int32(1) # 0 = loose (e.g. unstable collisions, minimum bias, requires a DT segment)
                                                          ## 1 = medium (e.g. cosmics, requires a stand alone muon)
                                                          ## 2 = tight (good collisions, requires a global muon)
#)


#process.load("UserCode/DTPhase2Trigger/DTPhase2Trigger_cfi")
#process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

#process.myDTNtuple.localDTmuons = cms.untracked.bool(False)
#process.myDTNtuple.outputFile = "DTNtuple.root"
#process.myDTNtuple.dtTrigSimDCCLabel = cms.InputTag("dtTriggerPrimitiveDigis")
#process.myDTNtuple.dtDigiLabel = cms.InputTag("dtunpacker")

#process.myDTNtuple.bmtfOutputDigis = cms.InputTag("bmtfDigis:BMTF")

### RPC unpacking
#process.load("EventFilter.RPCRawToDigi.rpcUnpackingModule_cfi")

### RPC recHit
#process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
#process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
#process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
#process.load("RecoLocalMuon.RPCRecHit.rpcRecHits_cfi")
#process.rpcRecHits.rpcDigiLabel = cms.InputTag('rpcUnpackingModule')

## process.p = cms.Path(process.DTMuonSelection * process.dtunpacker * process.twinMuxStage2Digis  * process.scalersRawToDigi * process.lumiProducer * process.dtTriggerPrimitiveDigis + process.BMTFStage2Digis + process.rpcUnpackingModule + process.rpcRecHits + process.myDTNtuple)
## process.p = cms.Path(process.DTMuonSelection * process.dtunpacker * process.twinMuxStage2Digis  * process.scalersRawToDigi * process.lumiProducer * process.dtTriggerPrimitiveDigis + process.bmtfDigis + process.rpcUnpackingModule + process.rpcRecHits + process.myDTNtuple)

#process.p = cms.Path(process.myDTNtuple)

### Output
##process.out = cms.OutputModule("PoolOutputModule"
                               ##, outputCommands = cms.untracked.vstring(
                               											##"keep *",
                                                                         ##"keep *_*_*_testRPCTwinMuxRawToDigi"
                                                                       ##, "keep *_*_*_DTNTandRPC"
																		##)
###                                , fileName = cms.untracked.string("file:cia.root")
                               ##, SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
##)





### FOR SIMULATION
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("DTNT",eras.Run2_2018)

process.load('Configuration/StandardSequences/Services_cff')

process.load('Configuration/StandardSequences/GeometryRecoDB_cff')  ##  solve STA problem
process.load('Configuration/EventContent/EventContent_cff')
process.load("Geometry.DTGeometryBuilder.dtGeometryDB_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

#for RAW
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")  #DB v2, at least since GR_E_V42

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic')

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",

  fileNames = cms.untracked.vstring
  (
    'file:/pool/cienciasrw/userstorage/vrbouza/ntuples/DT_ntuples/simulaciones/MuonGuns/CMSSW_10_1_1/src/UserCode/producers/Pt2_100_Eta1_24_B2B/step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM_5.root',
  ),
  secondaryFileNames = cms.untracked.vstring
  (
  )
)

process.load("UserCode/DTPhase2Trigger/DTPhase2Trigger_cfi")

process.myDTNtuple.localDTmuons = cms.untracked.bool(False)
process.myDTNtuple.outputFile = "DTNtuple.root"
process.myDTNtuple.bmtfOutputDigis = cms.InputTag("bmtfDigis:BMTF")
process.myDTNtuple.rpcLabel = cms.InputTag("muonRPCDigis")
process.myDTNtuple.runOnSimulation = cms.bool(True)


process.p = cms.Path(process.myDTNtuple)



