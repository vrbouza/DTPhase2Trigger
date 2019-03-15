import FWCore.ParameterSet.Config as cms

myDTNtuple = cms.EDAnalyzer('DTPhase2Trigger',
                            outputFile =cms.string("DTTree.root"),
                            #dtDigiLabel = cms.InputTag("muonDTDigis"),
                            dtDigiLabel = cms.InputTag("simMuonDTDigis"),
                            dtSegmentLabel = cms.InputTag("dt4DSegments"),
                            cscSegmentLabel = cms.InputTag("cscSegments"),
                            dtTrigTwinMuxInLabel= cms.InputTag("twinMuxStage2Digis","PhIn"),
                            dtTrigTwinMuxOutLabel = cms.InputTag("twinMuxStage2Digis","PhOut"),
                            dtTrigTwinMuxThLabel= cms.InputTag("twinMuxStage2Digis","ThIn"),
                            staMuLabel     = cms.InputTag("muons"),
                            # gmtLabel     = cms.InputTag("gtDigis"),  # legacy
                            gmtLabel     = cms.InputTag("gmtStage2Digis:Muon"),
                            gtLabel      = cms.InputTag("gtDigis"),  # legacy
                            rpcRecHitLabel = cms.InputTag("rpcRecHits", "", "RECO"),
                            dtDigiSize = cms.int32(1000),
                            dtSegmentSize = cms.int32(50),
                            cscSegmentSize = cms.int32(50),
                            dtTrigTwinMuxInSize = cms.int32(50),
                            dtTrigTwinMuxOutSize = cms.int32(50),
                            dtTrigTwinMuxThSize = cms.int32(50),
                            gmtSize       = cms.int32(50), # legacy
                            recoMuSize    = cms.int32(20),
                            rpcRecHitSize = cms.int32(300),
                            PrimaryVertexTag = cms.InputTag("offlinePrimaryVertices"),
                            TriggerTag = cms.InputTag("TriggerResults::HLT"),
                            trigFilterNames = cms.vstring("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q",
                                                          "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"),
                            TriggerEventTag = cms.InputTag("hltTriggerSummaryAOD::HLT"),
                            beamSpotTag  = cms.InputTag("offlineBeamSpot"),
                            puSummaryTag = cms.InputTag("addPileupInfo"),
                            scalersResults = cms.InputTag("scalersRawToDigi"),
                            # CB commented in code to avoid crashes
                            # lumiInputTag   = cms.InputTag("lumiProducer"),
                            runOnRaw = cms.bool(True),
                            runOnSimulation = cms.bool(False),
                            #runOnDigiSimLinks = cms.bool(False),
                            runOnDigiSimLinks = cms.bool(True),
                            localDTmuons    = cms.bool(False),
                            # bmtfInputPhDigis = cms.InputTag("BMTFStage2Digis"),
                            # bmtfInputThDigis = cms.InputTag("BMTFStage2Digis"),
                            # bmtfOutputDigis = cms.InputTag("BMTFStage2Digis"),
                            bmtfInputPhDigis = cms.InputTag("bmtfDigis"),
                            bmtfInputThDigis = cms.InputTag("bmtfDigis"),
                            bmtfOutputDigis = cms.InputTag("bmtfDigis"),
                            rpcLabel                  = cms.InputTag("rpcUnpackingModule"),
                            UnpackingRpcRecHitLabel = cms.InputTag("rpcRecHits", "", "DTNT"), # from unpacking RPC
                            OnlyBarrel = cms.bool(False), # Set False to save all rpc info (not only for Barrel)
                            #  UnpackingRpcRecHitLabel = cms.InputTag("rpcRecHits", "", "DTNTandRPC"), # from unpacking RPC
                            ## Parameters for retrieving the ttrig to correct the recHit times
                            tTrigModeConfig = cms.untracked.PSet(
                                  vPropWire = cms.double(24.4),
                                  doTOFCorrection = cms.bool(False),
                                  ##tofCorrType = cms.int32(2),  ## old
                                  tofCorrType = cms.int32(2),
                                  ##kFactor = cms.double(-1.3),
                                  wirePropCorrType = cms.int32(0),
                                  doWirePropCorrection = cms.bool(False),
                                  doT0Correction = cms.bool(True),
                                  tTrigLabel = cms.string(''),
                                  debug = cms.untracked.bool(False)
                              ),
                              tTrigMode = cms.untracked.string('DTTTrigSyncFromDB')
                            ## END Parameters for retrieving the ttrig to correct the recHit times
                              
)
