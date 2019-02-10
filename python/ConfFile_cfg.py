import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:/pool/ciencias/userstorage/vrbouza/ntuples/DT_ntuples/2018/Prompt/DTTree_SingleMuon_ZMuSkim315252-315510.root'
        #'file:/pool/ciencias/userstorage/vrbouza/ntuples/DT_ntuples/simulaciones/longevity_task_force/CMSSW_10_2_5/src/UserCode/DTDPGAnalysis/test/SingleMuPt_3_1000_Eta_1_24_N1M/SingleMuPt_3_1000_Eta_1_24_N1M_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM_0.root'
        'file:/pool/cienciasrw/userstorage/vrbouza/ntuples/DT_ntuples/simulaciones/MuonGuns/CMSSW_10_1_1/src/UserCode/producers/Pt2_100_Eta1_24_B2B/step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM_5.root'
    )
)

process.demo = cms.EDAnalyzer('DTPhase2Trigger',
                              tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)


process.p = cms.Path(process.demo)
