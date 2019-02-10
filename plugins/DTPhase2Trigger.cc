// -*- C++ -*-
//
// Package:    L1Trigger/DTPhase2Trigger
// Class:      DTPhase2Trigger
//
/**\class DTPhase2Trigger DTPhase2Trigger.cc L1Trigger/DTPhase2Trigger/plugins/DTPhase2Trigger.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Santiago Folgueras
//         Created:  Thu, 06 Dec 2018 09:44:33 GMT
//
//

// user include files (from CMSSW)
#include "UserCode/DTPhase2Trigger/interface/DTPhase2Trigger.h"



//
// constructors and destructor
//
DTPhase2Trigger::DTPhase2Trigger(const edm::ParameterSet& pset):
      rpcToken_(consumes<MuonDigiCollection<RPCDetId,RPCDigi> > (pset.getParameter<edm::InputTag>("rpcLabel"))),
      UnpackingRpcRecHitToken_(consumes<RPCRecHitCollection> (pset.getParameter<edm::InputTag>("UnpackingRpcRecHitLabel"))) {
   // now do whatever initialization is needed
//   my_debug = pset.getUntrackedParameter<bool>("debug");
//   my_DTTFnum = pset.getParameter<bool>("DTTFSectorNumbering");
//   my_params = pset;
    
//   digiLabel_ = pset.getParameter<edm::InputTag>("digiTag");
//   dt4DSegments = consumes<DTRecSegment4DCollection>(pset.getParameter < edm::InputTag > ("dt4DSegments"));
  
//   string outputfile = pset.getUntrackedParameter<string>("outputFileName");
//   my_rootfile = new TFile(outputfile.c_str(),"RECREATE");
//   my_tree = new TTree("tree", "L1T", 0);

  // get the tTrigDBInfo
//   theSync = DTTTrigSyncFactory::get()->create(pset.getUntrackedParameter<std::string>("tTrigMode"), pset.getUntrackedParameter<edm::ParameterSet>("tTrigModeConfig"));
  
  
  // Ttree gen. constructor     =================================================
  // get the tTrigDBInfo
  theSync = DTTTrigSyncFactory::get()->create(pset.getUntrackedParameter<std::string>("tTrigMode"),
                                pset.getUntrackedParameter<edm::ParameterSet>("tTrigModeConfig"));


  // get run configuration options
  runOnRaw_          = pset.getParameter<Bool_t>("runOnRaw");
  runOnSimulation_   = pset.getParameter<Bool_t>("runOnSimulation");
  runOnDigiSimLinks_ = pset.getParameter<Bool_t>("runOnDigiSimLinks");

  //get parameters from the configuration file
  //names of the different event collections
  dtDigiLabel_     = pset.getParameter<edm::InputTag>("dtDigiLabel");

  dtDigiToken_     = consumes<DTDigiCollection>(edm::InputTag(dtDigiLabel_));
  dtDigiTokenSim_  = consumes<DTDigiSimLinkCollection>(edm::InputTag(dtDigiLabel_));

  dtSegmentLabel_  = pset.getParameter<edm::InputTag>("dtSegmentLabel");
  dtSegmentToken_  = consumes<DTRecSegment4DCollection>(edm::InputTag(dtSegmentLabel_));

  cscSegmentLabel_ = pset.getParameter<edm::InputTag>("cscSegmentLabel");
  cscSegmentToken_ = consumes<CSCSegmentCollection>(edm::InputTag(cscSegmentLabel_));

  dtTrigTwinMuxInLabel_    = pset.getParameter<edm::InputTag>("dtTrigTwinMuxInLabel");
  dtTrigTwinMuxThLabel_    = pset.getParameter<edm::InputTag>("dtTrigTwinMuxThLabel");
  dtTrigTwinMuxInToken_    = consumes<L1MuDTChambPhContainer>(edm::InputTag(dtTrigTwinMuxInLabel_)); 
  dtTrigTwinMux_ThToken_ = consumes<L1MuDTChambThContainer>(edm::InputTag(dtTrigTwinMuxThLabel_));

  dtTrigTwinMuxOutLabel_  = pset.getParameter<edm::InputTag>("dtTrigTwinMuxOutLabel");
  dtTrigTwinMuxOutToken_  = consumes<L1MuDTChambPhContainer>(edm::InputTag(dtTrigTwinMuxOutLabel_));

  staMuLabel_      = pset.getParameter<edm::InputTag>("staMuLabel");
  staMuToken_      = consumes<reco::MuonCollection>(edm::InputTag(staMuLabel_));

  gmtLabel_        = pset.getParameter<edm::InputTag>("gmtLabel");
  gmtToken_        = consumes<l1t::MuonBxCollection>(edm::InputTag(gmtLabel_));

  triggerTag_      = pset.getParameter<edm::InputTag>("TriggerTag");
  triggerToken_    = consumes<edm::TriggerResults>(edm::InputTag(triggerTag_));

  triggerEventTag_   = pset.getParameter<edm::InputTag>("TriggerEventTag");
  triggerEventToken_ = consumes<trigger::TriggerEvent>(edm::InputTag(triggerEventTag_));
  trigFilterNames_    = pset.getParameter<std::vector<std::string>>("trigFilterNames");

  gtLabel_         = pset.getParameter<edm::InputTag>("gtLabel"); // legacy
  gtToken_         = consumes<L1GlobalTriggerReadoutRecord>(edm::InputTag(gtLabel_)); //legacy

  rpcRecHitLabel_  = pset.getParameter<edm::InputTag>("rpcRecHitLabel");
  rpcRecHitToken_  = consumes<RPCRecHitCollection>(edm::InputTag(rpcRecHitLabel_));

  //max size of the different saved objects (per event)
  digisSize_          = pset.getParameter<Int_t>("dtDigiSize");
  dtsegmentsSize_     = pset.getParameter<Int_t>("dtSegmentSize");
  cscsegmentsSize_    = pset.getParameter<Int_t>("cscSegmentSize");
  dtltTwinMuxOutSize_ = pset.getParameter<Int_t>("dtTrigTwinMuxOutSize");
  dtltTwinMuxInSize_  = pset.getParameter<Int_t>("dtTrigTwinMuxInSize");
  dtltTwinMuxThSize_  = pset.getParameter<Int_t>("dtTrigTwinMuxThSize");
  gmtSize_            = pset.getParameter<Int_t>("gmtSize"); 
  recoMuSize_         = pset.getParameter<Int_t>("recoMuSize");
  rpcRecHitSize_      = pset.getParameter<Int_t>("rpcRecHitSize"); 

  PrimaryVertexTag_   = pset.getParameter<edm::InputTag>("PrimaryVertexTag");
  PrimaryVertexToken_ = consumes<reco::VertexCollection>(edm::InputTag(PrimaryVertexTag_));

  beamSpotTag_        = pset.getParameter<edm::InputTag>("beamSpotTag");
  beamSpotToken_      = consumes<reco::BeamSpot>(edm::InputTag(beamSpotTag_));

  scalersSource_      = pset.getParameter<edm::InputTag>("scalersResults");
  scalersSourceToken_ = consumes<LumiScalersCollection>(edm::InputTag(scalersSource_));

  puSummaryTag_       = pset.getParameter<edm::InputTag>("puSummaryTag");
  puSummaryToken_     = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag(puSummaryTag_));
       
  // CB comment for now to avoid crash
  // lumiInputTag_      = pset.getParameter<edm::InputTag>("lumiInputTag");
  // lumiProducerToken_ = consumes<LumiDetails, edm::InLumi>(lumiInputTag_);
  
  bmtfPhInputTag_   = consumes<L1MuDTChambPhContainer>(pset.getParameter<edm::InputTag>("bmtfInputPhDigis"));
  bmtfThInputTag_   = consumes<L1MuDTChambThContainer>(pset.getParameter<edm::InputTag>("bmtfInputThDigis"));
  bmtfOutputTag_    = consumes<l1t::RegionalMuonCandBxCollection>(pset.getParameter<edm::InputTag>("bmtfOutputDigis"));

  OnlyBarrel_       = pset.getParameter<Bool_t>("OnlyBarrel");

  localDTmuons_     = pset.getUntrackedParameter<Bool_t>("localDTmuons", false);

  AnaTrackGlobalMu_ = pset.getUntrackedParameter<Bool_t>("AnaTrackGlobalMu", true);  // set to false when problems with tracks of global muons

  runLegacy_gmt_    = pset.getUntrackedParameter<Bool_t>("runLegacy_gmt", true); // Collection not available in cosmics 2017, at least on september runs
                                                                             // set to false in cosmics if needed
  outFile_          = pset.getParameter<std::string>("outputFile");

  initialize_Tree_variables();

  //counters
  idigis          = 0;
  idtsegments     = 0;
  icscsegments    = 0;
  idtltTwinMuxOut = 0;
  idtltTwinMuxIn  = 0;
  idtltTwinMux_th = 0;
  imuons          = 0;
  igmt            = 0;
  igtalgo         = 0; // legacy
  igttt           = 0; // legacy
  ihlt            = 0;
  // END CONSTRUCTOR =============================================
  
  hlpr = new TransformHelper();
  hlpr->doVarious = false;
  
}


DTPhase2Trigger::~DTPhase2Trigger() {
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
//   delete my_rootfile;
//   if (my_debug) cout << "[DTTrigTest] Destructor executed!!!" << endl;
  
  delete hlpr;
  
}


//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void
DTPhase2Trigger::beginJob()
{
//   h_digiTDC        = new TH1F("h_digiTDC","h_digiTDC",1601,-0.5,1600.5);
//   h_digiTDCPhase2  = new TH1F("h_digiTDCPhase2","h_digiTDCPhase2",3563*32,-0.5,3563*32+1);
// 
//   h_digiTime       = new TH1F("h_digiTime","h_digiTime",1275,-0.5,1274.5);
//   h_digiTimePhase2 = new TH1F("h_digiTimePhase2","h_digiTimePhase2",8907,-0.5,89075.5);
  
  
  initialize_Tree_variables();
  hlpr->initialise();
  
}


// ------------ method called for each event  ------------
void
DTPhase2Trigger::analyze(const edm::Event& event, const edm::EventSetup& context) {
  
  // Santi's code
//   Handle<DTDigiCollection> dtdigis;
//   iEvent.getByLabel(digiLabel_, dtdigis);
//   
//   Handle<DTRecSegment4DCollection> all4DSegments;
//   iEvent.getByToken(dt4DSegments, all4DSegments);
// 
// 
//   // The code (for now), will focus in one chamber and
//   // look at segments / digis there
//   DTChamberId selected_chamber_ID(0,1,6);//ciemat_chamber
// 
//   // Select event if there is a 4D segment in desired chamber
//   int DTSegmentCounterInChamber=0;
//   vector<const DTRecSegment4D*> my_segments4D;
//   for (DTRecSegment4DCollection::const_iterator segm = all4DSegments->begin(); segm!=all4DSegments->end(); ++segm){
//     if (segm->chamberId()!=selected_chamber_ID) continue;
//     if (!segm->hasPhi()) continue;
//     if (!segm->hasZed()) continue;
// 
//     DTSegmentCounterInChamber++;
//     my_segments4D.push_back(&(*segment1));
//   }
// 
//   //   DTRecSegment4DCollection::const_iterator segment;
//   //   for (segment = all4DSegments->begin();segment!=all4DSegments->end(); ++segment){}     
//   
//   /// FOCUS ON SELECTED CHAMBER:
//   DTDigiCollection::DigiRangeIterator dtLayerId_It;
//   for (dtLayerId_It=dtdigis->begin(); dtLayerId_It!=dtdigis->end(); ++dtLayerId_It){
//     for (DTDigiCollection::const_iterator digiIt = ((*dtLayerId_It).second).first;digiIt!=((*dtLayerId_It).second).second; ++digiIt){
//       const DTLayerId dtLId = (*dtLayerId_It).first;
//     
//       if (dtLId.wheel()!=selected_chamber_ID.wheel()) continue;
//       if (dtLId.sector()!=selected_chamber_ID.sector()) continue;
//       if (dtLId.station()!=selected_chamber_ID.station()) continue;
// 
//       int l  = dtLId.layer();
//       int sl = dtLId.superlayer();
//       
//       int digiTDC = (*digiIt).countsTDC();
//       int digiTDCPhase2 =  (*digiIt).countsTDC()+ iEvent.eventAuxiliary().bunchCrossing()*32;
// 
//       float ttrig = theSync->offset((*digiIt.wireId()));
//       float digiTIME = (*digiIt).time();
//       float digiTIMEPhase2 =  (*digiIt).time()+ iEvent.eventAuxiliary().bunchCrossing()*25-ttrig;//how to get the value of other station/chamber?
//       
//       
//       h_digiTDC->Fill(digiTDC); 
//       h_digiTDCPhase2->Fill(digiTDCPhase2);  
//       
//       h_digiTime->Fill(digiTIME); 
//       h_digiTimePhase2->Fill(digiTIMEPhase2);
//     }
//   } // end DTDigiCollection
//   
//   
//   // NOW build the L1
  
  // DTTREE PRODUCTION PRELIMINARIES ==========================================
  theSync->setES(context);

  edm::ESHandle<DTGeometry> dtGeomH;
  context.get<MuonGeometryRecord>().get(dtGeomH);
//   const DTGeometry* dtGeom_ = dtGeomH.product();

  //retrieve the beamspot info
  if (!localDTmuons_) {
    edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
    //event.getByLabel(beamSpotTag_, recoBeamSpotHandle);  // Doesn't work after 75X
    event.getByToken(beamSpotToken_, recoBeamSpotHandle);
    beamspot = *recoBeamSpotHandle;
    
    //retrieve the luminosity
    if (!runOnSimulation_) {                          // crashes with simulation trying to get the lumiperblock
      edm::Handle<LumiScalersCollection> lumiScalers; // The vector<LumiScalers> but not the <LumiScalersCollection>
//       event.getByLabel(scalersSource_, lumiScalers);  // Doesn't work after 75X
      event.getByToken(scalersSourceToken_, lumiScalers);
      LumiScalersCollection::const_iterator lumiIt = lumiScalers->begin();
      lumiperblock = lumiIt->instantLumi();
    }

    //retrieve gen PU info in MC
    if (runOnSimulation_) {
      edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
      event.getByToken(puSummaryToken_, puInfo);
      
      for (const auto & puInfoBx : (*puInfo)) {
        Int_t bx = puInfoBx.getBunchCrossing();
        if(bx == 0) {
          true_pileup   = puInfoBx.getTrueNumInteractions();
          actual_pileup = puInfoBx.getPU_NumInteractions();
          break;
        }
      }
    }
  }

  //retrieve the collections you are interested on in the event
  edm::Handle<DTDigiCollection> dtdigis;
  edm::Handle<MuonDigiCollection<DTLayerId, DTDigiSimLink>> dtdigisSim;
  if (runOnRaw_ && !runOnDigiSimLinks_)            event.getByToken(dtDigiToken_, dtdigis);
  else if (runOnSimulation_ && runOnDigiSimLinks_) event.getByToken(dtDigiTokenSim_, dtdigisSim); // Added to cope with simulation including Digis

  edm::Handle<DTRecSegment4DCollection> dtsegments4D;
  event.getByToken(dtSegmentToken_, dtsegments4D);

  context.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  edm::Handle<reco::VertexCollection> privtxs;
  if (!localDTmuons_) event.getByToken(PrimaryVertexToken_, privtxs);
  
  edm::Handle<CSCSegmentCollection> cscsegments;
  if (!localDTmuons_) event.getByToken(cscSegmentToken_, cscsegments);

//   edm::Handle<L1MuDTChambPhContainer> localTriggerTwinMuxOut;
//   Bool_t hasPhiTwinMuxOut = false;
//   if (runOnRaw_) hasPhiTwinMuxOut = event.getByToken(dtTrigTwinMuxOutToken_, localTriggerTwinMuxOut);
// 
//   edm::Handle<L1MuDTChambPhContainer> localTriggerTwinMuxIn;
//   Bool_t hasPhiTwinMuxIn = false;
//   if (runOnRaw_) hasPhiTwinMuxIn = event.getByToken(dtTrigTwinMuxInToken_, localTriggerTwinMuxIn);
// 
//   edm::Handle<L1MuDTChambThContainer> localTriggerTwinMux_Th;
//   Bool_t hasThetaTwinMux = false;
//   if (runOnRaw_) hasThetaTwinMux = event.getByToken(dtTrigTwinMux_ThToken_, localTriggerTwinMux_Th);

  edm::Handle<reco::MuonCollection> MuList;
  if (!localDTmuons_) event.getByToken(staMuToken_, MuList);

  edm::Handle<l1t::MuonBxCollection> gmt;   // legacy
  if (!localDTmuons_ && runLegacy_gmt_) event.getByToken(gmtToken_, gmt); // legacy

  edm::Handle< L1GlobalTriggerReadoutRecord > gtrc; // legacy
  if (runOnRaw_ && !localDTmuons_) event.getByToken(gtToken_, gtrc); // legacy

//   edm::ESHandle<L1GtTriggerMenu> menuRcd;
//   context.get<L1GtTriggerMenuRcd>().get(menuRcd);
//   const L1GtTriggerMenu* menu = menuRcd.product();

  edm::Handle<edm::TriggerResults> hltresults;
  if (!localDTmuons_) event.getByToken(triggerToken_, hltresults);

  edm::Handle<trigger::TriggerEvent> hltEvent;
  if (!localDTmuons_) event.getByToken(triggerEventToken_, hltEvent);

  edm::Handle<RPCRecHitCollection> rpcHits;
  if (!localDTmuons_) event.getByToken(rpcRecHitToken_, rpcHits);

  //clear the containers
  clear_Arrays();

  // Get the propagators
  context.get<TrackingComponentsRecord>().get("SmartPropagatorAny", propagatorAlong);
  context.get<TrackingComponentsRecord>().get("SmartPropagatorAnyOpposite", propagatorOpposite);

  //get the magnetic field
  context.get<IdealMagneticFieldRecord>().get(theBField);

  //Fill the event info block
  runnumber   = event.run();
  lumiblock   = event.getLuminosityBlock().luminosityBlock();
  eventNumber = event.eventAuxiliary().event();
  timestamp   = event.eventAuxiliary().time().value();
  bunchXing   = event.eventAuxiliary().bunchCrossing();
  orbitNum    = event.eventAuxiliary().orbitNumber();

  // CB comment for now to avoid crash
  // it's filling nothing at the moment lumiDetails->isValid() return false (30/04/2016  M.C.F)
  // if(!localDTmuons_ && !runOnSimulation_)  // Crashes when run in simulation no <LumiDetails> available
  // {
  //    edm::Handle<LumiDetails> lumiDetails;
  //    event.getLuminosityBlock().getByToken(lumiProducerToken_, lumiDetails);
  //    if(lumiDetails->isValid()){
  //      beam1Intensity = lumiDetails->lumiBeam1Intensity(bunchXing);
  //      beam2Intensity = lumiDetails->lumiBeam2Intensity(bunchXing);
  //    }
  // }

  //Primary vertex
  if (!localDTmuons_) {
    if ((*privtxs).size() != 0) {
      PV_x = (*privtxs)[0].position().x();
      PV_y = (*privtxs)[0].position().y();
      PV_z = (*privtxs)[0].position().z();
      
      PV_xxE = (*privtxs)[0].covariance(0,0);
      PV_yyE = (*privtxs)[0].covariance(1,1);
      PV_zzE = (*privtxs)[0].covariance(2,2);
      PV_xyE = (*privtxs)[0].covariance(0,1);
      PV_xzE = (*privtxs)[0].covariance(0,2);
      PV_yzE = (*privtxs)[0].covariance(1,2);
      
      PV_normchi2 = (*privtxs)[0].chi2() / (*privtxs)[0].ndof();
      
      PV_Nvtx = (*privtxs).size();
    }
    else {
      PV_x   = -999.;
      PV_y   = -999.;
      PV_z   = -999.;
      PV_xxE = -999.;
      PV_yyE = -999.;
      PV_zzE = -999.;
      PV_xyE = -999.;
      PV_xzE = -999.;
      PV_yzE = -999.;
      PV_normchi2 = -999.;
      PV_Nvtx = -999;
    }
  }
  // END OF DTTREE PRODUCTION PRELIMINARIES ==========================================
  
  hlpr->run(segm4D_wheel, segm4D_sector, segm4D_station, segm4D_phiHits_SuperLayer, segm4D_phiHits_Layer, segm4D_phiHits_Pos, segm4D_phiHits_Side, segm4D_phiHits_Wire, segm4D_phinhits);
  
}



// ------------ method called once each job just after ending the event loop  ------------
void
DTPhase2Trigger::endJob() {
//   my_rootfile->cd();

//   h_digiTDC->Write();
//   h_digiTDCPhase2->Write();
//   
//   h_digiTime->Write();
//   h_digiTimePhase2->Write();
  
  hlpr->finish();
  
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DTPhase2Trigger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}


void DTPhase2Trigger::clear_Arrays() {
  //digi variables
  digi_wheel.clear();
  digi_sector.clear();
  digi_station.clear();
  digi_sl.clear();
  digi_layer.clear();
  digi_wire.clear();
  digi_time.clear();

  //DT segment variables 
  segm4D_wheel.clear();
  segm4D_sector.clear();
  segm4D_station.clear();
  segm4D_hasPhi.clear();
  segm4D_hasZed.clear();
  segm4D_x_pos_loc.clear();
  segm4D_y_pos_loc.clear();
  segm4D_z_pos_loc.clear();
  segm4D_x_dir_loc.clear();
  segm4D_y_dir_loc.clear();
  segm4D_z_dir_loc.clear();
  segm4D_cosx.clear();
  segm4D_cosy.clear();
  segm4D_cosz.clear();
  segm4D_phi.clear();
  segm4D_theta.clear();
  segm4D_eta.clear();
  segm4D_t0.clear();
  segm4D_vdrift.clear();
  segm4D_phinormchi2.clear();
  segm4D_phinhits.clear();
  segm4D_znormchi2.clear();
  segm4D_znhits.clear();

  segm4D_hitsExpPos->Clear();
  segm4D_hitsExpWire->Clear();

  segm4D_phiHits_Pos->Clear();
  segm4D_phiHits_PosCh->Clear();
  segm4D_phiHits_PosErr->Clear();
  segm4D_phiHits_Side->Clear();
  segm4D_phiHits_Wire->Clear();
  segm4D_phiHits_Layer->Clear();
  segm4D_phiHits_SuperLayer->Clear();
  segm4D_phiHits_Time->Clear();
  segm4D_phiHits_TimeCali->Clear();
  segm4D_hitsExpPos->Clear();
  segm4D_hitsExpWire->Clear();

  segm4D_zHits_Pos->Clear();
  segm4D_zHits_PosCh->Clear();
  segm4D_zHits_PosErr->Clear();
  segm4D_zHits_Side->Clear();
  segm4D_zHits_Wire->Clear();
  segm4D_zHits_Layer->Clear();
  segm4D_zHits_Time->Clear();
  segm4D_zHits_TimeCali->Clear();

  //CSC segment variables 
  cscsegm_ring.clear();
  cscsegm_chamber.clear();
  cscsegm_station.clear();
  cscsegm_cosx.clear();
  cscsegm_cosy.clear();
  cscsegm_cosz.clear();
  cscsegm_phi.clear();
  cscsegm_eta.clear();
  cscsegm_normchi2.clear();
  cscsegm_nRecHits.clear();

  //TM Variables
  ltTwinMuxIn_wheel.clear();
  ltTwinMuxIn_sector.clear();
  ltTwinMuxIn_station.clear();
  ltTwinMuxIn_quality.clear();
  ltTwinMuxIn_bx.clear();
  ltTwinMuxIn_phi.clear();
  ltTwinMuxIn_phiB.clear();
  ltTwinMuxIn_is2nd.clear();

  ltTwinMuxOut_wheel.clear();
  ltTwinMuxOut_sector.clear();
  ltTwinMuxOut_station.clear();
  ltTwinMuxOut_quality.clear();
  ltTwinMuxOut_rpcbit.clear();
  ltTwinMuxOut_bx.clear();
  ltTwinMuxOut_phi.clear();
  ltTwinMuxOut_phiB.clear();
  ltTwinMuxOut_is2nd.clear();

  ltTwinMux_thWheel.clear();
  ltTwinMux_thSector.clear();
  ltTwinMux_thStation.clear();
  ltTwinMux_thBx.clear();
  ltTwinMux_thHits.clear();

  //muon variables
  Mu_isMuGlobal.clear();
  Mu_isMuTracker.clear();
  Mu_isMuTrackerArb.clear();
  Mu_isMuStandAlone.clear();
  Mu_isMuRPC.clear();
  
  Mu_px_mu.clear();
  Mu_py_mu.clear();
  Mu_pz_mu.clear();
  Mu_phi_mu.clear();
  Mu_eta_mu.clear();
  Mu_chargeMu.clear();

  Mu_nMatches.clear();
  Mu_numberOfChambers.clear();
  Mu_numberOfMatches.clear();
  Mu_numberOfMatchedStations.clear();
  Mu_stationMask.clear();

  STAMu_numberOfHits.clear();
  STAMu_segmIndex.clear();

  STAMu_recHitsSize.clear();
  STAMu_normchi2Mu.clear();
  STAMu_dxyMu.clear();
  STAMu_dzMu.clear();

  GLBMu_normchi2Mu.clear();
  GLBMu_dxyMu.clear();
  GLBMu_dzMu.clear();

  GLBMu_numberOfPixelHits.clear();
  GLBMu_numberOfTrackerHits.clear();

  GLBMu_tkIsoR03.clear();
  GLBMu_ntkIsoR03.clear();
  GLBMu_emIsoR03.clear();
  GLBMu_hadIsoR03.clear();

  TRKMu_normchi2Mu.clear();
  TRKMu_dxyMu.clear();
  TRKMu_dzMu.clear();

  TRKMu_numberOfPixelHits.clear();
  TRKMu_numberOfTrackerLayers.clear();

  TRKMu_tkIsoR03.clear();

  TRKMu_algo.clear();
  TRKMu_origAlgo.clear();

  STAMu_caloCompatibility.clear();
  STAMu_time.clear();
  STAMu_timeNDof.clear();

  RPCMu_numberOfRPCLayers.clear();

  STAMu_z_mb2.clear();
  STAMu_phi_mb2.clear();
  STAMu_pseta_mb2.clear();

  Mu_matches_Wh->Clear(); 
  Mu_matches_Sec->Clear();    
  Mu_matches_St->Clear();     

  Mu_matches_x->Clear();      
  Mu_matches_y->Clear();      

  Mu_matches_phi->Clear();    
  Mu_matches_eta->Clear();    

  Mu_matches_edgeX->Clear();  
  Mu_matches_edgeY->Clear();

  Mu_hlt_Dr->Clear();  

  //GMT
  gmt_bx.clear();
  gmt_phi.clear();
  gmt_eta.clear();
  gmt_pt.clear();
  gmt_qual.clear();
  gmt_charge.clear();
  gmt_tf_idx.clear();

  //GT
  gt_algo_bit.clear();
  gt_algo_bx.clear();
  gt_tt_bit.clear();
  gt_tt_bx.clear();

  //HLT
  hlt_path.clear();

  // RPC rec hits
  rpc_region.clear();
  rpc_clusterSize.clear();
  rpc_strip.clear();
  rpc_bx.clear();
  rpc_station.clear();
  rpc_sector.clear();
  rpc_layer.clear();
  rpc_subsector.clear();
  rpc_roll.clear();
  rpc_ring.clear();
  
        //Bmtf_Size.clear();
  Bmtf_Pt.clear();
  Bmtf_Eta.clear();
  Bmtf_Phi.clear();
  Bmtf_GlobalPhi.clear();
  Bmtf_qual.clear();
  Bmtf_ch.clear();
  Bmtf_bx.clear();
  Bmtf_processor.clear();
  Bmtf_trAddress.clear();
  Bmtf_wh.clear();
  Bmtf_FineBit.clear();

  Bmtf_phBx.clear();
  Bmtf_phWh.clear();
  Bmtf_phSe.clear();
  Bmtf_phSt.clear();
  Bmtf_phAng.clear();
  Bmtf_phBandAng.clear();
  Bmtf_phCode.clear();
  Bmtf_phTs2Tag.clear();

  Bmtf_thBx.clear();
  Bmtf_thWh.clear();
  Bmtf_thSe.clear();
  Bmtf_thSt.clear();
  Bmtf_thTheta.clear();
  Bmtf_thCode.clear();

  RpcDigi_TwinMux_bx.clear();
  RpcDigi_TwinMux_strip.clear();
  RpcDigi_TwinMux_region.clear();
  RpcDigi_TwinMux_ring.clear();
  RpcDigi_TwinMux_station.clear();
  RpcDigi_TwinMux_layer.clear();
  RpcDigi_TwinMux_sector.clear();
  RpcDigi_TwinMux_subsector.clear();
  RpcDigi_TwinMux_roll.clear();
  RpcDigi_TwinMux_trIndex.clear();
  RpcDigi_TwinMux_det.clear();
  RpcDigi_TwinMux_subdetId.clear();
  RpcDigi_TwinMux_rawId.clear();
  
  // Unpacking RPC rec hits
  RpcRechit_TwinMux_region.clear();
  RpcRechit_TwinMux_clusterSize.clear();
  RpcRechit_TwinMux_strip.clear();
  RpcRechit_TwinMux_bx.clear();
  RpcRechit_TwinMux_station.clear();
  RpcRechit_TwinMux_sector.clear();
  RpcRechit_TwinMux_layer.clear();
  RpcRechit_TwinMux_subsector.clear();
  RpcRechit_TwinMux_roll.clear();
  RpcRechit_TwinMux_ring.clear();
  
  return;
}


void DTPhase2Trigger::initialize_Tree_variables() {
  //Event variables
  runnumber   = 0;
  lumiblock   = 0;
  eventNumber = 0;
  timestamp   = 0.;
  bunchXing   = 0;
  orbitNum    = 0;
  true_pileup   = 0.;
  actual_pileup = 0;

  PV_x = 0.;
  PV_y = 0.;
  PV_z = 0.;
  PV_xxE = 0.;
  PV_yyE = 0.;
  PV_zzE = 0.;
  PV_xyE = 0.;
  PV_xzE = 0.;
  PV_yzE = 0.;
  PV_normchi2 = 0.;
  PV_Nvtx = 0.;

  lumiperblock = 0.;
  beam1Intensity = -1.;
  beam2Intensity = -1.;

  segm4D_phiHits_Pos    = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_PosCh  = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_PosErr = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_Side   = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_Wire   = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_Layer  = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_SuperLayer  = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_Time   = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_TimeCali   = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_hitsExpPos     = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_hitsExpWire    = new TClonesArray("TVectorF",dtsegmentsSize_);

  segm4D_zHits_Pos    = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_zHits_PosCh  = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_zHits_PosErr = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_zHits_Side   = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_zHits_Wire   = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_zHits_Layer  = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_zHits_Time   = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_zHits_TimeCali   = new TClonesArray("TVectorF",dtsegmentsSize_);

  Mu_matches_Wh  = new TClonesArray("TVectorF",recoMuSize_);
  Mu_matches_Sec = new TClonesArray("TVectorF",recoMuSize_);
  Mu_matches_St  = new TClonesArray("TVectorF",recoMuSize_);

  Mu_matches_x  = new TClonesArray("TVectorF",recoMuSize_);
  Mu_matches_y  = new TClonesArray("TVectorF",recoMuSize_);

  Mu_matches_phi  = new TClonesArray("TVectorF",recoMuSize_);
  Mu_matches_eta  = new TClonesArray("TVectorF",recoMuSize_);

  Mu_matches_edgeX  = new TClonesArray("TVectorF",recoMuSize_);
  Mu_matches_edgeY  = new TClonesArray("TVectorF",recoMuSize_);

  Mu_hlt_Dr  = new TClonesArray("TVectorF",recoMuSize_);

  return;
}

void DTPhase2Trigger::fill_dtsegments_variables(edm::Handle<DTRecSegment4DCollection> segments4D, const DTGeometry* dtGeom_) {
  idtsegments = 0;
  static TVectorF dummyfloat(1); dummyfloat(0) = -999.;
  for (DTRecSegment4DCollection::id_iterator chambIt = segments4D->id_begin(); chambIt != segments4D->id_end(); ++chambIt) {
    DTRecSegment4DCollection::range  range = segments4D->get(*chambIt);
    
    for (DTRecSegment4DCollection::const_iterator segment4D = range.first; segment4D != range.second; ++segment4D) {
      if (idtsegments >= dtsegmentsSize_) return;
      
      segm4D_wheel.push_back((*chambIt).wheel());
      segm4D_sector.push_back((*chambIt).sector());
      segm4D_station.push_back((*chambIt).station());
      
      const Bool_t hasPhi = segment4D->hasPhi();
      const Bool_t hasZed = segment4D->hasZed();
      segm4D_hasPhi.push_back(hasPhi);
      segm4D_hasZed.push_back(hasZed);
      
      segm4D_x_pos_loc.push_back(segment4D->localPosition().x());
      segm4D_y_pos_loc.push_back(segment4D->localPosition().y());
      segm4D_z_pos_loc.push_back(segment4D->localPosition().z());
      segm4D_x_dir_loc.push_back(segment4D->localDirection().x());
      segm4D_y_dir_loc.push_back(segment4D->localDirection().y());
      segm4D_z_dir_loc.push_back(segment4D->localDirection().z());
      
      if (hasPhi || hasZed) {
        TVectorF hitExpectedPos(12);
        TVectorF hitExpectedWire(12);
        std::vector<DTWireId> wireVector;
        
        for (Int_t kSL = 1; kSL < 4; kSL = kSL + 1) {
          if ((*chambIt).station() == 4 && kSL == 2) continue; //FRC 21-12-2016
          for (Int_t kL = 1; kL < 5; kL++) {
            wireVector.push_back(DTWireId((*chambIt).wheel(), (*chambIt).station(), (*chambIt).sector(), kSL, kL, 2));
          }
        }

        Int_t kkk = 0;
        const DTChamber* mych = dtGeom_->chamber(*chambIt);
        StraightLinePlaneCrossing segmentPlaneCrossing(((*mych).toGlobal(segment4D->localPosition())).basicVector(), ((*mych).toGlobal(segment4D->localDirection())).basicVector(), anyDirection);

        for (std::vector<DTWireId>::const_iterator wireIt = wireVector.begin(); wireIt != wireVector.end(); ++wireIt) {
          const DTLayer* layer = dtGeom_->layer(*wireIt);
          const DTChamber* chamber = dtGeom_->chamber(wireIt->layerId().chamberId());
          pair<Bool_t, Basic3DVector<Float_t> > ppt = segmentPlaneCrossing.position(layer->surface());
          Bool_t  success    = ppt.first; // check for failure
          Int_t   theExpWire = -999;
          Float_t theExpPos  = 999;

          if (success) {
            GlobalPoint segExrapolationToLayer(ppt.second);
            LocalPoint  segPosAtZWireLayer   = layer->toLocal(segExrapolationToLayer);
            LocalPoint  segPosAtZWireChamber = chamber->toLocal(segExrapolationToLayer);
            
            if ((kkk < 4 || kkk > 7) && hasPhi) {
              theExpPos  = segPosAtZWireChamber.x();
              theExpWire = layer->specificTopology().channel(segPosAtZWireLayer);
            }
            else if ((kkk >= 4 && kkk <= 7) && hasZed) {
              theExpPos = segPosAtZWireChamber.y();     //theExpPos = segPosAtZWire.x();
              //LocalPoint passPoint(-segPosAtZWire.y(),segPosAtZWire.x(),segPosAtZWire.z());
              theExpWire = layer->specificTopology().channel(segPosAtZWireLayer);
            }
          }
          hitExpectedWire[kkk] = theExpWire;
          hitExpectedPos[kkk]  = theExpPos;
          kkk++;
          if ((*chambIt).station() == 4 && kkk == 4) kkk += 4; //FRC 22-12-2016
        }
        
        new ((*segm4D_hitsExpPos)[idtsegments])  TVectorF(hitExpectedPos);
        new ((*segm4D_hitsExpWire)[idtsegments]) TVectorF(hitExpectedWire);
      }
      else {
        new ((*segm4D_hitsExpPos)[idtsegments])  TVectorF(dummyfloat);
        new ((*segm4D_hitsExpWire)[idtsegments]) TVectorF(dummyfloat);
      }
      
      const GeomDet*     geomDet   = theTrackingGeometry->idToDet(segment4D->geographicalId());
      const GlobalVector point_glb = geomDet->toGlobal(segment4D->localDirection());
      segm4D_cosx.push_back(point_glb.x());
      segm4D_cosy.push_back(point_glb.y());
      segm4D_cosz.push_back(point_glb.z());
      segm4D_phi.push_back(point_glb.phi());
      segm4D_theta.push_back(point_glb.theta());
      segm4D_eta.push_back(point_glb.eta());
      
      if (hasPhi) fill_dtphi_info(segment4D->phiSegment(), geomDet);
      else {
        segm4D_t0.push_back(-999.);
        segm4D_vdrift.push_back(-999.);
        segm4D_phinormchi2.push_back(-999.);
        segm4D_phinhits.push_back(-999);
        new ((*segm4D_phiHits_Pos)[idtsegments])        TVectorF(dummyfloat);
        new ((*segm4D_phiHits_PosCh)[idtsegments])      TVectorF(dummyfloat);
        new ((*segm4D_phiHits_PosErr)[idtsegments])     TVectorF(dummyfloat);
        new ((*segm4D_phiHits_Side)[idtsegments])       TVectorF(dummyfloat);
        new ((*segm4D_phiHits_Wire)[idtsegments])       TVectorF(dummyfloat);
        new ((*segm4D_phiHits_Layer)[idtsegments])      TVectorF(dummyfloat);
        new ((*segm4D_phiHits_SuperLayer)[idtsegments]) TVectorF(dummyfloat);
        new ((*segm4D_phiHits_Time)[idtsegments])       TVectorF(dummyfloat);
        new ((*segm4D_phiHits_TimeCali)[idtsegments])   TVectorF(dummyfloat);
      }
      
//       if (hasZed) fill_dtz_info(segment4D->zSegment(), geomDet);
//       else {
//         segm4D_znormchi2.push_back(-999.);
//         segm4D_znhits.push_back(-999);
//         new ((*segm4D_zHits_Pos)[idtsegments])      TVectorF(dummyfloat);
//         new ((*segm4D_zHits_PosCh)[idtsegments])    TVectorF(dummyfloat);
//         new ((*segm4D_zHits_PosErr)[idtsegments])   TVectorF(dummyfloat);
//         new ((*segm4D_zHits_Side)[idtsegments])     TVectorF(dummyfloat);
//         new ((*segm4D_zHits_Wire)[idtsegments])     TVectorF(dummyfloat);
//         new ((*segm4D_zHits_Layer)[idtsegments])    TVectorF(dummyfloat);
//         new ((*segm4D_zHits_Time)[idtsegments])     TVectorF(dummyfloat);
//         new ((*segm4D_zHits_TimeCali)[idtsegments]) TVectorF(dummyfloat);
//       }
      idtsegments++;
    }
  }
  return;
}


void DTPhase2Trigger::fill_dtphi_info(const DTChamberRecSegment2D* phiSeg, const GeomDet* chamb) {
  std::vector<DTRecHit1D> phirecHitslist = phiSeg->specificRecHits();
  
  //segment information
  segm4D_t0.push_back(phiSeg->t0());
  segm4D_vdrift.push_back(phiSeg->vDrift());
  segm4D_phinormchi2.push_back(phiSeg->chi2()/phiSeg->degreesOfFreedom());
  
  //rechits information
  const Int_t nphirecHits = phirecHitslist.size();
  segm4D_phinhits.push_back(nphirecHits);
  TVectorF phiPosRechits(nphirecHits);
  TVectorF phiPosChRechits(nphirecHits);
  TVectorF phiPosErrRechits(nphirecHits);
  TVectorF phiSideRechits(nphirecHits);
  TVectorF phiwireRechits(nphirecHits);
  TVectorF philayerRechits(nphirecHits);
  TVectorF phisuperlayerRechits(nphirecHits);
  TVectorF phiTimeRechits(nphirecHits);
  TVectorF phiTimeCaliRechits(nphirecHits);
  Int_t rechitscounter = 0;
  
  for (std::vector<DTRecHit1D>::const_iterator recHitsIt = phirecHitslist.begin(); recHitsIt!=phirecHitslist.end(); ++recHitsIt) {
    const GeomDet * layer = theTrackingGeometry->idToDet(recHitsIt->wireId().layerId());
    phiPosRechits(rechitscounter)        = recHitsIt->localPosition().x();
    phiPosChRechits(rechitscounter)      = chamb->toLocal(layer->toGlobal(recHitsIt->localPosition())).x();
    phiPosErrRechits(rechitscounter)     = recHitsIt->localPositionError().xx();
    phiSideRechits(rechitscounter)       = recHitsIt->lrSide();
    phiwireRechits(rechitscounter)       = recHitsIt->wireId().wire();
    philayerRechits(rechitscounter)      = recHitsIt->wireId().layerId().layer();
    phisuperlayerRechits(rechitscounter) = recHitsIt->wireId().layerId().superlayer();
    phiTimeRechits(rechitscounter)       = recHitsIt->digiTime();
    
    Float_t ttrig = theSync->offset(recHitsIt->wireId());
    phiTimeCaliRechits(rechitscounter)   = recHitsIt->digiTime() - ttrig;
    rechitscounter++;
  }
  
  new ((*segm4D_phiHits_Pos)[idtsegments])        TVectorF(phiPosRechits);
  new ((*segm4D_phiHits_PosCh)[idtsegments])      TVectorF(phiPosChRechits);
  new ((*segm4D_phiHits_PosErr)[idtsegments])     TVectorF(phiPosErrRechits);
  new ((*segm4D_phiHits_Side)[idtsegments])       TVectorF(phiSideRechits);
  new ((*segm4D_phiHits_Wire)[idtsegments])       TVectorF(phiwireRechits);
  new ((*segm4D_phiHits_Layer)[idtsegments])      TVectorF(philayerRechits);
  new ((*segm4D_phiHits_SuperLayer)[idtsegments]) TVectorF(phisuperlayerRechits);
  new ((*segm4D_phiHits_Time)[idtsegments])       TVectorF(phiTimeRechits);
  new ((*segm4D_phiHits_TimeCali)[idtsegments])   TVectorF(phiTimeCaliRechits);
  return;
}


// void DTPhase2Trigger::ttreegeneratorinit(const edm::ParameterSet& pset):
//       rpcToken_(consumes<MuonDigiCollection<RPCDetId,RPCDigi> > (pset.getParameter<edm::InputTag>("rpcLabel"))),
//       UnpackingRpcRecHitToken_(consumes<RPCRecHitCollection> (pset.getParameter<edm::InputTag>("UnpackingRpcRecHitLabel"))) {
//   // get the tTrigDBInfo
//   theSync = DTTTrigSyncFactory::get()->create(pset.getUntrackedParameter<std::string>("tTrigMode"),
//                                 pset.getUntrackedParameter<edm::ParameterSet>("tTrigModeConfig"));
// 
// 
//   // get run configuration options
//   runOnRaw_          = pset.getParameter<bool>("runOnRaw");
//   runOnSimulation_   = pset.getParameter<bool>("runOnSimulation");
//   runOnDigiSimLinks_ = pset.getParameter<bool>("runOnDigiSimLinks");
// 
//   //get parameters from the configuration file
//   //names of the different event collections
//   dtDigiLabel_     = pset.getParameter<edm::InputTag>("dtDigiLabel");
// 
//   dtDigiToken_     = consumes<DTDigiCollection>(edm::InputTag(dtDigiLabel_));
//   dtDigiTokenSim_  = consumes<DTDigiSimLinkCollection>(edm::InputTag(dtDigiLabel_));
// 
//   dtSegmentLabel_  = pset.getParameter<edm::InputTag>("dtSegmentLabel");
//   dtSegmentToken_  = consumes<DTRecSegment4DCollection>(edm::InputTag(dtSegmentLabel_));
// 
//   cscSegmentLabel_ = pset.getParameter<edm::InputTag>("cscSegmentLabel");
//   cscSegmentToken_ = consumes<CSCSegmentCollection>(edm::InputTag(cscSegmentLabel_));
// 
//   dtTrigTwinMuxInLabel_    = pset.getParameter<edm::InputTag>("dtTrigTwinMuxInLabel");
//   dtTrigTwinMuxThLabel_    = pset.getParameter<edm::InputTag>("dtTrigTwinMuxThLabel");
//   dtTrigTwinMuxInToken_    = consumes<L1MuDTChambPhContainer>(edm::InputTag(dtTrigTwinMuxInLabel_)); 
//   dtTrigTwinMux_ThToken_ = consumes<L1MuDTChambThContainer>(edm::InputTag(dtTrigTwinMuxThLabel_));
// 
//   dtTrigTwinMuxOutLabel_  = pset.getParameter<edm::InputTag>("dtTrigTwinMuxOutLabel");
//   dtTrigTwinMuxOutToken_  = consumes<L1MuDTChambPhContainer>(edm::InputTag(dtTrigTwinMuxOutLabel_));
// 
//   staMuLabel_      = pset.getParameter<edm::InputTag>("staMuLabel");
//   staMuToken_      = consumes<reco::MuonCollection>(edm::InputTag(staMuLabel_));
// 
//   gmtLabel_        = pset.getParameter<edm::InputTag>("gmtLabel");
//   gmtToken_        = consumes<l1t::MuonBxCollection>(edm::InputTag(gmtLabel_));
// 
//   triggerTag_      = pset.getParameter<edm::InputTag>("TriggerTag");
//   triggerToken_    = consumes<edm::TriggerResults>(edm::InputTag(triggerTag_));
// 
//   triggerEventTag_   = pset.getParameter<edm::InputTag>("TriggerEventTag");
//   triggerEventToken_ = consumes<trigger::TriggerEvent>(edm::InputTag(triggerEventTag_));
//   trigFilterNames_    = pset.getParameter<std::vector<std::string>>("trigFilterNames");
// 
//   gtLabel_         = pset.getParameter<edm::InputTag>("gtLabel"); // legacy
//   gtToken_         = consumes<L1GlobalTriggerReadoutRecord>(edm::InputTag(gtLabel_)); //legacy
// 
//   rpcRecHitLabel_  = pset.getParameter<edm::InputTag>("rpcRecHitLabel");
//   rpcRecHitToken_  = consumes<RPCRecHitCollection>(edm::InputTag(rpcRecHitLabel_));
// 
//   //max size of the different saved objects (per event)
//   digisSize_       = pset.getParameter<int>("dtDigiSize");
//   dtsegmentsSize_  = pset.getParameter<int>("dtSegmentSize");
//   cscsegmentsSize_ = pset.getParameter<int>("cscSegmentSize");
//   dtltTwinMuxOutSize_     = pset.getParameter<int>("dtTrigTwinMuxOutSize");
//   dtltTwinMuxInSize_ = pset.getParameter<int>("dtTrigTwinMuxInSize");
//   dtltTwinMuxThSize_ = pset.getParameter<int>("dtTrigTwinMuxThSize");
//   gmtSize_         = pset.getParameter<int>("gmtSize"); 
//   recoMuSize_       = pset.getParameter<int>("recoMuSize");
//   rpcRecHitSize_   = pset.getParameter<int>("rpcRecHitSize"); 
// 
//   PrimaryVertexTag_ = pset.getParameter<edm::InputTag>("PrimaryVertexTag");
//   PrimaryVertexToken_ =consumes<reco::VertexCollection>(edm::InputTag(PrimaryVertexTag_));
// 
//   beamSpotTag_      = pset.getParameter<edm::InputTag>("beamSpotTag");
//   beamSpotToken_    = consumes<reco::BeamSpot>(edm::InputTag(beamSpotTag_));
// 
//   scalersSource_    = pset.getParameter<edm::InputTag>("scalersResults");
//   scalersSourceToken_ = consumes<LumiScalersCollection>(edm::InputTag(scalersSource_));
// 
//   puSummaryTag_      = pset.getParameter<edm::InputTag>("puSummaryTag");
//   puSummaryToken_    = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag(puSummaryTag_));
//        
//   // CB comment for now to avoid crash
//   // lumiInputTag_      = pset.getParameter<edm::InputTag>("lumiInputTag");
//   // lumiProducerToken_ = consumes<LumiDetails, edm::InLumi>(lumiInputTag_);
//   
//   bmtfPhInputTag_ = consumes<L1MuDTChambPhContainer>(pset.getParameter<edm::InputTag>("bmtfInputPhDigis"));
//   bmtfThInputTag_ = consumes<L1MuDTChambThContainer>(pset.getParameter<edm::InputTag>("bmtfInputThDigis"));
//   bmtfOutputTag_ = consumes<l1t::RegionalMuonCandBxCollection>(pset.getParameter<edm::InputTag>("bmtfOutputDigis"));
// 
//   OnlyBarrel_ = pset.getParameter<bool>("OnlyBarrel");
// 
//   localDTmuons_    = pset.getUntrackedParameter<bool>("localDTmuons",false);
// 
//   AnaTrackGlobalMu_= pset.getUntrackedParameter<bool>("AnaTrackGlobalMu",true);  // set to false when problems with tracks of global muons
// 
// 
//   runLegacy_gmt_ = pset.getUntrackedParameter<bool>("runLegacy_gmt",true); // Collection not available in cosmics 2017, at least on september runs
//                                                                            // set to false in cosmics if needed
// 
//   outFile_         = pset.getParameter<std::string>("outputFile");
// 
//   initialize_Tree_variables();
// 
//   //counters
//   idigis       = 0;
//   idtsegments  = 0;
//   icscsegments = 0;
//   idtltTwinMuxOut     = 0;
//   idtltTwinMuxIn     = 0;
//   idtltTwinMux_th  = 0;
//   imuons       = 0;
//   igmt         = 0;
//   igtalgo      = 0; // legacy
//   igttt        = 0; // legacy
//   ihlt         = 0;
// }


// ================================
// ================================ Methods from other classes that are not the analyser ==============================
// ================================
TransformHelper::TransformHelper() {};
TransformHelper::~TransformHelper() {};


Float_t TransformHelper::getHitPosition(Int_t* sl, Int_t* l) {
  Int_t index = (*sl - 1) * 4 + (*l - 1); // Index for a given SL and layer
  return wirepos_z[index];
}


std::pair<Float_t, Float_t> TransformHelper::getLRHitPosition(Int_t* sl, Int_t* l, Float_t* pos, Int_t* lr) {
  Int_t index = (*sl - 1) * 4 + (*l - 1); // Index for a given SL and layer
  Float_t relativeChamberPosition = (*pos - wirepos_x[index]) + chamberLength + 1.35;
  Float_t relativeCellPosition    = fmod(relativeChamberPosition, cellLength);
  Float_t positionInCell          = relativeCellPosition - cellSemiLength;
  return {*pos, *pos + TMath::Power(-1, *lr) * 2 * abs(positionInCell)};
}


std::pair<Float_t, Float_t> TransformHelper::getWirePosition(Int_t sl, Int_t l, Int_t wire) {
  Int_t index = (sl - 1) * 4 + (l - 1); // Index for a given SL and layer
  Float_t x = wirepos_x[index] + wire * cellLength;
  return {x, wirepos_z[index]};
}


std::vector<Float_t> TransformHelper::getDigiPosition(Int_t sl, Int_t l, Int_t wire, Float_t dtime) {
  std::vector<Float_t> res;
  res.clear();
  Float_t posInCell = dtime * vdrift;
  std::pair<Float_t, Float_t> p = getWirePosition(sl, l, wire);
  res.push_back(std::get<0>(p) - posInCell);
  res.push_back(std::get<0>(p) + posInCell);
  res.push_back(std::get<1>(p));
  return res;
}


TH2D* TransformHelper::makeHoughTransform(TGraph* occupancy) {
  TH2D *linespace = new TH2D("linespace", "linespace", binsangle, 0, TMath::TwoPi(), binsrho, -250., 250.);
  Double_t x = 0, z = 0, rho = 0, phi = 0;
  
  // loop over the occupancy (left-right)
  for (Int_t i = 0; i < occupancy->GetN(); i++) {
    occupancy->GetPoint(i, x, z);
    phi = 0.; rho = 0.;
    while (phi < TMath::TwoPi()) {
        rho = x * TMath::Cos(phi) + z * TMath::Sin(phi);
        linespace->Fill(phi, rho);
        phi = phi + 0.001;
    }
  }
  return linespace;
}


// std::vector<std::pair<Float_t, Float_t>> TransformHelper::findLocalMaxima(TH2D* histo2D) {
//   std::vector<std::pair<Float_t, Float_t>> res;
//   res.clear();
//   TSpectrum2* spc;
//   spc = new TSpectrum2();
//   Int_t xtemp, ytemp, ztemp, npks = 0;
//   Float_t theta_pos = 0., rho_pos = 0, m = 0, n = 0;
//   
//   spc->Search(histo2D, 5, "noMarkov", 0.55);
//   cout << histo2D->ShowPeaks(1, "noMarkov", 0.8) << endl;
//   
//   npks = spc->GetNPeaks();
//   if (npks == 0) throw std::runtime_error("We have not found any peak.");
//   else           cout << "\n> Number of found peaks: " << npks << endl;
//   
//   Double_t* posx; Double_t* posy;
//   posx = spc->GetPositionX(); posy = spc->GetPositionY();
//   for (Int_t i = 0; i < npks; i++) {
//     xtemp = 0; ytemp = 0; ztemp = 0; 
//     Double_t zero = 0.;
//     histo2D->GetBinXYZ(histo2D->FindBin(posx[i], posy[i], zero), xtemp, ytemp, ztemp);
//     theta_pos = histo2D->GetXaxis()->GetBinCenter(xtemp);
//     rho_pos   = histo2D->GetYaxis()->GetBinCenter(ytemp);
//     
//     if ((abs(theta_pos - TMath::Pi()/2) < 0.2) || (abs(theta_pos - TMath::Pi()*3/2) < 0.2)) continue;
//     
//     m = -1. / TMath::Tan(theta_pos);
//     n = rho_pos / TMath::Sin(theta_pos);
//     res.push_back({m, n});
//   }
//   cout << "> Number of accepted peaks: " << res.size() << endl;
//   delete spc;
//   return res;
// }


std::pair<Float_t, Float_t> TransformHelper::findBestSegment(TH2D* linespace) {
  // This function goes back to point-space (for now only one)
  Int_t locx = 0;
  Int_t locy = 0;
  Int_t locz = 0;
  linespace->GetMaximumBin(locx, locy, locz);
  
//   Double_t theta_err = linespace->GetXaxis()->GetBinLowEdge(locx);
//   Double_t rho_err   = linespace->GetYaxis()->GetBinLowEdge(locy);
  Double_t theta_pos = linespace->GetXaxis()->GetBinCenter(locx);
  Double_t rho_pos   = linespace->GetYaxis()->GetBinCenter(locy);
  
  Float_t m = -1./TMath::Tan(theta_pos);
  Float_t n = rho_pos/TMath::Sin(theta_pos);
  cout << "\n> Maximum bin in the dual space. x: " << locx << " y: " << locy << " z: " << locz << endl;
  cout << "> Maximum bin in the dual space in terms of rho and theta. rho: " << rho_pos << " theta: " << theta_pos << endl;
  cout << "> Associated line in the direct space. m: " << m << " n: " << n << "\n" << endl;
  return {m, n};
}


void TransformHelper::print_canvas(TCanvas* canvas, TString output_name_without_ext) {
//   if (! fs::exists(path.Data())) {
//     system("mkdir -p " + path);
//   }
  canvas->Print(path + "/" + output_name_without_ext + ".png");
  canvas->Print(path + "/" + output_name_without_ext + ".pdf");
  canvas->Print(path + "/" + output_name_without_ext + ".root");
  return;
}


void TransformHelper::initialise(UInt_t neventsmax, Int_t mb, Int_t wh, Int_t se) {
  cout << "> Initialising things..." << endl;
  if ((mb <= 0) || (mb > 4))  chamber = -1;
  else                        chamber = mb;
  if ((wh < -2) || (wh > +2)) wheel   = -3;
  else                        wheel   = wh;
  if ((se <= 0) || (se > 14)) sector  = -1;
  else                        sector  = se;
  
  txtwh = ""; txtmb = ""; txtse = "";
  
  if (wheel == -3) txtwh = "ALL";
  else             txtwh += wheel;
  
  if (chamber == -1) txtmb = "ALL";
  else               txtmb += chamber;
  
  if (sector == -1) txtse = "ALL";
  else              txtse += sector;
  
  cout << "> We're going to process for wheel " << txtwh + "," << " sector " << txtse + "," << " and chamber " << txtmb + "." << endl;
  
  occupancy   = new TGraph();  // Geometrical plot: it will be the direct (or point) space plot.
  actualocc   = new TGraph();  // Geometrical plot: it will be the direct (or point) space plot.
  linespace   = new TH2D();    // Dual (or line) space plot.
  linespace1D = new TH1D();
  
  nevents  = 0;
  n_events_limit = neventsmax;
  nhitsmax = 8;
  filled   = false;
  path     = "/nfs/fanae/user/vrbouza/www/Proyectos/trigger_primitives/results/houghtrans/";
  
}


void TransformHelper::run(std::vector<Short_t> dtsegm4D_wheel, std::vector<Short_t> dtsegm4D_sector, std::vector<Short_t> dtsegm4D_station, TClonesArray *dtsegm4D_phi_hitsSuperLayer, TClonesArray *dtsegm4D_phi_hitsLayer, TClonesArray *dtsegm4D_phi_hitsPos, TClonesArray *dtsegm4D_phi_hitsSide, TClonesArray *dtsegm4D_phi_hitsWire, std::vector<Short_t> dtsegm4D_phinhits) {
  // We produce an artificial event taking segment info from various chambers and/or sectors and/or wheels until 20 hits are collected.
  if (filled) return;
  
  nevents += 1;
  if (nevents >= n_events_limit) return;
  if (fmod((nevents + 1), 1000) == 0) cout << "Event: " << nevents + 1 << endl;
  
  cout << nevents << endl;
  
  nhits = 0; nsegs = 0; actualhits = 0;
  Float_t z = 0, x1 = 0, x2 = 0;
  std::pair<Float_t, Float_t> p = {0, 0};
  
  for (Int_t idtsegments = 0; idtsegments < (Int_t)dtsegm4D_wheel.size(); idtsegments++) {                // i := segment index
    if ((dtsegm4D_wheel.at(idtsegments)   != wheel) &&   (wheel   != -3)) continue;
    if ((dtsegm4D_station.at(idtsegments) != chamber) && (chamber != -1)) continue;
    if ((dtsegm4D_sector.at(idtsegments)  != sector) &&  (sector  != -1)) continue;
    for (Int_t idthits = 0; idthits < dtsegm4D_phinhits.at(idtsegments); idthits++) {
      if ( ((TVectorF*)dtsegm4D_phi_hitsSuperLayer->At(idtsegments))[idthits] == 2) continue;

      z = getHitPosition( (Int_t*)(&(((TVectorF*)dtsegm4D_phi_hitsSuperLayer->At(idtsegments))[idthits])), (Int_t*)(&(((TVectorF*)dtsegm4D_phi_hitsLayer->At(idtsegments))[idthits])) );

      p = getLRHitPosition( (Int_t*)(&(((TVectorF*)dtsegm4D_phi_hitsSuperLayer->At(idtsegments))[idthits])),
                            (Int_t*)(&(((TVectorF*)dtsegm4D_phi_hitsLayer->At(idtsegments))[idthits])),
                            (Float_t*)(&(((TVectorF*)dtsegm4D_phi_hitsPos->At(idtsegments))[idthits])),
                            (Int_t*)(&(((TVectorF*)dtsegm4D_phi_hitsSide->At(idtsegments))[idthits])));
      x1 = p.first; x2 = p.second;
      
      occupancy->SetPoint(nhits,     x1, z);
      actualocc->SetPoint(actualhits,x1, z);
      occupancy->SetPoint(nhits + 1, x2, z);
      nhits += 2; actualhits += 1; nsegs += 1;
    }
  }
  cout << nhits << endl;
  
  // if enough
  if (nhits > nhitsmax) filled = true;
  else {
    for (UInt_t i = 0; i < nhits; i++)      occupancy->RemovePoint(i);
    for (UInt_t i = 0; i < actualhits; i++) actualocc->RemovePoint(i);
  }
}


void TransformHelper::finish() {
  if (! filled) cerr << "> Not filled properly!" << endl;
  cout << "> Number of hits plotted: " << nhits << endl;
  cout << "> Number of segments selected: " << nsegs << endl;
  
  // Now draw and save TGraph
  TCanvas* c1 = new TCanvas("c1", "Occupancy", 700, 700);
  c1->SetLeftMargin(0.15);
  
  occupancy->SetTitle("x-z hits");
  occupancy->SetMarkerColor(kBlack);
  occupancy->SetMinimum(-2);
  occupancy->SetMaximum(30);
  occupancy->SetMarkerStyle(5);
  occupancy->SetMarkerSize(1);
  occupancy->GetXaxis()->SetRange(0, 250);
  occupancy->GetXaxis()->SetLimits(-140, 80);
  occupancy->GetXaxis()->SetTitle("x coordinate (in cm)");
  occupancy->GetYaxis()->SetTitle("z coordinate (in cm)");
  occupancy->Draw("AP");
  print_canvas(c1, "Occupancy_MB" + txtmb + "_Wh" + txtwh + "_S" + txtse);
  delete c1;
  
  // now make HT
  TH2D* linespace; std::vector<std::pair<Float_t, Float_t>> maxima;
  linespace = makeHoughTransform(occupancy);
//   if (doVarious) maxima = findLocalMaxima(linespace);
  TCanvas* c2 = new TCanvas("c2", "HT", 700, 700);
  c2->SetLeftMargin(0.15);
  
  linespace->SetTitle("Linespace");
  linespace->GetXaxis()->SetTitle("#theta (rad)");
  linespace->GetYaxis()->SetTitle("#rho (cm)");
  linespace->Draw("COLZ");
  print_canvas(c2, "HoughTransform_MB" + txtmb + "_Wh" + txtwh + "_S" + txtse);
  delete c2;
  
  TCanvas* c3 = new TCanvas("c1", "Occupancy", 700, 700);
  c3->SetLeftMargin(0.15);
  occupancy->SetTitle("x-z hits");
  occupancy->SetMarkerColor(kBlack);
  occupancy->SetMinimum(-2);
  occupancy->SetMaximum(30);
  occupancy->SetMarkerStyle(5);
  occupancy->SetMarkerSize(1);
  occupancy->GetXaxis()->SetRange(0, 250);
  occupancy->GetXaxis()->SetLimits(-140, 80);
  occupancy->GetXaxis()->SetTitle("x coordinate (in cm)");
  occupancy->GetYaxis()->SetTitle("z coordinate (in cm)");
  occupancy->Draw("AP");
  
  actualocc->SetMarkerColor(kAzure);
  actualocc->SetMarkerStyle(5);
  actualocc->SetMarkerSize(1);
  actualocc->Draw("P");
  
  // Now try to find segments
  if (! doVarious) {
    Float_t m = 0, n = 0;
    std::pair<Float_t, Float_t> p = {0, 0};
    p = findBestSegment(linespace);
    m = p.first; n = p.second;
    TF1* segm = new TF1("seg", "[0]*x+[1]", -150, 150);
    segm->SetParameters(m, n);
    segm->SetLineColor(kRed);
    segm->SetLineWidth(2);
    segm->Draw("SAME");
    delete segm;
  }
  else {
    std::vector<TF1*> hmax;
    hmax.clear();
    for (UInt_t i = 0; i < maxima.size(); i++) {
      TString asdf = "segmax_";
      asdf += i;
      TF1* tmpsegm = new TF1(asdf, "[0]*x+[1]", -150, 150);
      cout << maxima.at(i).first << maxima.at(i).second << endl;
      tmpsegm->SetParameters(maxima.at(i).first, maxima.at(i).second);
      tmpsegm->SetLineColor(kRed);
      tmpsegm->SetLineWidth(2);
      hmax.push_back((TF1*)tmpsegm->Clone(asdf));
      delete tmpsegm;
      hmax.back()->Draw("SAME");
    }
  }
  print_canvas(c3, "Occupancy_MB" + txtmb + "_Wh" + txtwh + "_S" + txtse);
  delete c3;
  return;
}

void TransformHelper::getvalues(const edm::Event& event) {
//   event.getByLabel("dtsegm4D_wheel", dtsegm4D_wheel);
//   event.getByLabel("dtsegm4D_sector", dtsegm4D_sector);
//   event.getByLabel("dtsegm4D_station", dtsegm4D_station);
//   event.getByLabel("dtsegm4D_phinhits", dtsegm4D_phinhits);
//   event.getByLabel("dtsegm4D_phi_hitsSuperLayer", dtsegm4D_phi_hitsSuperLayer);
//   event.getByLabel("dtsegm4D_phiHits_Layer", dtsegm4D_phiHits_Layer);
//   event.getByLabel("dtsegm4D_phi_hitsPos", dtsegm4D_phi_hitsPos);
//   event.getByLabel("dtsegm4D_phi_hitsSide", dtsegm4D_phi_hitsSide);
//   event.getByLabel("dtsegm4D_phi_hitsWire", dtsegm4D_phi_hitsWire);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DTPhase2Trigger);
