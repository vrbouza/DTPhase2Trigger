#ifndef DTPHASE2TRIGGER_H
#define DTPHASE2TRIGGER_H
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


// System / std headers
#include <memory>
// #include <experimental/filesystem>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <vector>

// CMSSW headers
#include <FWCore/Framework/interface/ConsumesCollector.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/MuonData/interface/MuonDigiCollection.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DTDigi/interface/DTDigi.h"
#include "DataFormats/DTDigi/interface/DTDigiCollection.h"
#include "DataFormats/DTDigi/interface/DTLocalTriggerCollection.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h" // New trying to avoid crashes in the topology functions
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/GeomPropagators/interface/StraightLinePlaneCrossing.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "SimDataFormats/DigiSimLinks/interface/DTDigiSimLinkCollection.h"// To cope with Digis from simulation
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/DTObjects/interface/DTMtime.h"
#include "CondFormats/DataRecord/interface/DTMtimeRcd.h"

#include "CalibMuon/DTDigiSync/interface/DTTTrigBaseSync.h"
#include "CalibMuon/DTDigiSync/interface/DTTTrigSyncFactory.h"
// #include "CalibMuon/DTCalibration/plugins/DTCalibrationMap.h"

#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"  // New trying to avoid crashes in the topology functions
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h" // New trying to avoid crashes in the topology functions
#include "Geometry/DTGeometry/interface/DTTopology.h" // New trying to avoid crashes in the topology functions



// ROOT headers
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TVectorF.h"
#include "TClonesArray.h"
#include "TGaxis.h"


// Other headers
#include "UserCode/DTPhase2Trigger/interface/DefineTreeVariables.h"

// Namespaces
using namespace std;
using namespace edm;
using namespace cms;
// namespace fs = std::experimental::filesystem;

// Enums or/and const declarations & definitions
const Float_t wirepos_x[] = {-95.0142, -97.1103, -95.0094, -97.1046, 113.4, 115.5, 113.4, 115.5, -97.0756, -99.1991, -97.0935, -99.1933};
const Float_t wirepos_z[] = {0., 1.3, 2.6, 3.9, 99999., 99999., 99999., 99999., 23.8, 25.1, 26.4, 27.7};  // in cm

const Float_t cellLength  = 4.2;  // The length of a cell
const Float_t cellSemiLength = cellLength/2;  // length from the wire to the end of the cell
const Float_t cellQuarterLength = cellSemiLength/2;  // length from the wire to the middle of space between the border and the wire
const Float_t chamberLength  = 210.;  // The length of the whole chamber (MB1)

const Float_t vdrift = 0.00545;  // in cm/ns
// const Float_t vdrift = 6.5/386.74;  // in cm/ns ¿?¿??, de lo de Camilo

const Int_t   timeOffset = 0;
// const Int_t   flat_calib = 325;
const Int_t   flat_calib = 0;


// Other classes declarations
class TransformHelper {
  public:
    // Constructor(s) and destructor
    TransformHelper();
    ~TransformHelper();
    
    
    // Methods
    void initialise(UInt_t neventsmax = 20000, Int_t mb = 2, Int_t wh = 0, Int_t se = -1, Bool_t db = false);
    void run(std::vector<Short_t> dtsegm4D_wheel, std::vector<Short_t> dtsegm4D_sector, std::vector<Short_t> dtsegm4D_station, TClonesArray *dtsegm4D_phi_hitsSuperLayer, TClonesArray *dtsegm4D_phi_hitsLayer, TClonesArray *dtsegm4D_phi_hitsPos, TClonesArray *dtsegm4D_phi_hitsSide, TClonesArray *dtsegm4D_phi_hitsWire, std::vector<Short_t> dtsegm4D_phinhits);
    void run(std::vector<Short_t> digi_wheel, std::vector<Short_t> digi_sector, std::vector<Short_t> digi_station, std::vector<Short_t> digi_sl, std::vector<Short_t> digi_layer, std::vector<Short_t> digi_wire, std::vector<Float_t> digi_time);
    void run_one(edm::Handle<DTDigiSimLinkCollection> dtdigisSim, const DTGeometry* dtGeom_, const edm::Event& iEv, DTTTrigBaseSync *theSync, const DTMtime* mTimeMap);
    void finish();
    void finish_one();
    
    
    // Class' variables & objects
    TString path;
    
    UInt_t  binsangle;
    UInt_t  binsrho;
    UInt_t  nhits, nsegs, actualhits, nseldigis;
    UInt_t  nhitsmax;
    
    Bool_t  doVarious;
    
    Int_t  chamber, wheel, sector;
    TString txtwh, txtmb, txtse;
//     edm::Handle<std::vector<std::vector<Int_t>>>* dtsegm4D_wheel;
//     edm::Handle<std::vector<std::vector<Int_t>>>* dtsegm4D_sector;
//     edm::Handle<std::vector<std::vector<Int_t>>>* dtsegm4D_station;
//     edm::Handle<std::vector<std::vector<Int_t>>>* dtsegm4D_phinhits;
//     edm::Handle<std::vector<TClonesArray>>* dtsegm4D_phi_hitsSuperLayer;
//     edm::Handle<std::vector<TClonesArray>>* dtsegm4D_phiHits_Layer;
//     edm::Handle<std::vector<TClonesArray>>* dtsegm4D_phi_hitsPos;
//     edm::Handle<std::vector<TClonesArray>>* dtsegm4D_phi_hitsSide;
//     edm::Handle<std::vector<TClonesArray>>* dtsegm4D_phi_hitsWire;
    
  private:
    // Methods
    Float_t getHitPosition(Int_t sl, Int_t l);
    std::pair<Float_t, Float_t> getLRHitPosition(Int_t sl, Int_t l, Float_t pos, Int_t lr);
    std::pair<Float_t, Float_t> getWirePosition(Int_t sl, Int_t l, Int_t wire);
    std::vector<Float_t>        getDigiPosition(Int_t sl, Int_t l, Int_t wire, Float_t dtime);
    std::pair<Float_t, Float_t> getDigiPosition(Float_t xpos, Float_t dtime);
    std::pair<Float_t, Float_t> getDigiPosition(Float_t xpos, Float_t dtime, Float_t vd);
    TH2D* makeHoughTransform(TGraph* occupancy);
//     std::vector<std::pair<Float_t, Float_t>>   findLocalMaxima(TH2D* histo2D);
    std::pair<Float_t, Float_t> findBestSegment(TH2D* linespace);
    
    Float_t getPhase2Time(const edm::Event& iEv, Int_t w, Int_t c, Int_t s, Int_t sl, Int_t l, Int_t wi, Float_t digiTime, DTTTrigBaseSync *theSync);
    
    void DrawGraphWithReversedYAxis(TGraph *g);
    void print_canvas(TCanvas* canvas, TString output_name_without_ext);
    void getvalues(const edm::Event& event);
    
    // Class' variables & objects
    UInt_t nevents, n_events_limit;
    
    Float_t xlowlim, xhighlim, zlowlim, zhighlim;
    
    Bool_t filled;
    Bool_t debug;
    Bool_t doDigis;
    
    TGraph* occupancy; TGraph* actualocc;
    
    TH2D* linespace;
    
    TH1D* linespace1D;
    
    TF1* tmpsegm;
};


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class DTPhase2Trigger : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    // Essential methods
    explicit DTPhase2Trigger(const edm::ParameterSet&);
    ~DTPhase2Trigger();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    // Other methods
    
    
    // Other definitions
    UInt_t n_events_limit = 20000;
    Bool_t doFanae        = false;
    Bool_t doDigis        = true;
    Bool_t doOne          = true;

    TransformHelper* hlpr;
    

  private:
    // Essential methods
    virtual void beginJob() override;
    virtual void analyze(const edm::Event& event, const edm::EventSetup& context) override;
    virtual void endJob() override;
    
    void fill_dtsegments_variables(edm::Handle<DTRecSegment4DCollection> segments4D, const DTGeometry* dtGeom_);
    void fill_digi_variables(edm::Handle<DTDigiCollection> dtdigis);
    void fill_digi_variablesSim(edm::Handle<DTDigiSimLinkCollection> dtdigisSim);
    void fill_dtphi_info(const DTChamberRecSegment2D* phiSeg, const GeomDet* chamb);
    void fill_dtz_info(const DTSLRecSegment2D* zSeg, const GeomDet* geomDet);
    void clear_Arrays();
    void initialize_Tree_variables();
//     void ttreegeneratorinit(const edm::ParameterSet& pset);
    Float_t getPhase2Time(const edm::Event& iEv, Int_t w, Int_t c, Int_t s, Int_t sl, Int_t l, Int_t wi, Float_t digiTime);
    
    
    DTTTrigBaseSync *theSync;
    
    // Declarations from dttree production .h (not DefineTreeVariables.h)   ========================================
    edm::InputTag dtDigiLabel_;
    edm::EDGetTokenT<DTDigiCollection> dtDigiToken_ ;
    edm::EDGetTokenT<DTDigiSimLinkCollection> dtDigiTokenSim_ ;
    edm::InputTag dtSegmentLabel_;
    edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentToken_;
    edm::InputTag cscSegmentLabel_;
    edm::EDGetTokenT<CSCSegmentCollection> cscSegmentToken_;
    edm::InputTag dtTrigTwinMuxOutLabel_;
    edm::EDGetTokenT<L1MuDTChambPhContainer> dtTrigTwinMuxOutToken_ ;
    edm::InputTag dtTrigTwinMuxInLabel_;
    edm::InputTag dtTrigTwinMuxThLabel_;
    edm::EDGetTokenT<L1MuDTChambPhContainer> dtTrigTwinMuxInToken_ ;
    edm::EDGetTokenT<L1MuDTChambThContainer> dtTrigTwinMux_ThToken_ ;
    edm::InputTag staMuLabel_;
    edm::EDGetTokenT<reco::MuonCollection> staMuToken_;
    edm::InputTag gmtLabel_;
    edm::EDGetTokenT<l1t::MuonBxCollection> gmtToken_;
    edm::InputTag gtLabel_; // legacy
    edm::EDGetTokenT<L1GlobalTriggerReadoutRecord> gtToken_; // legacy
    edm::InputTag rpcRecHitLabel_;
    edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitToken_;
    edm::InputTag PrimaryVertexTag_;
    edm::EDGetTokenT<reco::VertexCollection> PrimaryVertexToken_ ;
    edm::InputTag beamSpotTag_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    edm::InputTag puSummaryTag_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puSummaryToken_;
    edm::InputTag scalersSource_;
    edm::EDGetTokenT<LumiScalersCollection> scalersSourceToken_;
    edm::InputTag triggerTag_;
    edm::EDGetTokenT<edm::TriggerResults> triggerToken_ ;
    edm::InputTag triggerEventTag_;
    edm::EDGetTokenT<trigger::TriggerEvent> triggerEventToken_ ;
    
    // CB comment for now to avoid crash
    // edm::InputTag lumiInputTag_;
    // edm::EDGetTokenT<LumiDetails> lumiProducerToken_ ;
    // edm::EDGetTokenT<LumiSummary> lumiSummaryToken_;
    edm::EDGetTokenT<L1MuDTChambPhContainer> bmtfPhInputTag_;
    edm::EDGetTokenT<L1MuDTChambThContainer> bmtfThInputTag_;
    edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection> bmtfOutputTag_; 
    edm::EDGetTokenT<MuonDigiCollection<RPCDetId,RPCDigi> > rpcToken_;   
    edm::InputTag UnpackingRpcRecHitLabel_;
    edm::EDGetTokenT<RPCRecHitCollection> UnpackingRpcRecHitToken_;
   
    Bool_t OnlyBarrel_;

    Bool_t debug;

    Bool_t runOnRaw_;
    Bool_t runOnSimulation_;
    Bool_t runOnDigiSimLinks_; // To use to read simulation digi sim links

    Bool_t localDTmuons_;
    Bool_t AnaTrackGlobalMu_;  // To avoid look to the global tracks (The muon collection: vector<reco::Muon> exit,  but not the global tracks:  vector<reco::Track> )
    Bool_t runLegacy_gmt_;
    std::string outFile_;

    std::vector<std::string> trigFilterNames_;

    const DTMtime* mTimeMap;
    edm::ESHandle<MagneticField> theBField;
    edm::ESHandle<Propagator> propagatorAlong;
    edm::ESHandle<Propagator> propagatorOpposite;
    edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;

    Int_t digisSize_;
    Int_t dtsegmentsSize_;
    Int_t cscsegmentsSize_;
    Int_t dtltTwinMuxOutSize_;
    Int_t dtltTwinMuxInSize_;
    Int_t dtltTwinMuxThSize_;
    Int_t gmtSize_;
    Int_t recoMuSize_;
    Int_t rpcRecHitSize_;

    //counters
    Short_t idigis;
    Short_t idtsegments;
    Short_t icscsegments;
    Short_t idtltTwinMuxOut;
    Short_t idtltTwinMuxIn;
    Short_t idtltTwinMux_th;
    Short_t imuons;
    Short_t igmt;
    Short_t igtalgo; // legacy
    Short_t igttt; // legacy
    Short_t ihlt;
    Short_t irpcrechits;
    Short_t irpcdigi_TwinMux;
    Short_t irpcrechits_TwinMux;
    Short_t bmtf_size;
    
    reco::BeamSpot beamspot;
    // End of declarations from dttree production   ========================================
    
};


#endif
