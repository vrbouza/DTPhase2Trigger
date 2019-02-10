//This is the list and types of the variables saved in the TTree;
//New variables must be declared here

#include "TString.h"

//event variables
Int_t runnumber;
Int_t lumiblock;
Long_t eventNumber;
ULong64_t timestamp;
Int_t bunchXing;
Long_t orbitNum;

Float_t true_pileup;
Int_t actual_pileup;

//primary vertex
Float_t PV_x;
Float_t PV_y;
Float_t PV_z;

Float_t PV_xxE;
Float_t PV_yyE;
Float_t PV_zzE;
Float_t PV_xyE;
Float_t PV_xzE;
Float_t PV_yzE;

Float_t PV_normchi2;
Float_t PV_Nvtx;

//luminosity
Float_t lumiperblock;
Float_t beam1Intensity;
Float_t beam2Intensity;

// HLT
std::vector<TString> hlt_path;

//digi variables
std::vector<Short_t> digi_wheel;
std::vector<Short_t> digi_sector;
std::vector<Short_t> digi_station;
std::vector<Short_t> digi_sl;
std::vector<Short_t> digi_layer;
std::vector<Short_t> digi_wire;
std::vector<Float_t> digi_time;

//DT segment variables
std::vector<Short_t> segm4D_wheel;
std::vector<Short_t> segm4D_sector;
std::vector<Short_t> segm4D_station;

std::vector<Short_t> segm4D_hasPhi;
std::vector<Short_t> segm4D_hasZed;

std::vector<Float_t> segm4D_x_pos_loc;
std::vector<Float_t> segm4D_y_pos_loc;
std::vector<Float_t> segm4D_z_pos_loc;
std::vector<Float_t> segm4D_x_dir_loc;
std::vector<Float_t> segm4D_y_dir_loc;
std::vector<Float_t> segm4D_z_dir_loc;

std::vector<Float_t> segm4D_cosx;
std::vector<Float_t> segm4D_cosy;
std::vector<Float_t> segm4D_cosz;
std::vector<Float_t> segm4D_phi;
std::vector<Float_t> segm4D_theta;
std::vector<Float_t> segm4D_eta;

std::vector<Float_t> segm4D_t0;
std::vector<Float_t> segm4D_vdrift;
std::vector<Float_t> segm4D_phinormchi2;
std::vector<Short_t> segm4D_phinhits;

std::vector<Float_t> segm4D_znormchi2;
std::vector<Short_t> segm4D_znhits;

TClonesArray *segm4D_phiHits_Pos;
TClonesArray *segm4D_phiHits_PosCh;
TClonesArray *segm4D_phiHits_PosErr;
TClonesArray *segm4D_phiHits_Side;
TClonesArray *segm4D_phiHits_Wire;
TClonesArray *segm4D_phiHits_Layer;
TClonesArray *segm4D_phiHits_SuperLayer;
TClonesArray *segm4D_phiHits_Time;
TClonesArray *segm4D_phiHits_TimeCali;

TClonesArray *segm4D_hitsExpPos;
TClonesArray *segm4D_hitsExpWire;

TClonesArray *segm4D_zHits_Pos;
TClonesArray *segm4D_zHits_PosCh;
TClonesArray *segm4D_zHits_PosErr;
TClonesArray *segm4D_zHits_Side;
TClonesArray *segm4D_zHits_Wire;
TClonesArray *segm4D_zHits_Layer;
TClonesArray *segm4D_zHits_Time;
TClonesArray *segm4D_zHits_TimeCali;

//CSC segment variables
std::vector<Short_t> cscsegm_ring;
std::vector<Short_t> cscsegm_chamber;
std::vector<Short_t> cscsegm_station;
std::vector<Float_t> cscsegm_cosx;
std::vector<Float_t> cscsegm_cosy;
std::vector<Float_t> cscsegm_cosz;
std::vector<Float_t> cscsegm_phi;
std::vector<Float_t> cscsegm_eta;
std::vector<Float_t> cscsegm_normchi2;
std::vector<Short_t> cscsegm_nRecHits;

//TM variables
std::vector<Short_t> ltTwinMuxIn_wheel;
std::vector<Short_t> ltTwinMuxIn_sector;
std::vector<Short_t> ltTwinMuxIn_station;
std::vector<Short_t> ltTwinMuxIn_quality;
std::vector<Short_t> ltTwinMuxIn_bx;
std::vector<Float_t> ltTwinMuxIn_phi;
std::vector<Float_t> ltTwinMuxIn_phiB;
std::vector<Short_t> ltTwinMuxIn_is2nd;

std::vector<Short_t> ltTwinMuxOut_wheel;
std::vector<Short_t> ltTwinMuxOut_sector;
std::vector<Short_t> ltTwinMuxOut_station;
std::vector<Short_t> ltTwinMuxOut_quality;
std::vector<Short_t> ltTwinMuxOut_rpcbit;
std::vector<Short_t> ltTwinMuxOut_bx;
std::vector<Float_t> ltTwinMuxOut_phi;
std::vector<Float_t> ltTwinMuxOut_phiB;
std::vector<Short_t> ltTwinMuxOut_is2nd;

std::vector<Short_t> ltTwinMux_thBx;
std::vector<Short_t> ltTwinMux_thWheel;
std::vector<Short_t> ltTwinMux_thSector;
std::vector<Short_t> ltTwinMux_thStation;
std::vector<Short_t> ltTwinMux_thHits;

//muon variables
std::vector<Short_t> Mu_isMuGlobal;
std::vector<Short_t> Mu_isMuTracker;
std::vector<Short_t> Mu_isMuTrackerArb;
std::vector<Short_t> Mu_isMuStandAlone;
std::vector<Short_t> Mu_isMuRPC;

std::vector<Int_t>   Mu_nMatches;
std::vector<Int_t>   Mu_numberOfChambers;
std::vector<Int_t>   Mu_numberOfMatches;
std::vector<Int_t>   Mu_numberOfMatchedStations;
std::vector<UInt_t>   Mu_stationMask;

std::vector<Float_t> Mu_px_mu;
std::vector<Float_t> Mu_py_mu;
std::vector<Float_t> Mu_pz_mu;
std::vector<Float_t> Mu_phi_mu;
std::vector<Float_t> Mu_eta_mu;
std::vector<Short_t> Mu_chargeMu;

std::vector<Int_t>   STAMu_numberOfHits;
std::vector<Int_t>   STAMu_segmIndex;

std::vector<Short_t> STAMu_recHitsSize;
std::vector<Float_t> STAMu_normchi2Mu;
std::vector<Float_t> STAMu_dxyMu;
std::vector<Float_t> STAMu_dzMu;

std::vector<Float_t> GLBMu_normchi2Mu;
std::vector<Float_t> GLBMu_dxyMu;
std::vector<Float_t> GLBMu_dzMu;

std::vector<Int_t> GLBMu_numberOfPixelHits;
std::vector<Int_t> GLBMu_numberOfTrackerHits;

std::vector<Float_t> GLBMu_tkIsoR03;
std::vector<Float_t> GLBMu_ntkIsoR03;
std::vector<Float_t> GLBMu_emIsoR03;
std::vector<Float_t> GLBMu_hadIsoR03;

std::vector<Float_t> TRKMu_normchi2Mu;
std::vector<Float_t> TRKMu_dxyMu;
std::vector<Float_t> TRKMu_dzMu;

std::vector<Int_t> TRKMu_numberOfPixelHits;
std::vector<Int_t> TRKMu_numberOfTrackerLayers;

std::vector<Float_t> TRKMu_tkIsoR03;

std::vector<Int_t> TRKMu_algo;
std::vector<Int_t> TRKMu_origAlgo;

std::vector<Float_t> STAMu_caloCompatibility;
std::vector<Float_t> STAMu_time;
std::vector<Float_t> STAMu_timeNDof;

std::vector<Int_t> RPCMu_numberOfRPCLayers;

TClonesArray *Mu_matches_Wh;  
TClonesArray *Mu_matches_Sec;
TClonesArray *Mu_matches_St;
TClonesArray *Mu_matches_x;
TClonesArray *Mu_matches_y;
TClonesArray *Mu_matches_phi;
TClonesArray *Mu_matches_eta;
TClonesArray *Mu_matches_edgeX;
TClonesArray *Mu_matches_edgeY;

TClonesArray *Mu_hlt_Dr;

std::vector<Float_t> STAMu_z_mb2;
std::vector<Float_t> STAMu_phi_mb2;
std::vector<Float_t> STAMu_pseta_mb2;

//GMT
std::vector<Short_t> gmt_bx;
std::vector<Float_t> gmt_phi;
std::vector<Float_t> gmt_eta;
std::vector<Float_t> gmt_pt;
std::vector<Short_t> gmt_qual;
std::vector<Short_t> gmt_charge;
std::vector<Int_t>   gmt_tf_idx;

//GT // legacy
std::vector<Short_t> gt_algo_bx;
std::vector<Short_t> gt_algo_bit;
std::vector<Short_t> gt_tt_bx;
std::vector<Short_t> gt_tt_bit;

//RPC
std::vector<Int_t>   rpc_region;
std::vector<Int_t>   rpc_clusterSize;
std::vector<Int_t>   rpc_strip;
std::vector<Int_t>   rpc_bx;
std::vector<Int_t>   rpc_station;
std::vector<Int_t>   rpc_sector;
std::vector<Int_t>   rpc_layer;
std::vector<Int_t>   rpc_subsector;
std::vector<Int_t>   rpc_roll;
std::vector<Int_t>   rpc_ring;

Int_t Bmtf_Size;
std::vector<Short_t> Bmtf_Pt;
std::vector<Short_t> Bmtf_Eta;
std::vector<Short_t> Bmtf_Phi;
std::vector<Short_t> Bmtf_GlobalPhi;
std::vector<Short_t> Bmtf_qual;
std::vector<Short_t> Bmtf_ch;
std::vector<Short_t> Bmtf_bx;
std::vector<Short_t> Bmtf_processor;
std::vector<Short_t> Bmtf_trAddress;
std::vector<Short_t> Bmtf_wh;
std::vector<Short_t> Bmtf_FineBit;

Int_t Bmtf_phSize;
std::vector<Int_t> Bmtf_phBx;
std::vector<Int_t> Bmtf_phWh;
std::vector<Int_t> Bmtf_phSe;
std::vector<Int_t> Bmtf_phSt;
std::vector<Float_t>  Bmtf_phAng;
std::vector<Float_t>  Bmtf_phBandAng;
std::vector<Int_t> Bmtf_phCode;
std::vector<Int_t> Bmtf_phTs2Tag;

Int_t Bmtf_thSize;
std::vector<Int_t>   Bmtf_thBx;
std::vector<Int_t>   Bmtf_thWh;
std::vector<Int_t>   Bmtf_thSe;
std::vector<Int_t>   Bmtf_thSt;
std::vector<Int_t> Bmtf_thTheta;
std::vector<Int_t> Bmtf_thCode;

std::vector<Int_t> RpcDigi_TwinMux_bx;
std::vector<Int_t> RpcDigi_TwinMux_strip;
std::vector<Int_t> RpcDigi_TwinMux_region;
std::vector<Int_t> RpcDigi_TwinMux_ring;
std::vector<Int_t> RpcDigi_TwinMux_station;
std::vector<Int_t> RpcDigi_TwinMux_layer;
std::vector<Int_t> RpcDigi_TwinMux_sector;
std::vector<Int_t> RpcDigi_TwinMux_subsector;
std::vector<Int_t> RpcDigi_TwinMux_roll;
std::vector<Int_t> RpcDigi_TwinMux_trIndex;
std::vector<Int_t> RpcDigi_TwinMux_det;
std::vector<Int_t> RpcDigi_TwinMux_subdetId;
std::vector<Int_t> RpcDigi_TwinMux_rawId;

std::vector<Int_t> RpcRechit_TwinMux_region;
std::vector<Int_t> RpcRechit_TwinMux_clusterSize;
std::vector<Int_t> RpcRechit_TwinMux_strip;
std::vector<Int_t> RpcRechit_TwinMux_bx;
std::vector<Int_t> RpcRechit_TwinMux_station;
std::vector<Int_t> RpcRechit_TwinMux_sector;
std::vector<Int_t> RpcRechit_TwinMux_layer;
std::vector<Int_t> RpcRechit_TwinMux_subsector;
std::vector<Int_t> RpcRechit_TwinMux_roll;
std::vector<Int_t> RpcRechit_TwinMux_ring;

