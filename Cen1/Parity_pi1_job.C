using namespace std;

#include "stdio.h"
#include "TString.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TObjArray.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "../MyRef.h"
#include "../src/StLorentzVectorD.h"

const int opt_TOF = 0;
const float PI = TMath::Pi();
const float mp = 0.938272;
const float mpi= 0.1396;
const float pt_trig_up = 1.6;
const float pt_trig_lo = .2;
const float pt_asso_up = 2.0;//1;
const float pt_asso_lo = 0.15;//0.15;
const float EtaCut = 1;
const float DcaCut = 2;
const Float_t Vz_offset = 0;
const Float_t Vz_cut = 30;    //40 for 19 GeV, 40 for 39 GeV, 75 for 11 GeV and 70 for 7GeV
const int cenDef[9] = {10, 22, 43, 76, 125, 193, 10, 200, 400}; // 200 GeV run11

const int bad_Ref_day2[12] = {13703,13800,13802,14005,14107,14301,14307,14600,15901,16304,16407,16705};
const int Nrun_MB1 = 134;
const int Nrun_MB2 = 60;
const int Nrun_MB5 = 107;
const int Nrun_MB6 = 42;
const int bad_Ref_day3_MB1[Nrun_MB1] = {133010,133011,133022,133027,133028,134065,127019,128021,133020,137004,127039,132024,132025,132026,132032,127020,127048,127049,128028,128029,128030,128031,128032,132052,133041,136039,127002,132061,134017,134063,135004,135039,137020,127030,128007,132062,133005,133039,133040,133053,133054,134006,134007,134018,134026,134028,134038,134041,135019,135020,135045,135057,136007,136032,136064,137010,127009,127046,128038,132019,132022,132044,132045,133002,133019,133038,133052,134005,134040,134055,135002,135033,135048,135049,136069,137003,127021,127022,127023,127024,128024,128025,132020,132021,132023,132033,132048,132051,132057,132063,132065,133021,134008,135012,135021,135024,135030,135054,136044,136081,136086,127003,127010,127011,127017,127018,127032,132009,132034,132043,132066,132069,133018,134023,134057,136005,136006,136014,136017,136022,136023,136024,136025,136027,136028,136029,136030,136031,136034,136053,136054,136070,136071,138017};//MB1
const int bad_Ref_day3_MB2[Nrun_MB2] = {139032,139043,139044,139045,142002,139042,140021,140029,142063,142064,142065,144004,138081,138082,138087,138088,138089,138090,138091,139002,139003,139006,139007,139008,139009,139010,139015,139016,139017,139018,139021,142016,142033,142061,142062,144051,138092,139019,139020,140016,140020,141003,141004,141026,141062,141065,142001,142013,142023,142034,142046,142068,142076,143009,143024,143058,144016,144028,144033,145003};  //MB2
const int bad_Ref_day3_MB5[Nrun_MB5] = {155050,155056,158010,165028,154043,154044,154045,155058,158069,158070,158072,158073,164067,154046,154047,155008,155009,156015,156062,156063,158074,162015,154067,155002,155012,155047,156008,156009,157023,157030,157052,158006,159023,160021,161006,161015,161060,162004,162028,162034,163024,163058,164009,164056,164066,164089,165013,154048,154066,155011,155021,155038,155051,155060,155062,155064,156004,156035,156056,157012,157014,157031,157038,157051,158015,158021,158026,158040,158041,158051,158054,158056,158057,158058,158061,159005,159021,159022,159024,160016,160025,161007,161014,161017,161020,161022,161053,162017,162030,162035,162055,162056,162057,162058,163006,163008,163015,164011,164037,164043,164086,165001,165003,165005,165007,165026,165031}; //MB5
const int bad_Ref_day3_MB6[Nrun_MB6] = {166051,170016,167014,170034,170050,170051,171009,171015,167049,168010,169028,169032,170007,165042,166052,166059,167002,167040,167048,169031,169059,170009,170012,170018,170020,170031,171004,171014,166002,166003,167015,167024,168009,168022,168060,168077,169033,169034,170044,170045,170054,170056};//MB6

//pion
const float PP0[9] = {8.81953e-01,8.69628e-01,8.79790e-01,8.64348e-01,8.44556e-01,8.16006e-01,7.67517e-01,7.20570e-01,6.79516e-01};
const float PP1[9] = {2.30849e-01,2.34070e-01,2.27089e-01,2.32444e-01,2.30446e-01,2.33877e-01,2.36398e-01,2.42266e-01,2.45421e-01};
const float PP2[9] = {1.40367e+00,1.47586e+00,1.29833e+00,1.34398e+00,1.39598e+00,1.36630e+00,1.43580e+00,1.43284e+00,1.41702e+00};

void Parity_pi1_job(int cen=1, int job=0){	//main_function

struct StRefMultCorr refmultCorrUtil  = StRefMultCorr("refmult") ;

delete gRandom;
gRandom = new TRandom3(0);
float Eweight = 1;
const int Phibin = 80;
char fname[200];
float PsiShiftE1=0,PsiShiftE2=0,PsiShiftE3=0,PsiShiftE4=0;
float PsiShiftW1=0,PsiShiftW2=0,PsiShiftW3=0,PsiShiftW4=0;
float PsiShiftf1=0,PsiShiftf2=0,PsiShiftf3=0,PsiShiftf4=0;
float PsiShiftb1=0,PsiShiftb2=0,PsiShiftb3=0,PsiShiftb4=0;
float PsiShiftF1=0,PsiShiftF2=0,PsiShiftF3=0,PsiShiftF4=0;
const int order = 1;
float PhiMean[4*order]={0,};
float PhiWgtFF[Phibin][4],PhiWgtRF[Phibin][4];
TProfile2D *TPCmean_FF, *TPCmean_RF;
TH2D *TPCPhi_FF, *TPCPhi_RF;
TH1D *TOF_eff;
TProfile2D *Read_TPC_EP_full, *Read_TPC_EP_east, *Read_TPC_EP_west, *Read_TPC_EP_for, *Read_TPC_EP_bac;
sprintf(fname,"cen%d.weight_pi_minbias5.root",cen);
TFile *fWgt=new TFile(fname,"READ");
if(!fWgt->IsOpen()) cout<<"no phi weight files!"<<endl;
if(fWgt->IsOpen()) {
	TOF_eff = (TH1D*)fWgt->Get("rc");
        TPCmean_FF = (TProfile2D*)fWgt->Get("TPCmeanPhi_FF");
        TPCmean_RF = (TProfile2D*)fWgt->Get("TPCmeanPhi_RF");
	TPCPhi_FF = (TH2D*)fWgt->Get("Hist_Phi_FF");
	TPCPhi_RF = (TH2D*)fWgt->Get("Hist_Phi_RF");
	Read_TPC_EP_full = (TProfile2D*)fWgt->Get("pTPC_EP_full");
	Read_TPC_EP_east = (TProfile2D*)fWgt->Get("pTPC_EP_east");
	Read_TPC_EP_west = (TProfile2D*)fWgt->Get("pTPC_EP_west");
        Read_TPC_EP_for = (TProfile2D*)fWgt->Get("pTPC_EP_for");
        Read_TPC_EP_bac = (TProfile2D*)fWgt->Get("pTPC_EP_bac");
	for(int j=0;j<4;j++) {
                float PhiMeanFF = TPCPhi_FF->ProjectionX("",j+1,j+1)->GetSum()/(float)Phibin;
                float PhiMeanRF = TPCPhi_RF->ProjectionX("",j+1,j+1)->GetSum()/(float)Phibin;
		for(int i=0;i<Phibin;i++) {
                	PhiWgtFF[i][j] = (TPCPhi_FF->GetBinContent(i+1,j+1)>0)? PhiMeanFF/(TPCPhi_FF->GetBinContent(i+1,j+1)):1;
                        PhiWgtRF[i][j] = (TPCPhi_RF->GetBinContent(i+1,j+1)>0)? PhiMeanRF/(TPCPhi_RF->GetBinContent(i+1,j+1)):1;
                	//cout<<" PhiWgt= "<<PhiWgt[i][j];
		}
	}
}  else {for(int j=0;j<4;j++) for(int i=0;i<Phibin;i++) { PhiWgtFF[i][j] = 1; PhiWgtRF[i][j] = 1; }

        }
/*
//prepare the files
  int nfile = 0;
  TChain* chain = new TChain("StrangenessDst");
  std::ifstream inlist(InputFileList);
  if(!inlist)cout<<" can't open file "<<endl;  //prepare the files
  char  *fileName=new char[200];
  while(!inlist.eof()) {
	inlist.getline(fileName,200,'\n');
	if(!inlist.eof()) nfile += chain->Add(fileName);
  }
*/
	TChain* chain = new TChain("StrangenessDst");
	int nfile = 0;
//        nfile += chain->Add("/media/Disk_Chen/Data/200GeV_run11/data_minbias5_new/*.flow.nt.root");
char filelist[200];
sprintf(filelist,"/media/Disk_Chen/Data/200GeV_run11/data_minbias5_new/*%d.flow.nt.root",job);
	nfile += chain->Add(filelist);
		cout <<"Added "<<nfile<<" files"<<endl;
		cout<<"# entries in chain: "<<chain->GetEntries()<<endl;

char fname_out[200];
sprintf(fname_out,"cen%d.ntuple_result_Parity_pi_DCA2_minbias5_eff_job%d.root",cen,job);
	TFile fout(fname_out,"RECREATE");
	//defining variables
	Float_t PVtxz, VPDvz, Bz, psi_E, psi_W, mod_E, mod_W, ZDCe, ZDCw, ZDCcoin, BBCco;			//run, event info
	Int_t   Run, Day, Day2, Day3, Event, Trigger, RefMult, TOFMult, Centrality, NPTracks;	//
	Float_t Charge, ndEdx, nSigma_p, nSigma_pi, DCAGlobal, Eta, Theta, Phi, Pt, P;		//track info	

	//defining histograms
	TH1D* hEventTally = new TH1D("EventTally","Event Tally",10,0,1);
	hEventTally->SetBit(TH1::kCanRebin);
	hEventTally->SetStats(0);

        TProfile *Ref_Day3 = new TProfile("Ref_Day3","RefMult vs Run",100000,80000,180000,-0.5,999.5);
        TProfile *TOF_Day3 = new TProfile("TOF_Day3","TOFMult vs Run",100000,80000,180000,0,5000);
        TProfile *NPT_Day3 = new TProfile("NPT_Day3","NPTracks vs Run",100000,80000,180000,0,2000);
        TProfile *TPC_Day3_cos2 = new TProfile("TPC_Day3_cos2","cos(2*psi) vs Run", 100000,80000,180000,-1,1);
        TProfile *TPC_Day3_sin2 = new TProfile("TPC_Day3_sin2","sin(2*psi) vs Run", 100000,80000,180000,-1,1);
	TH1D *Hist_RefMult = new TH1D("Hist_RefMult","Hist_RefMult",1000,-0.5,999.5);
	TH1D *Hist_TOFMult = new TH1D("Hist_TOFMult","Hist_TOFMult",5000,-0.5,4999.5);
        TProfile *p_RefMult = new TProfile("p_RefMult","p_RefMult",50,0.5,1000.5, 0,1000, "");
        TProfile *p_TOFMult = new TProfile("p_TOFMult","p_TOFMult",50,0.5,5000.5, 0,5000, "");

	TH1D *hBz      = new TH1D("hBz","magnetic field",10, -10, 10);
	TH1D *hTrigger = new TH1D("hTrigger","hTrigger",400, 0.5, 400.5);
	TH1D *hCentrality = new TH1D("hCentrality","hCentrality",10,0,10);
	TH1D *hVertexZ = new TH1D("hVertexZ","hVertexZ",100,-100,100);
	TH2D *hMult_Vz = new TH2D("hMult_Vz","hMult_Vz",1000,-0.5,999.5,100,-100,100);
        TH2D *hMult_Vz_new = new TH2D("hMult_Vz_new","hMult_Vz_new",1000,-0.5,999.5,100,-100,100);
	TH1D *hBBC_coin = new TH1D("hBBC_coin","hBBC_coin",120,0,60000);
        TProfile2D *pTPCmeanPhi_FF = new TProfile2D("TPCmeanPhi_FF","TPCmeanPhi_FF",
                                        8*order,0.5,8*order+0.5,10000,8000,18000,-1,1,"");
        TProfile2D *pTPCmeanPhi_RF = new TProfile2D("TPCmeanPhi_RF","TPCmeanPhi_RF",
                                        8*order,0.5,8*order+0.5,10000,8000,18000,-1,1,"");

        TH1D* Hist_proton = new TH1D("Hist_proton","Hist_proton",50,-0.5,49.5);
        TH1D* Hist_pbar = new TH1D("Hist_pbar","Hist_pbar",50,-0.5,49.5);
	TH1D* Hist_netP = new TH1D("Hist_netP","Hist_netP",99,-49.5,49.5);
        TH2D *hEtaPtDist = new TH2D("EtaPtDist","EtaPtDist",26, -1.3, 1.3,300,0,15);
	TH2D* Hist_Phi = new TH2D("Hist_Phi","Hist_Phi",Phibin,-PI,PI,4,0.5,4.5);
	TH2D* Hist_TPC_EP_east = new TH2D("Hist_TPC_EP_east","Hist_TPC_EP_east",36,0,PI,100,80,180);
        TH2D* Hist_TPC_EP_west = new TH2D("Hist_TPC_EP_west","Hist_TPC_EP_west",36,0,PI,100,80,180);
        TH2D* Hist_TPC_EP_for = new TH2D("Hist_TPC_EP_for","Hist_TPC_EP_for",36,0,PI,100,80,180);
        TH2D* Hist_TPC_EP_bac = new TH2D("Hist_TPC_EP_bac","Hist_TPC_EP_bac",36,0,PI,100,80,180);
        TH2D* Hist_TPC_EP_full = new TH2D("Hist_TPC_EP_full","Hist_TPC_EP_full",36,0,PI,100,80,180);
        TH2D* Hist_TPC_EP_east_flat = new TH2D("Hist_TPC_EP_east_flat","Hist_TPC_EP_east_flat",36,0,PI,100,80,180);
        TH2D* Hist_TPC_EP_west_flat = new TH2D("Hist_TPC_EP_west_flat","Hist_TPC_EP_west_flat",36,0,PI,100,80,180);
        TH2D* Hist_TPC_EP_for_flat = new TH2D("Hist_TPC_EP_for_flat","Hist_TPC_EP_for_flat",36,0,PI,100,80,180);
        TH2D* Hist_TPC_EP_bac_flat = new TH2D("Hist_TPC_EP_bac_flat","Hist_TPC_EP_bac_flat",36,0,PI,100,80,180);
        TH2D* Hist_TPC_EP_full_flat = new TH2D("Hist_TPC_EP_full_flat","Hist_TPC_EP_full_flat",36,0,PI,100,80,180);
	TProfile2D* pTPC_EP_east = new TProfile2D("pTPC_EP_east","pTPC_EP_east",4,0.5,4.5,100,80,180,-1,1,"");
        TProfile2D* pTPC_EP_west = new TProfile2D("pTPC_EP_west","pTPC_EP_west",4,0.5,4.5,100,80,180,-1,1,"");
        TProfile2D* pTPC_EP_for = new TProfile2D("pTPC_EP_for","pTPC_EP_for",4,0.5,4.5,100,80,180,-1,1,"");
        TProfile2D* pTPC_EP_bac = new TProfile2D("pTPC_EP_bac","pTPC_EP_bac",4,0.5,4.5,100,80,180,-1,1,"");
        TProfile2D* pTPC_EP_full = new TProfile2D("pTPC_EP_full","pTPC_EP_full",4,0.5,4.5,100,80,180,-1,1,"");
	TH1F* Hist_dif_count = new TH1F("Hist_dif_count","Hist_dif_count",500,-250,250);
        TH1F* Hist_ful_count = new TH1F("Hist_ful_count","Hist_ful_count",1000,0,1000);


	TH1D* Hist_Pt = new TH1D("Hist_Pt","Hist_Pt",300,0,15);
        TH2D* Hist_Phi_FF = new TH2D("Hist_Phi_FF","Hist_Phi_FF",Phibin,-PI,PI,4,0.5,4.5);
        TH2D* Hist_Phi_RF = new TH2D("Hist_Phi_RF","Hist_Phi_RF",Phibin,-PI,PI,4,0.5,4.5);
        TH2D* Hist_Phi_FF_new = new TH2D("Hist_Phi_FF_new","Hist_Phi_FF_new",Phibin,-PI,PI,4,0.5,4.5);
        TH2D* Hist_Phi_RF_new = new TH2D("Hist_Phi_RF_new","Hist_Phi_RF_new",Phibin,-PI,PI,4,0.5,4.5);
        TProfile *Hist_cos = new TProfile("Hist_cos","Hist_cos",2,0.5,2.5,-1,1,"");
        TProfile2D *Hist_cos_RefMult = new TProfile2D("Hist_cos_RefMult","Hist_cos_RefMult",2,0.5,2.5, 50,0.5,1000.5,-1,1,"");
	TProfile2D *Hist_cos_TOFMult = new TProfile2D("Hist_cos_TOFMult","Hist_cos_TOFMult",2,0.5,2.5, 50,0.5,5000.5,-1,1,"");
	TProfile2D *Hist_cos_ZDC = new TProfile2D("Hist_cos_ZDC","Hist_cos_ZDC",2,0.5,2.5, 50,0.5,20000.5,-1,1,"");
        TProfile2D *Hist_cos_BBC_RefMult = new TProfile2D("Hist_cos_BBC_RefMult","Hist_cos_BBC_RefMult",60,0,60000, 50,0.5,1000.5,-1,1,"");
       	TH1D* hDpt   = new TH1D("hDpt","hDpt",200,0,2);
	TH1D* hQinv  = new TH1D("hQinv","hQinv",1000,0,10);
	TH1D* hQinv2 = new TH1D("hQinv2","hQinv2",1100,-1,10);
        TH1D* Hist_InvM_os = new TH1D("Hist_InvM_os","Hist_invmass_os",1800,0.2,2);
        TH1D* Hist_InvM_ss = new TH1D("Hist_InvM_ss","Hist_invmass_ss",1800,0.2,2);

	TProfile *Hist_v2_pt_obs = new TProfile("Hist_v2_pt_obs","Hist_v2_pt_obs",300,0,15,-100,100,"");
        TProfile *Hist_v2_pt_pos_obs = new TProfile("Hist_v2_pt_pos_obs","Hist_v2_pt_pos_obs",300,0,15,-100,100,"");
        TProfile *Hist_v2_pt_neg_obs = new TProfile("Hist_v2_pt_neg_obs","Hist_v2_pt_neg_obs",300,0,15,-100,100,"");
	TProfile *p_v2_RefMult_obs = new TProfile("p_v2_RefMult_obs","p_v2_RefMult_obs",50,0.5, 1000.5,-100,100,"");
        TProfile *p_v2_TOFMult_obs = new TProfile("p_v2_TOFMult_obs","p_v2_TOFMult_obs",50,0.5,5000.5,-100,100,"");
        TProfile *p_v2_ZDC_obs = new TProfile("p_v2_ZDC_obs","p_v2_ZDC_obs",50,0.5,20000.5,-100,100,"");
        TProfile2D *p_v2_BBC_RefMult_obs = new TProfile2D("p_v2_BBC_RefMult_obs","p_v2_BBC_RefMult_obs",60,0,60000,50,0.5, 1000.5,-100,100,"");
        TProfile *p_RefMult_ZDC_obs = new TProfile("p_RefMult_ZDC_obs","p_RefMult_ZDC_obs",50,0.5,20000.5,0,1000,"");
        TProfile *p_TOFMult_ZDC_obs = new TProfile("p_TOFMult_ZDC_obs","p_TOFMult_ZDC_obs",50,0.5,20000.5,0,5000,"");

        TProfile *pParity_int_obs1 = new TProfile("Parity_int_obs1","Parity_int_obs1",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_obs2 = new TProfile("Parity_int_obs2","Parity_int_obs2",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_obs3 = new TProfile("Parity_int_obs3","Parity_int_obs3",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_r_obs1 = new TProfile("Parity_int_r_obs1","Parity_int_r_obs1",4,0.5,4.5,-1600,1600,"");
        TProfile *pParity_int_r_obs2 = new TProfile("Parity_int_r_obs2","Parity_int_r_obs2",4,0.5,4.5,-1600,1600,"");
        TProfile *pParity_int_r_obs3 = new TProfile("Parity_int_r_obs3","Parity_int_r_obs3",4,0.5,4.5,-1600,1600,"");
        TProfile *pParity_int_rr_obs1 = new TProfile("Parity_int_rr_obs1","Parity_int_rr_obs1",4,0.5,4.5,-1600,1600,"");
        TProfile *pParity_int_rr_obs2 = new TProfile("Parity_int_rr_obs2","Parity_int_rr_obs2",4,0.5,4.5,-1600,1600,"");
        TProfile *pParity_int_rr_obs3 = new TProfile("Parity_int_rr_obs3","Parity_int_rr_obs3",4,0.5,4.5,-1600,1600,"");
        TProfile *pParity_int_s_obs1 = new TProfile("Parity_int_s_obs1","Parity_int_s_obs1",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_s_obs2 = new TProfile("Parity_int_s_obs2","Parity_int_s_obs2",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_s_obs3 = new TProfile("Parity_int_s_obs3","Parity_int_s_obs3",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_ss_obs1 = new TProfile("Parity_int_ss_obs1","Parity_int_ss_obs1",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_ss_obs2 = new TProfile("Parity_int_ss_obs2","Parity_int_ss_obs2",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_ss_obs3 = new TProfile("Parity_int_ss_obs3","Parity_int_ss_obs3",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_w_obs1 = new TProfile("Parity_int_w_obs1","Parity_int_w_obs1",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_w_obs2 = new TProfile("Parity_int_w_obs2","Parity_int_w_obs2",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_w_obs3 = new TProfile("Parity_int_w_obs3","Parity_int_w_obs3",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_ww_obs1 = new TProfile("Parity_int_ww_obs1","Parity_int_ww_obs1",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_ww_obs2 = new TProfile("Parity_int_ww_obs2","Parity_int_ww_obs2",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_ww_obs3 = new TProfile("Parity_int_ww_obs3","Parity_int_ww_obs3",4,0.5,4.5,-100,100,"");

        TProfile *pParity_int_same_run = new TProfile("Parity_int_same_run","Parity_int_same_run",10000,8000,18000,-100,100,"");
        TProfile *pParity_int_rr_same_run = new TProfile("Parity_int_rr_same_run","Parity_int_rr_same_run",10000,8000,18000,-100,100,"");
        TProfile *pParity_int_ss_same_run = new TProfile("Parity_int_ss_same_run","Parity_int_ss_same_run",10000,8000,18000,-100,100,"");
        TProfile *pParity_int_ww_same_run = new TProfile("Parity_int_ww_same_run","Parity_int_ww_same_run",10000,8000,18000,-100,100,"");
        TProfile *pParity_int_oppo_run = new TProfile("Parity_int_oppo_run","Parity_int_oppo_run",10000,8000,18000,-100,100,"");
        TProfile *pParity_int_rr_oppo_run = new TProfile("Parity_int_rr_oppo_run","Parity_int_rr_oppo_run",10000,8000,18000,-100,100,"");
        TProfile *pParity_int_ss_oppo_run = new TProfile("Parity_int_ss_oppo_run","Parity_int_ss_oppo_run",10000,8000,18000,-100,100,"");
        TProfile *pParity_int_ww_oppo_run = new TProfile("Parity_int_ww_oppo_run","Parity_int_ww_oppo_run",10000,8000,18000,-100,100,"");
        TProfile *pParity_int_ss_ran_obs1 = new TProfile("Parity_int_ss_ran_obs1","Parity_int_ss_ran_obs1",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_ss_ran_obs2 = new TProfile("Parity_int_ss_ran_obs2","Parity_int_ss_ran_obs2",4,0.5,4.5,-100,100,"");
        TProfile *pParity_int_ss_ran_obs3 = new TProfile("Parity_int_ss_ran_obs3","Parity_int_ss_ran_obs3",4,0.5,4.5,-100,100,"");

        TProfile2D *pParity_eta_ss_obs1 = new TProfile2D("Parity_eta_ss_obs1","Parity_eta_ss_obs1",12,0.5,12.5,20,-1,1,-100,100,"");
        TProfile2D *pParity_eta_ss_obs2 = new TProfile2D("Parity_eta_ss_obs2","Parity_eta_ss_obs2",12,0.5,12.5,20,-1,1,-100,100,"");
        TProfile2D *pParity_eta_ss_obs3 = new TProfile2D("Parity_eta_ss_obs3","Parity_eta_ss_obs3",12,0.5,12.5,20,-1,1,-100,100,"");

        TProfile2D *pParity_Deta_ss_obs1 = new TProfile2D("Parity_Deta_ss_obs1","Parity_Deta_ss_obs1",12,0.5,12.5,20,0,2,-100,100,"");
        TProfile2D *pParity_Deta_ss_obs2 = new TProfile2D("Parity_Deta_ss_obs2","Parity_Deta_ss_obs2",12,0.5,12.5,20,0,2,-100,100,"");
        TProfile2D *pParity_Deta_ss_obs3 = new TProfile2D("Parity_Deta_ss_obs3","Parity_Deta_ss_obs3",12,0.5,12.5,20,0,2,-100,100,"");
        TProfile2D *pParity_pt_ss_obs1  = new TProfile2D("Parity_pt_ss_obs1","Parity_pt_ss_obs1",12,0.5,12.5,20,0,2.0,-100,100,"");
        TProfile2D *pParity_pt_ss_obs2  = new TProfile2D("Parity_pt_ss_obs2","Parity_pt_ss_obs2",12,0.5,12.5,20,0,2.0,-100,100,"");
        TProfile2D *pParity_pt_ss_obs3  = new TProfile2D("Parity_pt_ss_obs3","Parity_pt_ss_obs3",12,0.5,12.5,20,0,2.0,-100,100,"");
        TProfile2D *pParity_Dpt_ss_obs1 = new TProfile2D("Parity_Dpt_ss_obs1","Parity_Dpt_ss_obs1",12,0.5,12.5,200,0,2.0,-100,100,"");
        TProfile2D *pParity_Dpt_ss_obs2 = new TProfile2D("Parity_Dpt_ss_obs2","Parity_Dpt_ss_obs2",12,0.5,12.5,200,0,2.0,-100,100,"");
        TProfile2D *pParity_Dpt_ss_obs3 = new TProfile2D("Parity_Dpt_ss_obs3","Parity_Dpt_ss_obs3",12,0.5,12.5,200,0,2.0,-100,100,"");
        TProfile2D *pParity_Q_ss_obs1   = new TProfile2D("Parity_Q_ss_obs1","Parity_Q_ss_obs1",4,0.5,4.5,500,0,5.0,-100,100,"");
        TProfile2D *pParity_Q_ss_obs2   = new TProfile2D("Parity_Q_ss_obs2","Parity_Q_ss_obs2",4,0.5,4.5,500,0,5.0,-100,100,"");
        TProfile2D *pParity_Q_ss_obs3   = new TProfile2D("Parity_Q_ss_obs3","Parity_Q_ss_obs3",4,0.5,4.5,500,0,5.0,-100,100,"");

        TProfile *pParity_noHBT_ss_obs2   = new TProfile("Parity_noHBT_ss_obs2","Parity_noHBT_ss_obs2",4,0.5,4.5,-100,100,"");
        TProfile *pParity_noHBT_ss_obs3   = new TProfile("Parity_noHBT_ss_obs3","Parity_noHBT_ss_obs3",4,0.5,4.5,-100,100,"");

	Int_t nentries = chain->GetEntries();
//	nentries = 50000;
	//loop through events
	for(int i = 0; i < nentries; i++){

		if((i+1)%1000==0) cout<<"Processing entry == "<< i+1 <<" == out of "<<nentries<<".\n";
		chain->GetEntry(i);

		TLeaf* leaf_RunId   = chain->GetLeaf("mRunId");		
		TLeaf* leaf_EventId = chain->GetLeaf("mEventId");
		TLeaf* leaf_Trigger = chain->GetLeaf("mTrigger");
                TLeaf* leaf_Bz	    = chain->GetLeaf("mBfield");
                TLeaf* leaf_PrimaryVertexZ = chain->GetLeaf("mPrimaryVertexPositionZ");
		TLeaf* leaf_VPDVz   = chain->GetLeaf("mVPDEWdiff");
                TLeaf* leaf_RefMult = chain->GetLeaf("mRefMult");
                TLeaf* leaf_TOFMult = chain->GetLeaf("mTOFmult");
		TLeaf* leaf_Centrality = chain->GetLeaf("mCentrality");
		TLeaf* leaf_ZDC_EP  = chain->GetLeaf("mZDC_EP");
                TLeaf* leaf_ZDC	    = chain->GetLeaf("mZDC");
                TLeaf* leaf_ZDCcoin = chain->GetLeaf("mZDCCoin");
		TLeaf* leaf_BBCco   = chain->GetLeaf("mBBCCoin");
                TLeaf* leaf_NoTracks = chain->GetLeaf("mNumberOfV0s");

		Run	= (int)leaf_RunId->GetValue(0);
		Event	= (int)leaf_EventId->GetValue(0);
		Trigger = (int)leaf_Trigger->GetValue(0);
		Bz	= leaf_Bz->GetValue(0);
		PVtxz	= leaf_PrimaryVertexZ->GetValue(0);
		VPDvz   = leaf_VPDVz->GetValue(0);
		RefMult = (int)leaf_RefMult->GetValue(0);
                TOFMult = (int)leaf_TOFMult->GetValue(0);
		NPTracks= (int)leaf_NoTracks->GetValue(0);
		ZDCw	= leaf_ZDC->GetValue(1) + leaf_ZDC->GetValue(2) + leaf_ZDC->GetValue(3);
		ZDCe	= leaf_ZDC->GetValue(5) + leaf_ZDC->GetValue(6) + leaf_ZDC->GetValue(7);
                ZDCcoin = leaf_ZDCcoin->GetValue(0);
		BBCco   = leaf_BBCco->GetValue(0);
		psi_E   = leaf_ZDC_EP->GetValue(0);
		psi_W	= leaf_ZDC_EP->GetValue(2);
		mod_E	= leaf_ZDC_EP->GetValue(1);
                mod_W   = leaf_ZDC_EP->GetValue(3);
		Centrality = leaf_Centrality->GetValue(0);
		Day 	= (int)((Run-12000000)/1000); 
		Day2    = (int)((Run-12000000)/10);
                Day3    = (int)((Run-12000000)/1);

                int Bad =0;
                for(int jj=0;jj<12;jj++) if(Day2== bad_Ref_day2[jj]) {Bad = 1;break;}
                if(Day3<=138024) {for(int jj=0;jj<Nrun_MB1;jj++) if(Day3 == bad_Ref_day3_MB1[jj]) {Bad = 1;break;}}
                else if(Day3<=145020) {for(int jj=0;jj<Nrun_MB2;jj++) if(Day3 == bad_Ref_day3_MB2[jj]) {Bad = 1;break;}}
                else if(Day3<=154021) {}
                else if(Day3<=165031) {for(int jj=0;jj<Nrun_MB5;jj++) if(Day3 == bad_Ref_day3_MB5[jj]) {Bad = 1;break;}}
                else {for(int jj=0;jj<Nrun_MB6;jj++) if(Day3 == bad_Ref_day3_MB6[jj]) {Bad = 1;break;}}
                if(Bad) continue;  //bad run

		
		if(RefMult) {
			float gM = 0.9995 + 21.89/(4.191*RefMult-18.17) - 2.723e-5*(4.191*RefMult-18.17);
			Eweight = gM + 0.0009326*(gM-1)*PVtxz;
		}
		hBz->Fill(Bz);
                hTrigger->Fill(Trigger);
		if((Trigger%5)<2) continue;

                Centrality = 0;
                for(int j=0;j<9;j++) if(RefMult>cenDef[j]) Centrality = j+1;

                refmultCorrUtil.init(Run);
                if ( refmultCorrUtil.isBadRun(Run) ) continue;
                refmultCorrUtil.initEvent(RefMult, PVtxz, ZDCcoin) ;
                Int_t cent9  = 1 + refmultCorrUtil.getCentralityBin9() ;
                Eweight = refmultCorrUtil.getWeight();
                Centrality = cent9;

		hVertexZ->Fill(PVtxz);
		hEventTally->Fill("Total Event",1);

                if(TOFMult < -100+3.0*RefMult) continue;                   //remove pile-up
                if(TOFMult > 180 +5.2*RefMult) continue;
                if(TMath::Abs(PVtxz) > Vz_cut) continue;              //Z-vertex cut; track quality cut done in PkNtupleMaker
		if((PVtxz-VPDvz)>3 || (PVtxz-VPDvz)<-3) continue;

                Ref_Day3->Fill(Day3,RefMult);
                TOF_Day3->Fill(Day3,TOFMult);
                NPT_Day3->Fill(Day3,NPTracks);
		Hist_RefMult->Fill(RefMult);
		Hist_TOFMult->Fill(TOFMult);
                if(RefMult>10) {p_RefMult->Fill(RefMult,RefMult,Eweight);
                p_TOFMult->Fill(TOFMult,TOFMult,Eweight);}

                hCentrality->Fill(Centrality,Eweight);
		hMult_Vz->Fill(RefMult,PVtxz);
		hMult_Vz_new->Fill(RefMult,PVtxz,Eweight);
		if(cen && Centrality != cen) continue;

		hBBC_coin->Fill(BBCco);
                p_RefMult_ZDC_obs->Fill(ZDCe+ZDCw,RefMult,Eweight);
                p_TOFMult_ZDC_obs->Fill(ZDCe+ZDCw,TOFMult,Eweight);

                TLeaf* leaf_PtV0        = chain->GetLeaf("fV0s.mPtV0");
                TLeaf* leaf_EtaV0       = chain->GetLeaf("fV0s.mEtaV0");
                TLeaf* leaf_PhiV0       = chain->GetLeaf("fV0s.mPhiV0");
                TLeaf* leaf_ChargeV0    = chain->GetLeaf("fV0s.mChargeV0");
                TLeaf* leaf_DCAglobalV0 = chain->GetLeaf("fV0s.mDcaV0");
                TLeaf* leaf_ndEdxV0     = chain->GetLeaf("fV0s.mndEdxV0");
		TLeaf* leaf_TofFlagV0   = chain->GetLeaf("fV0s.mTofFlagV0");
                TLeaf* leaf_TofYLocalV0   = chain->GetLeaf("fV0s.mTofYLocalV0");
                TLeaf* leaf_TofM2V0   = chain->GetLeaf("fV0s.mTofM2V0");
                TLeaf* leaf_nSigmaPV0   = chain->GetLeaf("fV0s.mnSigma_pV0");
                TLeaf* leaf_nSigmaPiV0  = chain->GetLeaf("fV0s.mnSigma_piV0");
                TLeaf* leaf_nSigmaKV0   = chain->GetLeaf("fV0s.mnSigma_kV0");
		TLeaf* leaf_nSigmaEV0   = chain->GetLeaf("fV0s.mnSigma_eV0");

	  int Np = 0, Npbar = 0;
          int Ntof = 0;
//TPC EP reconstruction
          TVector2 mQ, mQ1, mQ2, mQ3, mQ4;
          Double_t mQx=0., mQy=0., mQx1=0., mQy1=0., mQx2=0., mQy2=0., mQx3=0., mQy3=0., mQx4=0., mQy4=0.;
          int Fcount = 0, Ecount = 0, Wcount =0;
          for(int trk = 0; trk < NPTracks; trk++) {
                float EtaAsso   = leaf_EtaV0->GetValue(trk);
                float PtAsso    = leaf_PtV0->GetValue(trk);
                float PhiAsso   = leaf_PhiV0->GetValue(trk);
                float DCAglAsso = leaf_DCAglobalV0->GetValue(trk);
                float ChargeAsso= leaf_ChargeV0->GetValue(trk);
		nSigma_p  = leaf_nSigmaPV0->GetValue(trk);
		float nSigma_K  = leaf_nSigmaKV0->GetValue(trk);
		float nSigma_e  = leaf_nSigmaEV0->GetValue(trk);
		float TOFflag   = leaf_TofFlagV0->GetValue(trk);
                if(TOFflag>0) Ntof++;

		// also nSigma_e>1
		if(PtAsso>0.4 && PtAsso<1 && DCAglAsso<1 && EtaAsso>-1 && EtaAsso<1 && nSigma_p>-2 && nSigma_p<2 && nSigma_K>2 && nSigma_e >1 && ChargeAsso >0) Np++;
                if(PtAsso>0.4 && PtAsso<1 && DCAglAsso<1 && EtaAsso>-1 && EtaAsso<1 && nSigma_p>-2 && nSigma_p<2 && nSigma_K>2 && nSigma_e >1 && ChargeAsso <0) Npbar++;

if(PtAsso > pt_asso_up || PtAsso < pt_asso_lo) continue;
if(DCAglAsso > DcaCut) continue;
if(EtaAsso>EtaCut || EtaAsso<-EtaCut) continue;
if(opt_TOF ==1 && TOFflag<1) continue;

                Fcount++;
                int n = (int)((PhiAsso+PI)/2./PI*Phibin);
                int fl = (PVtxz > 0)? ((ChargeAsso > 0)? 1:2):((ChargeAsso > 0)? 3:4);
                Hist_Phi->Fill(PhiAsso,fl,PtAsso);
                float W_phi = (EtaAsso>0)? PhiWgtFF[n][fl-1]:PhiWgtRF[n][fl-1];
                mQx += W_phi*PtAsso * cos(PhiAsso * 2.);
                mQy += W_phi*PtAsso * sin(PhiAsso * 2.);

                if(EtaAsso>0 && PVtxz>0 && ChargeAsso > 0) Hist_Phi_FF->Fill(PhiAsso,1,Eweight);
                if(EtaAsso>0 && PVtxz>0 && ChargeAsso < 0) Hist_Phi_FF->Fill(PhiAsso,2,Eweight);
                if(EtaAsso<0 && PVtxz>0 && ChargeAsso > 0) Hist_Phi_RF->Fill(PhiAsso,1,Eweight);
                if(EtaAsso<0 && PVtxz>0 && ChargeAsso < 0) Hist_Phi_RF->Fill(PhiAsso,2,Eweight);
                if(EtaAsso>0 && PVtxz<0 && ChargeAsso > 0) Hist_Phi_FF->Fill(PhiAsso,1+2,Eweight);
                if(EtaAsso>0 && PVtxz<0 && ChargeAsso < 0) Hist_Phi_FF->Fill(PhiAsso,2+2,Eweight);
                if(EtaAsso<0 && PVtxz<0 && ChargeAsso > 0) Hist_Phi_RF->Fill(PhiAsso,1+2,Eweight);
                if(EtaAsso<0 && PVtxz<0 && ChargeAsso < 0) Hist_Phi_RF->Fill(PhiAsso,2+2,Eweight);
          }
          if(Ntof<2) continue;

	  int net_Np = Np - Npbar;
	  Hist_proton->Fill(Np);
	  Hist_pbar->Fill(Npbar);
	  Hist_netP->Fill(net_Np);
          int iTrack[Fcount], Scount = Fcount/2 -1;
          for(int q=0;q<Fcount;q++) iTrack[q] = q;
          random_shuffle(iTrack,iTrack+Fcount);
          Fcount = 0;
          for(int trk = 0; trk < NPTracks; trk++) {
                float EtaAsso   = leaf_EtaV0->GetValue(trk);
                float PtAsso    = leaf_PtV0->GetValue(trk);
                float PhiAsso   = leaf_PhiV0->GetValue(trk);
                float DCAglAsso = leaf_DCAglobalV0->GetValue(trk);
                float ChargeAsso= leaf_ChargeV0->GetValue(trk);
                float TOFflag   = leaf_TofFlagV0->GetValue(trk);
if(PtAsso > pt_asso_up || PtAsso < pt_asso_lo) continue;
if(DCAglAsso > DcaCut) continue;
if(EtaAsso>EtaCut || EtaAsso<-EtaCut) continue;
if(opt_TOF ==1 && TOFflag<1) continue;

                int n = (int)((PhiAsso+PI)/2./PI*Phibin);
                int fl = (PVtxz > 0)? ((ChargeAsso > 0)? 1:2):((ChargeAsso > 0)? 3:4);
                float W_phi = (EtaAsso>0)? PhiWgtFF[n][fl-1]:PhiWgtRF[n][fl-1];
                if(iTrack[Fcount] > Scount) {mQx1 +=W_phi*PtAsso*cos(PhiAsso*2.); mQy1 +=W_phi*PtAsso*sin(PhiAsso * 2.); Ecount++;}
                else {mQx2 += W_phi*PtAsso * cos(PhiAsso * 2.); mQy2 += W_phi*PtAsso * sin(PhiAsso * 2.); Wcount++;}
		if(EtaAsso> 0.1) {mQx3 +=W_phi*PtAsso*cos(PhiAsso*2.); mQy3 +=W_phi*PtAsso*sin(PhiAsso * 2.);}
		if(EtaAsso<-0.1) {mQx4 +=W_phi*PtAsso*cos(PhiAsso*2.); mQy4 +=W_phi*PtAsso*sin(PhiAsso * 2.);}
                Fcount++;
          }
          if(mQx1==0 || mQy1==0 || mQx2==0 || mQy2==0) continue;
	  if(mQx3==0 || mQy3==0 || mQx4==0 || mQy4==0) continue;
          mQ.Set(mQx, mQy); mQ1.Set(mQx1, mQy1); mQ2.Set(mQx2, mQy2); mQ3.Set(mQx3, mQy3); mQ4.Set(mQx4, mQy4);
          float TPC_EP_full = 0.5*mQ.Phi();
          float TPC_EP_east = 0.5*mQ1.Phi();
          float TPC_EP_west = 0.5*mQ2.Phi();
	  float TPC_EP_for  = 0.5*mQ3.Phi();
 	  float TPC_EP_bac  = 0.5*mQ4.Phi();
          Hist_TPC_EP_full->Fill(TPC_EP_full,Day);
          Hist_TPC_EP_east->Fill(TPC_EP_east,Day);
          Hist_TPC_EP_west->Fill(TPC_EP_west,Day);
          Hist_TPC_EP_for->Fill(TPC_EP_for,Day);
          Hist_TPC_EP_bac->Fill(TPC_EP_bac,Day);
	  pTPC_EP_east->Fill(1,Day,cos(2*TPC_EP_east)); pTPC_EP_east->Fill(2,Day,sin(2*TPC_EP_east)); 
	  pTPC_EP_east->Fill(3,Day,cos(4*TPC_EP_east)); pTPC_EP_east->Fill(4,Day,sin(4*TPC_EP_east));
          pTPC_EP_west->Fill(1,Day,cos(2*TPC_EP_west)); pTPC_EP_west->Fill(2,Day,sin(2*TPC_EP_west));
          pTPC_EP_west->Fill(3,Day,cos(4*TPC_EP_west)); pTPC_EP_west->Fill(4,Day,sin(4*TPC_EP_west));
          pTPC_EP_for->Fill(1,Day,cos(2*TPC_EP_for)); pTPC_EP_for->Fill(2,Day,sin(2*TPC_EP_for));
          pTPC_EP_for->Fill(3,Day,cos(4*TPC_EP_for)); pTPC_EP_for->Fill(4,Day,sin(4*TPC_EP_for));
          pTPC_EP_bac->Fill(1,Day,cos(2*TPC_EP_bac)); pTPC_EP_bac->Fill(2,Day,sin(2*TPC_EP_bac));
          pTPC_EP_bac->Fill(3,Day,cos(4*TPC_EP_bac)); pTPC_EP_bac->Fill(4,Day,sin(4*TPC_EP_bac));
	  if(fWgt->IsOpen() && Read_TPC_EP_east->GetEntries()) {
		PsiShiftE1 = Read_TPC_EP_east->GetBinContent(1,Day-79); PsiShiftE2 = Read_TPC_EP_east->GetBinContent(2,Day-79);
                PsiShiftE3 = Read_TPC_EP_east->GetBinContent(3,Day-79); PsiShiftE4 = Read_TPC_EP_east->GetBinContent(4,Day-79);
	  }
          if(fWgt->IsOpen() && Read_TPC_EP_west->GetEntries()) {
                PsiShiftW1 = Read_TPC_EP_west->GetBinContent(1,Day-79); PsiShiftW2 = Read_TPC_EP_west->GetBinContent(2,Day-79);
                PsiShiftW3 = Read_TPC_EP_west->GetBinContent(3,Day-79); PsiShiftW4 = Read_TPC_EP_west->GetBinContent(4,Day-79);
          }
          if(fWgt->IsOpen() && Read_TPC_EP_for->GetEntries()) {
                PsiShiftf1 = Read_TPC_EP_for->GetBinContent(1,Day-79); PsiShiftf2 = Read_TPC_EP_for->GetBinContent(2,Day-79);
                PsiShiftf3 = Read_TPC_EP_for->GetBinContent(3,Day-79); PsiShiftf4 = Read_TPC_EP_for->GetBinContent(4,Day-79);
          }
          if(fWgt->IsOpen() && Read_TPC_EP_bac->GetEntries()) {
                PsiShiftb1 = Read_TPC_EP_bac->GetBinContent(1,Day-79); PsiShiftb2 = Read_TPC_EP_bac->GetBinContent(2,Day-79);
                PsiShiftb3 = Read_TPC_EP_bac->GetBinContent(3,Day-79); PsiShiftb4 = Read_TPC_EP_bac->GetBinContent(4,Day-79);
          }
	  float TPC_EP_east_new = TPC_EP_east, TPC_EP_west_new = TPC_EP_west;
	  float TPC_EP_for_new = TPC_EP_for, TPC_EP_bac_new = TPC_EP_bac;
	  TPC_EP_east_new += 2*(-PsiShiftE2*cos(2*TPC_EP_east)+PsiShiftE1*sin(2*TPC_EP_east))/(float)2 + 2*(-PsiShiftE4*cos(4*TPC_EP_east)+PsiShiftE3*sin(4*TPC_EP_east))/(float)4;
          TPC_EP_west_new += 2*(-PsiShiftW2*cos(2*TPC_EP_west)+PsiShiftW1*sin(2*TPC_EP_west))/(float)2 + 2*(-PsiShiftW4*cos(4*TPC_EP_west)+PsiShiftW3*sin(4*TPC_EP_west))/(float)4;
          TPC_EP_for_new += 2*(-PsiShiftf2*cos(2*TPC_EP_for)+PsiShiftf1*sin(2*TPC_EP_for))/(float)2 + 2*(-PsiShiftf4*cos(4*TPC_EP_for)+PsiShiftf3*sin(4*TPC_EP_for))/(float)4;
          TPC_EP_bac_new += 2*(-PsiShiftb2*cos(2*TPC_EP_bac)+PsiShiftb1*sin(2*TPC_EP_bac))/(float)2 + 2*(-PsiShiftb4*cos(4*TPC_EP_bac)+PsiShiftb3*sin(4*TPC_EP_bac))/(float)4;
	  if(TPC_EP_east_new>PI) TPC_EP_east_new -= PI;
	  if(TPC_EP_east_new< 0) TPC_EP_east_new += PI;
          if(TPC_EP_west_new>PI) TPC_EP_west_new -= PI;
          if(TPC_EP_west_new< 0) TPC_EP_west_new += PI;
          if(TPC_EP_for_new>PI) TPC_EP_for_new -= PI;
          if(TPC_EP_for_new< 0) TPC_EP_for_new += PI;
          if(TPC_EP_bac_new>PI) TPC_EP_bac_new -= PI;
          if(TPC_EP_bac_new< 0) TPC_EP_bac_new += PI;
          Hist_TPC_EP_east_flat->Fill(TPC_EP_east_new,Day);
          Hist_TPC_EP_west_flat->Fill(TPC_EP_west_new,Day);
          Hist_TPC_EP_for_flat->Fill(TPC_EP_for_new,Day);
          Hist_TPC_EP_bac_flat->Fill(TPC_EP_bac_new,Day);
	  mQx = mQ1.Mod()*cos(2*TPC_EP_east_new) + mQ2.Mod()*cos(2*TPC_EP_west_new);
          mQy = mQ1.Mod()*sin(2*TPC_EP_east_new) + mQ2.Mod()*sin(2*TPC_EP_west_new);
          if(mQx==0 || mQy==0) continue;
	  mQ.Set(mQx, mQy);
	  TPC_EP_full = 0.5*mQ.Phi();
	  pTPC_EP_full->Fill(1,Day,cos(2*TPC_EP_full)); pTPC_EP_full->Fill(2,Day,sin(2*TPC_EP_full));
          pTPC_EP_full->Fill(3,Day,cos(4*TPC_EP_full)); pTPC_EP_full->Fill(4,Day,sin(4*TPC_EP_full));
          if(fWgt->IsOpen() && Read_TPC_EP_full->GetEntries()) {
                PsiShiftF1 = Read_TPC_EP_full->GetBinContent(1,Day-79); PsiShiftF2 = Read_TPC_EP_full->GetBinContent(2,Day-79);
                PsiShiftF3 = Read_TPC_EP_full->GetBinContent(3,Day-79); PsiShiftF4 = Read_TPC_EP_full->GetBinContent(4,Day-79);
          }
          float TPC_EP_full_new = TPC_EP_full;
          TPC_EP_full_new += 2*(-PsiShiftF2*cos(2*TPC_EP_full)+PsiShiftF1*sin(2*TPC_EP_full))/(float)2 + 2*(-PsiShiftF4*cos(4*TPC_EP_full)+PsiShiftF3*sin(4*TPC_EP_full))/(float)4;
          if(TPC_EP_full_new>PI) TPC_EP_full_new -= PI;
          if(TPC_EP_full_new< 0) TPC_EP_full_new += PI;
          Hist_TPC_EP_full_flat->Fill(TPC_EP_full_new,Day);
	  mQx = mQ.Mod()*cos(2*TPC_EP_full_new);
	  mQy = mQ.Mod()*sin(2*TPC_EP_full_new);
          Hist_cos->Fill(1,cos(2.*TPC_EP_for_new-2.*TPC_EP_bac_new), Eweight);
	  Hist_cos->Fill(2,cos(2.*TPC_EP_east_new-2.*TPC_EP_west_new), Eweight);
          Hist_cos_RefMult->Fill(1,RefMult,cos(2.*TPC_EP_for_new-2.*TPC_EP_bac_new), Eweight);
          Hist_cos_RefMult->Fill(2,RefMult,cos(2.*TPC_EP_east_new-2.*TPC_EP_west_new), Eweight);
          Hist_cos_TOFMult->Fill(1,TOFMult,cos(2.*TPC_EP_for_new-2.*TPC_EP_bac_new), Eweight);
          Hist_cos_TOFMult->Fill(2,TOFMult,cos(2.*TPC_EP_east_new-2.*TPC_EP_west_new), Eweight);
          Hist_cos_ZDC->Fill(1,ZDCe+ZDCw,cos(2.*TPC_EP_for_new-2.*TPC_EP_bac_new), Eweight);
          Hist_cos_ZDC->Fill(2,ZDCe+ZDCw,cos(2.*TPC_EP_east_new-2.*TPC_EP_west_new), Eweight);
	  Hist_cos_BBC_RefMult->Fill(BBCco,RefMult,cos(2.*TPC_EP_for_new-2.*TPC_EP_bac_new), Eweight);
          Hist_dif_count->Fill(Ecount - Wcount);
          Hist_ful_count->Fill(Ecount + Wcount);
          TPC_Day3_cos2->Fill(Day3, cos(2*TPC_EP_full));
          TPC_Day3_sin2->Fill(Day3, sin(2*TPC_EP_full));
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//loop through matched primary tracks
		for(int trki = 0; trki < NPTracks; trki++){
			Pt	  = leaf_PtV0->GetValue(trki);
			Eta	  = leaf_EtaV0->GetValue(trki);
                        P         = Pt*cosh(Eta);
			Theta     = 2.*atan(exp(-Eta));
			Charge	  = leaf_ChargeV0->GetValue(trki);
			Phi	  = leaf_PhiV0->GetValue(trki);
			ndEdx	  = leaf_ndEdxV0->GetValue(trki);
			DCAGlobal = leaf_DCAglobalV0->GetValue(trki);
                        nSigma_p  = leaf_nSigmaPV0->GetValue(trki);
			nSigma_pi = leaf_nSigmaPiV0->GetValue(trki);		
                	float TOFflag   = leaf_TofFlagV0->GetValue(trki);
			float En  = sqrt(mpi*mpi+pow(Pt*cosh(Eta),2));
                        float eff = PP0[cen-1]*exp(-pow(PP1[cen-1]/Pt,PP2[cen-1]));
			float eff_tof = TOF_eff->GetBinContent(TOF_eff->FindBin(Pt));
			eff *= eff_tof;
                        float TofYLocal = leaf_TofYLocalV0->GetValue(trki);
                        float TofM2     = leaf_TofM2V0->GetValue(trki);

                        if(ndEdx<15||nSigma_pi<-2 || nSigma_pi>2) continue;
			if(TOFflag<1||TofYLocal<-1.8||TofYLocal>1.8||TofM2<-0.01 || TofM2>0.1) continue;

			float mQx_i = mQx, mQy_i = mQy;
	                int n = (int)((Phi+PI)/2./PI*Phibin);
        	        int fl = (PVtxz > 0)? ((Charge > 0)? 1:2):((Charge > 0)? 3:4);
                	float W_phi = (Eta>0)? PhiWgtFF[n][fl-1]:PhiWgtRF[n][fl-1];
                	if(Pt > pt_asso_lo && Pt < pt_asso_up && Eta < EtaCut && Eta > -EtaCut && DCAGlobal < DcaCut) {
				mQx_i -= W_phi*Pt * cos(Phi * 2.);
                		mQy_i -= W_phi*Pt * sin(Phi * 2.);
			}
			TVector2 mQ_i(mQx_i,mQy_i);
			float psi_F = 0.5*mQ_i.Phi();
			float psi_sub = TPC_EP_for_new;
			if(Eta>0) psi_sub = TPC_EP_bac_new;
			float v2 = cos(2*Phi - 2*psi_sub)*100;
			Hist_v2_pt_obs->Fill(Pt,v2, Eweight);
			if(Eta > -1 && Eta <1 && Pt > pt_trig_lo && Pt < pt_trig_up) {
				p_v2_RefMult_obs->Fill(RefMult,v2, Eweight);
				p_v2_BBC_RefMult_obs->Fill(BBCco,RefMult,v2, Eweight);
                                p_v2_TOFMult_obs->Fill(TOFMult,v2, Eweight);
                                p_v2_ZDC_obs->Fill(ZDCe+ZDCw,v2, Eweight);
			}
			if(Charge>0 && Eta > -1 && Eta <1 && ndEdx>10 && nSigma_pi>-2 && nSigma_pi<2) Hist_v2_pt_pos_obs->Fill(Pt,v2,Eweight);
                        if(Charge<0 && Eta > -1 && Eta <1 && ndEdx>10 && nSigma_pi>-2 && nSigma_pi<2) Hist_v2_pt_neg_obs->Fill(Pt,v2,Eweight);

                        Hist_Pt->Fill(Pt,Eweight);
                        if(DCAGlobal > DcaCut) continue;
                        hEtaPtDist->Fill(Eta,Pt,Eweight);

                        if(Pt < pt_trig_lo || Pt > pt_trig_up) continue;
                        if(Eta > EtaCut || Eta < -EtaCut) continue;
                        StLorentzVectorD FourP(Pt*cos(Phi), Pt*sin(Phi), Pt*sinh(Eta), sqrt(mpi*mpi+P*P));

//Corrections
	                if(fWgt->IsOpen() && (TPCmean_FF->GetEntries() || TPCmean_RF->GetEntries())) {
        	                for(int k=1;k<5;k++) {
					for(int kk=0;kk<order;kk++) {
                	                  if(Eta>0 && PVtxz>0) PhiMean[k-1+4*kk] = TPCmean_FF->GetBinContent(k+8*kk,Day2-7999);
                        	          if(Eta<0 && PVtxz>0) PhiMean[k-1+4*kk] = TPCmean_RF->GetBinContent(k+8*kk,Day2-7999);
                                          if(Eta>0 && PVtxz<0) PhiMean[k-1+4*kk] = TPCmean_FF->GetBinContent(k+4+8*kk,Day2-7999);
                                          if(Eta<0 && PVtxz<0) PhiMean[k-1+4*kk] = TPCmean_RF->GetBinContent(k+4+8*kk,Day2-7999);
					}
                        	}
	                }

			float cos1 =0, cos2=0, sin1=0, sin2=0, Phi_new = Phi;
			// Recentering parameters
                        if(Eta>0 && PVtxz>0) {
                        if(Charge > 0) {
					for(int kk=0;kk<order;kk++) {
						pTPCmeanPhi_FF->Fill(1+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
						pTPCmeanPhi_FF->Fill(2+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
					}
				        }
                        if(Charge < 0) {
                                        for(int kk=0;kk<order;kk++) {
                                                pTPCmeanPhi_FF->Fill(3+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
                                                pTPCmeanPhi_FF->Fill(4+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
                                        }
					}
                        }
                        if(Eta<0 && PVtxz>0) {
                        if(Charge > 0) {
                                        for(int kk=0;kk<order;kk++) {
                                                pTPCmeanPhi_RF->Fill(1+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
                                                pTPCmeanPhi_RF->Fill(2+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
                                        }
					}
                        if(Charge < 0) {
                                        for(int kk=0;kk<order;kk++) {
                                                pTPCmeanPhi_RF->Fill(3+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
                                                pTPCmeanPhi_RF->Fill(4+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
                                        }
					}
                        }
                        if(Eta>0 && PVtxz<0) {
                        if(Charge > 0) {
                                        for(int kk=0;kk<order;kk++) {
                                                pTPCmeanPhi_FF->Fill(1+4+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
                                                pTPCmeanPhi_FF->Fill(2+4+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
                                        }
					}
                        if(Charge < 0) {
                                        for(int kk=0;kk<order;kk++) {
                                                pTPCmeanPhi_FF->Fill(3+4+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
                                                pTPCmeanPhi_FF->Fill(4+4+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
                                        }
					}
                        }
                        if(Eta<0 && PVtxz<0) {
                        if(Charge > 0) {
                                        for(int kk=0;kk<order;kk++) {
                                                pTPCmeanPhi_RF->Fill(1+4+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
                                                pTPCmeanPhi_RF->Fill(2+4+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
                                        }
					}
                        if(Charge < 0) {
                                        for(int kk=0;kk<order;kk++) {
                                                pTPCmeanPhi_RF->Fill(3+4+8*kk,Day2,cos(kk*Phi+Phi),Eweight);
                                                pTPCmeanPhi_RF->Fill(4+4+8*kk,Day2,sin(kk*Phi+Phi),Eweight);
                                        }
					}
                        }

			if(Charge > 0) {cos1 = cos(Phi) - PhiMean[0]; sin1 = sin(Phi) - PhiMean[1];
                                        if(order>1) {
						double a2np = 1+PhiMean[0+4], a2nn = 1-PhiMean[0+4];
						double lambda_2nsp = PhiMean[1+4]/a2np;
						double lambda_2nsn = PhiMean[1+4]/a2nn;
						double cos11 = (cos1 - lambda_2nsn*sin1)/(1-lambda_2nsn*lambda_2nsp);
						double sin11 = (sin1 - lambda_2nsp*cos1)/(1-lambda_2nsn*lambda_2nsp);
						cos1 = cos11/a2np;
						sin1 = sin11/a2nn;
					}
					for(int jj=0;jj<order;jj++) Phi_new += -2*PhiMean[1+4*jj]*cos(jj*Phi+Phi)/int(jj+1) 
									       +2*PhiMean[0+4*jj]*sin(jj*Phi+Phi)/int(jj+1); 
					if(Phi_new> PI) Phi_new -= 2*PI;
                			if(Phi_new<-PI) Phi_new += 2*PI;
                                        if(Eta>0 && PVtxz>0) Hist_Phi_FF_new->Fill(Phi_new,1,Eweight);
                                        if(Eta<0 && PVtxz>0) Hist_Phi_RF_new->Fill(Phi_new,1,Eweight);
                                        if(Eta>0 && PVtxz<0) Hist_Phi_FF_new->Fill(Phi_new,1+2,Eweight);
                                        if(Eta<0 && PVtxz<0) Hist_Phi_RF_new->Fill(Phi_new,1+2,Eweight);}
                        if(Charge < 0) {cos1 = cos(Phi) - PhiMean[2]; sin1 = sin(Phi) - PhiMean[3];
                                        if(order>1) {
                                        	double a2np = 1+PhiMean[2+4], a2nn = 1-PhiMean[2+4];
                                        	double lambda_2nsp = PhiMean[3+4]/a2np;
                                        	double lambda_2nsn = PhiMean[3+4]/a2nn;
                                        	double cos11 = (cos1 - lambda_2nsn*sin1)/(1-lambda_2nsn*lambda_2nsp);
                                        	double sin11 = (sin1 - lambda_2nsp*cos1)/(1-lambda_2nsn*lambda_2nsp);
                                        	cos1 = cos11/a2np;
                                        	sin1 = sin11/a2nn;
					}
                                        for(int jj=0;jj<order;jj++) Phi_new += -2*PhiMean[3+4*jj]*cos(jj*Phi+Phi)/int(jj+1) 
                                                                               +2*PhiMean[2+4*jj]*sin(jj*Phi+Phi)/int(jj+1);
                                        if(Phi_new> PI) Phi_new -= 2*PI;
                                        if(Phi_new<-PI) Phi_new += 2*PI;
                                        if(Eta>0 && PVtxz>0) Hist_Phi_FF_new->Fill(Phi_new,2,Eweight);
                                        if(Eta<0 && PVtxz>0) Hist_Phi_RF_new->Fill(Phi_new,2,Eweight);
                                        if(Eta>0 && PVtxz<0) Hist_Phi_FF_new->Fill(Phi_new,2+2,Eweight);
                                        if(Eta<0 && PVtxz<0) Hist_Phi_RF_new->Fill(Phi_new,2+2,Eweight);}
								
			for(int trkj = trki+1; trkj < NPTracks; trkj++) {
                        	float Pt2        = leaf_PtV0->GetValue(trkj);
                        	float Eta2       = leaf_EtaV0->GetValue(trkj);
                                double P2        = Pt2*cosh(Eta2);
                        	int Charge2      = leaf_ChargeV0->GetValue(trkj);
                        	float Phi2       = leaf_PhiV0->GetValue(trkj);
				float DCAGlobal2 = leaf_DCAglobalV0->GetValue(trkj);
				float ndEdx2     = leaf_ndEdxV0->GetValue(trkj);
				float nSigma_pi2 = leaf_nSigmaPiV0->GetValue(trkj);
				float TOFflag2   = leaf_TofFlagV0->GetValue(trkj);
				float En2	 = sqrt(mpi*mpi+pow(Pt2*cosh(Eta2),2));
                                float eff2 = PP0[cen-1]*exp(-pow(PP1[cen-1]/Pt2,PP2[cen-1]));
                        	float eff_tof2 = TOF_eff->GetBinContent(TOF_eff->FindBin(Pt2));
                        	eff2 *= eff_tof2;
	                        float TofYLocal2 = leaf_TofYLocalV0->GetValue(trkj);
        	                float TofM22     = leaf_TofM2V0->GetValue(trkj);

				if(ndEdx2<15||nSigma_pi2<-2 || nSigma_pi2>2) continue;	
        	                if(TOFflag2<1||TofYLocal2<-1.8||TofYLocal2>1.8||TofM22<-0.01 || TofM22>0.1) continue;

				if(Pt2 < pt_trig_lo || Pt2 > pt_trig_up) continue;
				if(Eta2 > EtaCut || Eta2 < -EtaCut) continue;
				if(DCAGlobal2 > DcaCut) continue;

                                StLorentzVectorD FourP2(Pt2*cos(Phi2), Pt2*sin(Phi2), Pt2*sinh(Eta2), sqrt(mpi*mpi+P2*P2));
                                float InvM = (FourP+FourP2).m();
                                if(Charge*Charge2<0) Hist_InvM_os->Fill(InvM);
                                if(Charge*Charge2>0) Hist_InvM_ss->Fill(InvM);

				float mQx_j = mQx_i, mQy_j = mQy_i;			
	                        int n2 = (int)((Phi2+PI)/2./PI*Phibin);
        	                int fl2 = (PVtxz > 0)? ((Charge2 > 0)? 1:2):((Charge2 > 0)? 3:4);
                	        float W_phi2 = (Eta2>0)? PhiWgtFF[n2][fl2-1]:PhiWgtRF[n2][fl2-1];
                        	mQx_j -= W_phi2*Pt2 * cos(Phi2 * 2.);
                        	mQy_j -= W_phi2*Pt2 * sin(Phi2 * 2.);
				TVector2 mQ_j(mQx_j, mQy_j);
				float psi_F_new = 0.5*mQ_j.Phi();

	                        float Delt_phi1 = Phi_new - psi_F_new;
        	                if(Delt_phi1>PI) Delt_phi1 -= 2*PI;
                	        if(Delt_phi1<-PI) Delt_phi1 += 2*PI;
                        	Delt_phi1 = (Delt_phi1>0)? psi_F_new + (gRandom->Rndm())*PI:psi_F_new - (gRandom->Rndm())*PI;

				hDpt->Fill(fabs(Pt-Pt2),Eweight);
				float q_inv = pow(Pt*sin(Phi)-Pt2*sin(Phi2),2)+pow(Pt*cos(Phi)-Pt2*cos(Phi2),2)+pow(Pt*sinh(Eta)-Pt2*sinh(Eta2),2)-pow(En-En2,2);
				hQinv2->Fill(q_inv,Eweight);
//				if(q_inv<0) continue;
//				q_inv = (q_inv>0)? sqrt(q_inv):0;
				hQinv->Fill(q_inv,Eweight);
//Corrections
	                        if(fWgt->IsOpen() && (TPCmean_FF->GetEntries() || TPCmean_RF->GetEntries())) {
        	                        for(int k=1;k<5;k++) {
                                          for(int kk=0;kk<order;kk++) {
                                            if(Eta>0 && PVtxz>0) PhiMean[k-1+4*kk] = TPCmean_FF->GetBinContent(k+8*kk,Day2-7999);
                                            if(Eta<0 && PVtxz>0) PhiMean[k-1+4*kk] = TPCmean_RF->GetBinContent(k+8*kk,Day2-7999);
                                            if(Eta>0 && PVtxz<0) PhiMean[k-1+4*kk] = TPCmean_FF->GetBinContent(k+4+8*kk,Day2-7999);
                                            if(Eta<0 && PVtxz<0) PhiMean[k-1+4*kk] = TPCmean_RF->GetBinContent(k+4+8*kk,Day2-7999);
                                          }
                                	}
                        	}

				// Recentering parameters
				float Phi2_new = Phi2;
				if(Charge2 > 0) {cos2 = cos(Phi2) - PhiMean[0]; sin2 = sin(Phi2) - PhiMean[1];
                                        	if(order>1) {
                                       		 	double a2np = 1+PhiMean[0+4], a2nn = 1-PhiMean[0+4];
                                        	 	double lambda_2nsp = PhiMean[1+4]/a2np;
                                        	 	double lambda_2nsn = PhiMean[1+4]/a2nn;
                                        	 	double cos22 = (cos2 - lambda_2nsn*sin2)/(1-lambda_2nsn*lambda_2nsp);
                                        	 	double sin22 = (sin2 - lambda_2nsp*cos2)/(1-lambda_2nsn*lambda_2nsp);
                                        	 	cos2 = cos22/a2np;
                                        	 	sin2 = sin22/a2nn;
						}
                                        	for(int jj=0;jj<order;jj++) 
							Phi2_new += -2*PhiMean[1+4*jj]*cos(jj*Phi2+Phi2)/int(jj+1)
                                                                    +2*PhiMean[0+4*jj]*sin(jj*Phi2+Phi2)/int(jj+1);
						}
                                if(Charge2 < 0) {cos2 = cos(Phi2) - PhiMean[2]; sin2 = sin(Phi2) - PhiMean[3];
	                                        if(order>1) {
                                        		double a2np = 1+PhiMean[2+4], a2nn = 1-PhiMean[2+4];
                                        		double lambda_2nsp = PhiMean[3+4]/a2np;
                                        		double lambda_2nsn = PhiMean[3+4]/a2nn;
                                        		double cos22 = (cos2 - lambda_2nsn*sin2)/(1-lambda_2nsn*lambda_2nsp);
                                        		double sin22 = (sin2 - lambda_2nsp*cos2)/(1-lambda_2nsn*lambda_2nsp);
                                        		cos2 = cos22/a2np;
                                        		sin2 = sin22/a2nn;
						}
                                        	for(int jj=0;jj<order;jj++) 
							Phi2_new += -2*PhiMean[3+4*jj]*cos(jj*Phi2+Phi2)/int(jj+1)
                                                                    +2*PhiMean[2+4*jj]*sin(jj*Phi2+Phi2)/int(jj+1);
						}
                        	float Delt_phi2 = Phi2_new - psi_F_new;
                        	if(Delt_phi2>PI) Delt_phi2 -= 2*PI;
                        	if(Delt_phi2<-PI) Delt_phi2 += 2*PI;
                        	Delt_phi2 = (Delt_phi2>0)? psi_F_new + (gRandom->Rndm())*PI:psi_F_new - (gRandom->Rndm())*PI;
				float correlator0 = cos(Phi + Phi2 - 2*psi_F_new);
				float correlator1 = (cos1*cos(Phi2)-sin1*sin(Phi2))*cos(2*psi_F_new)
                                                        + (sin1*cos(Phi2)+cos1*sin(Phi2))*sin(2*psi_F_new);
				float correlator2 = (cos1*cos2-sin1*sin2)*cos(2*psi_F_new)
							+ (sin1*cos2+cos1*sin2)*sin(2*psi_F_new);
				float correlator3 = cos(Phi_new + Phi2 - 2*psi_F_new);
				float correlator4 = cos(Phi_new + Phi2_new - 2*psi_F_new);
				float correlator5 = cos(Delt_phi1 + Delt_phi2 +PI - 2*psi_F_new);
				float correlator6 = cos(Phi_new-psi_F_new)*cos(Phi2_new-psi_F_new);
                                float correlator7 = sin(Phi_new-psi_F_new)*sin(Phi2_new-psi_F_new);
				if(Charge>0 && Charge2>0) {
                                        pParity_int_obs1->Fill(1,100*correlator0);
                                        pParity_int_obs2->Fill(1,100*correlator0,Eweight);
                                        pParity_int_obs3->Fill(1,100*correlator0,Eweight/eff/eff2);
                                        pParity_int_w_obs1->Fill(1,100*correlator0,W_phi);
                                        pParity_int_w_obs2->Fill(1,100*correlator0,Eweight*W_phi);
                                        pParity_int_w_obs3->Fill(1,100*correlator0,Eweight*W_phi/eff/eff2);
                                        pParity_int_ww_obs1->Fill(1,100*correlator0,W_phi*W_phi2);
                                        pParity_int_ww_obs2->Fill(1,100*correlator0,Eweight*W_phi*W_phi2);
                                        pParity_int_ww_obs3->Fill(1,100*correlator0,Eweight*W_phi*W_phi2/eff/eff2);
                                        pParity_int_r_obs1->Fill(1,100*correlator1);
                                        pParity_int_r_obs2->Fill(1,100*correlator1,Eweight);
                                        pParity_int_r_obs3->Fill(1,100*correlator1,Eweight/eff/eff2);
                                        pParity_int_rr_obs1->Fill(1,100*correlator2);
                                        pParity_int_rr_obs2->Fill(1,100*correlator2,Eweight);
                                        pParity_int_rr_obs3->Fill(1,100*correlator2,Eweight/eff/eff2);
                                        pParity_int_s_obs1->Fill(1,100*correlator3);
                                        pParity_int_s_obs2->Fill(1,100*correlator3,Eweight);
                                        pParity_int_s_obs3->Fill(1,100*correlator3,Eweight/eff/eff2);
                                        pParity_int_ss_obs1->Fill(1,100*correlator4);
                                        pParity_int_ss_obs2->Fill(1,100*correlator4,Eweight);
                                        pParity_int_ss_obs3->Fill(1,100*correlator4,Eweight/eff/eff2);
                                        pParity_int_ss_ran_obs1->Fill(1,100*correlator5);
                                        pParity_int_ss_ran_obs2->Fill(1,100*correlator5,Eweight);
                                        pParity_int_ss_ran_obs3->Fill(1,100*correlator5,Eweight/eff/eff2);
                                        pParity_eta_ss_obs1->Fill(1,0.5*(Eta+Eta2),100*correlator4);
                                        pParity_eta_ss_obs2->Fill(1,0.5*(Eta+Eta2),100*correlator4,Eweight);
                                        pParity_eta_ss_obs3->Fill(1,0.5*(Eta+Eta2),100*correlator4,Eweight/eff/eff2);
                                        pParity_Deta_ss_obs1->Fill(1,fabs(Eta-Eta2),100*correlator4);
                                        pParity_Deta_ss_obs2->Fill(1,fabs(Eta-Eta2),100*correlator4,Eweight);
                                        pParity_Deta_ss_obs3->Fill(1,fabs(Eta-Eta2),100*correlator4,Eweight/eff/eff2);
                                        pParity_pt_ss_obs1->Fill(1,0.5*(Pt+Pt2),100*correlator4);
                                        pParity_pt_ss_obs2->Fill(1,0.5*(Pt+Pt2),100*correlator4,Eweight);
                                        pParity_pt_ss_obs3->Fill(1,0.5*(Pt+Pt2),100*correlator4,Eweight/eff/eff2);
                                        pParity_Dpt_ss_obs1->Fill(1,fabs(Pt-Pt2),100*correlator4);
                                        pParity_Dpt_ss_obs2->Fill(1,fabs(Pt-Pt2),100*correlator4,Eweight);
                                        pParity_Dpt_ss_obs3->Fill(1,fabs(Pt-Pt2),100*correlator4,Eweight/eff/eff2);
                                        pParity_eta_ss_obs1->Fill(1+4,0.5*(Eta+Eta2),100*correlator6);
                                        pParity_eta_ss_obs2->Fill(1+4,0.5*(Eta+Eta2),100*correlator6,Eweight);
                                        pParity_eta_ss_obs3->Fill(1+4,0.5*(Eta+Eta2),100*correlator6,Eweight/eff/eff2);
                                        pParity_Deta_ss_obs1->Fill(1+4,fabs(Eta-Eta2),100*correlator6);
                                        pParity_Deta_ss_obs2->Fill(1+4,fabs(Eta-Eta2),100*correlator6,Eweight);
                                        pParity_Deta_ss_obs3->Fill(1+4,fabs(Eta-Eta2),100*correlator6,Eweight/eff/eff2);
                                        pParity_pt_ss_obs1->Fill(1+4,0.5*(Pt+Pt2),100*correlator6);
                                        pParity_pt_ss_obs2->Fill(1+4,0.5*(Pt+Pt2),100*correlator6,Eweight);
                                        pParity_pt_ss_obs3->Fill(1+4,0.5*(Pt+Pt2),100*correlator6,Eweight/eff/eff2);
                                        pParity_Dpt_ss_obs1->Fill(1+4,fabs(Pt-Pt2),100*correlator6);
                                        pParity_Dpt_ss_obs2->Fill(1+4,fabs(Pt-Pt2),100*correlator6,Eweight);
                                        pParity_Dpt_ss_obs3->Fill(1+4,fabs(Pt-Pt2),100*correlator6,Eweight/eff/eff2);
                                        pParity_eta_ss_obs1->Fill(1+8,0.5*(Eta+Eta2),100*correlator7);
                                        pParity_eta_ss_obs2->Fill(1+8,0.5*(Eta+Eta2),100*correlator7,Eweight);
                                        pParity_eta_ss_obs3->Fill(1+8,0.5*(Eta+Eta2),100*correlator7,Eweight/eff/eff2);
                                        pParity_Deta_ss_obs1->Fill(1+8,fabs(Eta-Eta2),100*correlator7);
                                        pParity_Deta_ss_obs2->Fill(1+8,fabs(Eta-Eta2),100*correlator7,Eweight);
                                        pParity_Deta_ss_obs3->Fill(1+8,fabs(Eta-Eta2),100*correlator7,Eweight/eff/eff2);
                                        pParity_pt_ss_obs1->Fill(1+8,0.5*(Pt+Pt2),100*correlator7);
                                        pParity_pt_ss_obs2->Fill(1+8,0.5*(Pt+Pt2),100*correlator7,Eweight);
                                        pParity_pt_ss_obs3->Fill(1+8,0.5*(Pt+Pt2),100*correlator7,Eweight/eff/eff2);
                                        pParity_Dpt_ss_obs1->Fill(1+8,fabs(Pt-Pt2),100*correlator7);
                                        pParity_Dpt_ss_obs2->Fill(1+8,fabs(Pt-Pt2),100*correlator7,Eweight);
                                        pParity_Dpt_ss_obs3->Fill(1+8,fabs(Pt-Pt2),100*correlator7,Eweight/eff/eff2);
                                        pParity_Q_ss_obs1->Fill(1,q_inv,100*correlator4);
                                        pParity_Q_ss_obs2->Fill(1,q_inv,100*correlator4,Eweight);
                                        pParity_Q_ss_obs3->Fill(1,q_inv,100*correlator4,Eweight/eff/eff2);
                                        if(fabs(Pt-Pt2)>0.15 && fabs(Eta-Eta2)>0.15) pParity_noHBT_ss_obs2->Fill(1,100*correlator4,Eweight);
                                        if(fabs(Pt-Pt2)>0.15 && fabs(Eta-Eta2)>0.15) pParity_noHBT_ss_obs3->Fill(1,100*correlator4,Eweight/eff/eff2);
				}
                                if(Charge<0 && Charge2<0) {
                                        pParity_int_obs1->Fill(2,100*correlator0);
                                        pParity_int_obs2->Fill(2,100*correlator0,Eweight);
                                        pParity_int_obs3->Fill(2,100*correlator0,Eweight/eff/eff2);
                                        pParity_int_w_obs1->Fill(2,100*correlator0,W_phi);
                                        pParity_int_w_obs2->Fill(2,100*correlator0,Eweight*W_phi);
                                        pParity_int_w_obs3->Fill(2,100*correlator0,Eweight*W_phi/eff/eff2);
                                        pParity_int_ww_obs1->Fill(2,100*correlator0,W_phi*W_phi2);
                                        pParity_int_ww_obs2->Fill(2,100*correlator0,Eweight*W_phi*W_phi2);
                                        pParity_int_ww_obs3->Fill(2,100*correlator0,Eweight*W_phi*W_phi2/eff/eff2);
                                        pParity_int_r_obs1->Fill(2,100*correlator1);
                                        pParity_int_r_obs2->Fill(2,100*correlator1,Eweight);
                                        pParity_int_r_obs3->Fill(2,100*correlator1,Eweight/eff/eff2);
                                        pParity_int_rr_obs1->Fill(2,100*correlator2);
                                        pParity_int_rr_obs2->Fill(2,100*correlator2,Eweight);
                                        pParity_int_rr_obs3->Fill(2,100*correlator2,Eweight/eff/eff2);
                                        pParity_int_s_obs1->Fill(2,100*correlator3);
                                        pParity_int_s_obs2->Fill(2,100*correlator3,Eweight);
                                        pParity_int_s_obs3->Fill(2,100*correlator3,Eweight/eff/eff2);
                                        pParity_int_ss_obs1->Fill(2,100*correlator4);
                                        pParity_int_ss_obs2->Fill(2,100*correlator4,Eweight);
                                        pParity_int_ss_obs3->Fill(2,100*correlator4,Eweight/eff/eff2);
                                        pParity_int_ss_ran_obs1->Fill(2,100*correlator5);
                                        pParity_int_ss_ran_obs2->Fill(2,100*correlator5,Eweight);
                                        pParity_int_ss_ran_obs3->Fill(2,100*correlator5,Eweight/eff/eff2);
                                        pParity_eta_ss_obs1->Fill(2,0.5*(Eta+Eta2),100*correlator4);
                                        pParity_eta_ss_obs2->Fill(2,0.5*(Eta+Eta2),100*correlator4,Eweight);
                                        pParity_eta_ss_obs3->Fill(2,0.5*(Eta+Eta2),100*correlator4,Eweight/eff/eff2);
                                        pParity_Deta_ss_obs1->Fill(2,fabs(Eta-Eta2),100*correlator4);
                                        pParity_Deta_ss_obs2->Fill(2,fabs(Eta-Eta2),100*correlator4,Eweight);
                                        pParity_Deta_ss_obs3->Fill(2,fabs(Eta-Eta2),100*correlator4,Eweight/eff/eff2);
                                        pParity_pt_ss_obs1->Fill(2,0.5*(Pt+Pt2),100*correlator4);
                                        pParity_pt_ss_obs2->Fill(2,0.5*(Pt+Pt2),100*correlator4,Eweight);
                                        pParity_pt_ss_obs3->Fill(2,0.5*(Pt+Pt2),100*correlator4,Eweight/eff/eff2);
                                        pParity_Dpt_ss_obs1->Fill(2,fabs(Pt-Pt2),100*correlator4);
                                        pParity_Dpt_ss_obs2->Fill(2,fabs(Pt-Pt2),100*correlator4,Eweight);
                                        pParity_Dpt_ss_obs3->Fill(2,fabs(Pt-Pt2),100*correlator4,Eweight/eff/eff2);
                                        pParity_eta_ss_obs1->Fill(2+4,0.5*(Eta+Eta2),100*correlator6);
                                        pParity_eta_ss_obs2->Fill(2+4,0.5*(Eta+Eta2),100*correlator6,Eweight);
                                        pParity_eta_ss_obs3->Fill(2+4,0.5*(Eta+Eta2),100*correlator6,Eweight/eff/eff2);
                                        pParity_Deta_ss_obs1->Fill(2+4,fabs(Eta-Eta2),100*correlator6);
                                        pParity_Deta_ss_obs2->Fill(2+4,fabs(Eta-Eta2),100*correlator6,Eweight);
                                        pParity_Deta_ss_obs3->Fill(2+4,fabs(Eta-Eta2),100*correlator6,Eweight/eff/eff2);
                                        pParity_pt_ss_obs1->Fill(2+4,0.5*(Pt+Pt2),100*correlator6);
                                        pParity_pt_ss_obs2->Fill(2+4,0.5*(Pt+Pt2),100*correlator6,Eweight);
                                        pParity_pt_ss_obs3->Fill(2+4,0.5*(Pt+Pt2),100*correlator6,Eweight/eff/eff2);
                                        pParity_Dpt_ss_obs1->Fill(2+4,fabs(Pt-Pt2),100*correlator6);
                                        pParity_Dpt_ss_obs2->Fill(2+4,fabs(Pt-Pt2),100*correlator6,Eweight);
                                        pParity_Dpt_ss_obs3->Fill(2+4,fabs(Pt-Pt2),100*correlator6,Eweight/eff/eff2);
                                        pParity_eta_ss_obs1->Fill(2+8,0.5*(Eta+Eta2),100*correlator7);
                                        pParity_eta_ss_obs2->Fill(2+8,0.5*(Eta+Eta2),100*correlator7,Eweight);
                                        pParity_eta_ss_obs3->Fill(2+8,0.5*(Eta+Eta2),100*correlator7,Eweight/eff/eff2);
                                        pParity_Deta_ss_obs1->Fill(2+8,fabs(Eta-Eta2),100*correlator7);
                                        pParity_Deta_ss_obs2->Fill(2+8,fabs(Eta-Eta2),100*correlator7,Eweight);
                                        pParity_Deta_ss_obs3->Fill(2+8,fabs(Eta-Eta2),100*correlator7,Eweight/eff/eff2);
                                        pParity_pt_ss_obs1->Fill(2+8,0.5*(Pt+Pt2),100*correlator7);
                                        pParity_pt_ss_obs2->Fill(2+8,0.5*(Pt+Pt2),100*correlator7,Eweight);
                                        pParity_pt_ss_obs3->Fill(2+8,0.5*(Pt+Pt2),100*correlator7,Eweight/eff/eff2);
                                        pParity_Dpt_ss_obs1->Fill(2+8,fabs(Pt-Pt2),100*correlator7);
                                        pParity_Dpt_ss_obs2->Fill(2+8,fabs(Pt-Pt2),100*correlator7,Eweight);
                                        pParity_Dpt_ss_obs3->Fill(2+8,fabs(Pt-Pt2),100*correlator7,Eweight/eff/eff2);
                                        pParity_Q_ss_obs1->Fill(2,q_inv,100*correlator4);
                                        pParity_Q_ss_obs2->Fill(2,q_inv,100*correlator4,Eweight);
                                        pParity_Q_ss_obs3->Fill(2,q_inv,100*correlator4,Eweight/eff/eff2);
                                        if(fabs(Pt-Pt2)>0.15 && fabs(Eta-Eta2)>0.15) pParity_noHBT_ss_obs2->Fill(2,100*correlator4,Eweight);
                                        if(fabs(Pt-Pt2)>0.15 && fabs(Eta-Eta2)>0.15) pParity_noHBT_ss_obs3->Fill(2,100*correlator4,Eweight/eff/eff2);
                                }
                                if(Charge*Charge2>0) {
                                        pParity_int_same_run->Fill(Day2,100*correlator0,Eweight);
                                        pParity_int_rr_same_run->Fill(Day2,100*correlator2,Eweight);
                                        pParity_int_ss_same_run->Fill(Day2,100*correlator4,Eweight);
                                        pParity_int_ww_same_run->Fill(Day2,100*correlator0,Eweight*W_phi*W_phi2);
                                        pParity_int_obs1->Fill(3,100*correlator0);
                                        pParity_int_obs2->Fill(3,100*correlator0,Eweight);
                                        pParity_int_obs3->Fill(3,100*correlator0,Eweight/eff/eff2);
                                        pParity_int_w_obs1->Fill(3,100*correlator0,W_phi);
                                        pParity_int_w_obs2->Fill(3,100*correlator0,Eweight*W_phi);
                                        pParity_int_w_obs3->Fill(3,100*correlator0,Eweight*W_phi/eff/eff2);
                                        pParity_int_ww_obs1->Fill(3,100*correlator0,W_phi*W_phi2);
                                        pParity_int_ww_obs2->Fill(3,100*correlator0,Eweight*W_phi*W_phi2);
                                        pParity_int_ww_obs3->Fill(3,100*correlator0,Eweight*W_phi*W_phi2/eff/eff2);
                                        pParity_int_r_obs1->Fill(3,100*correlator1);
                                        pParity_int_r_obs2->Fill(3,100*correlator1,Eweight);
                                        pParity_int_r_obs3->Fill(3,100*correlator1,Eweight/eff/eff2);
                                        pParity_int_rr_obs1->Fill(3,100*correlator2);
                                        pParity_int_rr_obs2->Fill(3,100*correlator2,Eweight);
                                        pParity_int_rr_obs3->Fill(3,100*correlator2,Eweight/eff/eff2);
                                        pParity_int_s_obs1->Fill(3,100*correlator3);
                                        pParity_int_s_obs2->Fill(3,100*correlator3,Eweight);
                                        pParity_int_s_obs3->Fill(3,100*correlator3,Eweight/eff/eff2);
                                        pParity_int_ss_obs1->Fill(3,100*correlator4);
                                        pParity_int_ss_obs2->Fill(3,100*correlator4,Eweight);
                                        pParity_int_ss_obs3->Fill(3,100*correlator4,Eweight/eff/eff2);
                                        pParity_int_ss_ran_obs1->Fill(3,100*correlator5);
                                        pParity_int_ss_ran_obs2->Fill(3,100*correlator5,Eweight);
                                        pParity_int_ss_ran_obs3->Fill(3,100*correlator5,Eweight/eff/eff2);
                                        pParity_eta_ss_obs1->Fill(3,0.5*(Eta+Eta2),100*correlator4);
                                        pParity_eta_ss_obs2->Fill(3,0.5*(Eta+Eta2),100*correlator4,Eweight);
                                        pParity_eta_ss_obs3->Fill(3,0.5*(Eta+Eta2),100*correlator4,Eweight/eff/eff2);
                                        pParity_Deta_ss_obs1->Fill(3,fabs(Eta-Eta2),100*correlator4);
                                        pParity_Deta_ss_obs2->Fill(3,fabs(Eta-Eta2),100*correlator4,Eweight);
                                        pParity_Deta_ss_obs3->Fill(3,fabs(Eta-Eta2),100*correlator4,Eweight/eff/eff2);
                                        pParity_pt_ss_obs1->Fill(3,0.5*(Pt+Pt2),100*correlator4);
                                        pParity_pt_ss_obs2->Fill(3,0.5*(Pt+Pt2),100*correlator4,Eweight);
                                        pParity_pt_ss_obs3->Fill(3,0.5*(Pt+Pt2),100*correlator4,Eweight/eff/eff2);
                                        pParity_Dpt_ss_obs1->Fill(3,fabs(Pt-Pt2),100*correlator4);
                                        pParity_Dpt_ss_obs2->Fill(3,fabs(Pt-Pt2),100*correlator4,Eweight);
                                        pParity_Dpt_ss_obs3->Fill(3,fabs(Pt-Pt2),100*correlator4,Eweight/eff/eff2);
                                        pParity_eta_ss_obs1->Fill(3+4,0.5*(Eta+Eta2),100*correlator6);
                                        pParity_eta_ss_obs2->Fill(3+4,0.5*(Eta+Eta2),100*correlator6,Eweight);
                                        pParity_eta_ss_obs3->Fill(3+4,0.5*(Eta+Eta2),100*correlator6,Eweight/eff/eff2);
                                        pParity_Deta_ss_obs1->Fill(3+4,fabs(Eta-Eta2),100*correlator6);
                                        pParity_Deta_ss_obs2->Fill(3+4,fabs(Eta-Eta2),100*correlator6,Eweight);
                                        pParity_Deta_ss_obs3->Fill(3+4,fabs(Eta-Eta2),100*correlator6,Eweight/eff/eff2);
                                        pParity_pt_ss_obs1->Fill(3+4,0.5*(Pt+Pt2),100*correlator6);
                                        pParity_pt_ss_obs2->Fill(3+4,0.5*(Pt+Pt2),100*correlator6,Eweight);
                                        pParity_pt_ss_obs3->Fill(3+4,0.5*(Pt+Pt2),100*correlator6,Eweight/eff/eff2);
                                        pParity_Dpt_ss_obs1->Fill(3+4,fabs(Pt-Pt2),100*correlator6);
                                        pParity_Dpt_ss_obs2->Fill(3+4,fabs(Pt-Pt2),100*correlator6,Eweight);
                                        pParity_Dpt_ss_obs3->Fill(3+4,fabs(Pt-Pt2),100*correlator6,Eweight/eff/eff2);
                                        pParity_eta_ss_obs1->Fill(3+8,0.5*(Eta+Eta2),100*correlator7);
                                        pParity_eta_ss_obs2->Fill(3+8,0.5*(Eta+Eta2),100*correlator7,Eweight);
                                        pParity_eta_ss_obs3->Fill(3+8,0.5*(Eta+Eta2),100*correlator7,Eweight/eff/eff2);
                                        pParity_Deta_ss_obs1->Fill(3+8,fabs(Eta-Eta2),100*correlator7);
                                        pParity_Deta_ss_obs2->Fill(3+8,fabs(Eta-Eta2),100*correlator7,Eweight);
                                        pParity_Deta_ss_obs3->Fill(3+8,fabs(Eta-Eta2),100*correlator7,Eweight/eff/eff2);
                                        pParity_pt_ss_obs1->Fill(3+8,0.5*(Pt+Pt2),100*correlator7);
                                        pParity_pt_ss_obs2->Fill(3+8,0.5*(Pt+Pt2),100*correlator7,Eweight);
                                        pParity_pt_ss_obs3->Fill(3+8,0.5*(Pt+Pt2),100*correlator7,Eweight/eff/eff2);
                                        pParity_Dpt_ss_obs1->Fill(3+8,fabs(Pt-Pt2),100*correlator7);
                                        pParity_Dpt_ss_obs2->Fill(3+8,fabs(Pt-Pt2),100*correlator7,Eweight);
                                        pParity_Dpt_ss_obs3->Fill(3+8,fabs(Pt-Pt2),100*correlator7,Eweight/eff/eff2);
                                        pParity_Q_ss_obs1->Fill(3,q_inv,100*correlator4);
                                        pParity_Q_ss_obs2->Fill(3,q_inv,100*correlator4,Eweight);
                                        pParity_Q_ss_obs3->Fill(3,q_inv,100*correlator4,Eweight/eff/eff2);
                                        if(fabs(Pt-Pt2)>0.15 && fabs(Eta-Eta2)>0.15) pParity_noHBT_ss_obs2->Fill(3,100*correlator4,Eweight);
                                        if(fabs(Pt-Pt2)>0.15 && fabs(Eta-Eta2)>0.15) pParity_noHBT_ss_obs3->Fill(3,100*correlator4,Eweight/eff/eff2);
                                }
                                if(Charge*Charge2<0) {
                                        pParity_int_oppo_run->Fill(Day2,100*correlator0,Eweight);
                                        pParity_int_rr_oppo_run->Fill(Day2,100*correlator2,Eweight);
                                        pParity_int_ss_oppo_run->Fill(Day2,100*correlator4,Eweight);
                                        pParity_int_ww_oppo_run->Fill(Day2,100*correlator0,Eweight*W_phi*W_phi2);
                                        pParity_int_obs1->Fill(4,100*correlator0);
                                        pParity_int_obs2->Fill(4,100*correlator0,Eweight);
                                        pParity_int_obs3->Fill(4,100*correlator0,Eweight/eff/eff2);
                                        pParity_int_w_obs1->Fill(4,100*correlator0,W_phi);
                                        pParity_int_w_obs2->Fill(4,100*correlator0,Eweight*W_phi);
                                        pParity_int_w_obs3->Fill(4,100*correlator0,Eweight*W_phi/eff/eff2);
                                        pParity_int_ww_obs1->Fill(4,100*correlator0,W_phi*W_phi2);
                                        pParity_int_ww_obs2->Fill(4,100*correlator0,Eweight*W_phi*W_phi2);
                                        pParity_int_ww_obs3->Fill(4,100*correlator0,Eweight*W_phi*W_phi2/eff/eff2);
                                        pParity_int_r_obs1->Fill(4,100*correlator1);
                                        pParity_int_r_obs2->Fill(4,100*correlator1,Eweight);
                                        pParity_int_r_obs3->Fill(4,100*correlator1,Eweight/eff/eff2);
                                        pParity_int_rr_obs1->Fill(4,100*correlator2);
                                        pParity_int_rr_obs2->Fill(4,100*correlator2,Eweight);
                                        pParity_int_rr_obs3->Fill(4,100*correlator2,Eweight/eff/eff2);
                                        pParity_int_s_obs1->Fill(4,100*correlator3);
                                        pParity_int_s_obs2->Fill(4,100*correlator3,Eweight);
                                        pParity_int_s_obs3->Fill(4,100*correlator3,Eweight/eff/eff2);
                                        pParity_int_ss_obs1->Fill(4,100*correlator4);
                                        pParity_int_ss_obs2->Fill(4,100*correlator4,Eweight);
                                        pParity_int_ss_obs3->Fill(4,100*correlator4,Eweight/eff/eff2);
                                        pParity_int_ss_ran_obs1->Fill(4,100*correlator5);
                                        pParity_int_ss_ran_obs2->Fill(4,100*correlator5,Eweight);
                                        pParity_int_ss_ran_obs3->Fill(4,100*correlator5,Eweight/eff/eff2);
                                        pParity_eta_ss_obs1->Fill(4,0.5*(Eta+Eta2),100*correlator4);
                                        pParity_eta_ss_obs2->Fill(4,0.5*(Eta+Eta2),100*correlator4,Eweight);
                                        pParity_eta_ss_obs3->Fill(4,0.5*(Eta+Eta2),100*correlator4,Eweight/eff/eff2);
                                        pParity_Deta_ss_obs1->Fill(4,fabs(Eta-Eta2),100*correlator4);
                                        pParity_Deta_ss_obs2->Fill(4,fabs(Eta-Eta2),100*correlator4,Eweight);
                                        pParity_Deta_ss_obs3->Fill(4,fabs(Eta-Eta2),100*correlator4,Eweight/eff/eff2);
                                        pParity_pt_ss_obs1->Fill(4,0.5*(Pt+Pt2),100*correlator4);
                                        pParity_pt_ss_obs2->Fill(4,0.5*(Pt+Pt2),100*correlator4,Eweight);
                                        pParity_pt_ss_obs3->Fill(4,0.5*(Pt+Pt2),100*correlator4,Eweight/eff/eff2);
                                        pParity_Dpt_ss_obs1->Fill(4,fabs(Pt-Pt2),100*correlator4);
                                        pParity_Dpt_ss_obs2->Fill(4,fabs(Pt-Pt2),100*correlator4,Eweight);
                                        pParity_Dpt_ss_obs3->Fill(4,fabs(Pt-Pt2),100*correlator4,Eweight/eff/eff2);
                                        pParity_eta_ss_obs1->Fill(4+4,0.5*(Eta+Eta2),100*correlator6);
                                        pParity_eta_ss_obs2->Fill(4+4,0.5*(Eta+Eta2),100*correlator6,Eweight);
                                        pParity_eta_ss_obs3->Fill(4+4,0.5*(Eta+Eta2),100*correlator6,Eweight/eff/eff2);
                                        pParity_Deta_ss_obs1->Fill(4+4,fabs(Eta-Eta2),100*correlator6);
                                        pParity_Deta_ss_obs2->Fill(4+4,fabs(Eta-Eta2),100*correlator6,Eweight);
                                        pParity_Deta_ss_obs3->Fill(4+4,fabs(Eta-Eta2),100*correlator6,Eweight/eff/eff2);
                                        pParity_pt_ss_obs1->Fill(4+4,0.5*(Pt+Pt2),100*correlator6);
                                        pParity_pt_ss_obs2->Fill(4+4,0.5*(Pt+Pt2),100*correlator6,Eweight);
                                        pParity_pt_ss_obs3->Fill(4+4,0.5*(Pt+Pt2),100*correlator6,Eweight/eff/eff2);
                                        pParity_Dpt_ss_obs1->Fill(4+4,fabs(Pt-Pt2),100*correlator6);
                                        pParity_Dpt_ss_obs2->Fill(4+4,fabs(Pt-Pt2),100*correlator6,Eweight);
                                        pParity_Dpt_ss_obs3->Fill(4+4,fabs(Pt-Pt2),100*correlator6,Eweight/eff/eff2);
                                        pParity_eta_ss_obs1->Fill(4+8,0.5*(Eta+Eta2),100*correlator7);
                                        pParity_eta_ss_obs2->Fill(4+8,0.5*(Eta+Eta2),100*correlator7,Eweight);
                                        pParity_eta_ss_obs3->Fill(4+8,0.5*(Eta+Eta2),100*correlator7,Eweight/eff/eff2);
                                        pParity_Deta_ss_obs1->Fill(4+8,fabs(Eta-Eta2),100*correlator7);
                                        pParity_Deta_ss_obs2->Fill(4+8,fabs(Eta-Eta2),100*correlator7,Eweight);
                                        pParity_Deta_ss_obs3->Fill(4+8,fabs(Eta-Eta2),100*correlator7,Eweight/eff/eff2);
                                        pParity_pt_ss_obs1->Fill(4+8,0.5*(Pt+Pt2),100*correlator7);
                                        pParity_pt_ss_obs2->Fill(4+8,0.5*(Pt+Pt2),100*correlator7,Eweight);
                                        pParity_pt_ss_obs3->Fill(4+8,0.5*(Pt+Pt2),100*correlator7,Eweight/eff/eff2);
                                        pParity_Dpt_ss_obs1->Fill(4+8,fabs(Pt-Pt2),100*correlator7);
                                        pParity_Dpt_ss_obs2->Fill(4+8,fabs(Pt-Pt2),100*correlator7,Eweight);
                                        pParity_Dpt_ss_obs3->Fill(4+8,fabs(Pt-Pt2),100*correlator7,Eweight/eff/eff2);
                                        pParity_Q_ss_obs1->Fill(4,q_inv,100*correlator4);
                                        pParity_Q_ss_obs2->Fill(4,q_inv,100*correlator4,Eweight);
                                        pParity_Q_ss_obs3->Fill(4,q_inv,100*correlator4,Eweight/eff/eff2);
                                        if(fabs(Pt-Pt2)>0.15 && fabs(Eta-Eta2)>0.15) pParity_noHBT_ss_obs2->Fill(4,100*correlator4,Eweight);
                                        if(fabs(Pt-Pt2)>0.15 && fabs(Eta-Eta2)>0.15) pParity_noHBT_ss_obs3->Fill(4,100*correlator4,Eweight/eff/eff2);
                                }
			} // 2nd track

		}  //Track

	} // Event



	fout.Write();
	return;
}
