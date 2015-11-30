static Double_t chi(double res);
static Double_t resEventPlane(double chi);

void PlotParity() {
const int opt_tof =0;

char filename[100];
const int opt = 2;
const int reverseAxis=1;
const float xScale = 80.;
const int Ncen = 9;
float cent[9] = {75,65,55,45,35,25,15,7.5,2.5};
if(reverseAxis==1) for(int i=0;i<9;i++) cent[i] = xScale-cent[i];
float PP[9],NN[9],PN[9],SAM[9],Pos_parity[9],Neg_parity[9],resol[9];
float PP_err[9],NN_err[9],PN_err[9],SAM_err[9],Pos_parity_err[9],Neg_parity_err[9];
TF1 *fun = new TF1("fun","[0]");
for(int i=0;i<Ncen;i++) {
sprintf(filename,"./Cen%d/cen%d.ntuple_result_MSC_39_subEP.root",i+1,i+1);
//if(opt_tof) sprintf(filename,"cen%d.ntuple_result_Parity_UU_Vz30cm_P_TOF.root",i+1);
cout<<i+1<<"th file"<<endl;
TFile *f= new TFile(filename);

if(f->IsOpen()) {
	TProfile* Cos_ZZ_add = (TProfile*)f->Get("Hist_cos;1");
	TProfile* Parity_int_obs1_add = (TProfile*)f->Get("Parity_int_ss_obs1;1");
        TProfile* Parity_int_obs2_add = (TProfile*)f->Get("Parity_int_ss_obs2;1");

	if(Cos_ZZ_add->GetSum() > -99999) { 
		float res = Cos_ZZ_add->GetBinContent(2);
                float res_err = Cos_ZZ_add->GetBinError(2);
		float pp = Parity_int_obs2_add->GetBinContent(1)/100.;//lambda+pos
                float nn = Parity_int_obs2_add->GetBinContent(2)/100.;//lambda+neg
		float sm = Parity_int_obs2_add->GetBinContent(3)/100.;
                float pn = Parity_int_obs2_add->GetBinContent(4)/100.;
                float pp_err = Parity_int_obs1_add->GetBinError(1)/100.;
                float nn_err = Parity_int_obs1_add->GetBinError(2)/100.;
		float sm_err = Parity_int_obs1_add->GetBinError(3)/100.;
                float pn_err = Parity_int_obs1_add->GetBinError(4)/100.;
                float chi = chi(sqrt(res));
                float reso = resEventPlane(sqrt(2.)*chi);
                resol[i] = reso;
		PP[i] = pp/reso; NN[i] = nn/reso; PN[i] = pn/reso; SAM[i] = sm/reso;
		PP_err[i] = pp_err/reso;
                NN_err[i] = nn_err/reso;
                PN_err[i] = pn_err/reso;
                SAM_err[i] = sm_err/reso;

		cout<<PP[i]<<" "<<PP_err[i]<<" "<<NN[i]<<" "<<NN_err[i]<<endl;
	}
	Cos_ZZ_add->Reset();
	Parity_int_obs1_add->Reset();
        Parity_int_obs2_add->Reset();
}
f->Close();

}
/*
for(int i=0;i<Ncen;i++) cout<<resol[i]<<",";cout<<endl<<endl;
for(int i=0;i<Ncen;i++) cout<<SAM[i]<<",";cout<<endl;
for(int i=0;i<Ncen;i++) cout<<SAM_err[i]<<",";cout<<endl;
for(int i=0;i<Ncen;i++) cout<<PN[i]<<",";cout<<endl;
for(int i=0;i<Ncen;i++) cout<<PN_err[i]<<",";cout<<endl;
cout<<"the difference is "<<endl;
for(int i=0;i<Ncen;i++) cout<<PN[i]-SAM[i]<<",";cout<<endl;
*/

gStyle->SetOptStat(0);
gStyle->SetOptDate(0);
gStyle->SetOptFit(0);

  // make the graph page
  int canvasWidth = 700, canvasHeight = 500;             // landscape
  TCanvas* can = new TCanvas("Flow_v", "",canvasWidth, canvasHeight);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetFillColor(0);
  gPad->SetLeftMargin(0.17);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.02);
  gPad->SetBottomMargin(0.15);

  // make a histogram
  TString* histGraphName = new TString("Flow_v1");
  TH1F* histGraph = new TH1F(histGraphName->Data(), "", 24, 0, 80);
  histGraph->SetMaximum(0.01);
  histGraph->SetMinimum(-0.01);
  histGraph->SetLineColor(kBlack);
  histGraph->SetXTitle("% Most Central");
  histGraph->SetYTitle("#LT cos(#phi_{1} + #phi_{2} - 2 #psi_{TPC}) #GT");
  histGraph->GetYaxis()->SetTitleOffset(1.3);
  histGraph->GetYaxis()->SetTitleSize(0.055);
  histGraph->GetXaxis()->SetTitleSize(0.055);
  histGraph->GetXaxis()->SetTitleOffset(1.);
  histGraph->GetXaxis()->SetNdivisions(505);
  if (reverseAxis==1){
    double lsize=histGraph->GetLabelSize();
    histGraph->SetLabelSize(0.,"X");
  }
  histGraph->Draw();

  if (reverseAxis==1){
    double ymin=-0.001;
    TF1 *f1=new TF1("f1","80.-x",0.,80);
    TGaxis *A1=new TGaxis(0.,ymin,xScale,ymin,"f1",510,"+");
    A1->SetLabelSize(lsize);
    A1->SetLineColor(1);
    A1->Draw();

  }


  TGraphErrors* parity_pp = new TGraphErrors(Ncen,cent,PP,0,PP_err);
  parity_pp->SetMarkerStyle(30);
  parity_pp->SetMarkerSize(1.5);
  parity_pp->SetMarkerColor(4);
  parity_pp->SetLineColor(4);
if(opt == 2)  parity_pp->Draw("pe");
//  parity_pp->Print("all");
  TGraphErrors* parity_nn = new TGraphErrors(Ncen,cent,NN,0,NN_err);
  parity_nn->SetMarkerStyle(29);
  parity_nn->SetMarkerSize(1.5);
  parity_nn->SetMarkerColor(2);
  parity_nn->SetLineColor(2);
if(opt == 2)  parity_nn->Draw("pe");
/*
  TGraphErrors* parity_sam = new TGraphErrors(Ncen,cent,SAM,0,SAM_err);
  parity_sam->SetMarkerStyle(kFullCircle);
  parity_sam->SetMarkerSize(1.5);
  parity_sam->SetMarkerColor(6);
  parity_sam->SetLineColor(6);
if(opt == 1)  parity_sam->Draw("pe,C");
  TGraphErrors* parity_pn = new TGraphErrors(Ncen,cent,PN,0,PN_err);
  parity_pn->SetMarkerStyle(kOpenCircle);
  parity_pn->SetMarkerSize(1.5);
  parity_pn->SetMarkerColor(2);
  parity_pn->SetLineColor(1);
if(opt == 1)  parity_pn->Draw("pe.C");
*/
TLegend* legend = new TLegend(0.50, 0.20, 0.80, 0.40);
 legend->SetFillColor(0);
 legend->SetTextSize(0.07);
 legend->SetLineColor(0);
 legend->SetBorderSize(0.000001);
if(opt == 2)  legend->AddEntry(parity_pp,"#Lambda p","p");
if(opt == 2)  legend->AddEntry(parity_nn,"#Lambda #bar{p}","p");
//if(opt == 1)  legend->AddEntry(parity_sam,"(h^{+} h^{+} + h^{-} h^{-})/2","p");
//if(opt == 1)  legend->AddEntry(parity_pn,"(h^{+} h^{-} + h^{-} h^{+})/2","p");
  legend->Draw();

   TLatex *   tex = new TLatex(55,0.0008,"200 GeV Au+Au");
   tex->SetTextSize(0.05);
   tex->SetTextColor(1);
   tex->Draw();

   TLatex *   tex = new TLatex(1.5,.0008,"STAR preliminary");
   tex->SetTextSize(0.05);
   tex->SetTextColor(16);
//   tex->SetTextAngle(270);
   tex->Draw();

  // make the graph page
  int canvasWidth = 700, canvasHeight = 500;             // landscape
  TCanvas* can2 = new TCanvas("Flow_v2", "",canvasWidth, canvasHeight);
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetFillColor(0);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.02);
  gPad->SetBottomMargin(0.15);

  // make a histogram
  TString* histGraphName = new TString("Flow_v12");
  TH1F* histGraph = new TH1F(histGraphName->Data(), "", 24, 10, 80);
  histGraph->SetMaximum(1);
  histGraph->SetMinimum(0.);
  histGraph->SetLineColor(kBlack);
  histGraph->SetXTitle("% Most Central");
  histGraph->GetXaxis()->SetTitleOffset(0.9);
  histGraph->SetYTitle("EP resolution");
  histGraph->GetYaxis()->SetTitleOffset(1.1);
  histGraph->GetYaxis()->SetNdivisions(505);
    if (reverseAxis==1){
    double lsize=histGraph->GetLabelSize();
    histGraph->SetLabelSize(0.,"X");
  }
  histGraph->GetYaxis()->SetTitleSize(0.055);
  histGraph->GetXaxis()->SetTitleSize(0.055);
  histGraph->Draw();
  if (reverseAxis==1){
    double ymin=-0.;
    TF1 *f1=new TF1("f1","70.-x",0.,70);
    TGaxis *A1=new TGaxis(10.,ymin,80,ymin,"f1",510,"+");
    A1->SetLabelSize(lsize);
    A1->Draw();

  }


  TGraphErrors* parity_pos = new TGraphErrors(Ncen,cent,resol,0,0);
  parity_pos->SetMarkerStyle(kOpenStar);
  parity_pos->SetMarkerSize(1.5);
  parity_pos->SetMarkerColor(4);
  parity_pos->SetLineColor(4);
  parity_pos->Draw("pe");
 
TLegend* legend = new TLegend(0.40, 0.82, 0.80, 0.92);
 legend->SetFillColor(0);
 legend->SetTextSize(0.07);
 legend->SetLineColor(0);
 legend->SetBorderSize(0.000001);
 legend->AddEntry(parity_pos,"Run2011 200 GeV","p");
 legend->Draw();

   TLatex *   tex = new TLatex(15,.04,"STAR preliminary");
   tex->SetTextSize(0.045);
   tex->SetTextColor(1);
//   tex->SetTextAngle(270);
   tex->Draw();
}

//------------------------------------------------------------------------
static Double_t chi(double res) {
  // Calculates chi from the event plane resolution

  double chi   = 2.0;
  double delta = 1.0;

  for (int i = 0; i < 15; i++) {
//    chi   = (resEventPlane(chi) < res) ? chi + delta : chi - delta;
//    delta = delta / 2.;
      while(resEventPlane(chi) < res) {chi += delta;}
      delta = delta / 2.;
      while(resEventPlane(chi) > res) {chi -= delta;}
      delta = delta / 2.;
  }

  return chi;
}

//-----------------------------------------------------------------------

static Double_t resEventPlane(double chi) {
  // Calculates the event plane resolution as a function of chi

  double con = 0.626657;                   // sqrt(pi/2)/2
  double arg = chi * chi / 4.;

  Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) +
                                          TMath::BesselI1(arg));

  return res;
}
