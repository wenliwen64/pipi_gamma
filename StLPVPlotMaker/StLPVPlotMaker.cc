#include <iostream>
#include <string>
#include "TMath.h"
#include "TProfile.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "StLPVPlotMaker.hh"
ClassImp(StLPVPlotMaker)

StLPVPlotMaker::StLPVPlotMaker(){}

StLPVPlotMaker::~StLPVPlotMaker(){}

void StLPVPlotMaker::Init(int n_cent, int n_jobs){
    mNCent = n_cent;
    mNJobs = n_jobs;
}

void StLPVPlotMaker::Compute(){
    for(int i = 0; i < mNCent; i++){
	for(int j = 0; j < mNJobs; j++){
	    char filename[100]; 
	    sprintf(filename, "../Cen%d/cen%d.ntuple_result_Parity_pi_DCA2_minbias5_eff_job%d.root", i+1, i+1, j);
	    TFile* file_temp = new TFile(filename, "read");
            if(file_temp->IsOpen()){
                TProfile* Ep_diff_east_west = (TProfile*)file_temp->Get("Hist_cos;1");
                if(Ep_diff_east_west->GetSum() > -9999.){
                    double resolution_ini = Ep_diff_east_west->GetBinContent(2);// <cos(phi_east - phi_west)>
		    double resolution_ini_err = Ep_diff_east_west->GetBinError(2); // error

		    double chi_ini = chi(sqrt(resolution_ini)); // Find out the coresponding chi, and find out the full event plane chi by multipling sqrt(2)
		    double resolution_final = resEventPlane(sqrt(2.)*chi_ini);
		    std::cout << "Cent# " << i << " Job#" << j << " eventplane resolution == " << resolution_final << std::endl;

		    TProfile* Parity_int_ss_obs2 = (TProfile*)file_temp->Get("Parity_int_ss_obs2");
		    mGamma_ss[i][j] = Parity_int_ss_obs2->GetBinContent(3)/resolution_final/100;
		    mGamma_ss_err[i][j] = Parity_int_ss_obs2->GetBinError(3)/resolution_final/100;
		    mGamma_os[i][j] = Parity_int_ss_obs2->GetBinContent(4)/resolution_final/100;
		    mGamma_os_err[i][j] = Parity_int_ss_obs2->GetBinError(4)/resolution_final/100;
		}
	    }
            delete file_temp;
	}

        double weighted_gamma_os_final[9] = {0, };
        double weighted_gamma_ss_final[9] = {0, };
        double weights_os_final[9] = {0, };
        double weights_ss_final[9] = {0, };

        for(int j = 0; j < mNJobs; j++){
            weighted_gamma_ss_final[i] += mGamma_ss[i][j] / (mGamma_ss_err[i][j]*mGamma_ss_err[i][j]);
            weights_ss_final[i] += 1/(mGamma_ss_err[i][j]*mGamma_ss_err[i][j]);
            weighted_gamma_os_final[i] += mGamma_os[i][j] / (mGamma_os_err[i][j]*mGamma_os_err[i][j]);
            weights_os_final[i] += 1/(mGamma_os_err[i][j]*mGamma_os_err[i][j]);
	}
       
        mGamma_ss_final[i] = weighted_gamma_ss_final[i] / weights_ss_final[i];
        mGamma_ss_err_final[i] = sqrt(1 / weights_ss_final[i]);
        mGamma_os_final[i] = weighted_gamma_os_final[i] / weights_os_final[i];
        mGamma_os_err_final[i] = sqrt(1 / weights_os_final[i]);

        mGamma_diff_final[i] = mGamma_os_final[i] - mGamma_ss_final[i];
        mGamma_diff_err_final[i] = sqrt(mGamma_os_err_final[i]*mGamma_os_err_final[i] + mGamma_ss_err_final[i]*mGamma_ss_err_final[i]);

        std::cout << "Cent#" << i << "gamma => " << mGamma_diff_final[i] << " with err = " << mGamma_diff_err_final[i] << std::endl;
    }
}

double StLPVPlotMaker::resEventPlane(double chi){
  // Calculates the event plane resolution as a function of chi

    double con = 0.626657;                   // sqrt(pi/2)/2
    double arg = chi * chi / 4.;

    double res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));

    return res;
}

double StLPVPlotMaker::chi(double res){
// Calculates chi from the event plane resolution

    double chi = 2.0;
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

void StLPVPlotMaker::DrawSeparate(){
    double cent[9] = {75,65,55,45,35,25,15,7.5,2.5};
    const int reverseAxis = 1; 
    const float xScale = 80.;
    const int opt = 2;

    // Define the plot style
    gStyle->SetOptStat(0);
    gStyle->SetOptDate(0);
    gStyle->SetOptFit(0);

    if(reverseAxis==1) for(int i=0;i<9;i++) cent[i] = xScale-cent[i];

    // make the graph page
    int canvasWidth = 700, canvasHeight = 500;             // landscape
    TCanvas* can = new TCanvas("Flow_v", "", canvasWidth, canvasHeight);
    gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.02);
    gPad->SetBottomMargin(0.15);
    //gPad->SetTicks(1, 1);

    // make a histogram
    TString* histGraphName = new TString("Flow_v1");
    TH1F* histGraph = new TH1F(histGraphName->Data(), "", 24, 0, 80);
    histGraph->SetMaximum(0.002);
    histGraph->SetMinimum(-0.0012);
    histGraph->SetLineColor(kBlack);
    histGraph->SetXTitle("% Most Central");
    histGraph->SetYTitle("#LT cos(#phi_{1} + #phi_{2} - 2 #psi_{TPC}) #GT");
    histGraph->GetYaxis()->SetTitleOffset(1.3);
    histGraph->GetYaxis()->SetTitleSize(0.055);
    histGraph->GetXaxis()->SetTitleSize(0.055);
    histGraph->GetXaxis()->SetTitleOffset(1.);
    histGraph->GetXaxis()->SetNdivisions(505);
    double lsize = 0;
    if(reverseAxis==1){
	lsize = histGraph->GetLabelSize();
	histGraph->SetLabelSize(0.,"X");
    }
    histGraph->Draw();

    if(reverseAxis == 1){
	double ymin = -0.0012;
	TF1 *f1=new TF1("f1", "80.-x", 0., 80);
	TGaxis *A1=new TGaxis(0., ymin, xScale, ymin, "f1", 510, "+");
	A1->SetLabelSize(lsize);
	A1->SetLineColor(1);
	A1->Draw();
    }

    TGraphErrors* parity_pp = new TGraphErrors(mNCent, cent, mGamma_ss_final, 0, mGamma_ss_err_final);
    parity_pp->SetMarkerStyle(30);
    parity_pp->SetMarkerSize(1.5);
    parity_pp->SetMarkerColor(4);
    parity_pp->SetLineColor(4);
    if(opt == 2)  parity_pp->Draw("pe");

    TGraphErrors* parity_nn = new TGraphErrors(mNCent, cent, mGamma_os_final, 0, mGamma_os_err_final);
    parity_nn->SetMarkerStyle(29);
    parity_nn->SetMarkerSize(1.5);
    parity_nn->SetMarkerColor(2);
    parity_nn->SetLineColor(2);
    if(opt == 2)  parity_nn->Draw("pe");

    TLegend* legend = new TLegend(0.50, 0.70, 0.70, 0.80);
    legend->SetFillColor(0);
    legend->SetTextSize(0.07);
    legend->SetLineColor(0);
    legend->SetBorderSize(0);
    if(opt == 2)  legend->AddEntry(parity_pp, "#pi^{+}#pi^{+} and #pi^{-}#pi^{-}", "p");
    if(opt == 2)  legend->AddEntry(parity_nn, "#pi^{+}#pi^{-}", "p");
    legend->Draw();

    TLatex*   tex = new TLatex(35, 0.0008, "200 GeV Au+Au Period5");
    tex->SetTextSize(0.05);
    tex->SetTextColor(1);
    tex->Draw();
}

void StLPVPlotMaker::DrawDiff(){
    double cent[9] = {75,65,55,45,35,25,15,7.5,2.5};
    const int reverseAxis = 1; 
    const float xScale = 80.;
    const int opt = 2;

    // Define the plot style
    gStyle->SetOptStat(0);
    gStyle->SetOptDate(0);
    gStyle->SetOptFit(0);

    if(reverseAxis==1) for(int i=0;i<9;i++) cent[i] = xScale-cent[i];
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
    double lsize = 0;
    if (reverseAxis==1){
	lsize = histGraph->GetLabelSize();
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

    TGraphErrors* parity_diff = new TGraphErrors(mNCent, cent, 0, mGamma_diff_final, mGamma_diff_err_final);
    parity_diff->SetMarkerStyle(kOpenStar);
    parity_diff->SetMarkerSize(1.5);
    parity_diff->SetMarkerColor(4);
    parity_diff->SetLineColor(1);
    parity_diff->Draw("pe");

    TLegend* legend = new TLegend(0.40, 0.82, 0.80, 0.92);
    legend->SetFillColor(0);
    legend->SetTextSize(0.07);
    legend->SetLineColor(0);
    legend->SetBorderSize(0);
    legend->AddEntry(parity_diff, "Run2011 200 GeV Period5", "p");
    legend->Draw();

    TLatex *   tex = new TLatex(15,.04,"STAR preliminary");
    tex->SetTextSize(0.045);
    tex->SetTextColor(1);
    //tex->Draw();
}
