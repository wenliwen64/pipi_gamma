Double_t chi(double res);
Double_t resEventPlane(double chi);
void plot(){
    const int job_no = 10;
    float gamma_ss[job_no];
    float gamma_ss_err[job_no];
    float gamma_os[job_no];
    float gamma_os_err[job_no];
    for(int i = 0; i < 10; i++){
       char filename[100];
       sprintf(filename, "Cen4/cen4.ntuple_result_Parity_pi_DCA2_minbias5_eff_job%d.root", i);
       TFile* file_temp = new TFile(filename, "read");
       if(file_temp->IsOpen()){
	   TProfile* Ep_diff_east_west = (TProfile*)file_temp->Get("Hist_cos;1");

           if(Ep_diff_east_west->GetSum() > -9999.){
	       double resolution_ini = Ep_diff_east_west->GetBinContent(2);
	       double resolution_ini_err = Ep_diff_east_west->GetBinError(2);

	       double chi_ini = chi(sqrt(resolution_ini));
	       double resolution_final = resEventPlane(sqrt(2.)*chi_ini);
	       std::cout << "job# " << i <<" eventplane resolution == " << resolution_final << std::endl;

	       TProfile* Parity_int_ss_obs2 = (TProfile*)file_temp->Get("Parity_int_ss_obs2");
	       gamma_ss[i] = Parity_int_ss_obs2->GetBinContent(3)/resolution_final/100;
	       gamma_ss_err[i] = Parity_int_ss_obs2->GetBinError(3)/resolution_final/100;
	       gamma_os[i] = Parity_int_ss_obs2->GetBinContent(4)/resolution_final/100;
	       gamma_os_err[i] = Parity_int_ss_obs2->GetBinError(4)/resolution_final/100;
	   }
       }
    }
//Combine all of these 10 jobs results

    double cen1_pipi_gamma_ss = 0;
    double cen1_pipi_gamma_os = 0;
    double inverse_error_pipi_gamma_ss = 0;
    double inverse_error_pipi_gamma_os = 0;
    for(int i = 0; i < 10; i++){
        cen1_pipi_gamma_ss += gamma_ss[i] / (gamma_ss_err[i]*gamma_ss_err[i]);
        inverse_error_pipi_gamma_ss += 1/(gamma_ss_err[i]*gamma_ss_err[i]);
        cen1_pipi_gamma_os += gamma_os[i] / (gamma_os_err[i]*gamma_os_err[i]);
        inverse_error_pipi_gamma_os += 1/(gamma_os_err[i]*gamma_os_err[i]);
    }
   
    cout<< cen1_pipi_gamma_ss/inverse_error_pipi_gamma_ss << "--error = " << sqrt(1/inverse_error_pipi_gamma_ss) << " =========" << cen1_pipi_gamma_os/inverse_error_pipi_gamma_os << "--error = "<< sqrt(1/inverse_error_pipi_gamma_os) << endl;
/*
    TFile* file_la = new TFile("/Users/lwen/Documents/CVE_Project/root_file_example/cen1.ntuple_result_Parity_Lambda.eff.Period5.root", "read");
    cout<<"happy"<<endl;
    const int bin_no = 50;
    double x_la[bin_no] = {};
    double x_la_err[bin_no] = {};
    double gamma_la[bin_no] = {};  
    double gamma_la_err[bin_no] = {};  
    if(file_la->IsOpen()){
	TProfile* Ep_diff_east_west = (TProfile*)file_la->Get("Hist_cos;1");
	TProfile* Parity_int_ss_la = (TProfile*)file_la->Get("Parity_int_ss_la");
	if(Ep_diff_east_west->GetSum() > -99999.){
	    double resolution_ini = Ep_diff_east_west->GetBinContent(2);
	    double resolution_ini_err = Ep_diff_east_west->GetBinError(2);

	    double chi_ini = chi(sqrt(resolution_ini));
	    double resolution_final = resEventPlane(sqrt(2.)*chi_ini);
            std::cout << "eventplane resolution == " << resolution_final << std::endl;
	    for(int i = 0; i < bin_no; i++){
                x_la[i] = Parity_int_ss_la->GetXaxis().GetBinCenter(i);  //0.25+i*0.5;
                x_la_err[i] = 5./bin_no;//0.25;
		gamma_la[i] = Parity_int_ss_la->GetBinContent(i+1)/resolution_final/100;
		gamma_la_err[i] = Parity_int_ss_la->GetBinError(i+1)/resolution_final/100;
                std::cout << gamma_la[i] << "===" << gamma_la_err[i] << "===" << std::endl;
	    }
	}
    }
*/
    //TGraphErrors* gr_err = new TGraphErrors(bin_no,x_la, gamma_la, x_la_err, gamma_la_err);
    //gr_err->SetMarkerStyle(8);
    //gr_err->Draw("AP");
}

Double_t chi(double res) {
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

Double_t resEventPlane(double chi) {
  // Calculates the event plane resolution as a function of chi

  double con = 0.626657;                   // sqrt(pi/2)/2
  double arg = chi * chi / 4.;

  Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) +
                                          TMath::BesselI1(arg));

  return res;
}
