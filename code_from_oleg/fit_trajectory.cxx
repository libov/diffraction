#include <iostream>
using namespace std;
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#ifndef __CINT__
//#include "Math/SpecFuncMathCore.h"
//#include "Math/SpecFuncMathMore.h"

#include "Math/SpecFunc.h"
#endif 

//#include <gsl/gsl_sf_hyperg.h>
//#include "Math/SpecFunc.h"
#include <TH1F.h>
#include <TDatime.h>
#include <TF1.h>
#include <TF2.h>
#include <TH2D.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFitter.h>
#include <TMinuit.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include "plot_style_utils.h"

Double_t dsigmael_dt(Double_t t, Double_t s);
Double_t dsigmael_dminust(Double_t *t, Double_t *s);
Double_t elastic_vertex(Double_t t);
Double_t alpha_t(Double_t t);

Double_t Re_alpha(Double_t s, Double_t* par, Double_t* lambda);
Double_t Im_alpha(Double_t s, Double_t* par, Double_t* lambda, Double_t& im);
Double_t Im_Amplitude(Double_t s, Double_t t, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm_Roper, Int_t include_bg_jenk,
		      Double_t par_jenk[]); 
Double_t DD_cross(Double_t s, Double_t t, Double_t M2, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper,
		  Int_t include_bg_jenk, Double_t par_jenk[]);

Double_t DD_cross_bRisto(Double_t s, Double_t t, Double_t M2, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper,
		  Int_t include_bg_jenk, Double_t par_jenk[]);

Double_t DD_cross_mod(Double_t s, Double_t t, Double_t M2, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper,
		  Int_t include_bg_jenk, Double_t par_jenk[]);

Double_t sigma_t_PpN(Double_t t, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm_Roper, Int_t include_bg_jenk,
		     Double_t par_jenk[], Double_t M2);

Double_t Theta_function(Double_t s1, Double_t s2);

Double_t dsigma3_dtdM1dM2(Double_t s, Double_t t, Double_t M2_1, Double_t M2_2, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper, Int_t include_bg_jenk, Double_t par_jenk[]);
Double_t dsigma3_dt1dt2dM(Double_t s1, Double_t s2, Double_t t1, Double_t t2, Double_t M2, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper, Int_t include_bg_jenk, Double_t par_jenk[]);
Double_t kinematic_factor(Double_t x, Double_t t, Double_t M2);
Double_t elastic_form_factor(Double_t t);
Double_t Risto_form_factor(Double_t t);


Double_t double_integral_DD_cross(Double_t s, Double_t* par, Double_t* lambda_par,
                                  Double_t t_min, Double_t t_max, Double_t M2_min, Double_t M2_max, Int_t n_int,
                                  Bool_t include_Roper, Int_t bg_parameter, Double_t c_cross, Double_t c_bg, Double_t b, Double_t norm_Roper,
                                  Int_t include_bg_jenk, Double_t par_jenk[]);

Double_t triple_integral_dsigma3_dt1dt2dM(Double_t s, Double_t t_min, Double_t t_max, Double_t M2_min, Double_t M2_max, Int_t n_int, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper, Int_t include_bg_jenk, Double_t par_jenk[]);

Double_t triple_integral_dsigma3_dt1dt2dM_copy(Double_t s, Double_t t_min, Double_t t_max, Double_t M2_min, Double_t M2_max, Int_t n_int, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper, Int_t include_bg_jenk, Double_t par_jenk[]);
Double_t triple_integral_dsigma3_dtdM1dM2_copy(Double_t s, Double_t t_min, Double_t t_max, Double_t M2_min, Double_t M2_max, Int_t n_int, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper, Int_t include_bg_jenk, Double_t par_jenk[]);

Double_t double_integral_DD_cross_mod(Double_t s, Double_t* par, Double_t* lambda_par,
                                  Double_t t_min, Double_t t_max, Double_t M2_min, Double_t M2_max, Int_t n_int,
                                  Bool_t include_Roper, Int_t bg_parameter, Double_t c_cross, Double_t c_bg, Double_t b, Double_t norm_Roper,
				      Int_t include_bg_jenk, Double_t par_jenk[]);

Double_t integral_M2_DD_cross_mod(Double_t s, Double_t t, Double_t* par, Double_t* lambda_par,
				  Double_t M2_min, Double_t M2_max, Int_t n_int, Bool_t include_Roper, Int_t bg_parameter,
				  Double_t c_cross, Double_t c_bg, Double_t b, Double_t norm_Roper, Int_t include_bg_jenk,
				  Double_t par_jenk[]);

Double_t integral_M2_DD_cross(Double_t s, Double_t t, Double_t* par, Double_t* lambda_par,
                              Double_t M2_min, Double_t M2_max, Int_t n_int, Bool_t include_Roper, Int_t bg_parameter,
                              Double_t c_cross, Double_t c_bg, Double_t b, Double_t norm_Roper, Int_t include_bg_jenk,
                              Double_t par_jenk[]);

Double_t integral_M2_DD_cross_bRisto(Double_t s, Double_t t, Double_t* par, Double_t* lambda_par,
                              Double_t M2_min, Double_t M2_max, Int_t n_int, Bool_t include_Roper, Int_t bg_parameter,
                              Double_t c_cross, Double_t c_bg, Double_t b, Double_t norm_Roper, Int_t include_bg_jenk,
                              Double_t par_jenk[]);


const Double_t s0 = 1;//GeV
const Int_t n_sum = 3;
Bool_t g_b_include_bg_jenk = kFALSE;//kTRUE;
Int_t g_include_bg_jenk = 0;//1;
Bool_t g_include_Roper = kFALSE;
Double_t global_param[7] = {
  -0.410839,
  -0.460021,
  0.508309,
  4.01083,
  4558.44,
  2.44,
  11.7
};
Double_t global_lambda_param[n_sum] = {
  0.840083,
  2.09707,
  11.1778
};
Double_t g_par_jenk[5];
//Double_t g_c_cross = 506.388;//150.131;//77.8957;//0.444576;
Double_t DD_cross_array[1000][1000] = {0.};
Double_t DD_cross_array_mod[1000][1000] = {0.};

Double_t ds3_dt1dt2dM_cross_array[100][100][100] = {0.};
Double_t ds3_dtdM1dM2_cross_array[100][100][100] = {0.};


const Int_t n_exp_data = 33;

Double_t g_c_cross = 977.441;
Double_t g_c_bg = 27.1566;
Double_t g_b = 10.;
Double_t g_norm_Roper = 0.072;
const Int_t global_n_int = 100;
//Double_t sigmatppn[global_n_int] = {0.};

//void calculate_sigmatppn();


int main()
{
Double_t C_cross[3];
C_cross[0] = g_c_cross;
C_cross[1] = g_c_bg;
C_cross[2] = g_b;

  g_par_jenk[0] = 9.91689;//5.16134e+06;//0.5 * 8.46365e+06;
  g_par_jenk[1] = 9.85517;//1.20674;//1.02355;//0.964596;
  g_par_jenk[2] = 0.493714;//0.409537;//0.206181;
  g_par_jenk[3] = -0.000665368;//-1.10369;//-1.77793;
  g_par_jenk[4] = 506.386;//180.427;//50. * 4.47663;
  cerr << "Program is starting" << endl;
  TF1* f_elastic = new TF1("f_elastic", dsigmael_dminust, 0, 1, 1);
  f_elastic->SetParameter(0, 7000);
  cerr << f_elastic->Eval(-0.1) << endl;
  const Int_t n_re_straight = 200;
  Double_t x_cross[n_re_straight], y_cross[n_re_straight];
  Double_t y_modified_cross[n_re_straight];
  Double_t y_bg[n_re_straight] = {0};
  Double_t t_draw1 = -0.05;
  TString s_t_draw1; s_t_draw1.Form("%.2f", t_draw1);
  for(Int_t i=0; i<n_re_straight; i++)
    {
      x_cross[i] = 9.*(i+0.5)/(double)n_re_straight;
      if(g_include_bg_jenk) {
	y_bg[i] = DD_cross(7.*7.e6, t_draw1, x_cross[i], global_param, global_lambda_param, g_include_Roper, g_c_cross, 0.036, -1, g_par_jenk);
      }
      y_cross[i] = y_bg[i] + DD_cross(7.*7.e6, t_draw1, x_cross[i], global_param, global_lambda_param, g_include_Roper, g_c_cross, 0.036, g_include_bg_jenk, g_par_jenk);
      y_modified_cross[i] = DD_cross_mod(7.*7.e6, t_draw1, x_cross[i], global_param, global_lambda_param, g_include_Roper, g_c_cross, 0.036, g_include_bg_jenk, g_par_jenk) *9.;
      cerr << "x = " << x_cross[i] << ", y = " << y_cross[i] << ", y_mod = " << y_modified_cross[i] << endl;
    }
  TGraph* gCross = new TGraph(n_re_straight, x_cross, y_cross);
  gCross->SetLineWidth(2);
  TGraph* gCrossMod = new TGraph(n_re_straight, x_cross, y_modified_cross);
  gCrossMod->SetLineWidth(2);
  gCrossMod->SetLineColor(kRed);  
  TGraph* gBg = new TGraph(n_re_straight, x_cross, y_bg);
  gBg->SetLineStyle(8);
  gBg->SetLineWidth(2);
  /////////////////////////
  //
  //
  //        DRAW!
  //
  //
  ////////////////////////

  TH2D* h_windowCross = new TH2D("h_windowCross", "", 10, 2., 8., 100, 0., 250.);
  TLegend *legBg2 = new TLegend( 0.4, 0.5, 0.93, 0.89);
  legBg2->SetBorderSize(0);
  legBg2->SetFillColor(0);
  legBg2->AddEntry(gCross, "#sqrt{s} = 7 TeV, t = "+s_t_draw1+" GeV^{2}", "l");
  TCanvas* cCross = new TCanvas("cCross", "cCross", 800, 800);
  cCross->Divide(1, 1);
  make_clean_pads(cCross, 1, 0, 0);
  cCross->Draw();  
  cCross->GetPad(1)->cd();
  sign_window(cCross->GetPad(1), h_windowCross, "M^{2}_{X} (GeV^{2})", "#frac{d^{2}#sigma}{dtdM^{2}_{X}}  #(){#frac{mb}{GeV^{4}}}", "");

  h_windowCross->Draw();
  gCross->Draw("C SAME");
  gCrossMod->Draw("C SAME");
  gBg->Draw("C SAME");
  legBg2->Draw();
  cCross->Print("SDD.eps");


  f_elastic->SetLineWidth(2);
  f_elastic->SetLineColor(kBlue);
  TH2D* h_windowElastic = new TH2D("h_windowElastic", "", 10, 0., 0.4, 100, 1.e-10, 1.2);
  TLegend *legElastic = new TLegend( 0.6, 0.79, 0.93, 0.89);
  legElastic->SetBorderSize(0);
  legElastic->SetFillColor(0);
  legElastic->AddEntry(f_elastic, "#sqrt{s} = 7 TeV", "l");
  TCanvas* cElastic = new TCanvas("cElastic", "cElastic", 800, 800);
  cElastic->Divide(1, 1);
  make_clean_pads(cElastic, 1, 0, 0);
  cElastic->Draw();
  cElastic->GetPad(1)->cd();
  //  cElastic->GetPad(1)->SetLogy();
  sign_window(cElastic->GetPad(1), h_windowElastic, "|t| (GeV)", "const #upoint d#sigma_{el}/dt (arb. units)", "");
  h_windowElastic->Draw();
  f_elastic->Draw("SAME");
  legElastic->Draw();
  cElastic->Print("elastic.eps");



  Double_t x_int_cross_graph[n_exp_data+11] = {
    1.e-3, 2.e-3, 4.e-3, 6.e-3, 9.e-3,
    0.014, 0.0162, 0.0176, 0.0191, 0.02,
    0.02, 0.0233, 0.0234, 0.0238, 0.0269,
    0.0272, 0.0274, 0.0305, 0.0323, 0.0324,
    0.0352, 0.0355, 0.0383, 0.0385, 0.0447,
    0.0537, 0.0623, 0.2, 0.546, 0.546,
    0.546, 0.546, 0.9, 1.8, 1.8,
    1.8, 1.8, 1.8,
    3., 4., 6., 9., 19.,
    30.
  };

  Double_t CD_int_cross[n_exp_data+11] = {0.};
  Double_t SD_int_cross[n_exp_data+11] = {0.};
  Double_t DD_int_cross[n_exp_data+11] = {0.};
  for(Int_t i=0; i<n_exp_data+11; i++)
    {
      //      x_int_cross[i] = 1.e-4 + 50.*(TMath::Exp(i+0.5))/TMath::Exp((double)n_re_straight_tot);
      cerr << i+1 << " from " << n_exp_data+11 << "... be patient :)" << endl;
      //      CD_int_cross[i] = triple_integral_dsigma3_dt1dt2dM_copy(x_int_cross_graph[i]*x_int_cross_graph[i]*1.e6,
      //							      -0.2, 0., 1.152, 1.5*1.5, 99,
      //							      global_param, global_lambda_param, g_include_Roper, 1., 0.036, 0, g_par_jenk);

      //      DD_int_cross[i] = triple_integral_dsigma3_dtdM1dM2_copy(x_int_cross_graph[i]*x_int_cross_graph[i]*1.e6,
      //							      -0.2, 0., 1.152, 8., 99,
      //							      global_param, global_lambda_param, g_include_Roper, 1., 0.036, 0, g_par_jenk);

//       SD_int_cross[i] =  double_integral_DD_cross_mod(x_int_cross_graph[i]*x_int_cross_graph[i]*1.e6,
// 						      global_param, 
// 						      global_lambda_param,
// 						      -0.2, 0.,
// 						      1.152, 8.,
// 						      99,
// 						      g_include_Roper, 
// 						      0,
// 						      C_cross[0], C_cross[1], C_cross[2], 
// 						      0.036, 0, g_par_jenk);

//       cout << "i = " << i << ", x = " << x_int_cross_graph[i]
// 	   << ", cross = " << SD_int_cross[i] << endl;

      //      y_int_cross[i] = double_integral_DD_cross(x_int_cross_graph[i]*x_int_cross_graph[i]*1.e6, 
      //						param, lambda_param, -0.2, 0., 1.152, 8., n_int, g_include_Roper, 1, C_cross[0], C_cross[1], C_cross[2], g_norm_Roper, 0, g_par_jenk);
      //      y_int_cross2[i] = double_integral_DD_cross(x_int_cross_graph[i]*x_int_cross_graph[i]*1.e6, param, lambda_param, -0.2, 0., 1.152, 8., n_int, g_include_Roper, 0, C_cross[0], C_cross[1], C_cross[2], g_norm_Roper, 0, g_par_jenk);//vorsicht!!!!
      //      y_int_cross3[i] = double_integral_DD_cross(x_int_cross_graph[i]*x_int_cross_graph[i]*1.e6, param, lambda_param, -0.2, 0., 1.152, 8., n_int, g_include_Roper, -1, C_cross[0], C_cross[1], C_cross[2], g_norm_Roper, 0, g_par_jenk);
      //      cerr << y_int_cross[i] << " " << y_int_cross2[i] << " " << y_int_cross3[i] << endl;
 //      cout << "integrated " << i<< " " << x_int_cross_data[i] << " " << y_int_cross[i] << endl;
    }


  Double_t checkDD[500] = {0.};
  Double_t x_checkDD[500] = {0.};
  Double_t M12_checkDD_fix = 3.;
  //  Double_t sigmaDD_check[100] = {0}.
  Double_t t_checkDD_fix = -0.05;
  for(Int_t i=0; i<500; i++) {
    x_checkDD[i] = (1.152 + (8.-1.152)*(i+1.)/500.)/M12_checkDD_fix;
    checkDD[i] = dsigma3_dtdM1dM2(7.*7.*1.e6, t_checkDD_fix, M12_checkDD_fix, x_checkDD[i]*M12_checkDD_fix, 
				  global_param, global_lambda_param, g_include_Roper, 1., 
				  0.036, 0, g_par_jenk);
    //    cout << "check DD " << i << " " << x_checkDD[i]*M12_checkDD_fix << " f(" << 
    //      x_checkDD[i] << ") = " << checkDD[i] << endl;

    cout << x_checkDD[i] << " " << checkDD[i] << endl;
    
  }


  //
  // SD and DD as funcion of t
  //

  Double_t x_cross_t[n_re_straight];
  Double_t x_cross_t_err[n_re_straight] = {0.};
  Double_t sd_cross_t[n_re_straight];
  Double_t dd_cross_t[n_re_straight];
  Double_t M2_draw = 1.5;
  TString s_M2_draw; s_M2_draw.Form("%.1f", M2_draw);
  for(Int_t i=0; i<n_re_straight; i++)
    {
      x_cross_t[i] = 0. + 1.5*(i+1)/((double)n_re_straight);
      sd_cross_t[i] = DD_cross(7.*7.e6, -x_cross_t[i], M2_draw, global_param, global_lambda_param, g_include_Roper, g_c_cross, 0.036, g_include_bg_jenk, g_par_jenk);
      dd_cross_t[i] = dsigma3_dtdM1dM2(7.*7.*1.e6, -x_cross_t[i], M2_draw, M2_draw, 
					global_param, global_lambda_param, g_include_Roper, 1., 
					0.036, 0, g_par_jenk);
      cout << "t = " << -x_cross_t[i] << " GeV, sd = " << sd_cross_t[i] << " mb/GeV^4, dd = " << dd_cross_t[i] << " mb/GeV^6" << endl;
    }

  TGraphErrors *gr_sd_cross_t = new TGraphErrors(n_re_straight, x_cross_t, sd_cross_t, x_cross_t_err, x_cross_t_err);

  TH2D* h_windowCross_t = new TH2D("h_windowCross_t", "", 10, 0., 1., 100, 1e-10, 1.);
  TLegend *legCross_t = new TLegend( 0.4, 0.5, 0.93, 0.89);
  legCross_t->SetBorderSize(0);
  legCross_t->SetFillColor(0);
  legCross_t->AddEntry(gr_sd_cross_t, "#sqrt{s} = 7 TeV, M^{2} = "+s_M2_draw+" GeV^{2}", "l");
  TCanvas* cCross_t = new TCanvas("cCross_t", "cCross_t", 800, 800);
  cCross_t->Divide(1, 1);
  make_clean_pads(cCross_t, 1, 0, 0);
  cCross_t->Draw();  
  cCross_t->GetPad(1)->cd();
  cCross_t->GetPad(1)->SetLogy();
  sign_window(cCross_t->GetPad(1), h_windowCross_t, "t (GeV^{2})", "#frac{d^{2}#sigma}{dtdM^{2}_{X}}  #(){#frac{mb}{GeV^{4}}}", "");

  h_windowCross_t->Draw();
  gr_sd_cross_t->Draw("C SAME");

  gBg->Draw("C SAME");
  legCross_t->Draw();
  cCross_t->Print("SDD_t.eps");


  //
  // SD and DD integrated over M2 as funcion of t
  //

  Double_t x_cross_intM2_t[n_re_straight];
  Double_t x_cross_intM2_t_err[n_re_straight] = {0.};
  Double_t sd_cross_mod_intM2_t[n_re_straight];
  Double_t sd_cross_intM2_t[n_re_straight];
  Double_t sd_cross_bRisto_intM2_t[n_re_straight];
  //  Double_t sd_cross_intM2_t[n_re_straight];
  Double_t dd_cross_intM2_t[n_re_straight];

  Int_t n_int_M2 = 50;
  for(Int_t i=0; i<n_re_straight; i++)
    {
      x_cross_intM2_t[i] = 0. + 1.5*(i+1)/((double)n_re_straight);
      sd_cross_mod_intM2_t[i] = integral_M2_DD_cross_mod(7.*7.e6, -x_cross_intM2_t[i], global_param, global_lambda_param, 1.152, 8., n_int_M2, g_include_Roper, 0, g_c_cross, g_c_bg, g_b, g_norm_Roper, g_include_bg_jenk, g_par_jenk);
      sd_cross_intM2_t[i] = integral_M2_DD_cross(7.*7.e6, -x_cross_intM2_t[i], global_param, global_lambda_param, 1.152, 8., n_int_M2, g_include_Roper, 0, g_c_cross, g_c_bg, g_b, g_norm_Roper, g_include_bg_jenk, g_par_jenk);
      sd_cross_bRisto_intM2_t[i] = integral_M2_DD_cross_bRisto(7.*7.e6, -x_cross_intM2_t[i], global_param, global_lambda_param, 1.152, 8., n_int_M2, g_include_Roper, 0, g_c_cross, g_c_bg, g_b, g_norm_Roper, g_include_bg_jenk, g_par_jenk);
      //      dd_cross_intM2_t[i] = dsigma3_dtdM1dM2(7.*7.*1.e6, -x_cross_intM2_t[i], 1.5, 1.5, 
      //					global_param, global_lambda_param, g_include_Roper, 1., 
      //					0.036, 0, g_par_jenk);
      cerr << "t = " << -x_cross_intM2_t[i] << " GeV, sd = " << sd_cross_intM2_t[i] << ", sd_mod = " << sd_cross_mod_intM2_t[i] << ", sd_bRisto = " << sd_cross_bRisto_intM2_t[i] << " mb/GeV^4" << endl;
    }

  TGraphErrors *gr_sd_cross_intM2_t = new TGraphErrors(n_re_straight, x_cross_intM2_t, sd_cross_intM2_t, x_cross_intM2_t_err, x_cross_intM2_t_err);
  TGraphErrors *gr_sd_cross_mod_intM2_t = new TGraphErrors(n_re_straight, x_cross_intM2_t, sd_cross_mod_intM2_t, x_cross_intM2_t_err, x_cross_intM2_t_err);
  TGraphErrors *gr_sd_cross_bRisto_intM2_t = new TGraphErrors(n_re_straight, x_cross_intM2_t, sd_cross_bRisto_intM2_t, x_cross_intM2_t_err, x_cross_intM2_t_err);

  gr_sd_cross_mod_intM2_t->SetLineColor(kBlue);
  gr_sd_cross_bRisto_intM2_t->SetLineColor(kRed);

  TH2D* h_windowCross_intM2_t = new TH2D("h_windowCross_intM2_t", "", 10, 0., 1., 100, 1e-4, 1.e3);
  TLegend *legCross_intM2_t = new TLegend( 0.4, 0.5, 0.93, 0.89);
  legCross_intM2_t->SetBorderSize(0);
  legCross_intM2_t->SetFillColor(0);
  legCross_intM2_t->AddEntry(gr_sd_cross_intM2_t, "#sqrt{s} = 7 TeV", "l");
  legCross_intM2_t->AddEntry(gr_sd_cross_mod_intM2_t, "#sqrt{s} = 7 TeV", "l");
  legCross_intM2_t->AddEntry(gr_sd_cross_bRisto_intM2_t, "#sqrt{s} = 7 TeV", "l");
  TCanvas* cCross_intM2_t = new TCanvas("cCross_intM2_t", "cCross_intM2_t", 800, 800);
  cCross_intM2_t->Divide(1, 1);
  make_clean_pads(cCross_intM2_t, 1, 0, 0);
  cCross_intM2_t->Draw();  
  cCross_intM2_t->GetPad(1)->cd();
  cCross_intM2_t->GetPad(1)->SetLogy();
  sign_window(cCross_intM2_t->GetPad(1), h_windowCross_intM2_t, "-t (GeV^{2})", "#frac{d^{2}#sigma}{dt}  #(){#frac{mb}{GeV^{2}}}", "");

  h_windowCross_intM2_t->Draw();
  gr_sd_cross_intM2_t->Draw("C SAME");
  gr_sd_cross_mod_intM2_t->Draw("C SAME");
  gr_sd_cross_bRisto_intM2_t->Draw("C SAME");

  gBg->Draw("C SAME");
  legCross_intM2_t->Draw();
  cCross_intM2_t->Print("SDD_intM2_t.eps");


  cerr << "Done" << endl;
  return 0;
}

Double_t dsigmael_dminust(Double_t *t, Double_t *s) {
  return dsigmael_dt(-t[0], s[0]);
}

Double_t dsigmael_dt(Double_t t, Double_t s)
{
  cout << "s = " << s << ", t = " << t << endl;
  Double_t res;
  Double_t elastic_vtx = elastic_vertex(t);
  res = TMath::Power(elastic_vtx, 4);
  res *= TMath::Power(s/s0, 2*alpha_t(t)-2);
  return res;
}

Double_t elastic_vertex(Double_t t)
{
  return TMath::Exp(5*t);
}

Double_t alpha_t(Double_t t)
{
  return 1.01+0.25*t;
}


Double_t Re_alpha(Double_t s, Double_t* par, Double_t* lambda_par)
{
  Double_t s_par[n_sum] = {0.};
  Double_t c1_par[n_sum] = {0.};
  //  Double_t lambda_par[n_sum] = {0.};
  Double_t alpha0 = par[0];
  Double_t delta = par[1];
  c1_par[0] = par[2];
  c1_par[1] = par[3];
  c1_par[2] = par[4];
  s_par[0] = 1.e-6 * TMath::Power(938.27 + 134.98, 2.);
  s_par[1] = par[5];
  s_par[2] = par[6];

  if(delta > 1) { cout << "Error: delta > 1 in Re_alpha()" << endl; exit(-1); }

  Double_t ReAlpha = alpha0;

  for(Int_t i=0; i<n_sum; i++)
    {      
      Double_t th1 = Theta_function(s_par[i], s);
      Double_t GGG = 0.;
      Double_t f21 = 0.;
      if(th1 != 0) {
	//	cerr << "calling tgamma(" << 1.-delta <<") = " << ROOT::Math::tgamma(1.-delta) << endl;
	//	cerr << "calling tgamma(" << lambda_par[i]+1.  <<") = " << ROOT::Math::tgamma(lambda_par[i]+1.) << endl;
	GGG = ROOT::Math::tgamma(1.-delta) * ROOT::Math::tgamma(lambda_par[i]+1.);
	//	cerr << "GGG 1 = " << GGG << endl;
	GGG /= ROOT::Math::tgamma(lambda_par[i]-delta+2.) * TMath::Power(s_par[i], 1. - delta);
	//	cerr << "GGG 2 = " << GGG << endl;
	if(s/s_par[i] != 1.) {//See Abramovitz Stegun p. 556 
	  //	  cerr << "calling ROOT::Math::hyperg(" << 1. << ", " << 1.-delta << ", " << lambda_par[i]-delta+2. << ", " << s/s_par[i] << "); c - a - b = " <<
	  //	    lambda_par[i] << endl;
	  if(lambda_par[i]-delta+2.<10)
	    f21 = ROOT::Math::hyperg(1., 1.-delta, lambda_par[i]-delta+2., s/s_par[i]);
	  else
	    f21 = ROOT::Math::hyperg(1., 1.-delta, 9.999, s/s_par[i]);
	  //	  cerr << "result is " << f21 << endl;//ROOT::Math::hyperg(1., 1.-delta, lambda_par[i]-delta+2., s/s_par[i]) << endl;
	}
	else {
	  f21 = ROOT::Math::tgamma(lambda_par[i]-delta+2.) *  ROOT::Math::tgamma(lambda_par[i]) ;
	  f21 /= ROOT::Math::tgamma(lambda_par[i]-delta+1.) *  ROOT::Math::tgamma(lambda_par[i] + 1.);
	}
      }
      
      Double_t th2 = Theta_function(s, s_par[i]);
      Double_t piscot = 0.;
      Double_t GGGsf21 = 0.; 
      if(th2 != 0) {
	piscot = TMath::Pi() * TMath::Power(s, delta - 1.);
	piscot *= TMath::Power( (s-s_par[i])/s, lambda_par[i]);
	piscot /= TMath::Tan( TMath::Pi()*(1.-delta) );

	GGGsf21 = ROOT::Math::tgamma(-1.*delta) * ROOT::Math::tgamma(lambda_par[i] + 1.);
	GGGsf21 *= TMath::Power(s_par[i], delta);
	GGGsf21 /= s;
	GGGsf21 /= ROOT::Math::tgamma(lambda_par[i]-delta+1.);

	if(s_par[i]/s != 1.) {//See Abramovitz Stegun p. 556
	  //	  cerr << "calling ROOT::Math::hyperg(" << delta-lambda_par[i] << ", " << 1. << ", " << delta+1. << ", " << s_par[i]/s << "); c - a - b = " <<
	  //	    lambda_par[i] << endl;
	  if(delta-lambda_par[i]>-10)
	    GGGsf21 *= ROOT::Math::hyperg(delta-lambda_par[i], 1., delta+1., s_par[i]/s);
	  else
	    GGGsf21 *= ROOT::Math::hyperg(-9.9999, 1., delta+1., s_par[i]/s);
	  
	  //	  cerr << "result is " << GGGsf21 << endl;//ROOT::Math::hyperg(delta-lambda_par[i], 1., delta+1., s_par[i]/s) << endl;
	  
	}
	else {
	  GGGsf21 *= ROOT::Math::tgamma(delta+1.) * ROOT::Math::tgamma(lambda_par[i]);
	  GGGsf21 /= ROOT::Math::tgamma(lambda_par[i] + 1.) * ROOT::Math::tgamma(delta);
	}	
      }

      Double_t A = GGG * f21 * th1 + (piscot - GGGsf21) * th2;
      ReAlpha += s * c1_par[i] * A / TMath::Pi();
      //      cout << "in ReAlpha: " << i << " " << GGG << " " << f21 << " " << th1 << " " << piscot << " " << GGGsf21 << " " << th2 << " " << A << " " << s * c1_par[i] * A / TMath::Pi() << endl;
    }
  
  return ReAlpha;
}


Double_t Im_alpha(Double_t s, Double_t* par, Double_t* lambda_par, Double_t& im)
{
  Double_t ImAlpha = 0.;
  Double_t ReAlphaPrime = 0.;
  Double_t G = 0.;

  Double_t s_par[n_sum] = {0.};
  Double_t c1_par[n_sum] = {0.};
  //  Double_t lambda_par[n_sum] = {0.};
  Double_t alpha0 = par[0];
  Double_t delta = par[1];
  c1_par[0] = par[2];
  c1_par[1] = par[3];
  c1_par[2] = par[4];
  s_par[0] = 1.e-6 * TMath::Power(938.27 + 134.98, 2.);
  s_par[1] = par[5];
  s_par[2] = par[6];

  if(delta > 1) { cout << "Error: delta > 1 in Re_alpha()" << endl; exit(-1); }

  Double_t ReAlpha = alpha0;

  for(Int_t i=0; i<n_sum; i++)
    {
      Double_t th1 = Theta_function(s_par[i], s);
      Double_t GGG = 0.;
      Double_t f21 = 0.;
      if(th1 != 0) {
	GGG = ROOT::Math::tgamma(1.-delta) * ROOT::Math::tgamma(lambda_par[i]+1.);
	GGG /= ROOT::Math::tgamma(lambda_par[i]-delta+2.) * TMath::Power(s_par[i], 1. - delta);

	if(s/s_par[i] != 1.) {//See Abramovitz Stegun p. 556 
	  if(lambda_par[i]-delta+2. < 10.)
	    f21 = ROOT::Math::hyperg(2., 1.-delta, lambda_par[i]-delta+2., s/s_par[i]);
	  else
	    f21 = ROOT::Math::hyperg(2., 1.-delta, 9.999, s/s_par[i]);
	}
	else {
	  f21 = ROOT::Math::tgamma(lambda_par[i]-delta+2.) *  ROOT::Math::tgamma(lambda_par[i]-1.) ;
	  f21 /= ROOT::Math::tgamma(lambda_par[i]-delta) *  ROOT::Math::tgamma(lambda_par[i] + 1.);
	}
      }

      Double_t th2 = Theta_function(s, s_par[i]);
      Double_t piscot = 0.;
      Double_t GGGsf21 = 0.; 
      if(th2 != 0) {
	piscot = TMath::Pi() * TMath::Power(s, delta - 1.);
	piscot *= TMath::Power( (s-s_par[i])/s, lambda_par[i]);
	piscot /= TMath::Tan( TMath::Pi()*(1.-delta) );
	piscot *= delta + lambda_par[i]*s_par[i]/(s-s_par[i]);

	GGGsf21 = ROOT::Math::tgamma(-1.*delta) * ROOT::Math::tgamma(lambda_par[i] + 1.);
	GGGsf21 *= TMath::Power(s_par[i], delta+1);
	GGGsf21 /= s*s;
	GGGsf21 /= 1. + delta;
	GGGsf21 /= ROOT::Math::tgamma(lambda_par[i]-delta);

	if(s_par[i]/s != 1.) {//See Abramovitz Stegun p. 556
	  if(1.+delta-lambda_par[i]>-10.)
	    GGGsf21 *= ROOT::Math::hyperg(1.+delta-lambda_par[i], 2., delta+2., s_par[i]/s);//ROOT::Math::hyperg(1.+delta-lambda_par[i], 2., delta+2., s_par[i]/s);
	  else GGGsf21 *= ROOT::Math::hyperg(-9.999, 2., delta+2., s_par[i]/s);
	}
	else {
	  GGGsf21 *= ROOT::Math::tgamma(delta+2.) * ROOT::Math::tgamma(lambda_par[i]-1.);
	  GGGsf21 /= ROOT::Math::tgamma(lambda_par[i] + 1.) * ROOT::Math::tgamma(delta);
	}	
      }

      Double_t B = GGG * f21 * th1 + (piscot - GGGsf21) * th2;
      ReAlphaPrime += c1_par[i] * B / TMath::Pi();
      if(Theta_function(s, s_par[i]) != 0) 
	ImAlpha += TMath::Power(s, delta) * c1_par[i] * TMath::Power(1. - s_par[i]/s, Re_alpha(s_par[i], par, lambda_par)) * Theta_function(s, s_par[i]);
      //      cout << "in ImAlpha: " << i << " " << GGG << " " << f21 << " " << th1 << " " << piscot << " " << GGGsf21 << " " << th2 << " " << B << " " 
      //	   << ReAlphaPrime << " " << ImAlpha << endl;
    }  

  if(ReAlphaPrime != 0)
    G = ImAlpha / (TMath::Sqrt(s) * ReAlphaPrime);
  else
    {
      G = 0.;
    }
  
  im = ImAlpha;
  return G;  
}



Double_t Im_Amplitude(Double_t s, Double_t t, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm_Roper, Int_t include_bg_jenk, Double_t par_jenk[])
{
  
  Double_t ImA = 0.;
  Double_t Bg = 0.;
  Double_t ImA_clean = 0.;
  if(include_bg_jenk == 1 || include_bg_jenk == 0) {
    for(Int_t i=0; i<n_sum; i++)
      {
	Double_t im;
	Im_alpha(s, par, lambda, im);
	Double_t A = 0.;
	A = TMath::Power(1. - t/0.71, -4.*(i+2)) * im;
	A /= TMath::Power(2*i+2.5-Re_alpha(s, par, lambda), 2.) + TMath::Power(im, 2.);
	ImA_clean += A;
      }
  }
  if(include_bg_jenk == -1 || include_bg_jenk == 1) {
    Double_t c4;
    c4 = TMath::Power(par_jenk[1]/(par_jenk[1] - t), 4.);
    Double_t s0 = TMath::Power(0.93827 + 0.13498, 2.);
    Double_t n_minus_a0a1 = 0.5 - par_jenk[2] - par_jenk[3]*TMath::Sqrt(s0);
    Double_t ds = s0 - s;

    if(ds >= 0)
      Bg  = 0.;//par_jenk[0] * c4 / (n_minus_a0a1 + par_jenk[3]*TMath::Sqrt(ds));
    else
      Bg = - par_jenk[0] * c4 * par_jenk[3] * TMath::Sqrt(-ds) / (n_minus_a0a1*n_minus_a0a1 - ds*par_jenk[3]*par_jenk[3]);
  }
  if(include_bg_jenk == -1) {
    ImA = Bg;
  }
  if(include_bg_jenk == 0) {
    ImA = ImA_clean;
    if(g_b_include_bg_jenk)
      ImA =  par_jenk[4] * ImA_clean;
  }
  if(include_bg_jenk == 1) {
    ImA = Bg + par_jenk[4] * ImA_clean;
  }
  //include Roper resonance
  if(include_Roper) {
    Double_t im;
    Im_alpha(s, par, lambda, im);
    Double_t de2 = TMath::Power(TMath::Sqrt(s)-1.440, 2.);
    Double_t G2 = TMath::Power(.325, 2.);
    Double_t Norm = norm_Roper;
    ImA += Norm * TMath::Power(1. - t/0.71, -4.) / (de2 + G2/4.);
  }
  return ImA;
}

Double_t elastic_form_factor(Double_t t)
{
  return TMath::Power(1.-t/0.71, -2);
}

Double_t Risto_form_factor(Double_t t)
{
  Double_t b = 10.;
  return TMath::Exp(b*t);
}


Double_t DD_cross_mod(Double_t s, Double_t t, Double_t M2, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper,
		  Int_t include_bg_jenk, Double_t par_jenk[])
{
  Double_t res = 0;
  res = sigma_t_PpN(0, global_param, global_lambda_param, g_include_Roper, 0.036, g_include_bg_jenk, g_par_jenk, M2) *
	TMath::Power(elastic_form_factor(t), 2) * TMath::Power(s/M2, 2*(1.08+0.25*t)-2)/ M2 ;

  return res;
}

Double_t DD_cross_bRisto(Double_t s, Double_t t, Double_t M2, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper,
		  Int_t include_bg_jenk, Double_t par_jenk[])
{
  Double_t res = 0;
  res = sigma_t_PpN(0, global_param, global_lambda_param, g_include_Roper, 0.036, g_include_bg_jenk, g_par_jenk, M2) *
    TMath::Power(Risto_form_factor(t), 2) * TMath::Power(s/M2, 2*(1.08+0.25*t)-2)/ M2 ;

  return res;
}


Double_t DD_cross(Double_t s, Double_t t, Double_t M2, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper,
		  Int_t include_bg_jenk, Double_t par_jenk[])
{
  Double_t cross = 0.;
  cross = TMath::Power(1. - t/0.71, -4.);

  Double_t m_proton = 0.938272013;
  Double_t m_proton2 = m_proton*m_proton;
  cross *= TMath::Power(s/M2, 2*(1.08 + 0.25*t) - 2);

  cross *= sigma_t_PpN(t, par, lambda, include_Roper, norm_Roper, include_bg_jenk, par_jenk, M2);
  
  Double_t x_bjorken = 0.;
  if(t != M2 - m_proton2)
    x_bjorken = -t / (M2 - m_proton2 - t);

  cross *= kinematic_factor(x_bjorken, t, M2);  

  //  if(M2 - m_proton2 != 0.)
  //    cross /= M2 - m_proton2;
  
  //   Double_t x_bjorken = 0.;
  //   if(t != M2 - m_proton2)
  //     x_bjorken = -t / (M2 - m_proton2 - t);
  //   Double_t x_bjorken2 = x_bjorken*x_bjorken;
  //   if(1 - 4*m_proton2*x_bjorken2/t != 0.)
  //     cross *= x_bjorken * (1-x_bjorken) * (1-x_bjorken) / TMath::Power((1 - 4*m_proton2*x_bjorken2/t), 1.5);
  //   if(M2 - m_proton2 != 0.)
  //     cross /= M2 - m_proton2;
  //   cross *= Im_Amplitude(M2, t, par, lambda, include_Roper, norm_Roper, include_bg_jenk, par_jenk) / m_proton;

  if(!g_b_include_bg_jenk)
    cross *= norm;

  return cross;
}

Double_t Theta_function(Double_t s1, Double_t s2)
{
  Double_t res = 0.;
  if(s1 > s2)
    res = 1.;
  else if(s1 == s2)
    {
      //      cout << "s1 == s2 in Theta_function" << endl;
      res = 0.5;
    }
  
  return res;
}

Double_t sigma_t_PpN(Double_t t, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm_Roper, Int_t include_bg_jenk,
		     Double_t par_jenk[], Double_t M2)
{
  Double_t res = 0;
  Double_t m_proton = 0.938272013;
  Double_t m_proton2 = m_proton*m_proton;

  //   Double_t x_bjorken2 = x_bjorken*x_bjorken;
  //   if(1 - 4*m_proton2*x_bjorken2/t != 0.)
  //     res = x_bjorken * (1-x_bjorken) * (1-x_bjorken) / TMath::Power((1 - 4*m_proton2*x_bjorken2/t), 1.5);

  //   if(M2 - m_proton2 != 0.)
  //     res /= M2 - m_proton2;

  //  res *= kinematic_factor(x_bjorken, t);

  res = Im_Amplitude(M2, t, par, lambda, include_Roper, norm_Roper, include_bg_jenk, par_jenk) / m_proton;

  return res;
}

Double_t kinematic_factor(Double_t x, Double_t t, Double_t M2)
{
  Double_t res = 0;
  Double_t m_proton = 0.938272013;
  Double_t m_proton2 = m_proton*m_proton;

  Double_t x2 = x*x;
  if(1 - 4*m_proton2*x2/t != 0.)
    res = x * (1-x) * (1-x) / TMath::Power((1 - 4*m_proton2*x2/t), 1.5);

  if(M2 - m_proton2 != 0.)
    res /= M2 - m_proton2;

  return res;
}

Double_t dsigma3_dtdM1dM2(Double_t s, Double_t t, Double_t M2_1, Double_t M2_2, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper, Int_t include_bg_jenk, Double_t par_jenk[])
{
  Double_t cross = sigma_t_PpN(t, par, lambda, include_Roper, norm_Roper, include_bg_jenk, par_jenk, M2_1) *
    sigma_t_PpN(t, par, lambda, include_Roper, norm_Roper, include_bg_jenk, par_jenk, M2_2);
  cross *= TMath::Power(s/(TMath::Power(TMath::Sqrt(M2_1) + TMath::Sqrt(M2_2) ,2)), 2*(1.08 + 0.25*t) - 2);
  cross /= M2_1 * M2_2;
  return cross;
}

Double_t dsigma3_dt1dt2dM(Double_t s1, Double_t s2, Double_t t1, Double_t t2, Double_t M2, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper, Int_t include_bg_jenk, Double_t par_jenk[])
{
  Double_t cross = TMath::Power(elastic_form_factor(t1), 2);
  cross *= TMath::Power(elastic_form_factor(t2), 2);
  cross *= TMath::Power(s1/M2, 2*(1.08+0.25*t1)-2);
  cross *= TMath::Power(s2/M2, 2*(1.08+0.25*t2)-2);
  cross *= sigma_t_PpN(0, par, lambda, include_Roper, norm_Roper, include_bg_jenk, par_jenk, M2);
  return cross;
}


Double_t double_integral_DD_cross(Double_t s, Double_t* par, Double_t* lambda_par,
                                  Double_t t_min, Double_t t_max, Double_t M2_min, Double_t M2_max, Int_t n_int,
                                  Bool_t include_Roper, Int_t bg_parameter, Double_t c_cross, Double_t c_bg, Double_t b, Double_t norm_Roper,
                                  Int_t include_bg_jenk, Double_t par_jenk[])
{
  Double_t res = 0.;
  Double_t integral = 0.;
  Double_t add_bg = 0.;
  Double_t dt = (t_max-t_min)/(double)n_int;
  Double_t dM2 = (M2_max-M2_min)/(double)n_int;
  if(dt<=0 || dM2<=0) {
    cout << "Error in double_integral_DD_cross: dt = " << dt << ", dM2 = " << dM2 << endl;
    exit(-1);
  }
  Double_t bg = 0.;

  //calculate cross section
  for(Int_t i=0; i<n_int+1; i++)
    {
      Double_t t_i = t_min + (t_max-t_min)*(i+0.5)/(double)n_int;
      for(Int_t j=0; j<n_int+1; j++)
        {
          //      cout << "dd: " << i << " " << j << endl;
          Double_t M2_j = M2_min + (M2_max-M2_min)*(j+0.5)/(double)n_int;
          DD_cross_array[i][j] = DD_cross(s, t_i, M2_j, par, lambda_par, include_Roper, c_cross, norm_Roper, include_bg_jenk, par_jenk);
        }
    }



  for(Int_t i=0; i<n_int; i++)
    {
      for(Int_t j=0; j<n_int; j++)
	{
	  Double_t cross = 0., value = 0.;
	  Double_t bg_temp = 0.;
	  cross = 0.25*(DD_cross_array[i][j] + DD_cross_array[i+1][j] + DD_cross_array[i][j+1] + DD_cross_array[i+1][j+1]);
	  value = cross;
	  integral += value;
	}
    }


  integral = integral * dt * dM2;
  //CALCULATE background
  bg = 1./b;//(TMath::Log(s) - TMath::Log(M2_min)) * (1. - TMath::Exp(-b*s)) / b;
  bg *= c_bg;
  //  if(g_include_bg) {
  if(bg_parameter == -1)
    res = bg;
  if(bg_parameter == 0)
    res = integral;
  if(bg_parameter == 1)
    res = integral + bg;

  return res;
}


Double_t double_integral_DD_cross_mod(Double_t s, Double_t* par, Double_t* lambda_par,
                                  Double_t t_min, Double_t t_max, Double_t M2_min, Double_t M2_max, Int_t n_int,
                                  Bool_t include_Roper, Int_t bg_parameter, Double_t c_cross, Double_t c_bg, Double_t b, Double_t norm_Roper,
                                  Int_t include_bg_jenk, Double_t par_jenk[])
{
  Double_t res = 0.;
  Double_t integral = 0.;
  Double_t add_bg = 0.;
  Double_t dt = (t_max-t_min)/(double)n_int;
  Double_t dM2 = (M2_max-M2_min)/(double)n_int;
  if(dt<=0 || dM2<=0) {
    cout << "Error in double_integral_DD_cross: dt = " << dt << ", dM2 = " << dM2 << endl;
    exit(-1);
  }
  Double_t bg = 0.;

  Double_t sigmatppn = 1;
  for(Int_t j=0; j<n_int+1; j++)
    {
      //      cout << "dd: " << i << " " << j << endl;
      Double_t M2_j = M2_min + (M2_max-M2_min)*(j+0.5)/(double)n_int;
      sigmatppn = sigma_t_PpN(0, par, lambda_par, include_Roper, norm_Roper, include_bg_jenk, par_jenk, M2_j);
      for(Int_t i=0; i<n_int+1; i++)
	{
	  Double_t t_i = t_min + (t_max-t_min)*(i+0.5)/(double)n_int;
          DD_cross_array_mod[i][j] = sigmatppn*
	    TMath::Power(elastic_form_factor(t_i), 2) * TMath::Power(s/M2_j, 2*(1.08+0.25*t_i)-2)/ M2_j ;//DD_cross_mod(s, t_i, M2_j, par, lambda_par, include_Roper, c_cross, norm_Roper, include_bg_jenk, par_jenk);
        }
    }



  for(Int_t i=0; i<n_int; i++)
    {
      for(Int_t j=0; j<n_int; j++)
	{
	  Double_t cross = 0., value = 0.;
	  Double_t bg_temp = 0.;
	  cross = 0.25*(DD_cross_array_mod[i][j] + DD_cross_array_mod[i+1][j] + DD_cross_array_mod[i][j+1] + DD_cross_array_mod[i+1][j+1]);
	  value = cross;
	  integral += value;
	}
    }


  integral = integral * dt * dM2;
  //CALCULATE background
  bg = 1./b;//(TMath::Log(s) - TMath::Log(M2_min)) * (1. - TMath::Exp(-b*s)) / b;
  bg *= c_bg;
  //  if(g_include_bg) {
  if(bg_parameter == -1)
    res = bg;
  if(bg_parameter == 0)
    res = integral;
  if(bg_parameter == 1)
    res = integral + bg;

  return res;
}

Double_t triple_integral_dsigma3_dt1dt2dM(Double_t s, Double_t t_min, Double_t t_max, Double_t M2_min, Double_t M2_max, Int_t n_int, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper, Int_t include_bg_jenk, Double_t par_jenk[])
{
  Double_t res = 0.;
  Double_t integral = 0.;
  Double_t dt = (t_max-t_min)/(double)n_int;
  Double_t dM2 = (M2_max-M2_min)/(double)n_int;
  if(dt<=0 || dM2<=0) {
    cout << "Error in double_integral_DD_cross: dt = " << dt << ", dM2 = " << dM2 << endl;
    exit(-1);
  }

  for(Int_t i=0; i<n_int+1; i++)
    {
      Double_t t_i = t_min + (t_max-t_min)*(i+0.5)/(double)n_int;
      for(Int_t j=0; j<n_int+1; j++)
	{ 
	  Double_t t_j = t_min + (t_max-t_min)*(j+0.5)/(double)n_int;
	  for(Int_t k=0; k<n_int+1; k++)
	    {
	      //      cout << "dd: " << i << " " << j << endl;
	      Double_t M2_k = M2_min + (M2_max-M2_min)*(k+0.5)/(double)n_int;
	      ds3_dt1dt2dM_cross_array[i][j][k] = dsigma3_dt1dt2dM(s/2., s/2., t_i, t_j, M2_k, par, lambda, include_Roper, g_c_cross, norm_Roper, include_bg_jenk, par_jenk);
	    }
	}
    }

  for(Int_t i=0; i<n_int; i++)
    {
      for(Int_t j=0; j<n_int; j++)
	{
	  for(Int_t k=0; k<n_int; k++)
	    {
	      Double_t cross = 0., value = 0.;
	      Double_t bg_temp = 0.;
	      cross = 0.125*(ds3_dt1dt2dM_cross_array[i][j][k] + 
			     ds3_dt1dt2dM_cross_array[i+1][j][k] + 
			     ds3_dt1dt2dM_cross_array[i][j+1][k] + 
			     ds3_dt1dt2dM_cross_array[i][j][k+1] +
			     ds3_dt1dt2dM_cross_array[i+1][j+1][k] +
			     ds3_dt1dt2dM_cross_array[i][j+1][k+1] + 
			     ds3_dt1dt2dM_cross_array[i+1][j][k+1] + 
			     ds3_dt1dt2dM_cross_array[i+1][j+1][k+1]);
	      value = cross;
	      integral += value;
	    }
	}
    }
  
  integral *= dt*dt*dM2;
  res = integral;
  return res;
}




Double_t triple_integral_dsigma3_dt1dt2dM_copy(Double_t s, Double_t t_min, Double_t t_max, Double_t M2_min, Double_t M2_max, Int_t n_int, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper, Int_t include_bg_jenk, Double_t par_jenk[])
{
  Double_t res = 0.;
  Double_t integral = 0.;
  Double_t dt = (t_max-t_min)/(double)n_int;
  Double_t dM2 = (M2_max-M2_min)/(double)n_int;
  if(dt<=0 || dM2<=0) {
    cout << "Error in double_integral_DD_cross: dt = " << dt << ", dM2 = " << dM2 << endl;
    exit(-1);
  }

  Double_t sigmatppn = 1;
  for(Int_t k=0; k<n_int+1; k++)
    {
      Double_t M2_k = M2_min + (M2_max-M2_min)*(k+0.5)/(double)n_int;      
      sigmatppn = sigma_t_PpN(0, par, lambda, include_Roper, norm_Roper, include_bg_jenk, par_jenk, M2_k);
      for(Int_t i=0; i<n_int+1; i++)
	{
	  Double_t t_i = t_min + (t_max-t_min)*(i+0.5)/(double)n_int;
	  for(Int_t j=0; j<n_int+1; j++)
	    { 
	      Double_t t_j = t_min + (t_max-t_min)*(j+0.5)/(double)n_int;
	      ds3_dt1dt2dM_cross_array[i][j][k] = sigmatppn;
	      ds3_dt1dt2dM_cross_array[i][j][k] *= TMath::Power(elastic_form_factor(t_i), 2);
	      ds3_dt1dt2dM_cross_array[i][j][k] *= TMath::Power(elastic_form_factor(t_j), 2);
	      ds3_dt1dt2dM_cross_array[i][j][k] *= TMath::Power(0.5*s/M2_k, 2*(1.08+0.25*t_i)-2);
	      ds3_dt1dt2dM_cross_array[i][j][k] *= TMath::Power(0.5*s/M2_k, 2*(1.08+0.25*t_j)-2);
	    }
	}
    }

  for(Int_t i=0; i<n_int; i++)
    {
      for(Int_t j=0; j<n_int; j++)
	{
	  for(Int_t k=0; k<n_int; k++)
	    {
	      Double_t cross = 0., value = 0.;
	      Double_t bg_temp = 0.;
	      cross = 0.125*(ds3_dt1dt2dM_cross_array[i][j][k] + 
			     ds3_dt1dt2dM_cross_array[i+1][j][k] + 
			     ds3_dt1dt2dM_cross_array[i][j+1][k] + 
			     ds3_dt1dt2dM_cross_array[i][j][k+1] +
			     ds3_dt1dt2dM_cross_array[i+1][j+1][k] +
			     ds3_dt1dt2dM_cross_array[i][j+1][k+1] + 
			     ds3_dt1dt2dM_cross_array[i+1][j][k+1] + 
			     ds3_dt1dt2dM_cross_array[i+1][j+1][k+1]);
	      value = cross;
	      integral += value;
	    }
	}
    }
  
  integral *= dt*dt*dM2;
  res = integral;
  return res;
}



Double_t triple_integral_dsigma3_dtdM1dM2_copy(Double_t s, Double_t t_min, Double_t t_max, Double_t M2_min, Double_t M2_max, Int_t n_int, Double_t* par, Double_t* lambda, Bool_t include_Roper, Double_t norm, Double_t norm_Roper, Int_t include_bg_jenk, Double_t par_jenk[])
{
  
  Double_t res = 0.;
  Double_t integral = 0.;
  Double_t dt = (t_max-t_min)/(double)n_int;
  Double_t dM2 = (M2_max-M2_min)/(double)n_int;
  if(dt<=0 || dM2<=0) {
    cout << "Error in double_integral_DD_cross: dt = " << dt << ", dM2 = " << dM2 << endl;
    exit(-1);
  }

  vector<Double_t> sigmatppn;
  for(Int_t k=0; k<n_int+1; k++)
    {
      Double_t M2_k = M2_min + (M2_max-M2_min)*double(k)/(double)n_int;      
      sigmatppn.push_back(sigma_t_PpN(0, par, lambda, include_Roper, norm_Roper, include_bg_jenk, par_jenk, M2_k));
    }
  for(Int_t k=0; k<n_int+1; k++)
    {
      Double_t M2_k = M2_min + (M2_max-M2_min)*double(k)/(double)n_int;
      for(Int_t i=0; i<n_int+1; i++)
	{
	  Double_t M2_i = M2_min + (M2_max-M2_min)*double(i)/(double)n_int;
	  for(Int_t j=0; j<n_int+1; j++)
	    { 
	      Double_t t_j = t_min + (t_max-t_min)*double(j)/(double)n_int;
	      ds3_dtdM1dM2_cross_array[i][j][k] = sigmatppn[k] * sigmatppn[i] / (M2_k * M2_i);
	      ds3_dtdM1dM2_cross_array[i][j][k] *= TMath::Power(s/(TMath::Power(TMath::Sqrt(M2_k)+TMath::Sqrt(M2_i), 2)), 2*(1.08+0.25*t_j)-2);
	    }
	}
    }
    

  for(Int_t i=0; i<n_int; i++)
    {
      for(Int_t j=0; j<n_int; j++)
	{
	  for(Int_t k=0; k<n_int; k++)
	    {
	      Double_t cross = 0., value = 0.;
	      Double_t bg_temp = 0.;
	      cross = 0.125*(ds3_dtdM1dM2_cross_array[i][j][k] + 
			     ds3_dtdM1dM2_cross_array[i+1][j][k] + 
			     ds3_dtdM1dM2_cross_array[i][j+1][k] + 
			     ds3_dtdM1dM2_cross_array[i][j][k+1] +
			     ds3_dtdM1dM2_cross_array[i+1][j+1][k] +
			     ds3_dtdM1dM2_cross_array[i][j+1][k+1] + 
			     ds3_dtdM1dM2_cross_array[i+1][j][k+1] + 
			     ds3_dtdM1dM2_cross_array[i+1][j+1][k+1]);
	      value = cross;
	      integral += value;
	    }
	}
    }
  
  integral *= dt*dM2*dM2;
  res = integral;
  return res;

}


Double_t integral_M2_DD_cross_mod(Double_t s, Double_t t, Double_t* par, Double_t* lambda_par,
                              Double_t M2_min, Double_t M2_max, Int_t n_int, Bool_t include_Roper, Int_t bg_parameter, Double_t c_cross, Double_t c_bg,
                              Double_t b, Double_t norm_Roper, Int_t include_bg_jenk, Double_t par_jenk[])
{
  Double_t integral = 0., res = 0., bg = 0.;
  Double_t dM2 = (M2_max-M2_min)/(double)n_int;
  if(dM2<=0) {
    cout << "Error in integral_M2_DD_cross: dM2 = " << dM2 << endl;
    exit(-1);
  }
  if(-t > s) {
    cout << "Error in integral_M2_DD_cross: -t = " << -t << " < s = " << s << endl;
  }
  for(Int_t i=0; i<n_int; i++)
    {
      Double_t M2_i = M2_min + (M2_max-M2_min)*(i+0.5)/(double)n_int;
      //      Double_t cross = DD_cross(s, t, M2_i, par, lambda_par, include_Roper, c_cross);
      //      Double_t bg_temp = DD_bg(s, t, M2_i, par, lambda_par, c_bg, b);
      //      Double_t value = cross + bg_temp;
      Double_t cross = 0., bg_temp = 0., value = 0.;
      /*      if(bg_parameter == -1) { //bg only
        bg_temp = DD_bg(s, t, M2_i, par, lambda_par, c_bg, b);
        value = bg_temp;
      }
      else if(bg_parameter == 0) {
        cross = DD_cross(s, t, M2_i, par, lambda_par, include_Roper, c_cross, norm_Roper);
        value = cross;
      }
      else if(bg_parameter == 1) {
        bg_temp = DD_bg(s, t, M2_i, par, lambda_par, c_bg, b);
        cross = DD_cross(s, t, M2_i, par, lambda_par, include_Roper, c_cross, norm_Roper);
        value = cross + bg_temp;
        }*/
      cross = DD_cross_mod(s, t, M2_i, par, lambda_par, include_Roper, c_cross, norm_Roper, include_bg_jenk, par_jenk);
      value = cross;
      integral += value;
    }

  integral = integral * dM2;
  //CALCULATE background
  //  cross_bg = 1./M2;
  //  cross_bg *=  TMath::Power(1. - t/0.71, -8.);
  bg = TMath::Exp(b*t);
  bg /= (1. - TMath::Exp(-b*s));
  bg *= c_bg;

  if(bg_parameter == -1)
    res = bg;
  if(bg_parameter == 0)
    res = integral;
  if(bg_parameter == 1)
    res = integral + bg;
  //integral += bg;
  //res = integral;
  return res;
}



Double_t integral_M2_DD_cross(Double_t s, Double_t t, Double_t* par, Double_t* lambda_par,
                              Double_t M2_min, Double_t M2_max, Int_t n_int, Bool_t include_Roper, Int_t bg_parameter, Double_t c_cross, Double_t c_bg,
                              Double_t b, Double_t norm_Roper, Int_t include_bg_jenk, Double_t par_jenk[])
{
  Double_t integral = 0., res = 0., bg = 0.;
  Double_t dM2 = (M2_max-M2_min)/(double)n_int;
  if(dM2<=0) {
    cout << "Error in integral_M2_DD_cross: dM2 = " << dM2 << endl;
    exit(-1);
  }
  if(-t > s) {
    cout << "Error in integral_M2_DD_cross: -t = " << -t << " < s = " << s << endl;
  }
  for(Int_t i=0; i<n_int; i++)
    {
      Double_t M2_i = M2_min + (M2_max-M2_min)*(i+0.5)/(double)n_int;
      //      Double_t cross = DD_cross(s, t, M2_i, par, lambda_par, include_Roper, c_cross);
      //      Double_t bg_temp = DD_bg(s, t, M2_i, par, lambda_par, c_bg, b);
      //      Double_t value = cross + bg_temp;
      Double_t cross = 0., bg_temp = 0., value = 0.;
      /*      if(bg_parameter == -1) { //bg only
        bg_temp = DD_bg(s, t, M2_i, par, lambda_par, c_bg, b);
        value = bg_temp;
      }
      else if(bg_parameter == 0) {
        cross = DD_cross(s, t, M2_i, par, lambda_par, include_Roper, c_cross, norm_Roper);
        value = cross;
      }
      else if(bg_parameter == 1) {
        bg_temp = DD_bg(s, t, M2_i, par, lambda_par, c_bg, b);
        cross = DD_cross(s, t, M2_i, par, lambda_par, include_Roper, c_cross, norm_Roper);
        value = cross + bg_temp;
        }*/
      cross = DD_cross(s, t, M2_i, par, lambda_par, include_Roper, c_cross, norm_Roper, include_bg_jenk, par_jenk);
      value = cross;
      integral += value;
    }

  integral = integral * dM2;
  //CALCULATE background
  //  cross_bg = 1./M2;
  //  cross_bg *=  TMath::Power(1. - t/0.71, -8.);
  bg = TMath::Exp(b*t);
  bg /= (1. - TMath::Exp(-b*s));
  bg *= c_bg;

  if(bg_parameter == -1)
    res = bg;
  if(bg_parameter == 0)
    res = integral;
  if(bg_parameter == 1)
    res = integral + bg;
  //integral += bg;
  //res = integral;
  return res;
}

Double_t integral_M2_DD_cross_bRisto(Double_t s, Double_t t, Double_t* par, Double_t* lambda_par,
                              Double_t M2_min, Double_t M2_max, Int_t n_int, Bool_t include_Roper, Int_t bg_parameter, Double_t c_cross, Double_t c_bg,
                              Double_t b, Double_t norm_Roper, Int_t include_bg_jenk, Double_t par_jenk[])
{
  Double_t integral = 0., res = 0., bg = 0.;
  Double_t dM2 = (M2_max-M2_min)/(double)n_int;
  if(dM2<=0) {
    cout << "Error in integral_M2_DD_cross: dM2 = " << dM2 << endl;
    exit(-1);
  }
  if(-t > s) {
    cout << "Error in integral_M2_DD_cross: -t = " << -t << " < s = " << s << endl;
  }
  for(Int_t i=0; i<n_int; i++)
    {
      Double_t M2_i = M2_min + (M2_max-M2_min)*(i+0.5)/(double)n_int;
      //      Double_t cross = DD_cross(s, t, M2_i, par, lambda_par, include_Roper, c_cross);
      //      Double_t bg_temp = DD_bg(s, t, M2_i, par, lambda_par, c_bg, b);
      //      Double_t value = cross + bg_temp;
      Double_t cross = 0., bg_temp = 0., value = 0.;
      /*      if(bg_parameter == -1) { //bg only
        bg_temp = DD_bg(s, t, M2_i, par, lambda_par, c_bg, b);
        value = bg_temp;
      }
      else if(bg_parameter == 0) {
        cross = DD_cross(s, t, M2_i, par, lambda_par, include_Roper, c_cross, norm_Roper);
        value = cross;
      }
      else if(bg_parameter == 1) {
        bg_temp = DD_bg(s, t, M2_i, par, lambda_par, c_bg, b);
        cross = DD_cross(s, t, M2_i, par, lambda_par, include_Roper, c_cross, norm_Roper);
        value = cross + bg_temp;
        }*/
      cross = DD_cross_bRisto(s, t, M2_i, par, lambda_par, include_Roper, c_cross, norm_Roper, include_bg_jenk, par_jenk);
      value = cross;
      integral += value;
    }

  integral = integral * dM2;
  //CALCULATE background
  //  cross_bg = 1./M2;
  //  cross_bg *=  TMath::Power(1. - t/0.71, -8.);
  bg = TMath::Exp(b*t);
  bg /= (1. - TMath::Exp(-b*s));
  bg *= c_bg;

  if(bg_parameter == -1)
    res = bg;
  if(bg_parameter == 0)
    res = integral;
  if(bg_parameter == 1)
    res = integral + bg;
  //integral += bg;
  //res = integral;
  return res;
}


