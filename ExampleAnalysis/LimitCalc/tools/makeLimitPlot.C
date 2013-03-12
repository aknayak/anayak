#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>
#include <cstring>

#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TGraphPainter.h"
#include "TMultiGraph.h"
#include "TTree.h"

using namespace std;


void plot( TString channel = "",
	   Float_t xMax_   = 1.0){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  
  TLegend* leg = new TLegend(0.14,0.60,0.41,0.85,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);


  float X[]        = {80,100,120,140,150,155,160};

  float obsY[]     = {0.,0.,0.,0.,0.,0.,0.};
  float expY[]     = {0.,0.,0.,0.,0.,0.,0.};

  float expX1sL[]  = {0.,0.,0.,0.,0.,0.,0.};
  float expX1sH[]  = {0.,0.,0.,0.,0.,0.,0.};
  float expY1sL[]  = {0.,0.,0.,0.,0.,0.,0.};
  float expY1sH[]  = {0.,0.,0.,0.,0.,0.,0.};

  float expX2sL[]  = {0.,0.,0.,0.,0.,0.,0.};
  float expX2sH[]  = {0.,0.,0.,0.,0.,0.,0.};
  float expY2sL[]  = {0.,0.,0.,0.,0.,0.,0.};
  float expY2sH[]  = {0.,0.,0.,0.,0.,0.,0.};

  int nMassPoints = 7;

  for(unsigned int i = 0 ; i < nMassPoints; i++){

    TFile fobs(Form("datacard_csbar_mH%.0f.txt_Hybrid_freqObsLimit.root",X[i]),"READ");
    TFile fexp(Form("datacard_csbar_mH%.0f.txt_HybridHybrid_limitbands.root",X[i]),"READ");

    if(fobs.IsZombie() || fexp.IsZombie()){
      cout << "Cannot open file for mass " << X[i] << endl;
      continue;
    }

    Double_t r_obs;
    TTree* Tobs = (TTree*)fobs.Get("T");
    Tobs->SetBranchAddress("limit",&r_obs);

    for(int k = 0 ; k< Tobs->GetEntries() ; k++){
      Tobs->GetEntry(k);
      obsY[i]    = r_obs;
    }

    Double_t r_median, r_mean, r_m1s, r_m2s, r_p1s, r_p2s;
    TTree* Texp = (TTree*)fexp.Get("T"); 
    Texp->SetBranchAddress("rmean",&r_mean); 
    Texp->SetBranchAddress("rmedian",&r_median); 
    Texp->SetBranchAddress("rm1s",&r_m1s); 
    Texp->SetBranchAddress("rm2s",&r_m2s); 
    Texp->SetBranchAddress("rp1s",&r_p1s); 
    Texp->SetBranchAddress("rp2s",&r_p2s); 
    
    for(int k = 0 ; k< Texp->GetEntries() ; k++){ 
      Texp->GetEntry(k); 
      expY2sL[i] = r_m2s;
      expY1sL[i] = r_m1s;
      expY[i]    = r_median;
      expY1sH[i] = r_p1s;
      expY2sH[i] = r_p2s;

    } 

  }

  for(int i1 = 0 ; i1 < nMassPoints ; i1++){
    expY1sH[i1] = TMath::Abs(expY1sH[i1]-expY[i1]);
    expY1sL[i1] = TMath::Abs(expY1sL[i1]-expY[i1]);
    expY2sH[i1] = TMath::Abs(expY2sH[i1]-expY[i1]);
    expY2sL[i1] = TMath::Abs(expY2sL[i1]-expY[i1]);

    cout << "Mass " << X[i1] << " => " << expY2sL[i1] << ", " << expY1sL[i1] << ", " << expY1sH[i1] << ", "
	 << expY2sH[i1] << ", " << endl;
  }

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("");

  TGraphAsymmErrors* expected = new TGraphAsymmErrors(nMassPoints, X, expY, expX1sL ,expX1sL , expX1sL, expX1sL);
  TGraphAsymmErrors* oneSigma = new TGraphAsymmErrors(nMassPoints, X, expY, expX1sL, expX1sL,  expY1sL, expY1sH);
  TGraphAsymmErrors* twoSigma = new TGraphAsymmErrors(nMassPoints, X, expY, expX2sL, expX2sL,  expY2sL, expY2sH);
  TGraphAsymmErrors* observed = new TGraphAsymmErrors(nMassPoints, X, obsY, expX1sL ,expX1sL , expX1sL, expX1sL);

 
  oneSigma->SetMarkerColor(kBlack);
  oneSigma->SetMarkerStyle(kFullCircle);
  oneSigma->SetFillColor(kGreen);
  oneSigma->SetFillStyle(1001);

  twoSigma->SetMarkerColor(kBlack);
  twoSigma->SetMarkerStyle(kFullCircle);
  twoSigma->SetFillColor(kYellow);
  twoSigma->SetFillStyle(1001);

  expected->SetMarkerColor(kBlack);
  expected->SetMarkerStyle(kFullCircle);
  expected->SetMarkerSize(1.5);
  expected->SetLineColor(kBlack);
  expected->SetLineWidth(2);

  observed->SetMarkerColor(kBlue);
  observed->SetMarkerStyle(1);
  observed->SetLineColor(kBlue);
  observed->SetLineWidth(4);

  mg->Add(twoSigma);
  mg->Add(oneSigma);
  mg->Add(expected);
  mg->Add(observed);

  mg->Draw("ALP3");

  c1->cd();
  gPad->Modified();
  mg->GetXaxis()->SetLimits(80,160);
  mg->GetYaxis()->SetTitleOffset(1.02);
  mg->SetMinimum(0.);
  mg->SetMaximum(xMax_);
  mg->GetXaxis()->SetTitle("m_{H^{+}} (GeV)");
  mg->GetYaxis()->SetTitle("95% CL limit for BR(t#rightarrow bH^{#pm})");

  leg->SetHeader(Form("#splitline{CMS Preliminary #sqrt{s}=8 TeV}{17.7 fb^{-1}}"));

  leg->AddEntry(expected,"Expected","P");
  leg->AddEntry(observed,"Observed","L");

  leg->Draw();

  TPaveText *pl2 = new TPaveText(0.68,0.75,0.88,0.88, "brNDC");
  pl2->SetTextSize(0.032);
  pl2->SetFillColor(0);
  pl2->SetTextFont(132);
  pl2->SetBorderSize(0);
  pl2->SetTextAlign(11);
  pl2->AddText("t #rightarrow H^{#pm}b, H^{+} #rightarrow c#bar{s}");
  //pl2->AddText(label.c_str());
  pl2->AddText("BR(H^{+} #rightarrow c#bar{s}) = 1");
  pl2->Draw();

  gPad->RedrawAxis();
  gPad->SaveAs(Form("limit_csbar_%s.png",channel.Data()));
  gPad->SaveAs(Form("limit_csbar_%s.pdf",channel.Data()));

}


void plotAll(){
  
  plot( "",1);


}
