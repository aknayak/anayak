#include "MiniTree/Selection/interface/HistogramPlotter.hh" 
#include <iostream> 
#include <iomanip> 
#include <math.h>
 
ClassImp(HistogramPlotter) 

  void HistogramPlotter::InitHist(TString dirname, TString parentDir, TFile *file)
{
  std::string name(dirname);
  std::string fullname;
  if(parentDir.Length() != 0){
    TDirectory *d = file->GetDirectory(parentDir.Data());
    d->mkdir(name.c_str());
    fullname = std::string(parentDir+"/"+dirname);
  }
  else{file->mkdir(name.c_str()); fullname = name;}
  TDirectory *d = file->GetDirectory(fullname.c_str());
  file->cd(fullname.c_str());
  addHisto("pt_jet", fullname, 500, 0., 100.);
  TH1 *h1 = getHisto("pt_jet", fullname);
  h1->SetDirectory(d);
  addHisto("eta_jet", fullname, 50, -5.0, 5.0);
  h1 = getHisto("eta_jet", fullname);
  h1->SetDirectory(d);
  addHisto("phi_jet", fullname, 63, -M_PI, M_PI);
  h1 = getHisto("phi_jet", fullname);
  h1->SetDirectory(d);
  addHisto("multi_jet", fullname, 20, 0., 20.);
  h1 = getHisto("multi_jet", fullname);
  h1->SetDirectory(d);
  addHisto("btag_jet", fullname, 100, -10., 10.);
  h1 = getHisto("btag_jet", fullname);
  h1->SetDirectory(d);
  addHisto("btagmulti_jet", fullname, 10, 0., 10.);
  h1 = getHisto("btagmulti_jet", fullname);
  h1->SetDirectory(d);

  addHisto("pt_ele", fullname, 500, 0., 100.);
  h1 = getHisto("pt_ele", fullname);
  h1->SetDirectory(d);
  addHisto("eta_ele", fullname, 50, -5.0, 5.0);
  h1 = getHisto("eta_ele", fullname);
  h1->SetDirectory(d);
  addHisto("phi_ele", fullname, 63, -M_PI, M_PI);
  h1 = getHisto("phi_ele", fullname);
  h1->SetDirectory(d);

  addHisto("pt_mu", fullname, 500, 0., 100.);
  h1 = getHisto("pt_mu", fullname);
  h1->SetDirectory(d);
  addHisto("eta_mu", fullname, 50, -5.0, 5.0);
  h1 = getHisto("eta_mu", fullname);
  h1->SetDirectory(d);
  addHisto("phi_mu", fullname, 63, -M_PI, M_PI);
  h1 = getHisto("phi_mu", fullname);
  h1->SetDirectory(d);

  addHisto("pt_tau", fullname, 500, 0., 100.);
  h1 = getHisto("pt_tau", fullname);
  h1->SetDirectory(d);
  addHisto("eta_tau", fullname, 50, -5.0, 5.0);
  h1 = getHisto("eta_tau", fullname);
  h1->SetDirectory(d);
  addHisto("phi_tau", fullname, 63, -M_PI, M_PI);
  h1 = getHisto("phi_tau", fullname);
  h1->SetDirectory(d);

  addHisto("pt_met", fullname, 500, 0., 100.);
  h1 = getHisto("pt_met", fullname);
  h1->SetDirectory(d);
  addHisto("phi_met", fullname, 63, -M_PI, M_PI);
  h1 = getHisto("phi_met", fullname);
  h1->SetDirectory(d);

}  

void HistogramPlotter::addHisto(TString name, TString dirname, int range, double min, double max)
{
  //TString fullname = name+"_"+dirname; 
  TString fullname = dirname+"/"+name;
  std::string hname(fullname);
  histos1_[fullname] = new TH1D(name.Data(), hname.c_str(), range, min, max); 
}

void HistogramPlotter::add2DHisto(TString name, TString dirname, int range1, double min1, double max1, int range2, double min2, double max2)
{
  //TString fullname = name+"_"+dirname;
  TString fullname = dirname+"/"+name;
  std::string hname(fullname); 
  histos2_[fullname] = new TH2D(name.Data(), hname.c_str(), range1, min1, max1, range2, min2, max2);
}


void HistogramPlotter::fillHisto(TString name, TString dirname, double value, double weight)
{
  //TString fullname = name+"_"+dirname;
  //TString fullname = dirname+"/"+name;
  TH1* h = getHisto(name, dirname);
  if(h != 0) h->Fill(value, weight);
}

TH1* HistogramPlotter::getHisto(TString name, TString dirname)
{
  //TString fullname = name+"_"+dirname;
  TString fullname = dirname+"/"+name;
  TH1 * h = 0;
  if(histos1_.find(fullname) != histos1_.end())h = histos1_[fullname];
  else if(histos2_.find(fullname) != histos2_.end())h = histos2_[fullname];
  return h;
}

