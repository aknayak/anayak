
/*
Run me like
  .L TauFakeRateEstimator.C+;
  TauFakeRateEstimator t;
  t.ComputeFakeRate("file_list");

*/
#include <iomanip>
#include <iostream>
#include <fstream>

#include "TRandom2.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TH2.h"
#include "TTimeStamp.h"
#include "Math/VectorUtil.h"

#include "MiniTree/Selection/interface/Reader.h"
#include "MiniTree/Selection/interface/ObjectSelector.hh"

using namespace std;

class TauFakeRateEstimator : public ObjectSelector
{
public :
  TauFakeRateEstimator() : ObjectSelector()
  {
    DRMIN_JET = 0.5;
    DRMIN_ELE = 0.5;
    DRMIN_TAU = 0.5;
    METCUT_   = 50.0;

    TAUPRONGS_ = 0;
 };
  ~TauFakeRateEstimator() { };

  
  void ComputeFakeRateDiJetData(const char* file_list,  string myKey="PFlow", string tauIsoLabel = "Loose", TString histoFiles_="TFR_DiJet_data");
  void ComputeFakeRateDiJetMC(string myKey="PFlow", string tauIsoLabel = "Loose", TString histoFiles_="TFR_DiJet_mc");
  void ComputeFakeRateWJet(const char* file_list,  string myKey="PFlow", string tauIsoLabel = "Loose", bool isData = true, TString histoFiles_="TFR_WJet_data");

  void addFakeRate(TH1* hn, TH1* hd, TH1* hf);
  void addFakeRate2D(TH2* hn, TH2* hd, TH2* hf);
  double getErrorFraction( double a, double b );
  void DefineHistos();

  void processEvents();

private :
  int    TAUPRONGS_;
  double DRMIN_JET, DRMIN_ELE, DRMIN_TAU, METCUT_;
  Reader *evR;
  map<string, TH1D*> histos_;
  map<string, TH2D*> histos2_;
};

/////Fake Rate QCD Data ////
void TauFakeRateEstimator::ComputeFakeRateDiJetData(const char* file_list,  string myKey, string tauIsoLabel, TString histoFiles_){

  bool isPFlow = (myKey.find("PFlow") != string::npos) ? true : false;
  string eAlgo("PFlow"), mAlgo("PFlow"), jAlgo("PFlow"), tAlgo("PFlow");
  if(myKey.find("PFlow") != string::npos){jAlgo = "PFlow"; tAlgo = "PFlow";}
  else if(myKey.find("HPS") != string::npos){jAlgo = "Jets"; tAlgo = "Taus"; mAlgo = "Muons", eAlgo = "Electrons";}
  
  evR = new Reader();
  
  int nEntries = evR->AssignEventTreeFromList(file_list);
  cout << nEntries << " events are available" << endl;
  if( nEntries == 0) {return; }

  //create the file ////////////////////////////////
  TString Filename_ = histoFiles_ +"_"+myKey+"_"+tauIsoLabel+".root";
  TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
  outFile_->SetCompressionLevel( 9 );
  //////////////////////////////////////////////////

  //define histograms
  DefineHistos();

  MyEvent *ev;
  int nTriggEvent = 0;
  for(int i=0; i<nEntries; ++i){
    Long64_t ientry = evR->LoadTree(i);
    if (ientry < 0) break;

    ev = evR->GetNewEventFromList(i);
    if(ev==0) continue;

    //trigger
    bool passTrig = false;
    vector<string> trig = ev->hlt;
    for(size_t it = 0; it < trig.size(); it++){
      if(trig[it].find("Jet") != string::npos) passTrig = true;
    }
    if(!passTrig) continue;
    nTriggEvent++;
    
    //get objects //
    vector<MyVertex> Vertices = ev->PrimaryVtxs;
    if(Vertices.size() <= 0){
      cout<<" no vertexes , exit"<<endl;
      continue;
    }
    vector<MyMuon> pfMuons = evR->getMuons(ev, mAlgo);
    vector<MyElectron> pfElectrons = evR->getElectrons(ev, eAlgo);
    vector<MyTau> pfTaus = evR->getTaus(ev, tAlgo);
    vector<MyJet> pfJets = evR->getJets(ev, jAlgo);

    // preselect objects //
    vector<int> m_init; m_init.clear();
    preSelectMuons(&m_init, pfMuons, Vertices[0], isPFlow);
    vector<int> e_init; e_init.clear();
    preSelectElectrons(&e_init, pfElectrons, Vertices[0], isPFlow);
    vector<int> t_init; t_init.clear();
    preSelectTaus( &t_init, pfTaus, TAUPRONGS_, tauIsoLabel, Vertices[0]);
    vector<int> j_init; j_init.clear();
    preSelectJets(jAlgo, &j_init, pfJets);

    // clean objects //
    vector<int> e_final; e_final.clear();
    ElectronCleaning( pfElectrons, pfMuons, &e_init, &e_final, &m_init, DRMIN_ELE);
    vector<int> t_final; t_final.clear();
    TauCleaning(pfTaus, pfMuons, pfElectrons, &t_init, &t_final, &m_init, &e_final, DRMIN_TAU);
    vector<int> j_afterLeptonRemoval; j_afterLeptonRemoval.clear();
    JetCleaning(pfJets, pfMuons, pfElectrons, pfTaus, &j_init, &j_afterLeptonRemoval, &m_init, &e_final, &t_final, DRMIN_JET, false);
    
    // remove jets matched to trigger & apply dZ cut and antilepton discr
    vector<int> j_final; j_final.clear();
    int nTrigJet = 0;
    for(size_t ijet = 0; ijet < j_afterLeptonRemoval.size(); ijet++){
      int ind =  j_afterLeptonRemoval[ijet];
      if(pfJets[ind].triggerJet_pt > 40)nTrigJet++;
    }
    //if(nTrigJet)cout<<" nTriggered jet "<<nTrigJet<<endl;
    double zvertex = Vertices[0].XYZ.z();
    for(size_t ijet = 0; ijet < j_afterLeptonRemoval.size(); ijet++){
      int ind =  j_afterLeptonRemoval[ijet];
      if(nTrigJet == 1 && pfJets[ind].triggerJet_pt > 40) continue; //not to use trigger matched Jet
      if(fabs(zvertex - pfJets[ind].tau_vertex.z()) > 0.2 )continue; //dZ < 0.2
      //if(pfJets[ind].tau_againstElectronMedium < 0.5 || pfJets[ind].tau_againstMuonMedium < 0.5)continue;  //discriminator against lepton
      j_final.push_back(ind);
    }

    //loop over jets 
    double nvtxs = Vertices.size();
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      histos_["pt_alljet"]->Fill(pfJets[ind_jet].p4.pt());
      histos_["eta_alljet"]->Fill(pfJets[ind_jet].p4.eta());
      histos_["phi_alljet"]->Fill(pfJets[ind_jet].p4.phi());
      double Rjet = ((pfJets[ind_jet].etaetaMoment + pfJets[ind_jet].phiphiMoment) > 0) ? sqrt(pfJets[ind_jet].etaetaMoment + pfJets[ind_jet].phiphiMoment) : 0.;
      histos_["r_alljet"]->Fill(Rjet);
      histos_["nvtx_alljet"]->Fill(nvtxs);

      //fill for selected taus
      double mindr = 0.5; //int jtau = -1;
      for(size_t itau = 0; itau < t_final.size(); itau++){
	int ind_tau = t_final[itau];
	double idr= DeltaR(pfTaus[ind_tau].p4, pfJets[ind_jet].p4);
	if(idr < mindr){
	  mindr = idr;
	  //jtau = ind_tau;
	}
      }
      if(mindr < 0.4){
	histos_["pt_taujet"]->Fill(pfJets[ind_jet].p4.pt());
	histos_["eta_taujet"]->Fill(pfJets[ind_jet].p4.eta());
	histos_["phi_taujet"]->Fill(pfJets[ind_jet].p4.phi());
	histos_["r_taujet"]->Fill(Rjet);
	histos_["nvtx_taujet"]->Fill(nvtxs);
      }
    }
  }
  cout<<" no. of triggered Events "<<nTriggEvent<<endl;

  //get fakerates
  addFakeRate(histos_["pt_taujet"], histos_["pt_alljet"], histos_["pt_fakerate"]);
  addFakeRate(histos_["eta_taujet"], histos_["eta_alljet"], histos_["eta_fakerate"]);
  addFakeRate(histos_["phi_taujet"], histos_["phi_alljet"], histos_["phi_fakerate"]);
  addFakeRate(histos_["r_taujet"], histos_["r_alljet"], histos_["r_fakerate"]);
  addFakeRate(histos_["nvtx_taujet"], histos_["nvtx_alljet"], histos_["nvtx_fakerate"]);

  outFile_->Write();
  outFile_->Close();

}

/////Fake Rate QCD MC ////
void TauFakeRateEstimator::ComputeFakeRateDiJetMC(string myKey, string tauIsoLabel, TString histoFiles_){

  bool isPFlow = (myKey.find("PFlow") != string::npos) ? true : false;
  string eAlgo("PFlow"), mAlgo("PFlow"), jAlgo("PFlow"), tAlgo("PFlow");
  if(myKey.find("PFlow") != string::npos){jAlgo = "PFlow"; tAlgo = "PFlow";}
  else if(myKey.find("HPS") != string::npos){jAlgo = "Jets"; tAlgo = "Taus"; mAlgo = "Muons"; eAlgo = "Electrons";}
								 

  //create the file ////////////////////////////////
  TString Filename_ = histoFiles_ +"_"+myKey+"_"+tauIsoLabel+".root";
  TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
  outFile_->SetCompressionLevel( 9 );
  //////////////////////////////////////////////////

  //define histograms
  DefineHistos();

  evR = new Reader();

  double CrossBin[5] = {815900000., 53120000., 6359000., 784300., 115100.};  //cross section
  double EventsBin[5] = {8213600., 1500000., 1500000., 1500000., 1500000.};
  TString MCFileList[5] = {"file_qcd_P15", "file_qcd_Pt30", "file_qcd_Pt50", "file_qcd_Pt80", "file_qcd_Pt120"};
  int nFiles = 5;
  double Lumi = 0.007653;
  double ProcEvents[5] = {-1, -1, -1, -1, -1};
  
  for(int ifile = 1; ifile <= nFiles; ifile++){
    
    int nEntries = evR->AssignEventTreeFromList(MCFileList[ifile-1]);
    cout << nEntries << " events are available" << endl;
    if( nEntries == 0) {continue; }

    double weight = Lumi*CrossBin[ifile-1]/EventsBin[ifile-1];
    if(ProcEvents[ifile-1] != -1 && nEntries > ProcEvents[ifile-1]){
      weight = (weight*nEntries/ProcEvents[ifile-1]);
    }

    int nTriggEvent = 0;
    int nProcEvents = 0;
    
    for(int i=0; i<nEntries; ++i){
      Long64_t ientry = evR->LoadTree(i);
      if (ientry < 0) break;
      
      nProcEvents++;
      if(ProcEvents[ifile-1] != -1 && nProcEvents > ProcEvents[ifile-1])break;
      
      MyEvent *ev = evR->GetNewEventFromList(i);
      if(ev==0) continue;
      
      //trigger
      bool passTrig = false;
      vector<string> trig = ev->hlt;
      for(size_t it = 0; it < trig.size(); it++){
	if(trig[it].find("Jet") != string::npos) passTrig = true;
      }
      if(!passTrig) continue;
      nTriggEvent++; 
      
      //get MC weight for pileup
      double evtWeight = weight;
      vector<double>puweights = ev->sampleInfo.puWeights;
      if(puweights.size() > 0) evtWeight *= puweights[0];

      //get objects //
      vector<MyVertex> Vertices = ev->PrimaryVtxs;
      if(Vertices.size() <= 0){
	cout<<" no vertexes , exit"<<endl;
	continue;
      }
      vector<MyMuon> pfMuons = evR->getMuons(ev, mAlgo);
      vector<MyElectron> pfElectrons = evR->getElectrons(ev, eAlgo);
      vector<MyTau> pfTaus = evR->getTaus(ev, tAlgo);
      vector<MyJet> pfJets = evR->getJets(ev, jAlgo);
      
      // preselect objects //
      vector<int> m_init; m_init.clear();
      preSelectMuons(&m_init, pfMuons, Vertices[0], isPFlow);
      vector<int> e_init; e_init.clear();
      preSelectElectrons(&e_init, pfElectrons, Vertices[0], isPFlow);
      vector<int> t_init; t_init.clear();
      preSelectTaus( &t_init, pfTaus, TAUPRONGS_, tauIsoLabel, Vertices[0]);
      vector<int> j_init; j_init.clear();
      preSelectJets(jAlgo, &j_init, pfJets);
      
      // clean objects //
      vector<int> e_final; e_final.clear();
      ElectronCleaning( pfElectrons, pfMuons, &e_init, &e_final, &m_init, DRMIN_ELE);
      vector<int> t_final; t_final.clear();
      TauCleaning(pfTaus, pfMuons, pfElectrons, &t_init, &t_final, &m_init, &e_final, DRMIN_TAU);
      vector<int> j_afterLeptonRemoval; j_afterLeptonRemoval.clear();
      JetCleaning(pfJets, pfMuons, pfElectrons, pfTaus, &j_init, &j_afterLeptonRemoval, &m_init, &e_final, &t_final, DRMIN_JET, false);
      
      // remove jets matched to trigger & apply dZ cut and antilepton discr
      vector<int> j_final; j_final.clear();
      int nTrigJet = 0;
      for(size_t ijet = 0; ijet < j_afterLeptonRemoval.size(); ijet++){
	int ind =  j_afterLeptonRemoval[ijet];
	if(pfJets[ind].triggerJet_pt > 40)nTrigJet++;
      }
      double zvertex = Vertices[0].XYZ.z();
      for(size_t ijet = 0; ijet < j_afterLeptonRemoval.size(); ijet++){
	int ind =  j_afterLeptonRemoval[ijet];
	if(nTrigJet == 1 && pfJets[ind].triggerJet_pt > 40) continue; //not to use trigger matched Jet
	if(fabs(zvertex - pfJets[ind].tau_vertex.z()) > 0.2 )continue; //dZ < 0.2
	if(pfJets[ind].tau_againstElectronMedium < 0.5 || pfJets[ind].tau_againstMuonMedium < 0.5)continue;  //discriminator against lepton
	j_final.push_back(ind);
      }
      
      //loop over jets 
      double nvtxs = Vertices.size();
      for(size_t ijet = 0; ijet < j_final.size(); ijet++){
	int ind_jet = j_final[ijet];
	histos_["pt_alljet"]->Fill(pfJets[ind_jet].p4.pt(), evtWeight);
	histos_["eta_alljet"]->Fill(pfJets[ind_jet].p4.eta(), evtWeight);
	histos_["phi_alljet"]->Fill(pfJets[ind_jet].p4.phi(), evtWeight);
	double Rjet = ((pfJets[ind_jet].etaetaMoment + pfJets[ind_jet].phiphiMoment) > 0) ? sqrt(pfJets[ind_jet].etaetaMoment + pfJets[ind_jet].phiphiMoment) : 0.;
	histos_["r_alljet"]->Fill(Rjet, evtWeight);
	histos_["nvtx_alljet"]->Fill(nvtxs, evtWeight);
	
	//fill for selected taus
	double mindr = 0.5; //int jtau = -1;
	for(size_t itau = 0; itau < t_final.size(); itau++){
	  int ind_tau = t_final[itau];
	  double idr= DeltaR(pfTaus[ind_tau].p4, pfJets[ind_jet].p4);
	  if(idr < mindr){
	    mindr = idr;
	    //jtau = ind_tau;
	  }
	}
	if(mindr < 0.4){
	  histos_["pt_taujet"]->Fill(pfJets[ind_jet].p4.pt(), evtWeight);
	  histos_["eta_taujet"]->Fill(pfJets[ind_jet].p4.eta(), evtWeight);
	  histos_["phi_taujet"]->Fill(pfJets[ind_jet].p4.phi(), evtWeight);
	  histos_["r_taujet"]->Fill(Rjet, evtWeight);
	  histos_["nvtx_taujet"]->Fill(nvtxs, evtWeight);
	}
      }
    }
    cout<<" no. of triggered Events "<<nTriggEvent<<endl;
  }

  //get fakerates
  addFakeRate(histos_["pt_taujet"], histos_["pt_alljet"], histos_["pt_fakerate"]);
  addFakeRate(histos_["eta_taujet"], histos_["eta_alljet"], histos_["eta_fakerate"]);
  addFakeRate(histos_["phi_taujet"], histos_["phi_alljet"], histos_["phi_fakerate"]);
  addFakeRate(histos_["r_taujet"], histos_["r_alljet"], histos_["r_fakerate"]);
  addFakeRate(histos_["nvtx_taujet"], histos_["nvtx_alljet"], histos_["nvtx_fakerate"]);

  outFile_->Write();
  outFile_->Close();

}

///W+Jet Fake rate ///////
void TauFakeRateEstimator::ComputeFakeRateWJet(const char* file_list,  string myKey, string tauIsoLabel, bool isData, TString histoFiles_){

  bool isPFlow = (myKey.find("PFlow") != string::npos) ? true : false;
  string eAlgo("PFlow"), mAlgo("PFlow"), jAlgo("PFlow"), tAlgo("PFlow"), metAlgo("PFlow");
  if(myKey.find("PFlow") != string::npos){jAlgo = "PFlow"; tAlgo = "PFlow";}
  else if(myKey.find("HPS") != string::npos){jAlgo = "Jets"; tAlgo = "Taus"; metAlgo = "METsPF"; mAlgo = "Muons"; eAlgo = "Electrons";}
  
  
  evR = new Reader();
  
  int nEntries = evR->AssignEventTreeFromList(file_list);
  cout << nEntries << " events are available" << endl;
  if( nEntries == 0) {return; }
  
  //create the file ////////////////////////////////
  TString Filename_ = histoFiles_ +"_"+myKey+"_"+tauIsoLabel+".root";
  TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
  outFile_->SetCompressionLevel( 9 );
  //////////////////////////////////////////////////

  //define histograms
  DefineHistos();
  histos_["w_mt"] = new TH1D("w_mt", "mt (mu, met)", 200, 0., 200.);

  MyEvent *ev;
  int nTriggEvent = 0, nSelEvents = 0;
  for(int i=0; i<nEntries; ++i){
    Long64_t ientry = evR->LoadTree(i);
    if (ientry < 0) break;

    ev = evR->GetNewEventFromList(i);
    if(ev==0) continue;

    //trigger
    bool passTrig = false;
    vector<string> trig = ev->hlt;
    for(size_t it = 0; it < trig.size(); it++){
      if(trig[it].find("Mu") != string::npos) passTrig = true;
    }
    if(!passTrig) continue;
    nTriggEvent++;
    
    //Apply PU re-weighting
    double evtWeight = 1.0;
    if(!isData){
      vector<double>puweights = ev->sampleInfo.puWeights;
      if(puweights.size() > 0) evtWeight *= puweights[0];
    }

    //get objects //
    vector<MyVertex> Vertices = ev->PrimaryVtxs;
    if(Vertices.size() <= 0){
      cout<<" no vertexes , exit"<<endl;
      continue;
    }
    vector<MyMuon> pfMuons = evR->getMuons(ev, mAlgo);
    vector<MyElectron> pfElectrons = evR->getElectrons(ev, eAlgo);
    vector<MyTau> pfTaus = evR->getTaus(ev, tAlgo);
    vector<MyJet> pfJets = evR->getJets(ev, jAlgo);
    MyMET met = evR->getMET(ev, metAlgo);

    // preselect objects //
    vector<int> m_init; m_init.clear();
    preSelectMuons(&m_init, pfMuons, Vertices[0], isPFlow);
    vector<int> e_init; e_init.clear();
    preSelectElectrons(&e_init, pfElectrons, Vertices[0], isPFlow);
    vector<int> t_init; t_init.clear();
    preSelectTaus( &t_init, pfTaus, TAUPRONGS_, tauIsoLabel, Vertices[0]);
    vector<int> j_init; j_init.clear();
    preSelectJets(jAlgo, &j_init, pfJets);

    // clean objects //
    vector<int> e_final; e_final.clear();
    ElectronCleaning( pfElectrons, pfMuons, &e_init, &e_final, &m_init, DRMIN_ELE);
    vector<int> t_final; t_final.clear();
    TauCleaning(pfTaus, pfMuons, pfElectrons, &t_init, &t_final, &m_init, &e_final, DRMIN_TAU);
    vector<int> j_afterLeptonRemoval; j_afterLeptonRemoval.clear();
    JetCleaning(pfJets, pfMuons, pfElectrons, pfTaus, &j_init, &j_afterLeptonRemoval, &m_init, &e_final, &t_final, DRMIN_JET, false);
    
    // apply dZ cut and antilepton discr
    vector<int> j_final; j_final.clear();
    double zvertex = Vertices[0].XYZ.z();
    for(size_t ijet = 0; ijet < j_afterLeptonRemoval.size(); ijet++){
      int ind =  j_afterLeptonRemoval[ijet];
      if(fabs(zvertex - pfJets[ind].tau_vertex.z()) > 0.2 )continue; //dZ < 0.2
      //if(pfJets[ind].tau_againstElectronMedium < 0.5 || pfJets[ind].tau_againstMuonMedium < 0.5)continue;  //discriminator against lepton
      j_final.push_back(ind);
    }

    //Apply Lepton selection//////////////////////////////
    //int nLepton = e_init.size() + m_init.size();
    int nLepton = m_init.size();  //only muon +jet events
    if(nLepton != 1)continue;
    if(e_final.size() > 0)continue;

    //see if we have loose muons
    if( looseMuonVeto( m_init[0],pfMuons, isPFlow) ) continue;
    // see if we have loose electrons
    if( looseElectronVeto(-1,pfElectrons, isPFlow) ) continue;

    //Apply selection on mT(mu, MET)
    double metPt = met.p4.pt();
    double leptonPt(0), deltaPhi(0);
    if(metPt < 30) continue;
    if( m_init.size() == 1 ){ 
      int m_i = m_init[0]; leptonPt = TMath::Abs(pfMuons[m_i].p4.pt());     
      //deltaPhi = muons[m_i].DeltaPhi( met );     
      deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(pfMuons[m_i].p4, met.p4);
    }
    //if( e_init.size() == 1 ){ int e_i = e_init[0]; leptonPt = TMath::Abs(electrons[e_i].Pt()); deltaPhi = electrons[e_i].DeltaPhi( met);  }
    double mt = sqrt (  2*leptonPt*metPt*(1 - cos(deltaPhi) ) ) ;
    //mon_.fillHisto("w_mt", "fakerate", mt);
    histos_["w_mt"]->Fill(mt, evtWeight);
    if(mt < 50) continue;
    nSelEvents++;
      
    //loop over jets 
    double nvtxs = Vertices.size();
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];

      histos_["pt_alljet"]->Fill(pfJets[ind_jet].p4.pt(), evtWeight);
      histos_["eta_alljet"]->Fill(pfJets[ind_jet].p4.eta(), evtWeight);
      histos_["phi_alljet"]->Fill(pfJets[ind_jet].p4.phi(), evtWeight);
      double Rjet = ((pfJets[ind_jet].etaetaMoment + pfJets[ind_jet].phiphiMoment) > 0) ? sqrt(pfJets[ind_jet].etaetaMoment + pfJets[ind_jet].phiphiMoment) : 0.;
      histos_["r_alljet"]->Fill(Rjet, evtWeight);
      histos_["nvtx_alljet"]->Fill(nvtxs, evtWeight);

      //fill for selected taus
      double mindr = 0.5; //int jtau = -1;
      for(size_t itau = 0; itau < t_final.size(); itau++){
	int ind_tau = t_final[itau];
	double idr= DeltaR(pfTaus[ind_tau].p4, pfJets[ind_jet].p4);
	if(idr < mindr){
	  mindr = idr;
	  //jtau = ind_tau;
	}
      }
      if(mindr < 0.4){
	histos_["pt_taujet"]->Fill(pfJets[ind_jet].p4.pt(), evtWeight);
	histos_["eta_taujet"]->Fill(pfJets[ind_jet].p4.eta(), evtWeight);
	histos_["phi_taujet"]->Fill(pfJets[ind_jet].p4.phi(), evtWeight);
	histos_["r_taujet"]->Fill(Rjet, evtWeight);
	histos_["nvtx_taujet"]->Fill(nvtxs, evtWeight);
      }
    }
  }
  cout<<" no. of triggered Events "<<nTriggEvent<<endl;
  cout<<" no. of selected Events "<<nSelEvents<<endl;

  //get fakerates
  addFakeRate(histos_["pt_taujet"], histos_["pt_alljet"], histos_["pt_fakerate"]);
  addFakeRate(histos_["eta_taujet"], histos_["eta_alljet"], histos_["eta_fakerate"]);
  addFakeRate(histos_["phi_taujet"], histos_["phi_alljet"], histos_["phi_fakerate"]);
  addFakeRate(histos_["r_taujet"], histos_["r_alljet"], histos_["r_fakerate"]);
  addFakeRate(histos_["nvtx_taujet"], histos_["nvtx_alljet"], histos_["nvtx_fakerate"]);

  outFile_->Write();
  outFile_->Close();

}

void TauFakeRateEstimator::addFakeRate(TH1* hn, TH1* hd, TH1* hf){
  
  TH1D * n; TH1 *d;
  n = (TH1D*)hn->Clone();
  d = (TH1D*)hd->Clone();

  int nBinsX1 = n->GetNbinsX(); int nBinsX2 = d->GetNbinsX();

  if( nBinsX1 != nBinsX2 ){
    //cout<<endl<<endl<<"H1 N bins is = "<<nBinsX1<<" but H2 N bin is = "<<nBinsX2<<" ... quit! "<<endl;
    return;
  }

  for(int i=1; i<= nBinsX1;i++){
    double numContent = n->GetBinContent(i); 
    double denContent = d->GetBinContent(i);
    double err = getErrorFraction(numContent,denContent);
    //cout<<endl<<" i "<<i<<" num = "<<numContent<<" den = "<<denContent<<endl;

    if(denContent){ hf->SetBinContent(i, numContent/denContent); hf->SetBinError(i,err); } 
    
  }

  delete n; delete d;

}
  
void TauFakeRateEstimator::addFakeRate2D(TH2* hn, TH2* hd, TH2* hf){

  TH2D * n; TH2D *d;
  n = (TH2D*)hn->Clone();
  d = (TH2D*)hd->Clone();

  int nBinsX1 = n->GetNbinsX(); int nBinsX2 = d->GetNbinsX();
  int nBinsY1 = n->GetNbinsY(); int nBinsY2 = d->GetNbinsY();
  
  if( (nBinsX1 != nBinsX2) || (nBinsY1 != nBinsY2)){
    //cout<<endl<<endl<<"H1 N bins is = "<<nBinsX1<<" but H2 N bin is = "<<nBinsX2<<" ... quit! "<<endl;
    return;
  }

  for(int i=1; i<= nBinsX1;i++){
    for(int j=1; j<=nBinsY1;j++){
      double numContent = n->GetBinContent(i,j);
      double denContent = d->GetBinContent(i,j);
      double err = getErrorFraction(numContent,denContent);
      //cout<<endl<<" i "<<i<<" num = "<<numContent<<" den = "<<denContent<<endl;
      
      if(denContent){ hf->SetBinContent(i, j, numContent/denContent); hf->SetBinError(i,j,err); }
    }
  }

  delete n; delete d;

}

double TauFakeRateEstimator::getErrorFraction( double a,double b){
  double ret(0);
  if (b){
    double temp = fabs(a)/(b*b)+ fabs( (a*a)/(b*b*b)) ;
    ret = sqrt(temp);
  }
  return ret;
};

void TauFakeRateEstimator::DefineHistos()
{

  double ptBin[16] = {0, 10, 20, 22.5, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 150, 250};

  histos_["pt_taujet"] = new TH1D("pt_taujet", "pt_taujet", 15, ptBin);
  histos_["pt_alljet"] = new TH1D("pt_alljet", "pt_alljet", 15, ptBin);
  histos_["pt_fakerate"] = new TH1D("pt_fakerate", "pt_fakerate", 15, ptBin);
  histos_["eta_taujet"] = new TH1D("eta_taujet", "eta_taujet", 10, -2.5, 2.5);
  histos_["eta_alljet"] = new TH1D("eta_alljet", "eta_alljet", 10, -2.5, 2.5);
  histos_["eta_fakerate"] = new TH1D("eta_fakerate", "eta_fakerate", 10, -2.5, 2.5);
  histos_["phi_taujet"] = new TH1D("phi_taujet", "phi_taujet", 14, -3.5, 3.5);
  histos_["phi_alljet"] = new TH1D("phi_alljet", "phi_alljet", 14, -3.5, 3.5);
  histos_["phi_fakerate"] = new TH1D("phi_fakerate", "phi_fakerate", 14, -3.5, 3.5);
  histos_["r_taujet"] = new TH1D("r_taujet", "r_taujet", 8, 0., 0.32);
  histos_["r_alljet"] = new TH1D("r_alljet", "r_alljet", 8, 0., 0.32);
  histos_["r_fakerate"] = new TH1D("r_fakerate", "r_fakerate", 8, 0., 0.32);
  histos_["nvtx_taujet"] = new TH1D("nvtx_taujet", "nvtx_taujet", 30, 0, 30.);
  histos_["nvtx_alljet"] = new TH1D("nvtx_alljet", "nvtx_alljet", 30, 0, 30.);
  histos_["nvtx_fakerate"] = new TH1D("nvtx_fakerate", "nvtx_fakerate", 30, 0, 30.);
}


void TauFakeRateEstimator::processEvents()
{
  
  ComputeFakeRateDiJetData("file_jet_dcs", "HPS", "Loose", "TFR_DiJet_data2012_w1_dcs"); 
  ComputeFakeRateDiJetData("file_jet_dcs", "HPS", "Medium", "TFR_DiJet_data2012_w1_dcs"); 
  ComputeFakeRateDiJetData("file_jet_dcs", "HPS", "Tight", "TFR_DiJet_data2012_w1_dcs"); 
  ComputeFakeRateDiJetData("file_jet_dcs", "HPS", "LooseDB", "TFR_DiJet_data2012_w1_dcs"); 
  ComputeFakeRateDiJetData("file_jet_dcs", "HPS", "MediumDB", "TFR_DiJet_data2012_w1_dcs"); 
  ComputeFakeRateDiJetData("file_jet_dcs", "HPS", "TightDB", "TFR_DiJet_data2012_w1_dcs"); 
  ComputeFakeRateDiJetData("file_jet_dcs", "HPS", "LooseCombDB", "TFR_DiJet_data2012_w1_dcs"); 
  ComputeFakeRateDiJetData("file_jet_dcs", "HPS", "MediumCombDB", "TFR_DiJet_data2012_w1_dcs"); 
  ComputeFakeRateDiJetData("file_jet_dcs", "HPS", "TightCombDB", "TFR_DiJet_data2012_w1_dcs");

  /*
  ComputeFakeRateDiJetMC("HPS", "Loose", "TFR_DiJet_mc");
  ComputeFakeRateDiJetMC("HPS", "Medium", "TFR_DiJet_mc");
  ComputeFakeRateDiJetMC("HPS", "Tight", "TFR_DiJet_mc");
  ComputeFakeRateDiJetMC("HPS", "LooseDB", "TFR_DiJet_mc");
  ComputeFakeRateDiJetMC("HPS", "MediumDB", "TFR_DiJet_mc");
  ComputeFakeRateDiJetMC("HPS", "TightDB", "TFR_DiJet_mc");
  ComputeFakeRateDiJetMC("HPS", "LooseCombDB", "TFR_DiJet_mc");
  ComputeFakeRateDiJetMC("HPS", "MediumCombDB", "TFR_DiJet_mc");
  ComputeFakeRateDiJetMC("HPS", "TightCombDB", "TFR_DiJet_mc");
  */
  //ComputeFakeRateWJet("file_muon", "HPS", "Loose", true, "TFR_WJet_data");
  //ComputeFakeRateWJet("file_muon", "HPS", "Medium", true, "TFR_WJet_data");
  //ComputeFakeRateWJet("file_muon", "HPS", "Tight", true, "TFR_WJet_data");
  //ComputeFakeRateWJet("file_muon", "HPS", "LooseDB", true, "TFR_WJet_data");
  //ComputeFakeRateWJet("file_muon", "HPS", "MediumDB", true, "TFR_WJet_data");
  //ComputeFakeRateWJet("file_muon", "HPS", "TightDB", true, "TFR_WJet_data");
  //ComputeFakeRateWJet("file_isomuon", "HPS", "LooseCombDB", true, "TFR_WJet_data2012-w1");
  //ComputeFakeRateWJet("file_isomuon", "HPS", "MediumCombDB", true, "TFR_WJet_data2012-w1");
  //ComputeFakeRateWJet("file_isomuon", "HPS", "TightCombDB", true, "TFR_WJet_data2012-w1");
  
  ComputeFakeRateWJet("file_isomuon_dcs", "HPS", "Loose", true, "TFR_WJet_data2012-w1-dcs");
  ComputeFakeRateWJet("file_isomuon_dcs", "HPS", "Medium", true, "TFR_WJet_data2012-w1-dcs");
  ComputeFakeRateWJet("file_isomuon_dcs", "HPS", "Tight", true, "TFR_WJet_data2012-w1-dcs");
  ComputeFakeRateWJet("file_isomuon_dcs", "HPS", "LooseDB", true, "TFR_WJet_data2012-w1-dcs");
  ComputeFakeRateWJet("file_isomuon_dcs", "HPS", "MediumDB", true, "TFR_WJet_data2012-w1-dcs");
  ComputeFakeRateWJet("file_isomuon_dcs", "HPS", "TightDB", true, "TFR_WJet_data2012-w1-dcs");

  ComputeFakeRateWJet("file_isomuon_dcs", "HPS", "LooseCombDB", true, "TFR_WJet_data2012-w1-dcs");
  ComputeFakeRateWJet("file_isomuon_dcs", "HPS", "MediumCombDB", true, "TFR_WJet_data2012-w1-dcs");
  ComputeFakeRateWJet("file_isomuon_dcs", "HPS", "TightCombDB", true, "TFR_WJet_data2012-w1-dcs");
  
}
