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
#include "MiniTree/Selection/interface/MomentumVec.h"
#include "LumiReweightingStandAlone.h"
#include "MiniTree/Selection/interface/UncertaintyComputer.hh"
#include "MiniTree/Selection/interface/HistogramPlotter.hh"

using namespace std;

class HplusAnalyzer : public ObjectSelector, HistogramPlotter
{
public :
  HplusAnalyzer() : ObjectSelector(), HistogramPlotter()
  {
    DRMIN_JET = 0.5;
    DRMIN_ELE = 0.5;
    METCUT_   = 30.0;
    /*
    std::vector<float> MCPUDist, DataPUDist;
    float dataDist_2012A[60] = {1.82661e-05,0.000124476,0.000494428,0.00141418,0.00322671,0.0062338,0.0105956,0.0162651,0.0229752,0.0302751,0.0376034,0.0443794,0.0500918,0.0543697,0.0570214,0.0580386,0.057568,0.0558629,0.0532267,0.0499626,0.0463369,0.0425599,0.0387819,0.0351004,0.0315723,0.0282278,0.0250814,0.0221399,0.0194073,0.0168863,0.0145789,0.0124856,0.0106047,0.0089316,0.00745897,0.00617645,0.00507131,0.00412902,0.0033339,0.00266978,0.0021206,0.00167088,0.00130612,0.00101301,0.000779625,0.000595444,0.000451357,0.000339596,0.000253632,0.000188051,0,0,0,0,0,0,0,0,0,0};
    
    float mc_Dist[60] = {0, 0, 0.0005, 0.000375, 0.003, 0.004125, 0.006, 0.01, 0.013125, 0.01725, 0.021625, 0.025125, 0.02725, 0.03225, 0.03875, 0.042125, 0.04125, 0.045625, 0.046, 0.047375, 0.046, 0.0465, 0.041875, 0.042125, 0.040875, 0.03325, 0.03625, 0.0365, 0.03075, 0.027125, 0.02425, 0.019375, 0.02075, 0.018, 0.01925, 0.0135, 0.012125, 0.01275, 0.00925, 0.00775, 0.006, 0.0055, 0.005375, 0.004375, 0.0035, 0.00325, 0.0025, 0.001625, 0.0015, 0.000875, 0.000875, 0.0015, 0.00075, 0.000125, 0.00025, 0.000125, 0.000125, 0, 0.00025, 0.000375};
    //ttjet
    for(int i = 0; i < 60; i++){
      DataPUDist.push_back(dataDist_2012A[i]);
      MCPUDist.push_back(mc_Dist[i]);
    }
    LumiWeights_ = reweight::LumiReWeighting(MCPUDist, DataPUDist);
    */

    LumiWeights_ = reweight::LumiReWeighting("MC_Pileup_Summer2012_600bins.root","Data_Pileup_2012ABC_600bins.root", "pileup", "pileup"); //crashes, don't know why
    
    PShiftDown_ = reweight::PoissonMeanShifter(-0.5);
    PShiftUp_ = reweight::PoissonMeanShifter(0.5);

    //cross sections
    xss["WJETS"] = 36257.0; 
    xss["TTBAR"] = 225.2; 
    xss["ZJETS"] = 3504.0; 
    xss["QCD"]   = 134680; 
    xss["TOPT"]  = 11.1; 
  };
  ~HplusAnalyzer() {
    delete evR;
  };
  
  void CutFlowAnalysis(TString url,  string myKey="PFlow", bool isData = true, string evtType="data");
  void CutFlowProcessor(TString url,  string myKey="PFlow", TString cutflowType="base", bool isData = true, TFile *outFile_=0);
  void CreateAnalHistos(TString flowType, TFile* outFile_);
  void processEvents();

private :
  double DRMIN_JET, DRMIN_ELE, METCUT_;
  Reader *evR;
  
  reweight::LumiReWeighting LumiWeights_;
  reweight::PoissonMeanShifter PShiftUp_; //pileup syst up
  reweight::PoissonMeanShifter PShiftDown_; //pileup syst down
  std::map<string, double> xss;
  ofstream outfile_;
};

void HplusAnalyzer::CutFlowAnalysis(TString url, string myKey, bool isData, string evtType){

  //create the output file ////////////////////////////////
  TString Filename_ = evtType+"_selection.root";
  TFile *outFile_ = TFile::Open( Filename_, "RECREATE" );
  outFile_->SetCompressionLevel( 9 );
  
  TString debug_Filename_ = Filename_+"_debug.txt";
  string debug_file(debug_Filename_);
  outfile_.open(debug_file.c_str());

  CutFlowProcessor(url, myKey, "base", isData, outFile_);
  //CutFlowProcessor(url, myKey, "JESPlus", isData, outFile_);
  //CutFlowProcessor(url, myKey, "JESMinus", isData, outFile_);
  //CutFlowProcessor(url, myKey, "JERPlus", isData, outFile_);
  //CutFlowProcessor(url, myKey, "JERMinus", isData, outFile_);
  //CutFlowProcessor(url, myKey, "METUCPlus", isData, outFile_);
  //CutFlowProcessor(url, myKey, "METUCMinus", isData, outFile_);
  //CutFlowProcessor(url, myKey, "bTagPlus", isData, outFile_);
  //CutFlowProcessor(url, myKey, "bTagMinus", isData, outFile_);
  //CutFlowProcessor(url, myKey, "PUPlus", isData, outFile_); 
  //CutFlowProcessor(url, myKey, "PUMinus", isData, outFile_);
  
  outfile_.close();
  outFile_->Write(); 
  outFile_->Close(); 
}
  
void HplusAnalyzer::CutFlowProcessor(TString url,  string myKey, TString cutflowType, bool isData, TFile *outFile_){

  outfile_<<"///// Begin processing "<<cutflowType<<"  selection ///////"<<endl;

  bool isPFlow = (myKey.find("PFlow") != string::npos) ? true : false;

  string eAlgo("Electrons"), mAlgo("Muons"), jAlgo("Jets"), metAlgo("METsPF");
  //  if(myKey.find("PFlow") != string::npos)jAlgo = "PFlow";
  //  else if(myKey.find("HPS") != string::npos)jAlgo = "JetsAK5PF";
  
  //Uncertainty variations, JES, JER, MET unclustered, bTag
  int jes = 0, jer = 0, metuc = 0, bscale = 0;
  if(cutflowType.Contains("JESPlus"))jes = 1;
  else if (cutflowType.Contains("JESMinus"))jes = -1;
  else if (cutflowType.Contains("JERPlus"))jer = 1;
  else if (cutflowType.Contains("JERMinus"))jer = -1;
  else if (cutflowType.Contains("METUCPlus"))metuc = 1;
  else if (cutflowType.Contains("METUCMinus"))metuc = -1;
  else if (cutflowType.Contains("bTagPlus"))bscale = 1;
  else if (cutflowType.Contains("bTagMinus"))bscale = -1; 

  double Lumi = 1000.0;

  evR = new Reader();

  TFile *f = TFile::Open(url);
  if(f==0) return ;
  if(f->IsZombie()) { f->Close(); return; }

  int nEntries = evR->AssignEventTreeFrom(f);
  if( nEntries == 0) {return; }

  //get initial number of events
  TH1F* inputcf = (TH1F*)(f->Get("allEventsFilter/totalEvents"))->Clone("inputcf");
  double initialEvents = inputcf->GetBinContent(1);
  outfile_<<"Input file : "<<url<<endl;
  outfile_<<"Available input sample : "<<initialEvents<<endl;
  
  //define histograms 
  CreateAnalHistos(cutflowType, outFile_);

  double sampleWeight(1.0);

  MyEvent *ev;
  int nTriggEvent = 0, nSelEvents = 0, matchjetcount= 0, threepairjet = 0;
  double nVerticesFailCount = 0.0;
  for(int i=0; i<nEntries; ++i){
    Long64_t ientry = evR->LoadTree(i);
    if (ientry < 0) break;

    ev = evR->GetNewEvent(i);
    if(ev==0) continue;
    
    // apply PU re-weighting

    double evtWeight = 1.0;
    if(!isData){
      //get sample information
      if(i < 1){
        string sampleName = ev->sampleInfo.sampleName;
        sampleWeight = xss[sampleName] * Lumi / initialEvents;
        outfile_<<"Scale factor for lumi "<<Lumi<<" pb is "<< sampleWeight<<endl;
      }
      evtWeight *= sampleWeight; // upto this only sigma*lumi weight is applied
      
      //vector<double>puweights = ev->sampleInfo.puWeights;
      //if(puweights.size() > 0) evtWeight *= puweights[0];

      //vector<double>pu = ev->sampleInfo.pileup;
      vector<double>pu = ev->sampleInfo.truepileup;
      if(pu.size() > 0) {
	//int npu = pu[0];
	float npu = pu[0];
	//double weight = LumiWeights_.ITweight(npu);
	double weight = LumiWeights_.weight(npu);
	//	cout<<"pu weight"<<weight<<endl;
	if (cutflowType.Contains("PUPlus"))weight = weight*PShiftUp_.ShiftWeight( npu ); 
	else if (cutflowType.Contains("PUMinus"))weight = weight*PShiftDown_.ShiftWeight( npu );
	evtWeight *= weight;
	cout<<" pileup weight "<<weight<<endl;
      }
      
    }
    double nCutPass = 0.0;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);

    //trigger
    bool passTrig = false;
    vector<string> trig = ev->hlt;
    for(size_t it = 0; it < trig.size(); it++){
      if(trig[it].find("Mu") != string::npos) passTrig = true;
    }
    if(!passTrig){
      //      cout << "not satisfying trigger" << endl;
      continue;
    }
    nTriggEvent++;

    //get objects //
    vector<MyVertex> Vertices = ev->PrimaryVtxs;
    
    if(Vertices.size() <= 0){
      nVerticesFailCount+=evtWeight;
      cout<<" no vertexes , exit"<<endl;
      continue;
    }
    
    vector<MyMuon> pfMuons = evR->getMuons(ev, mAlgo);
    vector<MyElectron> pfElectrons = evR->getElectrons(ev, eAlgo);
    vector<MyJet> pfJets = evR->getJets(ev, jAlgo);
    MyMET met = evR->getMET(ev, metAlgo);
    vector<MyTau>pfTaus; pfTaus.clear(); 

    // preselect objects 
    vector<int> m_init; m_init.clear();
    preSelectMuons(&m_init, pfMuons, Vertices[0], isPFlow);
    vector<int> e_init; e_init.clear();
    preSelectElectrons(&e_init, pfElectrons, Vertices[0], isPFlow);
    vector<int> j_init; j_init.clear();
    preSelectJets(jAlgo, &j_init, pfJets, jes, jer);
    
    // clean objects //
    vector<int> e_final; e_final.clear();
    ElectronCleaning( pfElectrons, pfMuons, &e_init, &e_final, &m_init, DRMIN_ELE);
    vector<int> j_final; j_final.clear();
    vector<int> t_final; t_final.clear();
    JetCleaning(pfJets, pfMuons, pfElectrons, pfTaus, &j_init, &j_final, &m_init, &e_final, &t_final, DRMIN_JET, false);

    //Apply Lepton selection//////////////////////////////
    
    int nLepton = m_init.size();  // this condition proof that only muon + jet events
    if(nLepton != 1)continue;
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight); 

    if( looseMuonVeto( m_init[0],pfMuons, isPFlow) ) continue;
    if( looseElectronVeto(-1,pfElectrons, isPFlow) ) continue;
    int m_i = m_init[0];
    
    // Fill histogram after trigger and one offline isolated muon and applied 2nd lepton veto
    fillHisto("pt_mu", cutflowType, pfMuons[m_i].p4.pt(), evtWeight); 
    fillHisto("eta_mu", cutflowType, pfMuons[m_i].p4.eta(), evtWeight);
    fillHisto("phi_mu", cutflowType, pfMuons[m_i].p4.phi(), evtWeight);
    // vertex just after one lepton selection
    double pri_vtxs = Vertices.size();
    fillHisto("nvtx", cutflowType, pri_vtxs, evtWeight);

    // Apply Jet Selection
    int nJet = j_final.size();
    fillHisto("multi_jet", cutflowType, nJet, evtWeight);
    if(nJet < 4)continue;  // this condition implies event should contains at least 4 jets
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);

    int threejet = 0;
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      double jetPt = jetPtWithJESJER(pfJets[ind_jet], jes, jer); 
      fillHisto("pt_jet", cutflowType, jetPt, evtWeight);
      fillHisto("eta_jet", cutflowType, pfJets[ind_jet].p4.eta(), evtWeight);
      fillHisto("phi_jet", cutflowType, pfJets[ind_jet].p4.phi(), evtWeight);
      if(jetPt > 30 )threejet++;
    }
    if(threejet < 3) continue; // three jets should have pt > 30 GeV
    
    // Met distribution   
    double   leptonPt(0), deltaPhi(0);
    double metPt = 0; //met.p4.pt();
    if(!metuc)metPt = metWithJESJER(pfJets, &j_final, met, jes, jer);
    else metPt = metWithUncl(pfJets, &j_final, pfMuons, &m_init, pfElectrons, &e_final, met, metuc);
    fillHisto("pt_met", cutflowType, metPt, evtWeight);

    if(metPt < 30) continue;  // Missing transverse energy cut 30 GeV
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight);

    int m_j = m_init[0];
    leptonPt = TMath::Abs(pfMuons[m_j].p4.pt());
    deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(pfMuons[m_j].p4, met.p4);
    double mt = sqrt (  2*leptonPt*metPt*(1 - cos(deltaPhi) ) ) ;
    fillHisto("wmt", cutflowType, mt, evtWeight);

    // Here putting the criteria that event should contain at least one b-tagged jet  
    int count = 0; 
    for(size_t ijet = 0; ijet < j_final.size(); ijet++){
      int ind_jet = j_final[ijet];
      //if(pfJets[ind_jet].bDiscriminator["combinedSecondaryVertexBJetTags"] > 0.679 && abs(pfJets[ind_jet].p4.eta()) < 2.4) count++;
      bool isBtag = getBtagWithSF(pfJets[ind_jet], isData, bscale, true); 
      if(isBtag)count++;
      fillHisto("btag_jet", cutflowType, pfJets[ind_jet].bDiscriminator["combinedSecondaryVertexBJetTags"], evtWeight); 
    }
    if(count <= 0) continue;
    fillHisto("btagmulti_jet", cutflowType, count, evtWeight);
    nCutPass++;
    fillHisto("cutflow", cutflowType, nCutPass, evtWeight); 
    nSelEvents++;  // this is the counter for counting final number of events  
    
  }

  outfile_ << "Number of times HadP and HadQ matches with jets" << matchjetcount << endl;
  outfile_ << "total number of selected events " << nSelEvents <<endl; 
  outfile_ << "No of times three pair jet matched:    " << threepairjet << endl;

  f->Close(); 
  delete f;
}

void HplusAnalyzer::CreateAnalHistos(TString cutflowType, TFile* outFile_)
{

  //Define Histograms 
  InitHist(cutflowType, "", outFile_); 
  addHisto("cutflow", cutflowType, 10, 0., 10.); 
  addHisto("wmt", cutflowType, 100, 0., 200.);
  addHisto("nvtx", cutflowType, 50, 0., 50.);
  //InitHist("OneLep", cutflowType, outFile_); 
  //InitHist("OneLep2J", cutflowType, outFile_); 
  //InitHist("OneLep2J1B", cutflowType, outFile_); 

}

void HplusAnalyzer::processEvents(){ 
  //CutFlowAnalysis("/tmp/gkole/8TeV/53x/wjet/wjet_su12_all.root", "PF",false, "wjet"); 
  //CutFlowAnalysis("/tmp/gkole/8TeV/53x/ttjet/ttjet_su12_all.root", "PF",false, "ttbar");
  CutFlowAnalysis("ttjet_all.root", "PF",false, "ttbar");
  //CutFlowAnalysis("/afs/cern.ch/work/a/anayak/CMS/HPlusCSbar/CMSSW_5_3_5/src/MiniTree/Selection/test/mc_tau.root", "PF",false, "ttbar");
} 
