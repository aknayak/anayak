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
#include "TTimeStamp.h"

#include "MiniTree/Selection/interface/Reader.h"
//#include "Reader.C"

using namespace std;

class selectEvent 
{
public :
  selectEvent () { };
  ~selectEvent () { };

  void checkEvent(TString url);
  void checkEvent1(const char*);

private :
  Reader *evR;
};

void selectEvent::checkEvent(TString url)
{
  evR = new Reader();
  TFile *f = TFile::Open(url);
  if(f==0) return ;
  if(f->IsZombie()) { f->Close(); return; }  

  int nEntries = evR->AssignEventTreeFrom(f);
  cout << nEntries << " events are available" << endl;
  if( nEntries == 0) { f->Close(); return; }
  
  for(int i=0; i<nEntries; ++i){
    MyEvent *ev = evR->GetNewEvent(i);
    if(ev==0) continue;

    vector<MyJet> allJets = ev->Jets;
    cout<<" nJets : "<<allJets.size()<<endl;
    for(size_t ijet=0; ijet < allJets.size(); ++ijet){
      if(allJets[ijet].jetName.find("PFlow") == std::string::npos)continue;
      cout<<" jet pT "<<allJets[ijet].p4.pt()<<endl;
    }
    vector<MyElectron> allElectrons = ev->Electrons;
    cout<<" nElectrons : "<<allElectrons.size()<<endl;
    for(size_t iele=0; iele < allElectrons.size(); ++iele){
      if(allElectrons[iele].GetName().find("PFlow") == std::string::npos)continue;
      cout<<" electron pT "<<allElectrons[iele].p4.pt()<<endl;
    }
    vector<MyMuon> allMuons = ev->Muons;
    cout<<" nMuons : "<<allMuons.size()<<endl;
    for(size_t imu=0; imu < allMuons.size(); ++imu){
      if(allMuons[imu].GetName().find("PFlow") == std::string::npos)continue;
      cout<<" muon pT "<<allMuons[imu].p4.pt()<<endl;
    }
    
    if(!(ev->isData)){
      vector<MyMCParticle>allMCParticles = ev->mcParticles;
      cout<<" nMCParticles : "<<allMCParticles.size()<<endl;
      for(size_t imc=0; imc < allMCParticles.size(); ++imc){
	cout<<" Part Type "<<allMCParticles[imc].pid<<" status "<<allMCParticles[imc].status<<" Pt "<<allMCParticles[imc].p4Gen.pt()<<" motherline "<<allMCParticles[imc].motherLine<<endl;
	if(allMCParticles[imc].mother.size() > 0)cout<<" mother id "<<allMCParticles[imc].mother[0]<<endl;
      }
    }

  }
}

void selectEvent::checkEvent1(const char *file_list)
{
  evR = new Reader();
  
  int nEntries = evR->AssignEventTreeFromList(file_list);
  cout << nEntries << " events are available" << endl;
  if( nEntries == 0) {return; }

  MyEvent *ev;
  int npassed = 0;
  for(int i=0; i<nEntries; ++i){
    Long64_t ientry = evR->LoadTree(i);
    if (ientry < 0) break;
    
    ev = evR->GetNewEventFromList(i);
    if(ev==0) continue;

    vector<MyJet> allJets = ev->Jets;
    cout<<" nJets : "<<allJets.size()<<endl;
    for(size_t ijet=0; ijet < allJets.size(); ++ijet){
      if(allJets[ijet].jetName.find("PFlow") == std::string::npos)continue;
      cout<<" jet pT "<<allJets[ijet].p4.pt()<<endl;
    }
    vector<MyElectron> allElectrons = ev->Electrons;
    cout<<" nElectrons : "<<allElectrons.size()<<endl;
    for(size_t iele=0; iele < allElectrons.size(); ++iele){
      if(allElectrons[iele].GetName().find("PFlow") == std::string::npos)continue;
      cout<<" electron pT "<<allElectrons[iele].p4.pt()<<endl;
    }
    vector<MyMuon> allMuons = ev->Muons;
    cout<<" nMuons : "<<allMuons.size()<<endl;
    for(size_t imu=0; imu < allMuons.size(); ++imu){
      if(allMuons[imu].GetName().find("PFlow") == std::string::npos)continue;
      cout<<" muon pT "<<allMuons[imu].p4.pt()<<endl;
    }
    
    if(!(ev->isData)){
      vector<MyMCParticle>allMCParticles = ev->mcParticles;
      cout<<" nMCParticles : "<<allMCParticles.size()<<endl;
      for(size_t imc=0; imc < allMCParticles.size(); ++imc){
	cout<<" Part Type "<<allMCParticles[imc].pid<<" status "<<allMCParticles[imc].status<<" Pt "<<allMCParticles[imc].p4Gen.pt()<<" motherline "<<allMCParticles[imc].motherLine<<endl;
	if(allMCParticles[imc].mother.size() > 0)cout<<" mother id "<<allMCParticles[imc].mother[0]<<endl;
      }
    }
    
    npassed++;
  }
  delete ev;
  cout<<"Total events "<<npassed<<endl;
}
